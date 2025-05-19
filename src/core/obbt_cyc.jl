### Optimality-Based Bound Tightening for Optimal Power Flow Relaxations
# Run either with or without cycles.
### Adapted from the PowerModels.jl pacakge.
const _LOGGER = Memento.getlogger(@__MODULE__)

function run_cyc_obbt_opf!(data::Dict{String,<:Any}, model;
    max_iter::Int = 100, # 100
    time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf,
    upper_bound_constraint::Bool = false,
    rel_gap_tol::Float64 = Inf,
    min_bound_width::Float64 = 1e-2,
    improvement_tol::Float64 = 1e-3,
    precision::Int = 4,
    termination::Symbol = :avg,
    kwargs...)

    nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
    result = solve_ac_opf(data, nlp_solver)
    upper_bound = result["objective"]
    # println("upper bound: ", upper_bound)

    Memento.info(_LOGGER, "maximum OBBT iterations set to default value of $max_iter")
    Memento.info(_LOGGER, "maximum time limit for OBBT set to default value of $time_limit seconds")

    # model_relaxation = instantiate_model(data, model_type, PowerModels.build_opf)
    # (_IM.ismultinetwork(model_relaxation)) && (Memento.error(_LOGGER, "OBBT is not supported for multi-networks"))
    # (ismulticonductor(model_relaxation)) && (Memento.error(_LOGGER, "OBBT is not supported for multi-conductor networks"))
    #
    # # check for model_type compatability with OBBT
    # _check_variables(model_relaxation)

    # check for other keyword argument consistencies
    _check_obbt_options(upper_bound, rel_gap_tol, upper_bound_constraint)

    # check termination norm criteria for obbt
    (termination != :avg && termination != :max && termination != :gop) && (Memento.error(_LOGGER, "OBBT termination criteria can only be :max or :avg or :gop"))

    # pass status
    status_pass = [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]

    # compute initial relative gap between relaxation objective and upper_bound
    JuMP.optimize!(model)
    current_relaxation_objective = objective_value(model)
    if upper_bound < current_relaxation_objective
        Memento.error(_LOGGER, "the upper bound provided to OBBT is not a valid ACOPF upper bound")
    end
    # if !(result_relaxation["termination_status"] in status_pass)
    #     Memento.warn(_LOGGER, "initial relaxation solve status is $(result_relaxation["termination_status"])")
    #     if result_relaxation["termination_status"] == :SubOptimal
    #         Memento.warn(_LOGGER, "continuing with the bound-tightening algorithm")
    #     end
    # end
    current_rel_gap = Inf
    if !isinf(upper_bound)
        current_rel_gap = (upper_bound - current_relaxation_objective)/upper_bound
        Memento.info(_LOGGER, "Initial relaxation gap = $current_rel_gap")
    end

    model = _build_model(data, model)
    # print(model)
    if upper_bound_constraint
        @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] + gen["cost"][3] for (i,gen) in ref[:gen]) <= upper_bound)
    end
    # (upper_bound_constraint) && (_constraint_obj_bound(model_bt, upper_bound))

    stats = Dict{String,Any}()
    # stats["model_type"] = model_type
    stats["initial_relaxation_objective"] = current_relaxation_objective
    stats["initial_rel_gap_from_ub"] = current_rel_gap
    stats["upper_bound"] = upper_bound

    # vm = var(model_bt, :vm)
    # td = var(model_bt, :td)
    # buses = ids(model_bt, :bus)
    # buspairs = ids(model_bt, :buspairs)
    buses = keys(ref[:bus])
    buspairs = keys(ref[:buspairs])

    vm_lb = Dict{Any,Float64}( [bus => JuMP.lower_bound(vm[bus]) for bus in buses] )
    vm_ub = Dict{Any,Float64}( [bus => JuMP.upper_bound(vm[bus]) for bus in buses] )
    td_lb = Dict{Any,Float64}( [bp => JuMP.lower_bound(td[bp]) for bp in buspairs] )
    td_ub = Dict{Any,Float64}( [bp => JuMP.upper_bound(td[bp]) for bp in buspairs] )

    vm_range_init = sum([vm_ub[bus] - vm_lb[bus] for bus in buses])
    stats["vm_range_init"] = vm_range_init
    stats["avg_vm_range_init"] = vm_range_init/length(buses)

    td_range_init = sum([td_ub[bp] - td_lb[bp] for bp in buspairs])
    stats["td_range_init"] = td_range_init
    stats["avg_td_range_init"] = td_range_init/length(buspairs)

    vm_range_final = 0.0
    td_range_final = 0.0

    total_vm_reduction = Inf
    max_vm_reduction = Inf
    avg_vm_reduction = Inf

    total_td_reduction = Inf
    max_td_reduction = Inf
    avg_td_reduction = Inf

    final_relaxation_objective = NaN

    current_iteration = 0
    time_elapsed = 0.0
    parallel_time_elapsed = 0.0

    bool_td_bound_reduction = true
    bool_td_range = true

    bool_vm_bound_reduction = true
    bool_vm_range = true

    check_termination = true
    (termination == :avg) && (check_termination = (avg_vm_reduction > improvement_tol || avg_td_reduction > improvement_tol))
    (termination == :max) && (check_termination = (max_vm_reduction > improvement_tol || max_td_reduction > improvement_tol))
    (termination == :gop) && (check_termination = ((bool_td_bound_reduction || bool_vm_bound_reduction) && (bool_td_range || bool_vm_range))) # The termination conditions (11) and (12) from Gaopinath et. al. 2020.

    while check_termination
        iter_start_time = time()
        total_vm_reduction = 0.0
        avg_vm_reduction = 0.0
        max_vm_reduction = 0.0
        max_vm_iteration_time = 0.0
        bool_vm_bound_reduction = false
        bool_vm_range = false

        total_td_reduction = 0.0
        avg_td_reduction = 0.0
        max_td_reduction = 0.0
        max_td_iteration_time = 0.0
        bool_td_bound_reduction = false
        bool_td_range = false

        # bound-tightening for the vm variables
        for bus in buses
            (vm_ub[bus] - vm_lb[bus] < min_bound_width) && (continue)

            start_time = time()
            # vm lower bound solve
            lb = NaN
            @objective(model, Min, model[:vm][bus])
            JuMP.optimize!(model)
            if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                (nlb > vm_lb[bus]) && (lb = nlb)
            else
                Memento.warn(_LOGGER, "BT minimization problem for vm[$bus] errored - change tolerances.")
                continue
            end

            #vm upper bound solve
            ub = NaN
            @objective(model, Max, model[:vm][bus])
            JuMP.optimize!(model)
            if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                (nub < vm_ub[bus]) && (ub = nub)
            else
                Memento.warn(_LOGGER, "BT maximization problem for vm[$bus] errored - change tolerances.")
                continue
            end
            end_time = time() - start_time
            max_vm_iteration_time = max(end_time, max_vm_iteration_time)

            # sanity checks
            (lb > ub) && (Memento.warn(_LOGGER, "bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
            (!isnan(lb) && lb > vm_ub[bus]) && (lb = vm_lb[bus])
            (!isnan(ub) && ub < vm_lb[bus]) && (ub = vm_ub[bus])
            isnan(lb) && (lb = vm_lb[bus])
            isnan(ub) && (ub = vm_ub[bus])

            # vm bound-reduction computation
            vm_reduction = 0.0
            prev_ub = vm_ub[bus]
            prev_lb = vm_lb[bus]
            if (ub - lb >= min_bound_width)
                vm_reduction = (vm_ub[bus] - vm_lb[bus]) - (ub - lb)
                vm_lb[bus] = lb
                vm_ub[bus] = ub
            else
                mean = (ub + lb)/2.0
                if (mean - min_bound_width/2.0 < vm_lb[bus])
                    lb = vm_lb[bus]
                    ub = vm_lb[bus] + min_bound_width
                elseif (mean + min_bound_width/2.0 > vm_ub[bus])
                    ub = vm_ub[bus]
                    lb = vm_ub[bus] - min_bound_width
                else
                    lb = mean - (min_bound_width/2.0)
                    ub = mean + (min_bound_width/2.0)
                end
                vm_reduction = (vm_ub[bus] - vm_lb[bus]) - (ub - lb)
                vm_lb[bus] = lb
                vm_ub[bus] = ub
            end

            bool_vm_bound_reduction = bool_vm_bound_reduction || (prev_ub - vm_ub[bus] > improvement_tol) || (vm_lb[bus] - prev_lb > improvement_tol)
            bool_vm_range = bool_vm_range || (vm_ub[bus] - vm_lb[bus] > min_bound_width)
            # println("prev_ub - vm_ub[bus]: ", prev_ub - vm_ub[bus])
            # println("vm_lb[bus] - prev_lb: ", vm_lb[bus] - prev_lb)
            # println("vm_ub[bus] - vm_lb[bus]: ", vm_ub[bus] - vm_lb[bus])
            # println("bool_vm_bound_reduction: ", bool_vm_bound_reduction)
            # println("bool_vm_range: ", bool_vm_range)

            total_vm_reduction += (vm_reduction)
            max_vm_reduction = max(vm_reduction, max_vm_reduction)
        end
        avg_vm_reduction = total_vm_reduction/length(buses)

        vm_range_final = sum([vm_ub[bus] - vm_lb[bus] for bus in buses])

        # bound-tightening for the td variables
        for bp in buspairs
            (td_ub[bp] - td_lb[bp] < min_bound_width) && (continue)

            start_time = time()
            # td lower bound solve
            lb = NaN
            @objective(model, Min, model[:td][bp])
            JuMP.optimize!(model)
            if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                (nlb > td_lb[bp]) && (lb = nlb)
                # println("previous lb for td[$bp]: ", td_lb[bp])
                # println("BT lb for td[$bp]: ", lb)
            else
                Memento.warn(_LOGGER, "BT minimization problem for td[$bp] errored - change tolerances")
                continue
            end

            # td upper bound solve
            ub = NaN
            @objective(model, Max, model[:td][bp])
            JuMP.optimize!(model)
            if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                (nub < td_ub[bp]) && (ub = nub)
                # println("previous ub for td[$bp]: ", td_ub[bp])
                # println("BT ub for td[$bp]: ", ub)
            else
                Memento.warn(_LOGGER, "BT maximization problem for td[$bp] errored - change tolerances.")
                continue
            end
            end_time = time() - start_time
            max_td_iteration_time = max(end_time, max_td_iteration_time)

            # sanity checks
            (lb > ub) && (Memento.warn(_LOGGER, "bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
            (!isnan(lb) && lb > td_ub[bp]) && (lb = td_lb[bp])
            (!isnan(ub) && ub < td_lb[bp]) && (ub = td_ub[bp])
            isnan(lb) && (lb = td_lb[bp])
            isnan(ub) && (ub = td_ub[bp])

            # td bound-reduction computation
            td_reduction = 0.0
            prev_ub = td_ub[bp]
            prev_lb = td_lb[bp]
            if (ub - lb >= min_bound_width)
                td_reduction = (td_ub[bp] - td_lb[bp]) - (ub - lb)
                td_lb[bp] = lb
                td_ub[bp] = ub
            else
                mean = (lb + ub)/2.0
                if (mean - min_bound_width/2.0 < td_lb[bp])
                    lb = td_lb[bp]
                    ub = td_lb[bp] + min_bound_width
                elseif (mean + min_bound_width/2.0 > td_ub[bp])
                    ub = td_ub[bp]
                    lb = td_ub[bp] - min_bound_width
                else
                    lb = mean - (min_bound_width/2.0)
                    ub = mean + (min_bound_width/2.0)
                end
                td_reduction = (td_ub[bp] - td_lb[bp]) - (ub - lb)
                td_lb[bp] = lb
                td_ub[bp] = ub
            end

            bool_td_bound_reduction = bool_td_bound_reduction || (prev_ub - td_ub[bp] > improvement_tol) || (td_lb[bp] - prev_lb > improvement_tol)
            bool_td_range = bool_td_range || (td_ub[bp] - td_lb[bp] > min_bound_width)
            # println("final td[$bp] lb: ", td_lb[bp])
            # println("final td[$bp] ub: ", td_ub[bp])
            # println("prev_ub - td_ub[bp]: ", prev_ub - td_ub[bp])
            # println("td_lb[bp] - prev_lb: ", td_lb[bp] - prev_lb)
            # println("td_ub[bp] - td_lb[bp]", td_ub[bp] - td_lb[bp])
            # println("bool_td_bound_reduction: ", bool_td_bound_reduction)
            # println("bool_td_range: ", bool_td_range)

            total_td_reduction += (td_reduction)
            max_td_reduction = max(td_reduction, max_td_reduction)
        end
        avg_td_reduction = total_td_reduction/length(buspairs)

        td_range_final = sum([td_ub[bp] - td_lb[bp] for bp in buspairs])

        parallel_time_elapsed += max(max_vm_iteration_time, max_td_iteration_time)

        time_elapsed += (time() - iter_start_time)

        # populate the modifications, update the data
        modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
        PowerModels.update_data!(data, modifications)

        # run the qc relaxation for the updated bounds
        model = _build_model(data, model)
        JuMP.optimize!(model)

        if termination_status(model) in status_pass
            current_rel_gap = (upper_bound - objective_value(model))/upper_bound
            final_relaxation_objective = objective_value(model)
        else
            Memento.warn(_LOGGER, "relaxation solve failed in iteration $(current_iteration+1)")
            Memento.warn(_LOGGER, "using the previous iteration's gap to check relative gap stopping criteria")
        end

        # Rebuild the bound-tightening model
        if upper_bound_constraint
            @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] + gen["cost"][3] for (i,gen) in ref[:gen]) <= upper_bound)
        end

        Memento.info(_LOGGER, "iteration $(current_iteration+1), vm range: $vm_range_final, td range: $td_range_final, relaxation obj: $final_relaxation_objective")

        # termination criteria update
        (termination == :avg) && (check_termination = (avg_vm_reduction > improvement_tol || avg_td_reduction > improvement_tol))
        (termination == :max) && (check_termination = (max_vm_reduction > improvement_tol || max_td_reduction > improvement_tol))
        (termination == :gop) && (check_termination = ((bool_td_bound_reduction || bool_vm_bound_reduction) && (bool_td_range || bool_vm_range)))
        # interation counter update
        current_iteration += 1
        # check all the stopping criteria
        (current_iteration >= max_iter) && (Memento.info(_LOGGER, "maximum iteration limit reached"); break)
        (time_elapsed > time_limit) && (Memento.info(_LOGGER, "maximum time limit reached"); break)
        if (!isinf(rel_gap_tol)) && (current_rel_gap < rel_gap_tol)
            Memento.info(_LOGGER, "relative optimality gap < $rel_gap_tol")
            break
        end

    end

    branches_vad_same_sign_count = 0
    for (key, branch) in data["branch"]
        is_same_sign = (branch["angmax"] >=0 && branch["angmin"] >= 0) || (branch["angmax"] <=0 && branch["angmin"] <= 0)
        (is_same_sign) && (branches_vad_same_sign_count += 1)
    end

    stats["final_relaxation_objective"] = final_relaxation_objective
    stats["final_rel_gap_from_ub"] = isnan(upper_bound) ? Inf : current_rel_gap
    stats["vm_range_final"] = vm_range_final
    stats["avg_vm_range_final"] = vm_range_final/length(buses)

    stats["td_range_final"] = td_range_final
    stats["avg_td_range_final"] = td_range_final/length(buspairs)

    stats["run_time"] = time_elapsed
    stats["iteration_count"] = current_iteration
    stats["sim_parallel_run_time"] = parallel_time_elapsed

    stats["vad_sign_determined"] = branches_vad_same_sign_count

    modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
    PowerModels.update_data!(data, modifications)
    model = _build_model(data, model)

    return data, stats, model

end

function run_cyc_obbt_ots!(data::Dict{String,<:Any}, model;
    max_iter::Int = 100, # 100
    time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf,
    upper_bound_constraint::Bool = false,
    rel_gap_tol::Float64 = Inf,
    min_bound_width::Float64 = 1e-2,
    improvement_tol::Float64 = 1e-3,
    precision::Int = 4,
    termination::Symbol = :avg,
    kwargs...)

    # cycle_3_4 = []
    # if params["cycle_cuts"] && !params["separation"] && !params["callback"] 
    #     append!(cycle_3_4, cycle["cyc_3"])
    #     if params["cycle_max_bnd"] >= 4
    #         append!(cycle_3_4, cycle["cyc_4"])
    #     end
    #     for i in 1:length(cycle_3_4)
    #         cycle_3_4[i] = Tuple(cycle_3_4[i])
    #     end
    # end
    # buspairs = keys(ref[:buspairs])
    # z_lb = Dict{Any,Float64}([bp => 0 for bp in buspairs])
    # z_ub = Dict{Any,Float64}([bp => 1 for bp in buspairs])
    # zc_lb = Dict{Any,Float64}([cyc => 0 for cyc in cycle_3_4])
    # zc_ub = Dict{Any,Float64}([cyc => 1 for cyc in cycle_3_4])
    # model = _build_ots_model(data, model, false, z_lb, z_ub, zc_lb, zc_ub)

    # opf UB is also ots UB.
    nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
    result = solve_ac_opf(data, nlp_solver)
    upper_bound = result["objective"]

    Memento.info(_LOGGER, "maximum OBBT iterations set to default value of $max_iter")
    Memento.info(_LOGGER, "maximum time limit for OBBT set to default value of $time_limit seconds")

    # model_relaxation = instantiate_model(data, model_type, PowerModels.build_opf)
    # (_IM.ismultinetwork(model_relaxation)) && (Memento.error(_LOGGER, "OBBT is not supported for multi-networks"))
    # (ismulticonductor(model_relaxation)) && (Memento.error(_LOGGER, "OBBT is not supported for multi-conductor networks"))
    #
    # # check for model_type compatability with OBBT
    # _check_variables(model_relaxation)

    # check for other keyword argument consistencies
    _check_obbt_options(upper_bound, rel_gap_tol, upper_bound_constraint)

    # check termination norm criteria for obbt
    (termination != :avg && termination != :max && termination != :gop) && (Memento.error(_LOGGER, "OBBT termination criteria can only be :max or :avg or :gop"))

    # pass status
    status_pass = [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]

    # Relax binary constraints for z and zc.
    cycle_3_4 = []
    if params["cycle_cuts"] && !params["separation"] && !params["callback"] #&& !params["off_insights"]
        append!(cycle_3_4, cycle["cyc_3"])
        if params["cycle_max_bnd"] >= 4
            append!(cycle_3_4, cycle["cyc_4"])
        end
        for i in 1:length(cycle_3_4)
            cycle_3_4[i] = Tuple(cycle_3_4[i])
        end
        for cyc in cycle_3_4
            unset_binary(model[:zc][cyc])
            set_lower_bound(model[:zc][cyc], 0)
            set_upper_bound(model[:zc][cyc], 1)
        end
    end

    for bp in keys(ref[:buspairs])
        unset_binary(model[:z][bp])
        set_lower_bound(model[:z][bp], 0)
        set_upper_bound(model[:z][bp], 1)
    end
    # Write_model(model, "results/model2.txt") #@

    # compute initial relative gap between relaxation objective and upper_bound
    JuMP.optimize!(model)
    current_relaxation_objective = objective_value(model)
    if upper_bound < current_relaxation_objective
        Memento.error(_LOGGER, "the upper bound provided to OBBT is not a valid ACOTS upper bound")
    end
    # if !(result_relaxation["termination_status"] in status_pass)
    #     Memento.warn(_LOGGER, "initial relaxation solve status is $(result_relaxation["termination_status"])")
    #     if result_relaxation["termination_status"] == :SubOptimal
    #         Memento.warn(_LOGGER, "continuing with the bound-tightening algorithm")
    #     end
    # end
    current_rel_gap = Inf
    if !isinf(upper_bound)
        current_rel_gap = (upper_bound - current_relaxation_objective)/upper_bound
        Memento.info(_LOGGER, "Initial relaxation gap = $current_rel_gap")
    end

    buses = keys(ref[:bus])
    buspairs = keys(ref[:buspairs])

    z_lb = Dict{Any,Float64}([bp => 0 for bp in buspairs])
    z_ub = Dict{Any,Float64}([bp => 1 for bp in buspairs])
    zc_lb = Dict{Any,Float64}([cyc => 0 for cyc in cycle_3_4])
    zc_ub = Dict{Any,Float64}([cyc => 1 for cyc in cycle_3_4])

    model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)
    # print(model)
    if upper_bound_constraint
        @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
    end
    # (upper_bound_constraint) && (_constraint_obj_bound(model_bt, upper_bound))

    stats = Dict{String,Any}()
    # stats["model_type"] = model_type
    stats["initial_relaxation_objective"] = current_relaxation_objective
    stats["initial_rel_gap_from_ub"] = current_rel_gap
    stats["upper_bound"] = upper_bound

    vm_lb = Dict{Any,Float64}([bus => JuMP.lower_bound(model[:vm][bus]) for bus in buses])
    vm_ub = Dict{Any,Float64}([bus => JuMP.upper_bound(model[:vm][bus]) for bus in buses])
    td_lb = Dict{Any,Float64}( [bp => ref[:buspairs][bp]["angmin"] for bp in buspairs] )
    td_ub = Dict{Any,Float64}( [bp => ref[:buspairs][bp]["angmax"] for bp in buspairs] )

    vm_range_init = sum([vm_ub[bus] - vm_lb[bus] for bus in buses])
    stats["vm_range_init"] = vm_range_init
    stats["avg_vm_range_init"] = vm_range_init/length(buses)

    td_range_init = sum([td_ub[bp] - td_lb[bp] for bp in buspairs])
    stats["td_range_init"] = td_range_init
    stats["avg_td_range_init"] = td_range_init/length(buspairs)

    z_range_init = sum(z_ub[bp] - z_lb[bp] for bp in buspairs)
    stats["z_range_init"] = z_range_init
    stats["avg_z_range_init"] = z_range_init / length(buspairs)
    if params["cycle_cuts"] && !params["separation"] && !params["callback"] # && !params["off_insights"]
        zc_range_init = sum(zc_ub[cyc] - zc_lb[cyc] for cyc in cycle_3_4)
        stats["zc_range_init"] = zc_range_init
        stats["avg_zc_range_init"] = zc_range_init / length(cycle_3_4)
    end

    vm_range_final = 0.0
    td_range_final = 0.0
    z_range_final = 0.0
    zc_range_final = 0.0

    total_vm_reduction = Inf
    max_vm_reduction = Inf
    avg_vm_reduction = Inf

    total_td_reduction = Inf
    max_td_reduction = Inf
    avg_td_reduction = Inf

    total_z_reduction = Inf
    max_z_reduction = Inf
    avg_z_reduction = Inf

    total_zc_reduction = Inf
    max_zc_reduction = Inf
    avg_zc_reduction = Inf

    final_relaxation_objective = NaN

    current_iteration = 0
    time_elapsed = 0.0
    parallel_time_elapsed = 0.0

    bool_td_bound_reduction = true
    bool_vm_bound_reduction = true
    bool_td_range = true
    bool_vm_range = true
    bool_z_bound_reduction = true
    bool_z_range = true
    bool_zc_bound_reduction = true
    bool_zc_range = true

    check_termination = true
    (termination == :avg) && (check_termination = (avg_vm_reduction > improvement_tol || avg_td_reduction > improvement_tol || avg_z_reduction > improvement_tol || avg_zc_reduction > improvement_tol))
    (termination == :max) && (check_termination = (max_vm_reduction > improvement_tol || max_td_reduction > improvement_tol || max_z_reduction > improvement_tol || max_zc_reduction > improvement_tol))
    (termination == :gop) && (check_termination = ((bool_td_bound_reduction || bool_vm_bound_reduction || bool_z_bound_reduction || bool_zc_bound_reduction) && (bool_td_range || bool_vm_range || bool_z_range || bool_zc_range))) # The termination conditions (11) and (12) from Gaopinath et. al. 2020.

    # check_termination = false
    while check_termination
        # if cnt0 >= 1
        #     break
        # end
        # cnt0 += 1

        iter_start_time = time()
        total_vm_reduction = 0.0
        avg_vm_reduction = 0.0
        max_vm_reduction = 0.0
        max_vm_iteration_time = 0.0
        bool_vm_bound_reduction = false
        bool_vm_range = false

        total_td_reduction = 0.0
        avg_td_reduction = 0.0
        max_td_reduction = 0.0
        max_td_iteration_time = 0.0
        bool_td_bound_reduction = false
        bool_td_range = false

        total_z_reduction = 0.0
        avg_z_reduction = 0.0
        max_z_reduction = 0.0
        max_z_iteration_time = 0.0
        bool_z_bound_reduction = false
        bool_z_range = false

        total_zc_reduction = 0.0
        avg_zc_reduction = 0.0
        max_zc_reduction = 0.0
        max_zc_iteration_time = 0.0
        bool_zc_bound_reduction = false
        bool_zc_range = false

        # bound-tightening for the vm variables
        if params["obbt_vm"]
            for bus in buses
                (vm_ub[bus] - vm_lb[bus] < min_bound_width) && (continue)

                start_time = time()
                # vm lower bound solve
                lb = NaN
                @objective(model, Min, model[:vm][bus])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    (nlb > vm_lb[bus]) && (lb = nlb)
                    # lb = objective_value(model)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for vm[$bus] errored - change tolerances.")
                    continue
                end

                #vm upper bound solve
                ub = NaN
                @objective(model, Max, model[:vm][bus])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    (nub < vm_ub[bus]) && (ub = nub)
                    # ub = objective_value(model)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for vm[$bus] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_vm_iteration_time = max(end_time, max_vm_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > vm_ub[bus]) && (lb = vm_lb[bus])
                (!isnan(ub) && ub < vm_lb[bus]) && (ub = vm_ub[bus])
                isnan(lb) && (lb = vm_lb[bus])
                isnan(ub) && (ub = vm_ub[bus])

                # vm bound-reduction computation
                vm_reduction = 0.0
                prev_ub = vm_ub[bus]
                prev_lb = vm_lb[bus]
                if (ub - lb >= min_bound_width)
                    vm_reduction = (vm_ub[bus] - vm_lb[bus]) - (ub - lb)
                    vm_lb[bus] = lb
                    vm_ub[bus] = ub
                else
                    mean = (ub + lb)/2.0
                    if (mean - min_bound_width/2.0 < vm_lb[bus])
                        lb = vm_lb[bus]
                        ub = vm_lb[bus] + min_bound_width
                    elseif (mean + min_bound_width/2.0 > vm_ub[bus])
                        ub = vm_ub[bus]
                        lb = vm_ub[bus] - min_bound_width
                    else
                        lb = mean - (min_bound_width/2.0)
                        ub = mean + (min_bound_width/2.0)
                    end
                    vm_reduction = (vm_ub[bus] - vm_lb[bus]) - (ub - lb)
                    vm_lb[bus] = lb
                    vm_ub[bus] = ub
                end

                bool_vm_bound_reduction = bool_vm_bound_reduction || (prev_ub - vm_ub[bus] > improvement_tol) || (vm_lb[bus] - prev_lb > improvement_tol)
                bool_vm_range = bool_vm_range || (vm_ub[bus] - vm_lb[bus] > min_bound_width)
                # println("prev_ub - vm_ub[bus]: ", prev_ub - vm_ub[bus])
                # println("vm_lb[bus] - prev_lb: ", vm_lb[bus] - prev_lb)
                # println("vm_ub[bus] - vm_lb[bus]: ", vm_ub[bus] - vm_lb[bus])
                # println("bool_vm_bound_reduction: ", bool_vm_bound_reduction)
                # println("bool_vm_range: ", bool_vm_range)

                total_vm_reduction += (vm_reduction)
                max_vm_reduction = max(vm_reduction, max_vm_reduction)
            end
            avg_vm_reduction = total_vm_reduction/length(buses)

            vm_range_final = sum([vm_ub[bus] - vm_lb[bus] for bus in buses])
        end

        # bound-tightening for the td variables
        if params["obbt_td"]
            for bp in buspairs
                (td_ub[bp] - td_lb[bp] < min_bound_width) && (continue)

                # leaf_gens(ref)
                # @objective(model, Min, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]))
                # @objective(model, Min, 0)
                # Write_model(model, "results/model_debug.txt")#@
                # JuMP.optimize!(model)
                # println("termination status for debug: ", termination_status(model))

                start_time = time()
                # td lower bound solve
                lb = NaN
                # ub = NaN
                if params["ext_cons"]
                    @objective(model, Min, model[:td][bp])
                else
                    @objective(model, Min, model[:va][bp[1]] - model[:va][bp[2]])
                end
                # Write_model(model, "results/model_lb.txt")#@
                JuMP.optimize!(model)
                # println("termination status for min td[$bp]: ", termination_status(model))
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    (nlb > td_lb[bp]) && (lb = nlb)
                    # lb = objective_value(model)
                    # println("previous lb for td[$bp]: ", td_lb[bp])
                    # println("BT lb for td[$bp]: ", lb)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for td[$bp] errored - change tolerances")
                    continue
                end

                # td upper bound solve
                ub = NaN
                if params["ext_cons"]
                    @objective(model, Max, model[:td][bp])
                else
                    @objective(model, Max, model[:va][bp[1]] - model[:va][bp[2]])
                end
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    (nub < td_ub[bp]) && (ub = nub)
                    # ub = objective_value(model)
                    # println("previous ub for td[$bp]: ", td_ub[bp])
                    # println("BT ub for td[$bp]: ", ub)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for td[$bp] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_td_iteration_time = max(end_time, max_td_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > td_ub[bp]) && (lb = td_lb[bp])
                (!isnan(ub) && ub < td_lb[bp]) && (ub = td_ub[bp])
                isnan(lb) && (lb = td_lb[bp])
                isnan(ub) && (ub = td_ub[bp])

                # td bound-reduction computation
                td_reduction = 0.0
                prev_ub = td_ub[bp]
                prev_lb = td_lb[bp]
                if (ub - lb >= min_bound_width)
                    td_reduction = (td_ub[bp] - td_lb[bp]) - (ub - lb)
                    td_lb[bp] = lb
                    td_ub[bp] = ub
                    # if cnt >= 2
                    #     break
                    # end
                    # cnt += 1
                else
                    mean = (lb + ub)/2.0
                    if (mean - min_bound_width/2.0 < td_lb[bp])
                        lb = td_lb[bp]
                        ub = td_lb[bp] + min_bound_width
                    elseif (mean + min_bound_width/2.0 > td_ub[bp])
                        ub = td_ub[bp]
                        lb = td_ub[bp] - min_bound_width
                    else
                        lb = mean - (min_bound_width/2.0)
                        ub = mean + (min_bound_width/2.0)
                    end
                    td_reduction = (td_ub[bp] - td_lb[bp]) - (ub - lb)
                    td_lb[bp] = lb
                    td_ub[bp] = ub
                end

                bool_td_bound_reduction = bool_td_bound_reduction || (prev_ub - td_ub[bp] > improvement_tol) || (td_lb[bp] - prev_lb > improvement_tol)
                bool_td_range = bool_td_range || (td_ub[bp] - td_lb[bp] > min_bound_width)
                # println("final td[$bp] lb: ", td_lb[bp])
                # println("final td[$bp] ub: ", td_ub[bp])
                # println("prev_ub - td_ub[bp]: ", prev_ub - td_ub[bp])
                # println("td_lb[bp] - prev_lb: ", td_lb[bp] - prev_lb)
                # println("td_ub[bp] - td_lb[bp]", td_ub[bp] - td_lb[bp])
                # println("bool_td_bound_reduction: ", bool_td_bound_reduction)
                # println("bool_td_range: ", bool_td_range)

                total_td_reduction += (td_reduction)
                max_td_reduction = max(td_reduction, max_td_reduction)
                # break
            end
            avg_td_reduction = total_td_reduction/length(buspairs)
            td_range_final = sum([td_ub[bp] - td_lb[bp] for bp in buspairs])
        end

        # bound-tightening for the z variables
        if params["obbt_z"]
            for bp in buspairs
                (z_ub[bp] - z_lb[bp] < min_bound_width) && (continue)

                start_time = time()
                # z lower bound solve
                lb = NaN
                @objective(model, Min, model[:z][bp])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    # nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    nlb = objective_value(model)
                    (nlb > 10.0^(-precision + 1)) && (lb = 1)
                    # println("previous lb for z[$bp]: ", z_lb[bp])
                    # println("BT lb for z[$bp]: ", lb)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for z[$bp] errored - change tolerances")
                    continue
                end

                # z upper bound solve
                ub = NaN
                @objective(model, Max, model[:z][bp])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    # nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    nub = objective_value(model)
                    (nub < 1.0 - 10.0^(-precision + 1)) && (ub = 0)
                    # println("previous ub for z[$bp]: ", z_ub[bp])
                    # println("BT ub for z[$bp]: ", ub)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for z[$bp] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_z_iteration_time = max(end_time, max_z_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "z bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > z_ub[bp]) && (lb = z_lb[bp])
                (!isnan(ub) && ub < z_lb[bp]) && (ub = z_ub[bp])
                isnan(lb) && (lb = z_lb[bp])
                isnan(ub) && (ub = z_ub[bp])

                # z bound-reduction computation
                z_lb[bp] = lb
                z_ub[bp] = ub
            end
        end

        # bound-tightening for the z variables
        if params["obbt_zc"]
            for cyc in cycle_3_4
                (zc_ub[cyc] - zc_lb[cyc] < min_bound_width) && (continue)

                start_time = time()
                # z lower bound solve
                lb = NaN
                @objective(model, Min, model[:zc][cyc])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    # nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    nlb = objective_value(model)
                    (nlb > 10.0^(-precision + 1)) && (lb = 1)
                    # println("previous lb for z[$cyc]: ", z_lb[cyc])
                    # println("BT lb for z[$cyc]: ", lb)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for zc[$cyc] errored - change tolerances")
                    continue
                end

                # z upper bound solve
                ub = NaN
                @objective(model, Max, model[:zc][cyc])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    # nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    nub = objective_value(model)
                    (nub < 1.0 - 10.0^(-precision + 1)) && (ub = 0)
                    # println("previous ub for z[$bp]: ", z_ub[bp])
                    # println("BT ub for z[$bp]: ", ub)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for zc[$cyc] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_zc_iteration_time = max(end_time, max_zc_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "zc bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > zc_ub[cyc]) && (lb = zc_lb[cyc])
                (!isnan(ub) && ub < zc_lb[cyc]) && (ub = zc_ub[cyc])
                isnan(lb) && (lb = zc_lb[cyc])
                isnan(ub) && (ub = zc_ub[cyc])
            end
        end

        parallel_time_elapsed += max(max_vm_iteration_time, max_td_iteration_time, max_z_iteration_time, max_zc_iteration_time)

        time_elapsed += (time() - iter_start_time)

        # populate the modifications, update the data
        modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
        PowerModels.update_data!(data, modifications)

        # run the qc relaxation for the updated bounds, with binary variables relaxed.
        model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)
        JuMP.optimize!(model)
        # Write_model(model, "results/model_iterEnd.txt")

        # println(termination_status(model)) #@

        if termination_status(model) in status_pass
            current_rel_gap = (upper_bound - objective_value(model)) / upper_bound
            final_relaxation_objective = objective_value(model)
        else
            Memento.warn(_LOGGER, "relaxation solve failed in iteration $(current_iteration+1)")
            Memento.warn(_LOGGER, "using the previous iteration's gap to check relative gap stopping criteria")
        end

        # Rebuild the bound-tightening model
        if upper_bound_constraint
            @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
        end

        Memento.info(_LOGGER, "iteration $(current_iteration+1), vm range: $vm_range_final, td range: $td_range_final, relaxation obj: $final_relaxation_objective")

        # termination criteria update
        (termination == :avg) && (check_termination = (avg_vm_reduction > improvement_tol || avg_td_reduction > improvement_tol))
        (termination == :max) && (check_termination = (max_vm_reduction > improvement_tol || max_td_reduction > improvement_tol))
        (termination == :gop) && (check_termination = ((bool_td_bound_reduction || bool_vm_bound_reduction) && (bool_td_range || bool_vm_range)))
        # interation counter update
        current_iteration += 1
        # check all the stopping criteria
        (current_iteration >= max_iter) && (Memento.info(_LOGGER, "maximum iteration limit reached"); break)
        (time_elapsed > time_limit) && (Memento.info(_LOGGER, "maximum time limit reached"); break)

        # This one doesn't seem necessary as the lower bound is from ACOTS LP relaxations.
        if (!isinf(rel_gap_tol)) && (current_rel_gap < rel_gap_tol)
            Memento.info(_LOGGER, "relative optimality gap < $rel_gap_tol")
            break
        end

    end

    branches_vad_same_sign_count = 0
    for (key, branch) in data["branch"]
        is_same_sign = (branch["angmax"] >=0 && branch["angmin"] >= 0) || (branch["angmax"] <=0 && branch["angmin"] <= 0)
        (is_same_sign) && (branches_vad_same_sign_count += 1)
    end

    stats["final_relaxation_objective"] = final_relaxation_objective
    stats["final_rel_gap_from_ub"] = isnan(upper_bound) ? Inf : current_rel_gap

    stats["vm_range_final"] = vm_range_final
    stats["avg_vm_range_final"] = vm_range_final/length(buses)
    stats["percent_vm_range_reduction"] = 1 - vm_range_final / vm_range_init

    stats["td_range_final"] = td_range_final
    stats["avg_td_range_final"] = td_range_final/length(buspairs)
    stats["percent_td_range_reduction"] = 1 - td_range_final / td_range_init

    stats["run_time"] = time_elapsed
    stats["iteration_count"] = current_iteration
    stats["sim_parallel_run_time"] = parallel_time_elapsed

    stats["vad_sign_determined"] = branches_vad_same_sign_count

    if params["obbt_z"]
        num_freeLines = 0
        num_zero_lines = 0
        for bp in keys(ref[:buspairs])
            if z_ub[bp] - z_lb[bp] > 0.5
                num_freeLines += 1
            end
            if z_ub[bp] == 0
                num_zero_lines += 1
            end
        end
        stats["ratio_free_lines"] = num_freeLines/length(keys(ref[:buspairs]))
        stats["num_zero_lines"] = num_zero_lines
    end

    modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
    PowerModels.update_data!(data, modifications)
    model = _build_ots_model(data, model, false, z_lb, z_ub, zc_lb, zc_ub)

    return data, stats, model

end

function run_cyc_obbt_ots_perRound!(data::Dict{String,<:Any}, model;
    max_iter::Int = 100, # 100
    time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf,
    upper_bound_constraint::Bool = false,
    rel_gap_tol::Float64 = Inf,
    min_bound_width::Float64 = 1e-2,
    improvement_tol::Float64 = 1e-3,
    precision::Int = 4,
    termination::Symbol = :avg,
    kwargs...)

    time1_init = time()

    # opf UB is also ots UB.
    nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
    result = solve_ac_opf(data, nlp_solver)
    upper_bound = result["objective"]

    Memento.info(_LOGGER, "maximum OBBT iterations set to default value of $max_iter")
    Memento.info(_LOGGER, "maximum time limit for OBBT set to default value of $time_limit seconds")

    println("time1: $(time() - time1_init)")

    # model_relaxation = instantiate_model(data, model_type, PowerModels.build_opf)
    # (_IM.ismultinetwork(model_relaxation)) && (Memento.error(_LOGGER, "OBBT is not supported for multi-networks"))
    # (ismulticonductor(model_relaxation)) && (Memento.error(_LOGGER, "OBBT is not supported for multi-conductor networks"))
    #
    # # check for model_type compatability with OBBT
    # _check_variables(model_relaxation)

    time2_init = time()
    # check for other keyword argument consistencies
    _check_obbt_options(upper_bound, rel_gap_tol, upper_bound_constraint)

    # check termination norm criteria for obbt
    (termination != :avg && termination != :max && termination != :gop) && (Memento.error(_LOGGER, "OBBT termination criteria can only be :max or :avg or :gop"))

    # pass status
    status_pass = [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]

    # Relax binary constraints for z and zc.
    cycle_3_4 = []
    if params["cycle_cuts"] && !params["separation"] && !params["callback"] # && !params["off_insights"]
        append!(cycle_3_4, cycle["cyc_3"])
        if params["cycle_max_bnd"] >= 4
            append!(cycle_3_4, cycle["cyc_4"])
        end
        for i in 1:length(cycle_3_4)
            cycle_3_4[i] = Tuple(cycle_3_4[i])
        end
        for cyc in cycle_3_4
            unset_binary(model[:zc][cyc])
            set_lower_bound(model[:zc][cyc], 0)
            set_upper_bound(model[:zc][cyc], 1)
        end
    end
    for bp in keys(ref[:buspairs])
        unset_binary(model[:z][bp])
        set_lower_bound(model[:z][bp], 0)
        set_upper_bound(model[:z][bp], 1)
    end
    # Write_model(model, "results/model2.txt") #@

    # compute initial relative gap between relaxation objective and upper_bound
    JuMP.optimize!(model)
    current_relaxation_objective = objective_value(model)
    if upper_bound < current_relaxation_objective
        Memento.error(_LOGGER, "the upper bound provided to OBBT is not a valid ACOTS upper bound")
    end
    # if !(result_relaxation["termination_status"] in status_pass)
    #     Memento.warn(_LOGGER, "initial relaxation solve status is $(result_relaxation["termination_status"])")
    #     if result_relaxation["termination_status"] == :SubOptimal
    #         Memento.warn(_LOGGER, "continuing with the bound-tightening algorithm")
    #     end
    # end
    current_rel_gap = Inf
    if !isinf(upper_bound)
        current_rel_gap = (upper_bound - current_relaxation_objective)/upper_bound
        Memento.info(_LOGGER, "Initial relaxation gap = $current_rel_gap")
    end

    buses = keys(ref[:bus])
    buspairs = keys(ref[:buspairs])

    z_lb = Dict{Any,Float64}([bp => 0 for bp in buspairs])
    z_ub = Dict{Any,Float64}([bp => 1 for bp in buspairs])
    zc_lb = Dict{Any,Float64}([cyc => 0 for cyc in cycle_3_4])
    zc_ub = Dict{Any,Float64}([cyc => 1 for cyc in cycle_3_4])

    model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)
    # print(model)
    if upper_bound_constraint
        @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
    end
    # (upper_bound_constraint) && (_constraint_obj_bound(model_bt, upper_bound))

    stats = Dict{String,Any}()
    # stats["model_type"] = model_type
    stats["initial_relaxation_objective"] = current_relaxation_objective
    stats["initial_rel_gap_from_ub"] = current_rel_gap
    stats["upper_bound"] = upper_bound

    vm_lb = Dict{Any,Float64}([bus => JuMP.lower_bound(model[:vm][bus]) for bus in buses])
    vm_ub = Dict{Any,Float64}([bus => JuMP.upper_bound(model[:vm][bus]) for bus in buses])
    td_lb = Dict{Any,Float64}( [bp => ref[:buspairs][bp]["angmin"] for bp in buspairs] )
    td_ub = Dict{Any,Float64}( [bp => ref[:buspairs][bp]["angmax"] for bp in buspairs] )

    vm_range_init = sum([vm_ub[bus] - vm_lb[bus] for bus in buses])
    stats["vm_range_init"] = vm_range_init
    stats["avg_vm_range_init"] = vm_range_init/length(buses)

    td_range_init = sum([td_ub[bp] - td_lb[bp] for bp in buspairs])
    stats["td_range_init"] = td_range_init
    stats["avg_td_range_init"] = td_range_init/length(buspairs)

    z_range_init = sum(z_ub[bp] - z_lb[bp] for bp in buspairs)
    stats["z_range_init"] = z_range_init
    stats["avg_z_range_init"] = z_range_init / length(buspairs)
    if params["cycle_cuts"] && !params["separation"] && !params["callback"] #  && !params["off_insights"]
        zc_range_init = sum(zc_ub[cyc] - zc_lb[cyc] for cyc in cycle_3_4)
        stats["zc_range_init"] = zc_range_init
        stats["avg_zc_range_init"] = zc_range_init / length(cycle_3_4)
    end

    vm_range_final = 0.0
    td_range_final = 0.0
    z_range_final = 0.0
    zc_range_final = 0.0

    total_vm_reduction = Inf
    max_vm_reduction = Inf
    avg_vm_reduction = Inf

    total_td_reduction = Inf
    max_td_reduction = Inf
    avg_td_reduction = Inf

    total_z_reduction = Inf
    max_z_reduction = Inf
    avg_z_reduction = Inf

    total_zc_reduction = Inf
    max_zc_reduction = Inf
    avg_zc_reduction = Inf

    final_relaxation_objective = NaN

    current_iteration = 0
    time_elapsed = 0.0
    parallel_time_elapsed = 0.0

    bool_td_bound_reduction = true
    bool_vm_bound_reduction = true
    bool_td_range = true
    bool_vm_range = true
    bool_z_bound_reduction = true
    bool_z_range = true
    bool_zc_bound_reduction = true
    bool_zc_range = true

    check_termination = true
    (termination == :avg) && (check_termination = (avg_vm_reduction > improvement_tol || avg_td_reduction > improvement_tol || avg_z_reduction > improvement_tol || avg_zc_reduction > improvement_tol))
    (termination == :max) && (check_termination = (max_vm_reduction > improvement_tol || max_td_reduction > improvement_tol || max_z_reduction > improvement_tol || max_zc_reduction > improvement_tol))
    (termination == :gop) && (check_termination = ((bool_td_bound_reduction || bool_vm_bound_reduction || bool_z_bound_reduction || bool_zc_bound_reduction) && (bool_td_range || bool_vm_range || bool_z_range || bool_zc_range))) # The termination conditions (11) and (12) from Gaopinath et. al. 2020.

    println("time2: $(time() - time2_init)")

    time3_init = time()
    # check_termination = false
    while check_termination
        # if cnt0 >= 1
        #     break
        # end
        # cnt0 += 1

        iter_start_time = time()
        total_vm_reduction = 0.0
        avg_vm_reduction = 0.0
        max_vm_reduction = 0.0
        max_vm_iteration_time = 0.0
        bool_vm_bound_reduction = false
        bool_vm_range = false

        total_td_reduction = 0.0
        avg_td_reduction = 0.0
        max_td_reduction = 0.0
        max_td_iteration_time = 0.0
        bool_td_bound_reduction = false
        bool_td_range = false

        total_z_reduction = 0.0
        avg_z_reduction = 0.0
        max_z_reduction = 0.0
        max_z_iteration_time = 0.0
        bool_z_bound_reduction = false
        bool_z_range = false

        total_zc_reduction = 0.0
        avg_zc_reduction = 0.0
        max_zc_reduction = 0.0
        max_zc_iteration_time = 0.0
        bool_zc_bound_reduction = false
        bool_zc_range = false

        # bound-tightening for the vm variables
        if params["obbt_vm"]
            for bus in buses
                (vm_ub[bus] - vm_lb[bus] < min_bound_width) && (continue)

                start_time = time()
                # vm lower bound solve
                lb = NaN
                @objective(model, Min, model[:vm][bus])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    (nlb > vm_lb[bus]) && (lb = nlb)
                    # lb = objective_value(model)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for vm[$bus] errored - change tolerances.")
                    continue
                end

                #vm upper bound solve
                ub = NaN
                @objective(model, Max, model[:vm][bus])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    (nub < vm_ub[bus]) && (ub = nub)
                    # ub = objective_value(model)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for vm[$bus] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_vm_iteration_time = max(end_time, max_vm_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > vm_ub[bus]) && (lb = vm_lb[bus])
                (!isnan(ub) && ub < vm_lb[bus]) && (ub = vm_ub[bus])
                isnan(lb) && (lb = vm_lb[bus])
                isnan(ub) && (ub = vm_ub[bus])

                # vm bound-reduction computation
                vm_reduction = 0.0
                prev_ub = vm_ub[bus]
                prev_lb = vm_lb[bus]
                if (ub - lb >= min_bound_width)
                    vm_reduction = (vm_ub[bus] - vm_lb[bus]) - (ub - lb)
                    vm_lb[bus] = lb
                    vm_ub[bus] = ub
                else
                    mean = (ub + lb)/2.0
                    if (mean - min_bound_width/2.0 < vm_lb[bus])
                        lb = vm_lb[bus]
                        ub = vm_lb[bus] + min_bound_width
                    elseif (mean + min_bound_width/2.0 > vm_ub[bus])
                        ub = vm_ub[bus]
                        lb = vm_ub[bus] - min_bound_width
                    else
                        lb = mean - (min_bound_width/2.0)
                        ub = mean + (min_bound_width/2.0)
                    end
                    vm_reduction = (vm_ub[bus] - vm_lb[bus]) - (ub - lb)
                    vm_lb[bus] = lb
                    vm_ub[bus] = ub
                end

                bool_vm_bound_reduction = bool_vm_bound_reduction || (prev_ub - vm_ub[bus] > improvement_tol) || (vm_lb[bus] - prev_lb > improvement_tol)
                bool_vm_range = bool_vm_range || (vm_ub[bus] - vm_lb[bus] > min_bound_width)
                # println("prev_ub - vm_ub[bus]: ", prev_ub - vm_ub[bus])
                # println("vm_lb[bus] - prev_lb: ", vm_lb[bus] - prev_lb)
                # println("vm_ub[bus] - vm_lb[bus]: ", vm_ub[bus] - vm_lb[bus])
                # println("bool_vm_bound_reduction: ", bool_vm_bound_reduction)
                # println("bool_vm_range: ", bool_vm_range)

                total_vm_reduction += (vm_reduction)
                max_vm_reduction = max(vm_reduction, max_vm_reduction)

                # Modify the problem every time a pair of variable bounds is updated.
                modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
                PowerModels.update_data!(data, modifications)

                # run the qc relaxation for the updated bounds
                model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)

                # Rebuild the bound-tightening model
                if upper_bound_constraint
                    @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
                end

            end
            avg_vm_reduction = total_vm_reduction/length(buses)

            vm_range_final = sum([vm_ub[bus] - vm_lb[bus] for bus in buses])
        end

        # bound-tightening for the td variables
        if params["obbt_td"]
            for bp in buspairs
                (td_ub[bp] - td_lb[bp] < min_bound_width) && (continue)

                # leaf_gens(ref)
                # @objective(model, Min, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]))
                # @objective(model, Min, 0)
                # Write_model(model, "results/model_debug.txt")#@
                # JuMP.optimize!(model)
                # println("termination status for debug: ", termination_status(model))

                start_time = time()
                # td lower bound solve
                lb = NaN
                # ub = NaN
                if params["ext_cons"]
                    @objective(model, Min, model[:td][bp])
                else
                    @objective(model, Min, model[:va][bp[1]] - model[:va][bp[2]])
                end
                # Write_model(model, "results/model_lb.txt")#@
                JuMP.optimize!(model)
                # println("termination status for min td[$bp]: ", termination_status(model))
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    (nlb > td_lb[bp]) && (lb = nlb)
                    # lb = objective_value(model)
                    # println("previous lb for td[$bp]: ", td_lb[bp])
                    # println("BT lb for td[$bp]: ", lb)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for td[$bp] errored - change tolerances")
                    continue
                end

                # td upper bound solve
                ub = NaN
                if params["ext_cons"]
                    @objective(model, Max, model[:td][bp])
                else
                    @objective(model, Max, model[:va][bp[1]] - model[:va][bp[2]])
                end
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    (nub < td_ub[bp]) && (ub = nub)
                    # ub = objective_value(model)
                    # println("previous ub for td[$bp]: ", td_ub[bp])
                    # println("BT ub for td[$bp]: ", ub)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for td[$bp] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_td_iteration_time = max(end_time, max_td_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > td_ub[bp]) && (lb = td_lb[bp])
                (!isnan(ub) && ub < td_lb[bp]) && (ub = td_ub[bp])
                isnan(lb) && (lb = td_lb[bp])
                isnan(ub) && (ub = td_ub[bp])

                # td bound-reduction computation
                td_reduction = 0.0
                prev_ub = td_ub[bp]
                prev_lb = td_lb[bp]
                if (ub - lb >= min_bound_width)
                    td_reduction = (td_ub[bp] - td_lb[bp]) - (ub - lb)
                    td_lb[bp] = lb
                    td_ub[bp] = ub
                    # if cnt >= 2
                    #     break
                    # end
                    # cnt += 1
                else
                    mean = (lb + ub)/2.0
                    if (mean - min_bound_width/2.0 < td_lb[bp])
                        lb = td_lb[bp]
                        ub = td_lb[bp] + min_bound_width
                    elseif (mean + min_bound_width/2.0 > td_ub[bp])
                        ub = td_ub[bp]
                        lb = td_ub[bp] - min_bound_width
                    else
                        lb = mean - (min_bound_width/2.0)
                        ub = mean + (min_bound_width/2.0)
                    end
                    td_reduction = (td_ub[bp] - td_lb[bp]) - (ub - lb)
                    td_lb[bp] = lb
                    td_ub[bp] = ub
                end

                bool_td_bound_reduction = bool_td_bound_reduction || (prev_ub - td_ub[bp] > improvement_tol) || (td_lb[bp] - prev_lb > improvement_tol)
                bool_td_range = bool_td_range || (td_ub[bp] - td_lb[bp] > min_bound_width)
                # println("final td[$bp] lb: ", td_lb[bp])
                # println("final td[$bp] ub: ", td_ub[bp])
                # println("prev_ub - td_ub[bp]: ", prev_ub - td_ub[bp])
                # println("td_lb[bp] - prev_lb: ", td_lb[bp] - prev_lb)
                # println("td_ub[bp] - td_lb[bp]", td_ub[bp] - td_lb[bp])
                # println("bool_td_bound_reduction: ", bool_td_bound_reduction)
                # println("bool_td_range: ", bool_td_range)

                total_td_reduction += (td_reduction)
                max_td_reduction = max(td_reduction, max_td_reduction)

                modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
                PowerModels.update_data!(data, modifications)

                # run the qc relaxation for the updated bounds
                model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)

                # Rebuild the bound-tightening model
                if upper_bound_constraint
                    @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
                end
            end
            avg_td_reduction = total_td_reduction/length(buspairs)
            td_range_final = sum([td_ub[bp] - td_lb[bp] for bp in buspairs])
        end

        # bound-tightening for the z variables
        if params["obbt_z"]
            for bp in buspairs
                (z_ub[bp] - z_lb[bp] < min_bound_width) && (continue)

                start_time = time()
                # z lower bound solve
                lb = NaN
                @objective(model, Min, model[:z][bp])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    (nlb > 10.0^(-precision + 1)) && (lb = 1)
                    # println("previous lb for z[$bp]: ", z_lb[bp])
                    # println("BT lb for z[$bp]: ", lb)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for z[$bp] errored - change tolerances")
                    continue
                end

                # z upper bound solve
                ub = NaN
                @objective(model, Max, model[:z][bp])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    (nub < 1.0 - 10.0^(-precision + 1)) && (ub = 0)
                    # println("previous ub for z[$bp]: ", z_ub[bp])
                    # println("BT ub for z[$bp]: ", ub)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for z[$bp] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_z_iteration_time = max(end_time, max_z_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "z bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > z_ub[bp]) && (lb = z_lb[bp])
                (!isnan(ub) && ub < z_lb[bp]) && (ub = z_ub[bp])
                isnan(lb) && (lb = z_lb[bp])
                isnan(ub) && (ub = z_ub[bp])

                # z bound-reduction computation
                z_reduction = 0.0
                prev_ub = z_ub[bp]
                prev_lb = z_lb[bp]
                if (ub - lb >= min_bound_width)
                    z_reduction = (z_ub[bp] - z_lb[bp]) - (ub - lb)
                    z_lb[bp] = lb
                    z_ub[bp] = ub
                else
                    mean = (lb + ub)/2.0
                    if (mean - min_bound_width/2.0 < z_lb[bp])
                        lb = z_lb[bp]
                        ub = z_lb[bp] + min_bound_width
                    elseif (mean + min_bound_width/2.0 > z_ub[bp])
                        ub = z_ub[bp]
                        lb = z_ub[bp] - min_bound_width
                    else
                        lb = mean - (min_bound_width/2.0)
                        ub = mean + (min_bound_width/2.0)
                    end
                    z_reduction = (z_ub[bp] - z_lb[bp]) - (ub - lb)
                    z_lb[bp] = lb
                    z_ub[bp] = ub
                end

                bool_z_bound_reduction = bool_z_bound_reduction || (prev_ub - z_ub[bp] > improvement_tol) || (z_lb[bp] - prev_lb > improvement_tol)
                bool_z_range = bool_z_range || (z_ub[bp] - z_lb[bp] > min_bound_width)
                # println("final z[$bp] lb: ", z_lb[bp])
                # println("final z[$bp] ub: ", z_ub[bp])
                # println("prev_ub - z_ub[bp]: ", prev_ub - z_ub[bp])
                # println("z_lb[bp] - prev_lb: ", z_lb[bp] - prev_lb)
                # println("z_ub[bp] - z_lb[bp]", z_ub[bp] - z_lb[bp])
                # println("bool_z_bound_reduction: ", bool_z_bound_reduction)
                # println("bool_z_range: ", bool_z_range)

                total_z_reduction += (z_reduction)
                max_z_reduction = max(z_reduction, max_z_reduction)

                modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
                PowerModels.update_data!(data, modifications)

                # run the qc relaxation for the updated bounds
                model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)

                # Rebuild the bound-tightening model
                if upper_bound_constraint
                    @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
                end
            end
            avg_z_reduction = total_z_reduction/length(buspairs)
            z_range_final = sum([z_ub[bp] - z_lb[bp] for bp in buspairs])
        end

        # bound-tightening for the z variables
        if params["obbt_zc"]
            for cyc in cycle_3_4
                (zc_ub[cyc] - zc_lb[cyc] < min_bound_width) && (continue)

                start_time = time()
                # z lower bound solve
                lb = NaN
                @objective(model, Min, model[:zc][cyc])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nlb = floor(10.0^precision * objective_value(model))/(10.0^precision)
                    (nlb > 10.0^(-precision + 1)) && (lb = 1)
                    # println("previous lb for z[$cyc]: ", z_lb[cyc])
                    # println("BT lb for z[$cyc]: ", lb)
                else
                    Memento.warn(_LOGGER, "BT minimization problem for zc[$cyc] errored - change tolerances")
                    continue
                end

                # z upper bound solve
                ub = NaN
                @objective(model, Max, model[:zc][cyc])
                JuMP.optimize!(model)
                if (termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL)
                    nub = ceil(10.0^precision * objective_value(model))/(10.0^precision)
                    (nub < 1.0 - 10.0^(-precision + 1)) && (ub = 0)
                    # println("previous ub for z[$bp]: ", z_ub[bp])
                    # println("BT ub for z[$bp]: ", ub)
                else
                    Memento.warn(_LOGGER, "BT maximization problem for zc[$cyc] errored - change tolerances.")
                    continue
                end
                end_time = time() - start_time
                max_zc_iteration_time = max(end_time, max_zc_iteration_time)

                # sanity checks
                (lb > ub) && (Memento.warn(_LOGGER, "zc bt lb > ub - adjust tolerances in optimizer to avoid issue"); continue)
                (!isnan(lb) && lb > zc_ub[cyc]) && (lb = zc_lb[cyc])
                (!isnan(ub) && ub < zc_lb[cyc]) && (ub = zc_ub[cyc])
                isnan(lb) && (lb = zc_lb[cyc])
                isnan(ub) && (ub = zc_ub[cyc])

                # z bound-reduction computation
                zc_reduction = 0.0
                prev_ub = zc_ub[cyc]
                prev_lb = zc_lb[cyc]
                if (ub - lb >= min_bound_width)
                    zc_reduction = (zc_ub[cyc] - zc_lb[cyc]) - (ub - lb)
                    zc_lb[cyc] = lb
                    zc_ub[cyc] = ub
                else
                    mean = (lb + ub)/2.0
                    if (mean - min_bound_width/2.0 < zc_lb[cyc])
                        lb = zc_lb[cyc]
                        ub = zc_lb[cyc] + min_bound_width
                    elseif (mean + min_bound_width/2.0 > zc_ub[cyc])
                        ub = zc_ub[cyc]
                        lb = zc_ub[cyc] - min_bound_width
                    else
                        lb = mean - (min_bound_width/2.0)
                        ub = mean + (min_bound_width/2.0)
                    end
                    zc_reduction = (zc_ub[cyc] - zc_lb[cyc]) - (ub - lb)
                    zc_lb[cyc] = lb
                    zc_ub[cyc] = ub
                end

                bool_zc_bound_reduction = bool_zc_bound_reduction || (prev_ub - zc_ub[cyc] > improvement_tol) || (zc_lb[cyc] - prev_lb > improvement_tol)
                bool_zc_range = bool_zc_range || (zc_ub[cyc] - zc_lb[cyc] > min_bound_width)
                # println("final z[$bp] lb: ", z_lb[bp])
                # println("final z[$bp] ub: ", z_ub[bp])
                # println("prev_ub - z_ub[bp]: ", prev_ub - z_ub[bp])
                # println("z_lb[bp] - prev_lb: ", z_lb[bp] - prev_lb)
                # println("z_ub[bp] - z_lb[bp]", z_ub[bp] - z_lb[bp])
                # println("bool_z_bound_reduction: ", bool_z_bound_reduction)
                # println("bool_z_range: ", bool_z_range)

                total_zc_reduction += (zc_reduction)
                max_zc_reduction = max(zc_reduction, max_zc_reduction)

                modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
                PowerModels.update_data!(data, modifications)

                # run the qc relaxation for the updated bounds
                model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)

                # Rebuild the bound-tightening model
                if upper_bound_constraint
                    @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
                end
            end
            avg_zc_reduction = total_zc_reduction/length(buspairs)
            zc_range_final = sum([zc_ub[cyc] - zc_lb[cyc] for cyc in cycle_3_4])
        end

        parallel_time_elapsed += max(max_vm_iteration_time, max_td_iteration_time, max_z_iteration_time, max_zc_iteration_time)

        time_elapsed += (time() - iter_start_time)

        # populate the modifications, update the data
        # modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
        # PowerModels.update_data!(data, modifications)

        # run the qc relaxation for the updated bounds
        model = _build_ots_model(data, model, true, z_lb, z_ub, zc_lb, zc_ub)
        JuMP.optimize!(model)

        if termination_status(model) in status_pass
            current_rel_gap = (upper_bound - objective_value(model)) / upper_bound
            final_relaxation_objective = objective_value(model)
        else
            Memento.warn(_LOGGER, "relaxation solve failed in iteration $(current_iteration+1)")
            Memento.warn(_LOGGER, "using the previous iteration's gap to check relative gap stopping criteria")
        end

        # Rebuild the bound-tightening model
        if upper_bound_constraint
            @constraint(model, ub_con, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * model[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]) <= upper_bound)
        end

        Memento.info(_LOGGER, "iteration $(current_iteration+1), vm range: $vm_range_final, td range: $td_range_final, relaxation obj: $final_relaxation_objective")

        # termination criteria update
        (termination == :avg) && (check_termination = (avg_vm_reduction > improvement_tol || avg_td_reduction > improvement_tol))
        (termination == :max) && (check_termination = (max_vm_reduction > improvement_tol || max_td_reduction > improvement_tol))
        (termination == :gop) && (check_termination = ((bool_td_bound_reduction || bool_vm_bound_reduction) && (bool_td_range || bool_vm_range)))
        # interation counter update
        current_iteration += 1
        # check all the stopping criteria
        (current_iteration >= max_iter) && (Memento.info(_LOGGER, "maximum iteration limit reached"); break)
        (time_elapsed > time_limit) && (Memento.info(_LOGGER, "maximum time limit reached"); break)
        if (!isinf(rel_gap_tol)) && (current_rel_gap < rel_gap_tol)
            Memento.info(_LOGGER, "relative optimality gap < $rel_gap_tol")
            break
        end

    end

    println("time3: $(time() - time3_init)")

    time4_init = time()

    branches_vad_same_sign_count = 0
    for (key, branch) in data["branch"]
        is_same_sign = (branch["angmax"] >=0 && branch["angmin"] >= 0) || (branch["angmax"] <=0 && branch["angmin"] <= 0)
        (is_same_sign) && (branches_vad_same_sign_count += 1)
    end

    stats["final_relaxation_objective"] = final_relaxation_objective
    stats["final_rel_gap_from_ub"] = isnan(upper_bound) ? Inf : current_rel_gap
    stats["vm_range_final"] = vm_range_final
    stats["avg_vm_range_final"] = vm_range_final/length(buses)

    stats["td_range_final"] = td_range_final
    stats["avg_td_range_final"] = td_range_final/length(buspairs)

    stats["run_time"] = time_elapsed
    stats["iteration_count"] = current_iteration
    stats["sim_parallel_run_time"] = parallel_time_elapsed

    stats["vad_sign_determined"] = branches_vad_same_sign_count

    modifications = _create_modifications(data, vm_lb, vm_ub, td_lb, td_ub)
    PowerModels.update_data!(data, modifications)
    model = _build_ots_model(data, model, false, z_lb, z_ub, zc_lb, zc_ub)
    # set_optimizer_attributes(model, "Method" => 1) # Set method to primal simplex.

    println("time4: $(time() - time4_init)")
    println("time1234:  $(time() - time1_init)")

    return data, stats, model

end

function _check_variables(pm::AbstractPowerModel)
    try
        vm = var(pm, :vm)
    catch err
        (isa(error, KeyError)) && (Memento.error(_LOGGER, "OBBT is not supported for models without explicit voltage magnitude variables"))
    end

    try
        td = var(pm, :td)
    catch err
        (isa(error, KeyError)) && (Memento.error(_LOGGER, "OBBT is not supported for models without explicit voltage angle difference variables"))
    end
end

function _check_obbt_options(ub::Float64, rel_gap::Float64, ub_constraint::Bool)
    if ub_constraint && isinf(ub)
        Memento.error(_LOGGER, "the option upper_bound_constraint cannot be set to true without specifying an upper bound")
    end

    if !isinf(rel_gap) && isinf(ub)
        Memento.error(_LOGGER, "rel_gap_tol is specified without providing an upper bound")
    end
end

# function _constraint_obj_bound(model_bt, bound)
#     # model = PowerModels.check_cost_models(pm)
#     # if model != 2
#     #     Memento.error(_LOGGER, "Only cost models of type 2 is supported at this time, given cost model type $(model)")
#     # end
#     #
#     # cost_index = PowerModels.calc_max_cost_index(pm.data)
#     # if cost_index > 3
#     #     Memento.error(_LOGGER, "Only quadratic generator cost models are supported at this time, given cost model of order $(cost_index-1)")
#     # end
#     #
#     # PowerModels.standardize_cost_terms!(pm.data, order=2)
#     #
#     # from_idx = Dict(arc[1] => arc for arc in ref(pm, :arcs_from_dc))
#
#     JuMP.@constraint(model, sum(gen["cost"][1] * model[:pg][i]^2  + gen["cost"][2] * model[:pg][i] + gen["cost"][3] for (i,gen) in ref[:gen]) <= bound)
# end

function _create_modifications(data,
    vm_lb::Dict{Any,Float64}, vm_ub::Dict{Any,Float64},
    td_lb::Dict{Any,Float64}, td_ub::Dict{Any,Float64})

    modifications = Dict{String,Any}()

    modifications["per_unit"] = true
    modifications["bus"] = Dict{String,Any}()
    modifications["branch"] = Dict{String,Any}()

    for bus in keys(ref[:bus])
        # index = string(ref(pm, :bus, bus, "index"))
        index = string(bus)
        modifications["bus"][index] = Dict("vmin" => vm_lb[bus], "vmax" => vm_ub[bus])
    end

    for branch in keys(ref[:branch]) #ids(pm, :branch)
        # index = string(ref(pm, :branch, branch, "index"))
        index = string(branch)
        # f_bus = ref(pm, :branch, branch, "f_bus")
        # t_bus = ref(pm, :branch, branch, "t_bus")
        f_bus = data["branch"][index]["f_bus"]
        t_bus = data["branch"][index]["t_bus"]
        bp = (f_bus, t_bus)
        modifications["branch"][index] = Dict("angmin" => td_lb[bp], "angmax" => td_ub[bp])
    end

    return modifications
end

# Build OPF model.
function _build_model(data, _m)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    _m = Get_solver(params)

    (params["obbt_solver_log"]) ? (print_level = 1) : (print_level = 0)
    set_optimizer_attributes(_m, "OutputFlag" => print_level)
   #---------------;
   #  Variables    ;
   #---------------;
   vars = variables(_m, ref, params)
   include("core/Get_variables.jl")

   #---------------;
   #  Objective    ;
   #---------------;
   @objective(_m, Min, sum(gen["cost"][1]* _m[:pg][i]^2  + gen["cost"][2]*_m[:pg][i] + gen["cost"][3] for (i,gen) in ref[:gen]))

   #---------------;
   # Constraints   ;
   #---------------;
   # include("core/constraints.jl")
   create_constraints(_m, ref, params, bd_data)

   # Cycle cuts
   if params["cycle_cuts"] && !params["callback"]
      if params["separation"]
         include("core/util_sep.jl")
      else
         # include("core/Get_cycle_constraints.jl")
         create_cyc_constr(_m, ref, params)
      end
   end

   # if params["print_model"]
   #    Write_model(m, "results/model.txt")
   # end
   return _m
end

# Build OTS model.
function _build_ots_model(data, _m, relax, z_lb, z_ub, zc_lb, zc_ub)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    _m = Get_solver(params)

    (params["obbt_solver_log"]) ? (print_level = 1) : (print_level = 0)
    set_optimizer_attributes(_m, "OutputFlag" => print_level)
   #---------------;
   #  Variables    ;
   #---------------;
   if relax
       params["model"] = "ots_relax"
       vars, bd_data = variables(_m, ref, params)
       params["model"] = "ots"
   else
       vars, bd_data = variables(_m, ref, params)
   end
   include("core/Get_variables.jl")
   for bp in keys(ref[:buspairs])
       if z_lb[bp] == 1
           set_lower_bound(_m[:z][bp], 1)
       end
       if z_ub[bp] == 0
           set_upper_bound(_m[:z][bp], 0)
       end
   end
   if params["cycle_cuts"] && !params["separation"] && !params["callback"] # && !params["off_insights"]
       cycle_3_4 = []
       append!(cycle_3_4, cycle["cyc_3"])
       if params["cycle_max_bnd"] >= 4
           append!(cycle_3_4, cycle["cyc_4"])
       end
       for i in 1:length(cycle_3_4)
           cycle_3_4[i] = Tuple(cycle_3_4[i])
       end
       for cyc in cycle_3_4
           if zc_lb[cyc] == 1
               set_lower_bound(_m[:zc][cyc], 1)
           end
           if zc_ub[cyc] == 0
               set_upper_bound(_m[:zc][cyc], 0)
           end
       end
   end

   #---------------;
   #  Objective    ;
   #---------------;
   leaf_gens(ref)
   @objective(_m, Min, sum(gen["cost"][1]*_m[:pg][i]^2  + gen["cost"][2]*_m[:pg][i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * _m[:z][bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]))

   #---------------;
   # Constraints   ;
   #---------------;
   # include("core/constraints.jl")
   create_constraints(_m, ref, params, bd_data)

   # Cycle cuts
   if params["cycle_cuts"] && !params["separation"] && !params["callback"] # && !params["off_insights"]
      create_cyc_constr(_m, ref, params)
   end
   if params["cycle_cuts"] && params["callback"] && params["run_obbt"] && relax
        create_cyc_constr(_m, ref, params)
    end

   # if params["print_model"]
   # Write_model(_m, "results/model.txt")
   # end
   return _m
end
