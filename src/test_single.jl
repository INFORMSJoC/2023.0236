include("core/solver.jl")
include("core/utility.jl")
include("core/relaxations.jl")
include("core/network_data.jl")
include("core/variables.jl")
include("core/constraints.jl")
include("core/Get_cycle_constraints.jl")
include("core/obbt_cyc.jl")
instance = params["instance"]
data = network_data(instance, params)
PowerModels.standardize_cost_terms!(data, order=2)
bd_data = Dict()

# if params["run_obbt"]
#   println("Running OBBT ...")
#   optimizer = with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes");
#   data, stats = run_obbt_opf!(data, optimizer, upper_bound = 4.0e+06, min_bound_width = 1E-3, improvement_tol = 1E-4) #,upper_bound = 3.33e+06
#   ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
# else
#   ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
# end

ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

if params["compare_cyc_constr"] != 0
   for bp in keys(ref[:buspairs])
      ref[:buspairs][bp]["angmin"] = - params["compare_cyc_constr"]
      ref[:buspairs][bp]["angmax"] = params["compare_cyc_constr"]
   end
end

cycle = Dict()
if params["cycle_cuts"]
   cycle = Get_cycles(ref, params)
   # cycle = Get_cycle_basis(ref, params)
   println(">>>>>> Enumerated cycles with maximum size of $(params["cycle_max_bnd"]) <<<<<<")
   # @show cycle
end

#------------------------;
#  Optimization model    ;
#------------------------;
params["optimize"] && (params["formulate_problem"] = true)

if params["formulate_problem"]

  #----------------------#
  # Model initialization
  #----------------------#
  m = Get_solver(params)

  #---------------;
  #  Variables    ;
  #---------------;
  vars, bd_data = variables(m, ref, params)
  include("core/Get_variables.jl")
  if params["naive_bt"]
     ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
     include("naive_bound.jl")
  end

  #---------------;
  #  Objective    ;
  #---------------;
  if params["model"] == "opf"
     JuMP.@objective(m, Min, sum(gen["cost"][1]*pg[i]^2  + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]))
  elseif params["model"] == "ots" || params["model"] == "ots_relax"
     # JuMP.@objective(m, Min, sum(gen["cost"][1]*pg[i]^2  + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]))
     leaf_gens(ref)
     JuMP.@objective(m, Min, sum(gen["cost"][1]*pg[i]^2  + gen["cost"][2]*pg[i] for (i,gen) in ref[:gen]) + sum(ref[:gen][g]["cost"][3] * z[bp] for (g, bp) in ref[:gens_on_leaf]) + sum(ref[:gen][g]["cost"][3] for g in ref[:gens_not_on_leaf]))
  end

  #---------------;
  # Constraints   ;
  #---------------;
  create_constraints(m, ref, params, bd_data)

  # Cycle cuts
  single_run_time = 0
  if params["cycle_cuts"] && !params["separation"] && !params["callback"]
      create_cyc_constr(m, ref, params)
   end
   if params["cycle_cuts"] && params["callback"] && params["run_obbt"]
      create_cyc_constr(m, ref, params)
   end
#   if params["cycle_cuts"]
#      if params["separation"]
#         # include("core/separation2.jl")
#         # single_run_time += @elapsed separation()
#         include("core/util_sep.jl")
#         single_run_time += @elapsed include("core/separation.jl") #include("core/separation.jl")
#      elseif params["callback"]
#         include("core/util_sep.jl")
#         include("core/util_cb.jl")
#         single_run_time += include("core/callback.jl")
#      else
#         create_cyc_constr(m, ref, params)
#      end
#   end

   # OBBT
  stats = nothing
  if params["run_obbt"] #&& params["cycle_cuts"]
     if params["model"] == "opf"
        data, stats, m = run_cyc_obbt_opf!(data, m, termination = params["obbt_termination"], upper_bound_constraint=true, rel_gap_tol = 1e-5, min_bound_width = 1e-5, improvement_tol = 1e-4)
     elseif params["model"] == "ots"
        if params["obbt_perVar"]
           obbt_func_time = @elapsed data, stats, m = run_cyc_obbt_ots_perRound!(data, m, termination = params["obbt_termination"], upper_bound_constraint=false, rel_gap_tol = 1e-5, min_bound_width = 1e-5, improvement_tol = 1e-4)
        else
            if params["get_root_relax"]
               MOI.set(m, MOI.RawOptimizerAttribute("NodeLimit"), Inf)
            end # Set back to default before OBBT
            obbt_func_time = @elapsed data, stats, m = run_cyc_obbt_ots!(data, m, termination = params["obbt_termination"], upper_bound_constraint=false, rel_gap_tol = 1e-5, min_bound_width = 1e-5, improvement_tol = 1e-4)
            if params["get_root_relax"]
               MOI.set(m, MOI.RawOptimizerAttribute("NodeLimit"), 0)
            end
            if params["solver_log"]
               set_optimizer_attributes(m, "OutputFlag" => 1)
            end
         end
     end
  end

  # Cycle cuts via separation or callback.
  if params["cycle_cuts"]
     if params["separation"]
        # include("core/separation2.jl")
        # single_run_time += @elapsed separation()
        include("core/util_sep.jl")
        single_run_time += @elapsed include("core/separation.jl") #include("core/separation.jl")
     elseif params["callback"]
        include("core/util_sep.jl")
        include("core/util_cb.jl")
        single_run_time += include("core/callback.jl")
   #   else
   #      create_cyc_constr(m, ref, params)
     end
  end

  #----------;
  # Solve    ;
  #----------;
  if params["optimize"]
     if params["separation"] == false && params["callback"] != true
        single_run_time = @elapsed JuMP.optimize!(m)
     end
     println("***************************************")
     println("Case: $instance")
     println("Primal status: ", primal_status(m))

     if params["get_root_relax"]
         println("Root relaxation for $(params["model"]): ", objective_bound(m))
     elseif primal_status(m) != NO_SOLUTION
         println("QC-$(params["model"]) objective: ", JuMP.objective_value(m))
     else
        println("No feasible solution found.")
     end
     println("Single_run_time: ", single_run_time)
     println("***************************************")
  end

  if params["get_root_relax"]
      append!(obj_value, objective_bound(m))
  elseif primal_status(m) != NO_SOLUTION
      append!(obj_value, JuMP.objective_value(m))
  else
     push!(obj_value, "NaN")
  end
  push!(solution_status, string(termination_status(m)))
  append!(run_time_li, single_run_time)
  if params["separation"]
     append!(master_time_li, specs["mp_time"])
     append!(sep_time_li, specs["sp_time"])
     append!(num_cuts_li, specs["num_cuts"])
     append!(num_cs_cuts_li, specs["num_cs_cuts"])
     append!(num_w_cuts_li, specs["num_w_cuts"])
     append!(num_iter_li, specs["num_iterations"])
   end
   if params["run_obbt"]
      append!(obbt_time_li, stats["run_time"])
      append!(obbt_func_time_li, obbt_func_time)
      append!(obbt_iter_li, stats["iteration_count"])
      percent_vm_reduction = stats["percent_vm_range_reduction"]
      percent_td_reduction = stats["percent_td_range_reduction"]
      if params["obbt_z"]
         append!(obbt_freeLine_li, stats["ratio_free_lines"])
         append!(obbt_zeroLine_li, stats["num_zero_lines"])
      end
   end
   if params["model"] == "ots" || params["model"] == "ots_relax"
      if primal_status(m) != NO_SOLUTION
         z_bar = JuMP.value.(z)
         switched_off = []
         switched_off_leaf = []
         for bp in keys(ref[:buspairs])
            if z_bar[bp] < 0.5
               push!(switched_off, bp)
            end
         end
         for (g, bp) in ref[:gens_on_leaf]
            if (z_bar[bp] < 0.5) && ((bp in switched_off_leaf) == false)
               push!(switched_off_leaf, bp)
            end
         end
         opt_gap = MOI.get(m, MOI.RelativeGap())
         num_off_lines = length(switched_off)
         off_lines_li = switched_off
         num_off_lines_leaf = length(switched_off_leaf)
         off_lines_leaf_li = switched_off_leaf
      else
         opt_gap = " "
         num_off_lines = " "
         off_lines_li = " "
         num_off_lines_leaf = " "
         off_lines_leaf_li = " "
      end
      if params["callback"]
         append!(total_time_li, specs["total_time"])
         append!(sep_time_li, specs["sp_time"])
         append!(sep_build_time_li, specs["sp_total_time"] - specs["sp_time"])
         append!(num_cuts_li, specs["num_cuts"])
         append!(num_cs_cuts_li, specs["num_cs_cuts"])
         append!(num_w_cuts_li, specs["num_w_cuts"])
         append!(num_iter_li, specs["num_iterations"])
      end
   end
   if params["write_sol"] && params["model"] == "opf"
      string_cycle = ""
      if params["cycle_cuts"]
         string_cycle = "cycle_"
      end
      string_load = ""
      if params["load_multiplier"] != 1
         string_load = string("_load", params["load_multiplier"])
      end
      string_circ2 = ""
      if params["circ_cons"]
         string_circ2 = "circ_"
      end
      string_rmc2 = ""
      if params["relax"] == "rmc"
         string_rmc2 = "rmc_"
      end
      output_path = string("results/opf_sol/", params["instance"], "_", string_cycle, string_load, string_circ2, string_rmc2, "sol.txt")
      # li = []
      dict1 = Dict()
      for item in keys(vars)
         the_var = vars[item]
         temp_dict = Dict()
         if params["cycle_cuts"] == false
            arr = JuMP.value.(vars[item])
            for ind in keys(the_var)
               if length(size(the_var)) == 1
                  temp_dict[ind[1]] = arr[ind]
               elseif length(size(the_var)) == 2
                  temp_dict[ind[1], ind[2]] = arr[ind]
               end
            end
         else
            for item2 in the_var
               temp_dict[string(item2)] = JuMP.value(item2)
            end
         end
         dict1[item] = temp_dict
      end
      open(output_path, "w") do f
         write(f, string(dict1))
      end
   end
end
