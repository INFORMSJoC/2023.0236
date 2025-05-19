# Directly solving problems using PowerModels.jl. For ACOPF: exact acopf obejctive. For ACOTS: QCRMPowerModel.
include("core/solver.jl")
include("core/relaxations.jl")
include("core/network_data.jl")
include("core/variables.jl")
include("core/utility.jl")
instance = params["instance"]
data = network_data(instance, params)
PowerModels.standardize_cost_terms!(data, order=2)
pm_model = nothing

if params["run_obbt"]
  println("Running OBBT ...")
  optimizer = with_optimizer(Ipopt.Optimizer, print_level=0, sb="yes");
  data, stats = run_obbt_opf!(data, optimizer, upper_bound = 4.0e+06, min_bound_width = 1E-3, improvement_tol = 1E-4) #,upper_bound = 3.33e+06
  ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
else
  ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
end

if params["opf_exact_pm"]
  if params["solver_log"]
    nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer)
  else
    nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
  end
  # set_optimizer_attributes(m_ac, "print_level" => 0)
  single_run_time = @elapsed result = solve_ac_opf(data, nlp_solver)
elseif params["opf_qc_pm"]
  minlp_solver = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => params["time_limit"])
  single_run_time = @elapsed result = solve_opf(data, QCLSPowerModel, minlp_solver)
  # single_run_time = @elapsed result = run_opf(data, QCRMPowerModel, minlp_solver)
elseif params["ots_exact_pm"]
  if params["minlp_solver"] == "juniper"
    minlp_solver = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-4, "print_level"=>0), "time_limit" => params["time_limit"]) #, "log_levels"=>[]
  elseif params["minlp_solver"] == "knitro"
    minlp_solver = JuMP.optimizer_with_attributes(KNITRO.Optimizer, MOI.TimeLimitSec() => params["time_limit"]) #infeastol => 1e-9 (default: 1e-8), "outlev" => 0, , "maxtime_real" => params["time_limit"]
  end

  single_run_time = @elapsed result = solve_ots(data, ACPPowerModel, minlp_solver)
elseif params["ots_qc_pm"]
  if params["get_root_relax"]
    minlp_solver = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => params["time_limit"], JuMP.MOI.RawOptimizerAttribute("NodeLimit") => 0)
  elseif params["presolve_thread"] == true
    minlp_solver = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => params["time_limit"], "FeasibilityTol" => params["grb_feas_tol"]) #, "print_level"=>0, "sb"=>"yes" "MIPGap" => 0.0005
    println("minlp_solver1: ", minlp_solver)
  else
    minlp_solver = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => params["time_limit"], "Threads" => 1, "Presolve" => 0)
    println("minlp_solver2: ", minlp_solver)
  end
  # single_run_time = @elapsed result, pm_model = run_ots(data, QCRMPowerModel, minlp_solver)
  single_run_time = @elapsed result = solve_ots(data, QCRMPowerModel, minlp_solver)
end
println("***************************************")
println("Case: $instance")
if params["get_root_relax"]
  println("Objective lb: ", result["objective_lb"])
else
  println("Objective: ", result["objective"])
end
println("Single_run_time: ", single_run_time)
println("***************************************")
append!(num_nodes, length(data["bus"]))
append!(num_edges, length(data["branch"]))
if params["get_root_relax"]
  append!(obj_value, result["objective_lb"])
else
  append!(obj_value, result["objective"])
end
push!(solution_status, string(result["termination_status"]))
push!(prim_status, string(result["primal_status"]))
append!(run_time_li, single_run_time)
if params["model"] == "ots" || params["model"] == "ots_relax"
  if result["termination_status"] == NO_SOLUTION || result["primal_status"] == INFEASIBLE_POINT || result["primal_status"] == NO_SOLUTION # Juniper could stop at an infeasible point and still provide a solution.
    num_off_lines = " "
    off_lines_li = " "
    opt_gap = " "
    if !params["get_root_relax"]
      obj_value[end] = NaN
    end
  else
    switched_off = []
    for bus in keys(data["bus"])
        on_off = result["solution"]["branch"][string(bus)]["br_status"]
        if on_off == 0.0
            i = data["branch"][string(bus)]["f_bus"]
            j = data["branch"][string(bus)]["t_bus"]
            bp = (i, j)
            push!(switched_off, bp)
        end
    end
    num_off_lines = length(switched_off)
    off_lines_li = switched_off
    opt_gap = (result["objective"] - result["objective_lb"]) / result["objective"]

    # check_pm_mcc(result)
  end
 end
