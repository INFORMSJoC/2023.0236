using JuMP
# using CPLEX
using Gurobi
using Ipopt, Juniper, KNITRO
using PowerModels
using CSV, DataFrames
using LightGraphs
using Memento
using Printf
setlevel!(getlogger(PowerModels), "error")
const GRB_ENV = Gurobi.Env()

#--------------------------;
#  Load network data       ;
#--------------------------;
inst_li = CSV.read("data/instance_names/pglib_sad.csv", DataFrame, header = false)[!,"Column1"][1:4] #[1:34] [1:17] [1:24]
append!(inst_li, CSV.read("data/instance_names/pglib_api.csv", DataFrame, header = false)[!,"Column1"][1:4])#[1:17] [1:24]
append!(inst_li, CSV.read("data/instance_names/pglib_typ.csv", DataFrame, header = false)[!,"Column1"][1:4])#[1:17] [1:24]
# inst_li = CSV.read("data/instance_names/nesta_sad.csv", DataFrame, header = false)[!,"Column1"]
# append!(inst_li, CSV.read("data/instance_names/nesta_sad.csv", DataFrame, header = false)[!,"Column1"])
# append!(inst_li, CSV.read("data/instance_names/nesta_api.csv", DataFrame, header = false)[!,"Column1"])
# append!(inst_li, CSV.read("data/instance_names/nesta_typ.csv", DataFrame, header = false)[!,"Column1"])
# append!(inst_li, CSV.read("data/instance_names/nesta_extra.csv", DataFrame, header = false)[!,"Column1"])
# inst_li = CSV.read("data/instance_names/useful_inst.csv", DataFrame, header = false)[!,"Column1"]
# inst_li = CSV.read("data/instance_names/pglib_api.csv", DataFrame, header = false)[!,"Column1"]
# append!(inst_li, CSV.read("data/instance_names/pglib_sad.csv", DataFrame, header = false)[!,"Column1"])
# inst_li = CSV.read("data/instance_names/nesta_extra.csv", DataFrame, header = false)[!,"Column1"]

# Terminal arguments for testing
# empty!(ARGS)
function parse_args(li)
   args = []
   if length(ARGS) == 0
      args = li
   else
      push!(args, eval(Meta.parse(ARGS[1])))
      push!(args, eval(Meta.parse(ARGS[2])))
      push!(args, eval(Meta.parse(ARGS[3])))
      push!(args, ARGS[4])
      push!(args, eval(Meta.parse(ARGS[5])))
      push!(args, eval(Meta.parse(ARGS[6])))
      push!(args, eval(Meta.parse(ARGS[7])))
      push!(args, eval(Meta.parse(ARGS[8])))
      push!(args, eval(Meta.parse(ARGS[9])))
      push!(args, eval(Meta.parse(ARGS[10])))
      # for item in ARGS
      #    push!(args, eval(Meta.parse(item)))
      # end
   end

   return args
end
# Input: [index of start instances, index of finish instance, grb_feas_tol, relax, warm_start, repeat_first_inst, time_limit(hour), run_obbt, cycle_cuts, callback]
args = parse_args([1, 12, 1e-6, "tri", false, false, 2, false, false, false]) # tri, rmc # 1, 17, 24*3 = 72
# length(inst_li)

#-----------------;
#  Functions      ;
#-----------------;
include("core/params.jl")
# Remove unlisted instances.
# if params["get_root_relax"]
   # inst_li
# filter!(x -> x != "pglib_opf_case162_ieee_dtc__sad", inst_li) # 12
# filter!(x -> x != "pglib_opf_case162_ieee_dtc", inst_li) # 44
# filter!(x -> x != "pglib_opf_case300_ieee__sad", inst_li) # 16
# filter!(x -> x != "pglib_opf_case300_ieee", inst_li) # 48
# end
# Make sure the right model is chosen.
if params["opf_exact_pm"]
   @assert(params["model"] == "opf")
   @assert(params["ots_exact_pm"] == false)
   @assert(params["ots_qc_pm"] == false)
end
if params["ots_exact_pm"]
   @assert(params["model"] == "ots" || params["model"] == "ots_relax")
   @assert(params["opf_exact_pm"] == false)
   @assert(params["ots_qc_pm"] == false)
end
if params["ots_qc_pm"]
   @assert(params["model"] == "ots" || arams["model"] == "ots_relax")
   @assert(params["ots_exact_pm"] == false)
   @assert(params["opf_exact_pm"] == false)
end
# @assert(params["turn_off_lines"] == false)

num_nodes = []
num_edges = []
obj_value = []
run_time_li = []
solution_status = []
prim_status = []
master_time_li = []
sep_time_li = []
sep_build_time_li = []
num_cuts_li = []
num_cs_cuts_li = []
num_w_cuts_li = []
num_iter_li = []
obbt_time_li = []
obbt_func_time_li = []
obbt_iter_li = []
obbt_freeLine_li = []
obbt_zeroLine_li = []
num_off_lines = 0
off_lines_li = []
num_off_lines_leaf = 0
off_lines_leaf_li = []
opt_gap = 0.0
total_time_li = []
percent_vm_reduction = 0.0
percent_td_reduction = 0.0

density_li = []
num_3_cycle_li = []
num_4_cycle_li = []

# start = eval(Meta.parse(ARGS[1]))
# finish = eval(Meta.parse(ARGS[2]))
start = args[1]
finish = args[2]
# start = 1 # 1 9
# finish = 1 #length(inst_li) #34, length(inst_li) 10
inst_li = inst_li[start:finish]
if params["repeat_first_inst"]
   insert!(inst_li, 1, inst_li[1]) # repeat the first instance for more accurate computation time.
end

string_cycle_type = ""
if params["cycle_c_s_cuts"] && params["cycle_wr_wi_cuts"]
   string_cycle_type = "both"
elseif params["cycle_c_s_cuts"] && params["cycle_wr_wi_cuts"] == false
   string_cycle_type = "cs"
else
   string_cycle_type = "w"
end
string_bt = ""
if params["run_obbt"]
   string_bt = "_obbt"
elseif params["naive_bt"]
   string_bt = "_naivebt"
end
string_rotate = ""
if params["rotate"] == false
   string_rotate = "_noRotation"
elseif params["rotate_4w"] == false && params["cycle_max_bnd"] >= 4
   string_rotate = "_noRotation4wepr"
end
string_compare = ""
if params["compare_cyc_constr"] != 0
   degree = 120
   # degree = Int64(params["compare_cyc_constr"] * 180 / pi)
   string_compare = "_td$(degree)"
end
string_solver = ""
if params["solver"] == "gurobi"
   string_solver = "-grb"
end
string_separation = ""
if params["separation"]
   string_separation = "_sep"
end
string_ots = ""
if params["model"] == "ots"
   string_ots = "_ots"
end
string_cb = ""
if params["callback"]
   string_cb = "ots_cb_"
end
string_relax = ""
if params["relax"] == "rmc"
   string_relax = "_rmc"
end
string_load_multiplier = ""
if params["load_multiplier"] != 1
   string_load_multiplier = string("_load", params["load_multiplier"])
end
string_circ = ""
if params["circ_cons"]
   string_circ = "_circ"
end
string_spanTree = ""
if params["span_tree"]
   if params["cal_weight"]
      string_spanTree = "_spanTreeWeight"
   else
      string_spanTree = "_spanTree"
   end
end
string_root_relax = ""
if params["get_root_relax"]
   string_root_relax = "_rootRelax"
end
# output_path = string("results/", params["cycle_max_bnd"], "cycle_results_", params["cycle_relax"], "_", string_cycle_type, string_rotate, string_compare, string_solver, string_separation, ".csv")
# inst_li = ["pglib_opf_case162_ieee_dtc__sad", "pglib_opf_case179_goc__sad", "pglib_opf_case200_activ__sad", "pglib_opf_case500_goc__sad"]
inst_ind = start - 1
cnt = 0
for instance in inst_li
   global inst_ind += 1
   global cnt += 1
   if params["repeat_first_inst"] && cnt == 2
      inst_ind -= 1
   end
   params["instance"] = instance
   if params["opf_exact_pm"] || params["opf_qc_pm"]
      include("test_single_PM.jl")
      if params["opf_exact_pm"]
         output_path = string("results/exact", string_bt, "_results" , "_", params["correct_bus_numbering"], "_", params["correct_branch_directions"],".csv")
      elseif params["opf_qc_pm"]
         output_path = string("results/opf_qc_pm", string_bt, "_results" , "_", params["correct_bus_numbering"], "_", params["correct_branch_directions"],".csv")
      end
      CSV.write(output_path, DataFrame(Instance = instance, Nodes = num_nodes[end], Edges = num_edges[end], Objective_value = obj_value[end], Run_time = run_time_li[end], Solution_status = solution_status[end], inst_ind = inst_ind); append = true, header = params["header"])
   elseif params["ots_exact_pm"] || params["ots_qc_pm"]
      include("test_single_PM.jl")
      if params["ots_exact_pm"]
         output_path = string("results/ots_exact_pm", string_bt, "_results" , "_", params["correct_bus_numbering"], "_", params["correct_branch_directions"],".csv")
      elseif params["ots_qc_pm"]
         output_path = string("results/ots_qc_pm", string_bt, "_results" , "_", params["correct_bus_numbering"], "_", params["correct_branch_directions"], string_root_relax, ".csv")
      end
      CSV.write(output_path, DataFrame(Instance = instance, Nodes = num_nodes[end], Edges = num_edges[end], Objective_value = obj_value[end], Run_time = run_time_li[end], Solution_status = solution_status[end], Primal_status = prim_status[end], Opt_gap = opt_gap, Num_off_lines = num_off_lines, Off_lines = string(off_lines_li), inst_ind = inst_ind); append = true, header = params["header"])
   elseif params["plot_measures"]
      include("test_single_plot.jl")
      output_path = string("results/plot_measure.csv")
      CSV.write(output_path, DataFrame(Instance = instance, Density = density_li[end], Num_3_cycle = num_3_cycle_li[end], Num_4_cycle = num_4_cycle_li[end], Num_3_4_cycle = num_3_cycle_li[end] + num_4_cycle_li[end], inst_ind = inst_ind); append = true, header = params["header"])
   elseif params["cycle_cuts"] == false # No cycle cuts
      include("test_single.jl")
      output_path = string("results/noCycle", string_ots, string_bt, string_relax, "_results" , "_", params["correct_bus_numbering"], "_", params["correct_branch_directions"], string_solver, string_spanTree, string_root_relax, ".csv")
      df = DataFrame(Instance = instance, Objective_value = obj_value[end], Run_time = run_time_li[end], Solution_status = solution_status[end])
      if params["run_obbt"]
         df[!,"Total_time"] = [obbt_time_li[end] + run_time_li[end]]
         df[!,"Optimization_time"] = [run_time_li[end]]
         df[!,"OBBT_time"] = [obbt_time_li[end]]
         df[!,"OBBT_function_time"] = [obbt_func_time_li[end]]
         df[!,"OBBT_iter"] = [obbt_iter_li[end]]
         df[!, "Precent_vm_reduce"] = [percent_vm_reduction]
         df[!, "Precent_td_reduce"] = [percent_td_reduction]
         if params["obbt_z"]
            df[!, "Free_lines"] = [obbt_freeLine_li[end]]
            df[!, "Zero_lines"] = [obbt_zeroLine_li[end]]
         end
      end
      if params["model"] == "ots"
         df[!, "Opt_gap"] = [opt_gap]
         df[!, "num_off_lines"] = [num_off_lines]
         df[!, "off_lines"] = [string(off_lines_li)]
         df[!, "num_off_lines_leaf"] = [num_off_lines_leaf]
         df[!, "off_lines_leaf"] = [string(off_lines_leaf_li)]
      end
      df[!, "inst_ind"] = [inst_ind]
      CSV.write(output_path, df; append = true, header = params["header"])
   elseif params["separation"]
      include("test_single.jl")
      output_path = string("results/", params["cycle_max_bnd"], "cycle_results_", params["cycle_relax"], "_", string_cycle_type, string_bt, string_rotate, string_compare, string_solver, string_separation, ".csv")
      CSV.write(output_path, DataFrame(Instance = instance, Objective_value = obj_value[end], Total_solve_time = master_time_li[end] + sep_time_li[end], Solution_status = solution_status[end], Run_time = run_time_li[end],  Master_time = master_time_li[end], Sep_time = sep_time_li[end], Extra_time = run_time_li[end] - master_time_li[end] - sep_time_li[end], Num_iter = num_iter_li[end], Num_cuts = num_cuts_li[end], Num_cs_cuts = num_cs_cuts_li[end], Num_w_cuts = num_w_cuts_li[end], inst_ind = inst_ind); append = true, header = params["header"])
   elseif params["callback"]
      include("test_single.jl")
      output_path = string("results/", string_cb, params["cycle_max_bnd"], "cycle_results_", params["cycle_relax"], "_", string_cycle_type, string_bt, string_rotate, string_compare, string_solver, string_circ, string_spanTree, string_root_relax, ".csv")
      df = DataFrame(Instance = instance, Objective_value = obj_value[end], Total_time = total_time_li[end], total_time_no_sub_build = total_time_li[end] - sep_build_time_li[end], Solution_status = solution_status[end], Master_time = total_time_li[end] - sep_build_time_li[end] - sep_time_li[end], Sep_build_time = sep_build_time_li[end], Sep_solve_time = sep_time_li[end], Num_iter = num_iter_li[end], Num_cuts = num_cuts_li[end])
      if params["run_obbt"]
         df[!,"Run_time"] = [obbt_time_li[end] + run_time_li[end]]
         df[!,"Optimization_time"] = [run_time_li[end]]
         df[!,"OBBT_time"] = [obbt_time_li[end]]
         df[!,"OBBT_function_time"] = [obbt_func_time_li[end]]
         df[!,"OBBT_iter"] = [obbt_iter_li[end]]
         df[!, "Precent_vm_reduce"] = [percent_vm_reduction]
         df[!, "Precent_td_reduce"] = [percent_td_reduction]
         if params["obbt_z"]
            df[!, "Free_lines"] = [obbt_freeLine_li[end]]
            df[!, "Zero_lines"] = [obbt_zeroLine_li[end]]
         end
      end
      if params["model"] == "ots" || params["model"] == "ots_relax"
         df[!, "Opt_gap"] = [opt_gap]
         df[!, "num_off_lines"] = [num_off_lines]
         df[!, "off_lines"] = [string(off_lines_li)]
         df[!, "num_off_lines_leaf"] = [num_off_lines_leaf]
         df[!, "off_lines_leaf"] = [string(off_lines_leaf_li)]
      end
      df[!, "inst_ind"] = [inst_ind]
      CSV.write(output_path, df; append = true, header = params["header"])
      # CSV.write(output_path, DataFrame(Instance = instance, Objective_value = obj_value[end], Opt_gap = opt_gap, Total_time = total_time_li[end], total_time_no_sub_build = total_time_li[end] - sep_build_time_li[end], Solution_status = solution_status[end], Master_time = total_time_li[end] - sep_build_time_li[end] - sep_time_li[end], Sep_build_time = sep_build_time_li[end], Sep_solve_time = sep_time_li[end], Num_iter = num_iter_li[end], Num_cuts = num_cuts_li[end], Num_cs_cuts = num_cs_cuts_li[end], Num_w_cuts = num_w_cuts_li[end], Num_off_lines = num_off_lines, off_lines = string(off_lines_li), Num_off_lines_leaf = num_off_lines_leaf, off_lines_leaf = string(off_lines_leaf_li), inst_ind = inst_ind); append = true, header = params["header"])
   else # Cycle cuts
      include("test_single.jl")
      if params["model"] == "opf"
         output_path = string("results/", params["cycle_max_bnd"], "cycle_results_", params["cycle_relax"], "_", string_cycle_type, string_bt, string_relax, string_rotate, string_compare, string_solver, string_separation, string_circ, ".csv")
      elseif params["model"] == "ots" || params["model"] == "ots_relax"
         output_path = string("results/ots_", params["cycle_max_bnd"], "cycle_results_", params["cycle_relax"], "_", string_cycle_type, string_bt, string_relax, string_rotate, string_compare, string_solver, string_separation, string_load_multiplier, string_circ, string_spanTree, string_root_relax, ".csv")
      end
      df = DataFrame(Instance = instance, Objective_value = obj_value[end], Run_time = run_time_li[end], Solution_status = solution_status[end])
      if params["run_obbt"]
         df[!,"Run_time"] = [obbt_time_li[end] + run_time_li[end]]
         df[!,"Optimization_time"] = [run_time_li[end]]
         df[!,"OBBT_time"] = [obbt_time_li[end]]
         df[!,"OBBT_function_time"] = [obbt_func_time_li[end]]
         df[!,"OBBT_iter"] = [obbt_iter_li[end]]
         df[!, "Precent_vm_reduce"] = [percent_vm_reduction]
         df[!, "Precent_td_reduce"] = [percent_td_reduction]
         if params["obbt_z"]
            df[!, "Free_lines"] = [obbt_freeLine_li[end]]
            df[!, "Zero_lines"] = [obbt_zeroLine_li[end]]
         end
      end
      if params["model"] == "ots" || params["model"] == "ots_relax"
         df[!, "Opt_gap"] = [opt_gap]
         df[!, "num_off_lines"] = [num_off_lines]
         df[!, "off_lines"] = [string(off_lines_li)]
         df[!, "num_off_lines_leaf"] = [num_off_lines_leaf]
         df[!, "off_lines_leaf"] = [string(off_lines_leaf_li)]
      end
      df[!, "inst_ind"] = [inst_ind]
      CSV.write(output_path, df; append = true, header = params["header"])
   end
end
