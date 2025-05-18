###########################################################################
# All parameters for the optimization problem will be populated here
###########################################################################

# model options: ots, ots_relax, opf
params = Dict("model" => "ots", #opf: plain AC-OPF with QC relaxations, ots: optimal transission switching problem, ots_relax: relax binary variables.
            "relax" => args[4], #tri: lambda formulation (ext), rmc: recursive McCormick (pm) (utility.jl L487, fix lb of vv and cs to 0)
            "lnc" => true, #lifted non-linear cuts, keep these active always
            "rlt_cuts" => false,
            "run_obbt" => args[8],
            "obbt_termination" => :avg,
            "obbt_perVar" => false, # update the OBBT model every time a variable's bounds are updated.
            "warm_start" => args[5],
            "compare_cyc_constr" => 0, # change the td to compare mc and epr. pi/2
            "opf_exact_pm" => false, # solving the acopf problem locally directly with PM.jl. Obtains a local feasible solution. Used for upper bound.
            "ots_exact_pm" => false, # solving the exact acots formulation locally with PM.jl. Obtains a local feasible solution. Used for upper bound.
            "opf_qc_pm" => false,
            "ots_qc_pm" => false, # Solving the OTS QC via PowerModels.jl
            "ots_card_lb" => 0, # Number of lines must turn off in ots; default: 0.
            "ots_card_ub" => Inf, # UB for number of lines to turn on in ots; default: Inf.
            "grb_feas_tol" => args[3], # 1e-4. Useful for warm start. Default: 1e-6, min 1e-2, max 1e-9.
            "time_limit" => 3600 * args[7], # 3600 * 2
            "repeat_first_inst" => args[6], # repeat the first instance for more accurate computation time.
            "load_multiplier" => 1.0, # 1.16
            "ext_cons" => true, # Turn on/off extra constraints in Ext.
            "presolve_thread" => true, # false: turn off presolve, use single thread; true: use gurobi default setting.
            "obbt_vm" => true,
            "obbt_td" => true,
            "obbt_z" => true,
            "obbt_zc" => true,
            "valid_ieq" => false, # Valid inequality for sum(z) >=1 is load - capacity > 0.
            "preprocess_tightening" => false, # preprocssing to tighten the bounds of p/q and z.

            "minlp_solver" => "juniper", # knitro, juniper
            "span_tree" => false, # Find a spanning tree and keep it on.
            "cal_weight" => true, # Calculate the weight of the spanning tree using ACOPF solution from PM.jl exact.
            "partial_start" => false, # Provide partial warm start with all lines on.
            "get_root_relax" => false, # Get the root relaxation of the problem.
            "off_insights" => false, # To obtain insights on off lines for E vs ECB.

            # Cycles cuts
            "cycle_cuts" => args[9],
            "separation" => false,
            "callback" => args[10],
            "iter_limit" => Inf64, # Inf64, 20. For separation and callback. Limit on the number of iterations or callback nodes.
            "cut_limit" => 200, # Inf64, 20. For callback. Limit on the total number of cuts.
            "cycle_relax" => "epr", # epr: extreme point representation; mc: McCormick; none: bilinear
            "cycle_c_s_cuts" => true,
            "cycle_wr_wi_cuts" => true,
            "cycle_sw_cuts" => false, # adding Harsha's cycle constraints for s and w.
            "rotate_4w" => false,
            "rotate" => true, # whether include all rotations of the cycle constraints.
            "cycle_max_bnd" => 4,

            # Instance
            # "instance_type" => "pglib", #pglib, nesta
            "instance" => "pglib_opf_case14_ieee__sad",# pglib_opf_case3_lmbd__sad pglib_opf_case14_ieee__sad 
            "header" => true,
            "solver" => "gurobi", #"ipopt" "cplex" "gurobi"
            "separation_solver" => "gurobi", # "cplex"
            "print_model" => false,
            "solver_log" => true,
            "sp_solver_log" => false,
            "obbt_solver_log" => false,
            "formulate_problem" => false,
            "write_sol" => false, # Wrting OPF solution for warm start.

            # Parameters that are seldomly changed.
            "extra_rlt" => false, # Extra RLT cuts from Proving Global Optimality paper.
            "circ_cons" => false, # Constraint cos^2 + sin^2 = 1.
            "run_auto_obbt" => false, # running the built-in OBBT.
            "naive_bt" => false, # tightening the bounds with exact solutions.
            "turn_on_lines" => false, # turn on all lines in QCOTS, to compare with QCOPF
            "turn_off_lines" => false, # Turn off lines in QCOPF, to check the QCOTS solution.
            "plot_measures" => false, # printing out measures for plotting the relationship.
            "sep_tol" => 1e-4, # Numerical tolerance for cuts generated in separation problems.
            "debug" => false, # Debugging by fixing variable values.
            "optimize" => true,
            "LightGraphs" => true,
            "correct_bus_numbering" => true,
            "correct_branch_directions" => false
            )
