function Get_solver(params)
    if params["solver"] == "ipopt" && params["model"] != "ots"
        m =  Model(Ipopt.Optimizer)
        (params["solver_log"]) ? (print_level = 5) : (print_level = 0)
        set_optimizer_attributes(m, "print_level" => print_level, "sb" => "yes")
    end

    if params["solver"] == "cplex"
        m =  Model(CPLEX.Optimizer)
        (params["solver_log"]) ? (print_level = 1) : (print_level = 0)
        set_optimizer_attributes(m, "CPX_PARAM_SCRIND" => print_level)
    end

    if params["solver"] == "gurobi"
        m =  Model(() -> Gurobi.Optimizer(GRB_ENV))
        (params["solver_log"]) ? (print_level = 1) : (print_level = 0)
        set_optimizer_attributes(m, "OutputFlag" => print_level)
        set_optimizer_attributes(m, "TimeLimit" => params["time_limit"])
        if params["model"] == "ots" && params["warm_start"]
            set_optimizer_attributes(m, "FeasibilityTol" => params["grb_feas_tol"])
        end
        if params["cycle_cuts"] && params["cycle_relax"] == "none"
            set_optimizer_attributes(m, "NonConvex" => 2)
        end
        if params["presolve_thread"] == false
            set_optimizer_attributes(m, "Threads" => 1)
            set_optimizer_attributes(m, "Presolve" => 0)
        end
        if params["get_root_relax"]
            MOI.set(m, MOI.RawOptimizerAttribute("NodeLimit"), 0)
            # set_optimizer_attributes(m, " NodeLimit" => 0)
        end
    end

    return m
end

function Get_sp_solver(params)
    if params["separation_solver"] == "gurobi"
        m =  Model(() -> Gurobi.Optimizer(GRB_ENV))
        (params["sp_solver_log"]) ? (print_level = 1) : (print_level = 0)
        set_optimizer_attributes(m, "OutputFlag" => print_level)
        set_optimizer_attributes(m, "InfUnbdInfo" => 1)
    elseif params["separation_solver"] == "cplex"
        m =  Model(CPLEX.Optimizer)
        (params["sp_solver_log"]) ? (print_level = 1) : (print_level = 0)
        set_optimizer_attributes(m, "CPX_PARAM_SCRIND" => print_level)
    end

    return m
end
