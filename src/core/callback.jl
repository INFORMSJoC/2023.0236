# Adding cycle constraints through callback. For QCOTS.

# function acots_cb()
cb_calls = Cint[]
x_dict = Dict()
if params["cycle_max_bnd"] >= 3
    if params["cycle_c_s_cuts"]
        m_sp3_cs = Get_sp_solver(params)
        x_dict[(3, "cs")] = cb_build_sp(3, m_sp3_cs, params, "cs")
    end
    if params["cycle_wr_wi_cuts"]
        m_sp3_w = Get_sp_solver(params)
        x_dict[(3, "w")] = cb_build_sp(3, m_sp3_w, params, "w")
    end
end
if params["cycle_max_bnd"] >= 4
    if params["cycle_c_s_cuts"]
        m_sp4_cs = Get_sp_solver(params)
        x_dict[(4, "cs")] = cb_build_sp(4, m_sp4_cs, params, "cs")
    end
    if params["cycle_wr_wi_cuts"]
        m_sp4_w = Get_sp_solver(params)
        x_dict[(4, "w")] = cb_build_sp(4, m_sp4_w, params, "w")
    end
end

# New callback function using general callback.
#=
function my_callback_function(cb_data)
    # Define variables and constants
    init_time = time()
    specs["num_iterations"] += 1
    cb_where = callback_node_status(cb_data, m)  # New way to determine the callback node status

    # Log callback calls
    # push!(cb_calls, cb_where)

    # Decide when to run the callback
    if cb_where != MOI.CALLBACK_NODE_STATUS_INTEGER
        return
    end
    if specs["num_iterations"] > params["iter_limit"]
        specs["num_iterations"] -= 1
        return
    end
    if specs["num_cuts"] > params["cut_limit"]
        return
    end

    # Query callback attributes using the new MOI interface
    if cb_where == MOI.CALLBACK_NODE_STATUS_INTEGER
        if callback_node_status(cb_data, m) != GRB_OPTIMAL
            return  # Solution is something other than optimal.
        end
    end

    # Load primal solution before querying variable values
    # callback_load_variable_primal(cb_data)  # Updated function call

    if params["cycle_max_bnd"] >= 3
        if params["cycle_c_s_cuts"]
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(3, m_sp3_cs, x_dict, params, "cs", cb_data)
        end
        if params["cycle_wr_wi_cuts"]
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(3, m_sp3_w, x_dict, params, "w", cb_data)
        end
    end
    if params["cycle_max_bnd"] >= 4
        if params["cycle_c_s_cuts"]
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(4, m_sp4_cs, x_dict, params, "cs", cb_data)
        end
        if params["cycle_wr_wi_cuts"]
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(4, m_sp4_w, x_dict, params, "w", cb_util)
        end
    end
    specs["cb_func_time"] += time() - init_time
    return
end
=#

# Old solver-specific callback function.
function my_callback_function(cb_data, cb_where::Cint)
    init_time = time()
    specs["num_iterations"] += 1
    # You can reference variables outside the function as normal
    push!(cb_calls, cb_where)
    # You can select where the callback is run
    if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
        return
    end
    if specs["num_iterations"] > params["iter_limit"]
        specs["num_iterations"] -= 1
        return
    end
    if specs["num_cuts"] > params["cut_limit"]
        return
    end
    # You can query a callback attribute using GRBcbget
    if cb_where == GRB_CB_MIPNODE
        resultP = Ref{Cint}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
        if resultP[] != GRB_OPTIMAL
            return  # Solution is something other than optimal.
        end
    end
    # Before querying `callback_value`, you must call:
    Gurobi.load_callback_variable_primal(cb_data, cb_where)
    if params["cycle_max_bnd"] >= 3
        if params["cycle_c_s_cuts"]
            # m_sp3_cs = Get_sp_solver(params)
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(3, m_sp3_cs, x_dict, params, "cs", cb_data)
        end
        if params["cycle_wr_wi_cuts"]
            # m_sp3_w = Get_sp_solver(params)
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(3, m_sp3_w, x_dict, params, "w", cb_data)
        end
    end
    if params["cycle_max_bnd"] >= 4
        if params["cycle_c_s_cuts"]
            # m_sp4_cs = Get_sp_solver(params)
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(4, m_sp4_cs, x_dict, params, "cs", cb_data)
        end
        if params["cycle_wr_wi_cuts"]
            # m_sp4_w = Get_sp_solver(params)
            specs["sp_total_time"] += @elapsed cb_sp_per_cyc(4, m_sp4_w, x_dict, params, "w", cb_data)
        end
    end
    specs["cb_func_time"] += time() - init_time
    return
end

# You _must_ set this parameter if using lazy constraints.
# For old solver-specific callback
MOI.set(m, MOI.RawOptimizerAttribute("LazyConstraints"), 1)
MOI.set(m, Gurobi.CallbackFunction(), my_callback_function)

# For new callback
# set_attribute(m, MOI.LazyConstraintCallback(), my_callback_function)

specs["total_time"] = @elapsed optimize!(m)
# end
