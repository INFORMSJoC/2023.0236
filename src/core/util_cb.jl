# Mostly using functions from util_sep.jl.

# Create the big-M for disjunctive Benders cuts for cycles in ACOTS.
function cycle_big_M(gamma_bar, x, cutvareq_lb, cutvareq_ub)
  M = 0.0
  for i in 1:length(gamma_bar)
    if typeof(x[i]) == VariableRef
      ub = cutvareq_ub[i]
      lb = cutvareq_lb[i]
      if gamma_bar[i] > 0
        M += gamma_bar[i] * max(0, ub)
      else
        M += gamma_bar[i] * min(0, lb)
      end
    end
  end
  return M
end

# Getting the RHS of subproblem equality constraints, i.e. the cs_bar and si_bar, from cb_data.
function cb_create_beq(params, cyc, cs_w, cb_data)
  b_eq = []
  # if params["cycle_relax"] == "mc" && length(cyc) == 3
  #   if cs_w == "cs"
  #     cs_bar = round_densearray(JuMP.value.(cs))
  #     si_bar = round_densearray(JuMP.value.(si))
  #     # McCormick constraints
  #     append!(b_eq, [cs_bar[(cyc[1], cyc[3])], si_bar[(cyc[1], cyc[3])], cs_bar[(cyc[1], cyc[2])], si_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], si_bar[(cyc[2], cyc[3])]])
  #   elseif cs_w == "w"
  #     append!(b_eq, [0 for i in 1:6])
  #   end
  # end
  # if params["cycle_relax"] == "mc" && length(cyc) == 4
  #   if cs_w == "cs"
  #     append!(b_eq, [0 for i in 1:6])
  #   end
  # end
  if params["cycle_relax"] == "epr" && length(cyc) == 3
    if cs_w == "cs"
      expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
      # cs_bar = round_densearray(callback_value.(Ref(cb_data), cs))
      # si_bar = round_densearray(callback_value.(Ref(cb_data), si))
      cs_bar = round_densearray(callback_value.(cb_data, m[:cs]))
      si_bar = round_densearray(callback_value.(cb_data, m[:si]))
      master_vars_bar = [cs_bar[(cyc[1], cyc[3])], si_bar[(cyc[1], cyc[3])], cs_bar[(cyc[1], cyc[2])], si_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], si_bar[(cyc[2], cyc[3])]]
      master_vars_bar2 = [cs_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], cs_bar[(cyc[1], cyc[3])], si_bar[(cyc[1], cyc[2])], si_bar[(cyc[2], cyc[3])], si_bar[(cyc[1], cyc[3])]]
      append!(b_eq, master_vars_bar)
      append!(b_eq, 1)
      append!(b_eq, master_vars_bar2)
      append!(b_eq, [0 for i in 1:length(expairs_3)])
    elseif cs_w == "w"
      expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
      # wr_bar = round_densearray(callback_value.(Ref(cb_data), wr))
      # wi_bar = round_densearray(callback_value.(Ref(cb_data), wi))
      # w_bar = round_densearray(callback_value.(Ref(cb_data), w))
      wr_bar = round_densearray(callback_value.(cb_data, m[:wr]))
      wi_bar = round_densearray(callback_value.(cb_data, m[:wi]))
      w_bar = round_densearray(callback_value.(cb_data, m[:w]))
      master_vars_bar2 = [wr_bar[(cyc[1], cyc[2])], wr_bar[(cyc[2], cyc[3])], wr_bar[(cyc[1], cyc[3])], wi_bar[(cyc[1], cyc[2])], wi_bar[(cyc[2], cyc[3])], wi_bar[(cyc[1], cyc[3])], w_bar[cyc[1]], w_bar[cyc[2]], w_bar[cyc[3]]]
      append!(b_eq, [0 for i in 1:6])
      append!(b_eq, 1)
      append!(b_eq, master_vars_bar2)
      append!(b_eq, [0 for i in 1:length(expairs_w_3)])
    end
  end
  if params["cycle_relax"] == "epr" && length(cyc) == 4
    if cs_w == "cs"
      expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
      # cs_bar = round_densearray(callback_value.(Ref(cb_data), cs))
      # si_bar = round_densearray(callback_value.(Ref(cb_data), si))
      cs_bar = round_densearray(callback_value.(cb_data, m[:cs]))
      si_bar = round_densearray(callback_value.(cb_data, m[:si]))
      master_vars_bar2 = [cs_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], cs_bar[(cyc[3], cyc[4])], cs_bar[(cyc[1], cyc[4])], si_bar[(cyc[1], cyc[2])], si_bar[(cyc[2], cyc[3])], si_bar[(cyc[3], cyc[4])], si_bar[(cyc[1], cyc[4])]]
      # println("cyc: ", cyc)
      # println("master_var2: ", master_vars_bar2)
      append!(b_eq, [0 for i in 1:6])
      append!(b_eq, 1)
      append!(b_eq, master_vars_bar2)
      append!(b_eq, [0 for i in 1:length(expairs_4)])
    elseif cs_w == "w"
      expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
      # wr_bar = round_densearray(callback_value.(Ref(cb_data), wr))
      # wi_bar = round_densearray(callback_value.(Ref(cb_data), wi))
      # w_bar = round_densearray(callback_value.(Ref(cb_data), w))
      wr_bar = round_densearray(callback_value.(cb_data, m[:wr]))
      wi_bar = round_densearray(callback_value.(cb_data, m[:wi]))
      w_bar = round_densearray(callback_value.(cb_data, m[:w]))
      num_cyc_con = 2
      master_vars_bar2 = [wr_bar[(cyc[1], cyc[2])], wr_bar[(cyc[2], cyc[3])], wr_bar[(cyc[3], cyc[4])], wr_bar[(cyc[1], cyc[4])], wi_bar[(cyc[1], cyc[2])], wi_bar[(cyc[2], cyc[3])], wi_bar[(cyc[3], cyc[4])], wi_bar[(cyc[1], cyc[4])]]
      if params["rotate_4w"]
        append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
        num_cyc_con = 6
        append!(master_vars_bar2, [w_bar[cyc[1]], w_bar[cyc[2]], w_bar[cyc[3]], w_bar[cyc[4]]])
      end
      append!(b_eq, [0 for i in 1:num_cyc_con])
      append!(b_eq, 1)
      append!(b_eq, master_vars_bar2)
      append!(b_eq, [0 for i in 1:length(expairs_w_4)])
    end
  end
  return b_eq
end

# Initialize the subproblems.
function cb_build_sp(cyc_size, m_sp, params, cs_w)
  x = create_x(params, m_sp, cyc_size, cs_w)
  A_eq = create_Aeq(params, cyc_size, x, cs_w) # conefficient matrix for equality constraints
  num_constr_eq = size(A_eq)[1]
  b_eq = zeros(num_constr_eq)
  @constraint(m_sp, gamma_eq[i = 1:num_constr_eq], sum(A_eq[i][j] * x[j] for j in 1:length(x)) == b_eq[i])

  A_ieq = nothing
  b_ieq = nothing
  num_constr_ieq = 0
  if params["cycle_relax"] == "mc"
    # A_ieq = create_Aieq(params, cyc_size, x, cs_w) # coefficient matrix for inequality constraints.
    # num_constr_ieq = size(A_ieq)[1]
    # b_ieq = zeros(num_constr_ieq)
    # @constraint(m_sp, gamma_ieq[i = 1:num_constr_ieq], sum(A_ieq[i][j] * x[j] for j in 1:length(x)) >= b_ieq[i])
  end
  @objective(m_sp, Min, 0)
  return x
end

# Generate all Benders feasibility cuts, one per cycle (if needed).
function cb_sp_per_cyc(cyc_size, m_sp, x_dict, params, cs_w, cb_data)
  x = x_dict[(cyc_size, cs_w)]
  # println("gamma: ", m_sp[:gamma_eq])
  num_constr_eq = size(m_sp[:gamma_eq])[1]
  for cyc in cycle[string("cyc_", string(cyc_size))]
    # println(callback_value.(Ref(cb_data), z))
    # println(cb_data)
    # println("cb_data.model: ", cb_data.model)
    # println(cb_data.tree)
    # println(Ref(cb_data))
    # println("z: ", z)
    # z_bar = callback_value.(cb_data, z)
    # println("Here!")
    # println("m[:z]: ", m[:z])
    # println("z 1 3: ", callback_value(cb_data, m[:z][(1,3)]))
    # println("z 1 3: ", callback_value(cb_data, z[(1,3)]))
    # z_bar = callback_value.(Ref(cb_data), m[:z])
    # z_bar = callback_value.(Ref(cb_data), z)
    buspairs = cyc_to_buspair(cyc)
    # z_bar = Dict()
    # for bp in buspairs
    #   # z_bar = callback_value.(cb_data, m[:z])
      
    #   z_bar[bp] = callback_value(cb_data, m[:z][bp])
    # end
    z_bar = callback_value.(cb_data, m[:z])
    num_used_lines = sum(z_bar[bp] for bp in buspairs)
    if (num_used_lines <= cyc_size - 0.5)
      continue
    end
    b_eq = cb_create_beq(params, cyc, cs_w, cb_data) # ****** cb_data used!
    # println("cs_w: ", cs_w)
    # println("b_eq: ", b_eq)
    for i in 1:num_constr_eq
      set_normalized_rhs(m_sp[:gamma_eq][i], b_eq[i])
    end
    if params["cycle_relax"] == "mc"
      # b_ieq = create_bieq(params, cyc, cs_w)
      # for i in 1:num_constr_ieq
      #   set_normalized_rhs(gamma_ieq[i], b_ieq[i])
      # end
    elseif params["cycle_relax"] == "epr"
      modify_Aeq_epr_polytope(cyc, m_sp[:gamma_eq], x, cs_w)
    end
    if params["print_model"] #&& cyc_size == 4
      Write_model(m_sp, "results/separation_model.txt")
    end
    specs["sp_time"] += @elapsed JuMP.optimize!(m_sp)
    sp_status = termination_status(m_sp)
    if sp_status != OPTIMAL # when subproblem is infeasible.
      # println("Type of problem: ", cs_w)
      # println("Termination status of sub: ", sp_status)
      # Write_model(m_sp, "results/separation_model.txt")
      # specs["terminate"] = false
      gamma_eq_bar = round_vector(JuMP.dual.(m_sp[:gamma_eq]))
      cut_var_eq = create_cutvareq(params, cyc, cs_w)
      gamma_ieq_bar = 0
      cut_var_ieq = 0
      if params["cycle_relax"] == "mc"
        # gamma_ieq_bar = JuMP.dual.(gamma_ieq)
        # cut_var_ieq = create_cutvarieq(params, cyc, cs_w)
      end

      if params["cycle_relax"] == "mc"
        # @constraint(m, sum(gamma_eq_bar[i] * cut_var_eq[i] for i in 1 : num_constr_eq) + sum(gamma_ieq_bar[i] * cut_var_ieq[i] for i in 1 : num_constr_ieq) <= - params["sep_tol"])
      elseif params["cycle_relax"] == "epr"
        cutvareq_lb, cutvareq_ub = cutvareq_bds(cyc, cs_w)
        M = cycle_big_M(gamma_eq_bar, cut_var_eq, cutvareq_lb, cutvareq_ub)
        con = []
        cyc_tp = Tuple(cyc)
        push!(con, @build_constraint(zc[cyc_tp] >= 1 - sum((1 - z[bp]) for bp in buspairs) - 1e-3))
        push!(con, @build_constraint(zc[cyc_tp] <= 1 / cyc_size * sum(z[bp] for bp in buspairs) + 1e-3))
        # push!(con, @build_constraint(zc[cyc_tp] == 1))
        ind_var = []
        ind_cnst = []
        for i in 1:num_constr_eq
          if typeof(cut_var_eq[i]) == Int64
            append!(ind_cnst, i)
          else
            append!(ind_var, i)
          end
        end
        push!(con, @build_constraint(sum(gamma_eq_bar[i] * cut_var_eq[i] for i in ind_var) <= - sum(gamma_eq_bar[i] * cut_var_eq[i] for i in ind_cnst) * zc[cyc_tp] + M * (1 - zc[cyc_tp]))) #- params["sep_tol"]
        # println("con2: ", @build_constraint(sum(gamma_eq_bar[i] * cut_var_eq[i] for i in 1:num_constr_eq) <= 0))
        # println("M: ", M)
        # println("con: ", con)
        for i in 1:length(con)
          MOI.submit(m, MOI.LazyConstraint(cb_data), con[i])
        end
      end
      specs["num_cuts"] += length(con)
      if cs_w == "cs"
        specs["num_cs_cuts"] += length(con)
      elseif cs_w == "w"
        specs["num_w_cuts"] += length(con)
      end
    end
  end
end

# # Generate all cycle cuts for cycle length = cyc_size.
# function cb_solve_sp(cyc_size, cb_data)
#   if params["cycle_c_s_cuts"]
#     m_sp = Get_sp_solver(params)
#     specs["sp_total_time"] += @elapsed cb_sp_per_cyc(cyc_size, m_sp, params, "cs", cb_data)
#   end
#   if params["cycle_wr_wi_cuts"]
#     m_sp = Get_sp_solver(params)
#     specs["sp_total_time"] += @elapsed cb_sp_per_cyc(cyc_size, m_sp, params, "w", cb_data)
#   end
# end
