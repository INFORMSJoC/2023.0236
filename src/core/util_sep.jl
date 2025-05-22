# Utility variables and functions for separation.
specs = Dict("sp_time" => 0.0,
      "sp_total_time" => 0.0,
      "mp_time" => 0.0,
      "total_time" => 0.0,
      "num_iterations" => 0,
      "num_cuts" => 0,
      "num_cs_cuts" => 0,
      "num_w_cuts" => 0,
      "terminate" => false,
      "num_vars_li" => [0, 0],  # The number of variables from cs and w.
      "num_cons_li" => [0, 0],
      "cb_func_time" => 0.0
      )
# Get bounds for master_vars2
function mastervar2_bds(cyc, cs_w)
  master_vars2_lb = []
  master_vars2_ub = []
  if length(cyc) == 3
    if cs_w == "cs"
      master_vars2_lb = [bd_data["cos_min"][(cyc[1], cyc[2])], bd_data["cos_min"][(cyc[2], cyc[3])], bd_data["cos_min"][(cyc[1], cyc[3])], bd_data["sin_min"][(cyc[1], cyc[2])], bd_data["sin_min"][(cyc[2], cyc[3])], bd_data["sin_min"][(cyc[1], cyc[3])]]
      master_vars2_ub = [bd_data["cos_max"][(cyc[1], cyc[2])], bd_data["cos_max"][(cyc[2], cyc[3])], bd_data["cos_max"][(cyc[1], cyc[3])], bd_data["sin_max"][(cyc[1], cyc[2])], bd_data["sin_max"][(cyc[2], cyc[3])], bd_data["sin_max"][(cyc[1], cyc[3])]]
    elseif cs_w == "w"
      master_vars2_lb = [bd_data["wr_min"][(cyc[1], cyc[2])], bd_data["wr_min"][(cyc[2], cyc[3])], bd_data["wr_min"][(cyc[1], cyc[3])], bd_data["wi_min"][(cyc[1], cyc[2])], bd_data["wi_min"][(cyc[2], cyc[3])], bd_data["wi_min"][(cyc[1], cyc[3])], bd_data["w_min"][cyc[1]], bd_data["w_min"][cyc[2]], bd_data["w_min"][cyc[3]]]
      master_vars2_ub = [bd_data["wr_max"][(cyc[1], cyc[2])], bd_data["wr_max"][(cyc[2], cyc[3])], bd_data["wr_max"][(cyc[1], cyc[3])], bd_data["wi_max"][(cyc[1], cyc[2])], bd_data["wi_max"][(cyc[2], cyc[3])], bd_data["wi_max"][(cyc[1], cyc[3])], bd_data["w_max"][cyc[1]], bd_data["w_max"][cyc[2]], bd_data["w_max"][cyc[3]]]
    end
  end
  if length(cyc) == 4
    if cs_w == "cs"
      master_vars2_lb = [bd_data["cos_min"][(cyc[1], cyc[2])], bd_data["cos_min"][(cyc[2], cyc[3])], bd_data["cos_min"][(cyc[3], cyc[4])], bd_data["cos_min"][(cyc[1], cyc[4])], bd_data["sin_min"][(cyc[1], cyc[2])], bd_data["sin_min"][(cyc[2], cyc[3])], bd_data["sin_min"][(cyc[3], cyc[4])], bd_data["sin_min"][(cyc[1], cyc[4])]]
      master_vars2_ub = [bd_data["cos_max"][(cyc[1], cyc[2])], bd_data["cos_max"][(cyc[2], cyc[3])], bd_data["cos_max"][(cyc[3], cyc[4])], bd_data["cos_max"][(cyc[1], cyc[4])], bd_data["sin_max"][(cyc[1], cyc[2])], bd_data["sin_max"][(cyc[2], cyc[3])], bd_data["sin_max"][(cyc[3], cyc[4])], bd_data["sin_max"][(cyc[1], cyc[4])]]
    elseif cs_w == "w"
      master_vars2_lb = [bd_data["wr_min"][(cyc[1], cyc[2])], bd_data["wr_min"][(cyc[2], cyc[3])], bd_data["wr_min"][(cyc[3], cyc[4])], bd_data["wr_min"][(cyc[1], cyc[4])], bd_data["wi_min"][(cyc[1], cyc[2])], bd_data["wi_min"][(cyc[2], cyc[3])], bd_data["wi_min"][(cyc[3], cyc[4])], bd_data["wi_min"][(cyc[1], cyc[4])]]
      master_vars2_ub = [bd_data["wr_max"][(cyc[1], cyc[2])], bd_data["wr_max"][(cyc[2], cyc[3])], bd_data["wr_max"][(cyc[3], cyc[4])], bd_data["wr_max"][(cyc[1], cyc[4])], bd_data["wi_max"][(cyc[1], cyc[2])], bd_data["wi_max"][(cyc[2], cyc[3])], bd_data["wi_max"][(cyc[3], cyc[4])], bd_data["wi_max"][(cyc[1], cyc[4])]]
      if params["rotate_4w"]
        append!(master_vars2_lb, [bd_data["w_min"][cyc[1]], bd_data["w_min"][cyc[2]], bd_data["w_min"][cyc[3]], bd_data["w_min"][cyc[4]]])
        append!(master_vars2_ub, [bd_data["w_max"][cyc[1]], bd_data["w_max"][cyc[2]], bd_data["w_max"][cyc[3]], bd_data["w_max"][cyc[4]]])
      end
    end
  end

  return master_vars2_lb, master_vars2_ub
end

# Get bounds for variables in cut_var_eq.
function cutvareq_bds(cyc, cs_w)
  cutvareq_lb = []
  cutvareq_ub = []
  master_vars2_lb, master_vars2_ub = mastervar2_bds(cyc, cs_w)
  if length(cyc) == 3
    if cs_w == "cs"
      expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
      master_vars_lb = [bd_data["cos_min"][(cyc[1], cyc[3])], bd_data["sin_min"][(cyc[1], cyc[3])], bd_data["cos_min"][(cyc[1], cyc[2])], bd_data["sin_min"][(cyc[1], cyc[2])], bd_data["cos_min"][(cyc[2], cyc[3])], bd_data["sin_min"][(cyc[2], cyc[3])]]
      master_vars_ub = [bd_data["cos_max"][(cyc[1], cyc[3])], bd_data["sin_max"][(cyc[1], cyc[3])], bd_data["cos_max"][(cyc[1], cyc[2])], bd_data["sin_max"][(cyc[1], cyc[2])], bd_data["cos_max"][(cyc[2], cyc[3])], bd_data["sin_max"][(cyc[2], cyc[3])]]
      append!(cutvareq_lb, master_vars_lb)
      append!(cutvareq_lb, [nothing])
      append!(cutvareq_lb, master_vars2_lb)
      append!(cutvareq_lb, [nothing for i in 1:length(expairs_3)])
      append!(cutvareq_ub, master_vars_ub)
      append!(cutvareq_ub, [nothing])
      append!(cutvareq_ub, master_vars2_ub)
      append!(cutvareq_ub, [nothing for i in 1:length(expairs_3)])
    elseif cs_w == "w"
      expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
      append!(cutvareq_lb, [nothing for i in 1:6])
      append!(cutvareq_lb, [nothing])
      append!(cutvareq_lb, master_vars2_lb)
      append!(cutvareq_lb, [nothing for i in 1:length(expairs_w_3)])
      append!(cutvareq_ub, [nothing for i in 1:6])
      append!(cutvareq_ub, [nothing])
      append!(cutvareq_ub, master_vars2_ub)
      append!(cutvareq_ub, [nothing for i in 1:length(expairs_w_3)])
    end
  end
  if length(cyc) == 4
    if cs_w == "cs"
      expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
      append!(cutvareq_lb, [nothing for i in 1:6])
      append!(cutvareq_lb, [nothing])
      append!(cutvareq_lb, master_vars2_lb)
      append!(cutvareq_lb, [nothing for i in 1:length(expairs_4)])
      append!(cutvareq_ub, [nothing for i in 1:6])
      append!(cutvareq_ub, [nothing])
      append!(cutvareq_ub, master_vars2_ub)
      append!(cutvareq_ub, [nothing for i in 1:length(expairs_4)])
    elseif cs_w == "w"
      expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
      num_cyc_con = 2
      if params["rotate_4w"]
        append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
        num_cyc_con = 6
      end
      append!(cutvareq_lb, [nothing for i in 1:num_cyc_con])
      append!(cutvareq_lb, [nothing])
      append!(cutvareq_lb, master_vars2_lb)
      append!(cutvareq_lb, [nothing for i in 1:length(expairs_w_4)])
      append!(cutvareq_ub, [nothing for i in 1:num_cyc_con])
      append!(cutvareq_ub, [nothing])
      append!(cutvareq_ub, master_vars2_ub)
      append!(cutvareq_ub, [nothing for i in 1:length(expairs_w_4)])
    end
  end
  return cutvareq_lb, cutvareq_ub
end

function create_x(params, m_sp, cyc_size, cs_w)
  x = []
  if params["cycle_relax"] == "mc" && cyc_size == 3
    if cs_w == "cs"
      arcpairs = []
      push!(arcpairs, ((1, 2), (2, 3)))
      push!(arcpairs, ((2, 3), (1, 3)))
      push!(arcpairs, ((1, 3), (1, 2)))
      JuMP.@variable(m_sp, hcc[ap in arcpairs]) # lifted variable for cs * cs, for McCormick
      JuMP.@variable(m_sp, hss[ap in arcpairs])
      JuMP.@variable(m_sp, hcs[ap in arcpairs])
      JuMP.@variable(m_sp, hsc[ap in arcpairs])
      for ap in arcpairs
        push!(x, hcc[ap])
        push!(x, hss[ap])
        push!(x, hcs[ap])
        push!(x, hsc[ap])
      end
      specs["num_vars_li"][1] = length(x)
    elseif cs_w == "w"
      arcpairs = []
      busarcpairs = []
      push!(arcpairs, ((1, 2), (2, 3)))
      push!(arcpairs, ((2, 3), (1, 3)))
      push!(arcpairs, ((1, 3), (1, 2)))
      push!(busarcpairs, (2, 1, 3))
      push!(busarcpairs, (3, 1, 2))
      push!(busarcpairs, (1, 2, 3))
      JuMP.@variables(m_sp, begin
        wwr[ba in busarcpairs]
        wwi[ba in busarcpairs]
        wrr[ap in arcpairs]
        wii[ap in arcpairs]
        wri[ap in arcpairs]
        wir[ap in arcpairs]
      end)
      for ap in arcpairs
        push!(x, wrr[ap])
        push!(x, wii[ap])
        push!(x, wri[ap])
        push!(x, wir[ap])
      end
      for ba in busarcpairs
        push!(x, wwr[ba])
        push!(x, wwi[ba])
      end
      specs["num_vars_li"][2] = length(x)
    end
  end
  if params["cycle_relax"] == "mc" && cyc_size == 4
    if cs_w == "cs"
      arcpairs = []
      push!(arcpairs, ((1, 2), (3, 4)))
      push!(arcpairs, ((1, 4), (2, 3)))
      push!(arcpairs, ((1, 2), (2, 3)))
      push!(arcpairs, ((1, 4), (3, 4)))
      push!(arcpairs, ((2, 3), (3, 4)))
      push!(arcpairs, ((1, 4), (1, 2)))
      JuMP.@variables(m_sp, begin
        hcc[ap in arcpairs]
        hss[ap in arcpairs]
        hcs[ap in arcpairs]
        hsc[ap in arcpairs]
      end)
      for i in 1:(length(arcpairs)/2)
        ind = Int64(2 * i - 1)
        push!(x, hcc[arcpairs[ind]])
        push!(x, hss[arcpairs[ind]])
        ind = Int64(2 * i)
        push!(x, hcc[arcpairs[ind]])
        push!(x, hss[arcpairs[ind]])
        ind = Int64(2 * i - 1)
        push!(x, hcs[arcpairs[ind]])
        push!(x, hsc[arcpairs[ind]])
        ind = Int64(2 * i)
        push!(x, hcs[arcpairs[ind]])
        push!(x, hsc[arcpairs[ind]])
      end
    end
    specs["num_vars_li"][1] = length(x)
  end
  if params["cycle_relax"] == "epr" && cyc_size == 3
    if cs_w == "cs"
      expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
      JuMP.@variable(m_sp, λ_c[1:(2^6)] >= 0)
      JuMP.@variable(m_sp, x_c[ep in expairs_3])
      for ep in expairs_3
        push!(x, x_c[ep])
      end
      for i in 1:2^6
        push!(x, λ_c[i])
      end
      specs["num_vars_li"][1] = length(x)
    elseif cs_w == "w"
      expairs_w_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6), (8,3), (8,6), (9,1), (9,4), (7,2), (7,5)]
      JuMP.@variable(m_sp, λ_w[1:(2^9)] >= 0)
      JuMP.@variable(m_sp, x_w[ep in expairs_w_3])
      for ep in expairs_w_3
        push!(x, x_w[ep])
      end
      for i in 1:(2^9)
        push!(x, λ_w[i])
      end
      specs["num_vars_li"][2] = length(x)
    end
  end
  if params["cycle_relax"] == "epr" && cyc_size == 4
    if cs_w == "cs"
      expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
      JuMP.@variable(m_sp, λ_c[1:(2^8)] >= 0)
      JuMP.@variable(m_sp, x_c[ep in expairs_4])
      for ep in expairs_4
        push!(x, x_c[ep])
      end
      for i in 1:(2^8)
        push!(x, λ_c[i])
      end
      specs["num_vars_li"][1] = length(x)
    elseif cs_w == "w"
      expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
      num_master_var = 8
      if params["rotate_4w"]
          append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
          num_master_var = 12
      end
      JuMP.@variable(m_sp, λ_w[1:(2^num_master_var)] >= 0)
      JuMP.@variable(m_sp, x_w[ep in expairs_w_4])
      for ep in expairs_w_4
        push!(x, x_w[ep])
      end
      for i in 1:(2^num_master_var)
        push!(x, λ_w[i])
      end
    end
    specs["num_vars_li"][2] = length(x)
  end
  return x
end

function create_Aeq(params, cyc_size, x, cs_w)
  num_col = length(x)
  A_eq = []
  if params["cycle_relax"] == "mc" && cyc_size == 3
    if cs_w == "cs"
      # cycle constraints.
      temp = [1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1]
      for i in 1:6
        a = zeros(num_col)
        a[2 * i - 1] = temp[2 * i - 1]
        a[2 * i] = temp[2 * i]
        push!(A_eq, a)
      end
    elseif cs_w == "w"
      # start_x = specs["num_vars_li"][1]
      start_x = 0
      extra_x = 3 * 4
      temp = [1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1]
      for i in 1:6
        a = zeros(num_col)
        a[start_x + 2 * i - 1] = temp[2 * i - 1]
        a[start_x + 2 * i] = temp[2 * i]
        a[start_x + extra_x + i] = -1
        push!(A_eq, a)
      end
    end
  end
  if params["cycle_relax"] == "mc" && cyc_size == 4
    if cs_w == "cs"
      temp = [1, -1, -1, -1, 1, 1, 1, -1]
      temp2 = [temp[i] for j in 1:3 for i in 1:length(temp)]
      for i in 1:6
        a = zeros(num_col)
        for j in 1:4
          a[4 * (i - 1) + j] = temp2[4 * (i - 1) + j]
        end
        push!(A_eq, a)
      end
    end
  end
  if params["cycle_relax"] == "epr" && cyc_size == 3
    if cs_w == "cs"
      temp = [1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1]
      expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
      for i in 1:6
        a = zeros(num_col)
        a[2 * i - 1] = temp[2 * i - 1]
        a[2 * i] = temp[2 * i]
        push!(A_eq, a)
      end
      start_λ = length(expairs_3)
      a = zeros(num_col)
      for i in 1:64
        a[start_λ + i] = 1
      end
      push!(A_eq, a)
      dim = 6
      for j in 1:dim
        a = zeros(num_col)
        push!(A_eq, a) # Those coefficients will be set later by modify_Aeq_epr_polytope.
      end
      for ep in expairs_3
        a = zeros(num_col)
        push!(A_eq, a)
      end
      specs["num_cons_li"][1] = length(A_eq)
    elseif cs_w == "w"
      temp = [1, -1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1]
      expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
      # start_x = specs["num_vars_li"][1]
      start_x = 0
      for i in 1:6
        a = zeros(num_col)
        a[start_x + 2 * i - 1] = temp[2 * i - 1]
        a[start_x + 2 * i] = temp[2 * i]
        a[start_x + 12 + i] = -1
        push!(A_eq, a)
      end
      start_λ = length(expairs_w_3)
      a = zeros(num_col)
      for i in 1:(2^9)
        a[start_x + start_λ + i] = 1
      end
      push!(A_eq, a)
      dim = 9
      for j in 1:dim
        a = zeros(num_col)
        push!(A_eq, a)
      end
      for ep in expairs_w_3
        a = zeros(num_col)
        push!(A_eq, a)
      end
      specs["num_cons_li"][2] = length(A_eq)
    end
  end
  if params["cycle_relax"] == "epr" && cyc_size == 4
    if cs_w == "cs"
      temp1 = [1, -1, -1, -1, 1, 1, 1, -1]
      temp = []
      for i in 1:3
        append!(temp, temp1)
      end
      expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
      for i in 1:6
        a = zeros(num_col)
        a[4 * i - 3] = temp[4 * i - 3]
        a[4 * i - 2] = temp[4 * i - 2]
        a[4 * i - 1] = temp[4 * i - 1]
        a[4 * i] = temp[4 * i]
        push!(A_eq, a)
      end
      start_λ = length(expairs_4)
      a = zeros(num_col)
      for i in 1:(2^8)
        a[start_λ + i] = 1
      end
      push!(A_eq, a)
      dim = 8
      for j in 1:dim
        a = zeros(num_col)
        push!(A_eq, a)
      end
      for ep in expairs_4
        a = zeros(num_col)
        push!(A_eq, a)
      end
      specs["num_cons_li"][1] = length(A_eq)
    elseif cs_w == "w"
      temp1 = [1, -1, -1, -1, 1, 1, 1, -1]
      temp = [1, -1, -1, -1, 1, 1, 1, -1]
      expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
      num_master_var = 8
      num_cyc_con = 2
      if params["rotate_4w"]
        for i in 1:2
          append!(temp, temp1)
        end
        append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
        num_master_var = 12
        num_cyc_con = 6
      end
      # start_x = specs["num_vars_li"][1]
      start_x = 0
      for i in 1:num_cyc_con
        a = zeros(num_col)
        a[start_x + 4 * i - 3] = temp[4 * i - 3]
        a[start_x + 4 * i - 2] = temp[4 * i - 2]
        a[start_x + 4 * i - 1] = temp[4 * i - 1]
        a[start_x + 4 * i] = temp[4 * i]
        push!(A_eq, a)
      end
      start_λ = length(expairs_w_4)
      a = zeros(num_col)
      for i in 1:(2^num_master_var)
        a[start_x + start_λ + i] = 1
      end
      push!(A_eq, a)
      for j in 1:num_master_var
        a = zeros(num_col)
        push!(A_eq, a)
      end
      for ep in expairs_w_4
        a = zeros(num_col)
        push!(A_eq, a)
      end
      specs["num_cons_li"][2] = length(A_eq)
    end
  end
  return A_eq
end

function create_Aieq(params, cyc_size, x, cs_w)
  num_col = length(x)
  A_ieq = []
  if params["cycle_relax"] == "mc" && cyc_size == 3
    if cs_w == "cs"
      # McCormick constraints
      a = zeros(num_col)
      temp = [1, 1, -1, -1]
      for i in 1:12
        for j in 1:4
          a = zeros(num_col)
          a[i] = temp[j]
          push!(A_ieq, a)
        end
      end
    elseif cs_w == "w"
      # start_x = specs["num_vars_li"][1]
      start_x = 0
      extra_x = 3 * 4
      temp = [1, 1, -1, -1]
      for i in 1:(12 + 6)
        for j in 1:4
          a = zeros(num_col)
          a[start_x + i] = temp[j]
          push!(A_ieq, a)
        end
      end
    end
  end
  if params["cycle_relax"] == "mc" && cyc_size == 4
    if cs_w == "cs"
      a = zeros(num_col)
      temp = [1, 1, -1, -1]
      for i in 1:(4 * 2 * 3)
        for j in 1:4
          a = zeros(num_col)
          a[i] = temp[j]
          push!(A_ieq, a)
        end
      end
    end
  end
  return A_ieq
end

function create_beq(params, cyc, cs_w)
  b_eq = []
  if params["cycle_relax"] == "mc" && length(cyc) == 3
    if cs_w == "cs"
      cs_bar = round_densearray(JuMP.value.(cs))
      si_bar = round_densearray(JuMP.value.(si))
      # McCormick constraints
      append!(b_eq, [cs_bar[(cyc[1], cyc[3])], si_bar[(cyc[1], cyc[3])], cs_bar[(cyc[1], cyc[2])], si_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], si_bar[(cyc[2], cyc[3])]])
    elseif cs_w == "w"
      append!(b_eq, [0 for i in 1:6])
    end
  end
  if params["cycle_relax"] == "mc" && length(cyc) == 4
    if cs_w == "cs"
      append!(b_eq, [0 for i in 1:6])
    end
  end
  if params["cycle_relax"] == "epr" && length(cyc) == 3
    if cs_w == "cs"
      expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
      cs_bar = round_densearray(JuMP.value.(cs))
      si_bar = round_densearray(JuMP.value.(si))
      master_vars_bar = [cs_bar[(cyc[1], cyc[3])], si_bar[(cyc[1], cyc[3])], cs_bar[(cyc[1], cyc[2])], si_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], si_bar[(cyc[2], cyc[3])]]
      master_vars_bar2 = [cs_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], cs_bar[(cyc[1], cyc[3])], si_bar[(cyc[1], cyc[2])], si_bar[(cyc[2], cyc[3])], si_bar[(cyc[1], cyc[3])]]
      append!(b_eq, master_vars_bar)
      append!(b_eq, 1)
      append!(b_eq, master_vars_bar2)
      append!(b_eq, [0 for i in 1:length(expairs_3)])
    elseif cs_w == "w"
      expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
      wr_bar = round_densearray(JuMP.value.(wr))
      wi_bar = round_densearray(JuMP.value.(wi))
      w_bar = round_densearray(JuMP.value.(w))
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
      cs_bar = round_densearray(JuMP.value.(cs))
      si_bar = round_densearray(JuMP.value.(si))
      master_vars_bar2 = [cs_bar[(cyc[1], cyc[2])], cs_bar[(cyc[2], cyc[3])], cs_bar[(cyc[3], cyc[4])], cs_bar[(cyc[1], cyc[4])], si_bar[(cyc[1], cyc[2])], si_bar[(cyc[2], cyc[3])], si_bar[(cyc[3], cyc[4])], si_bar[(cyc[1], cyc[4])]]
      # println("cyc: ", cyc)
      # println("master_var2: ", master_vars_bar2)
      append!(b_eq, [0 for i in 1:6])
      append!(b_eq, 1)
      append!(b_eq, master_vars_bar2)
      append!(b_eq, [0 for i in 1:length(expairs_4)])
    elseif cs_w == "w"
      expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
      wr_bar = round_densearray(JuMP.value.(wr))
      wi_bar = round_densearray(JuMP.value.(wi))
      w_bar = round_densearray(JuMP.value.(w))
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

function create_bieq(params, cyc, cs_w)
  b_ieq = []
  if params["cycle_relax"] == "mc" && length(cyc) == 3
    if cs_w == "cs"
      cs_bar = round_densearray(JuMP.value.(cs))
      si_bar = round_densearray(JuMP.value.(si))
      arcpairs = []
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[2], cyc[3]), (cyc[1], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[3]), (cyc[1], cyc[2])))
      for ap in arcpairs
        arc1 = ap[1]
        arc2 = ap[2]
        temp = [(cs[arc1], cs[arc2]), (si[arc1], si[arc2]), (cs[arc1], si[arc2]), (si[arc1], cs[arc2])]
        temp2 = [(cs_bar[arc1], cs_bar[arc2]), (si_bar[arc1], si_bar[arc2]), (cs_bar[arc1], si_bar[arc2]), (si_bar[arc1], cs_bar[arc2])]
        for j in 1:4
          y1 = temp[j][1]
          y2 = temp[j][2]
          lb = [lower_bound(y1), lower_bound(y2)]
          ub = [upper_bound(y1), upper_bound(y2)]
          y1_bar = temp2[j][1]
          y2_bar = temp2[j][2]
          push!(b_ieq, lb[1]*y2_bar + lb[2]*y1_bar - lb[1]*lb[2])
          push!(b_ieq, ub[1]*y2_bar + ub[2]*y1_bar - ub[1]*ub[2])
          push!(b_ieq, - lb[1]*y2_bar - ub[2]*y1_bar + lb[1]*ub[2])
          push!(b_ieq, - ub[1]*y2_bar - lb[2]*y1_bar + ub[1]*lb[2])
        end
      end
    elseif cs_w == "w"
      wr_bar = round_densearray(JuMP.value.(wr))
      wi_bar = round_densearray(JuMP.value.(wi))
      w_bar = round_densearray(JuMP.value.(w))
      arcpairs = []
      busarcpairs = []
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[2], cyc[3]), (cyc[1], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[3]), (cyc[1], cyc[2])))
      push!(busarcpairs, (cyc[2], cyc[1], cyc[3]))
      push!(busarcpairs, (cyc[3], cyc[1], cyc[2]))
      push!(busarcpairs, (cyc[1], cyc[2], cyc[3]))
      for ap in arcpairs
        arc1 = ap[1]
        arc2 = ap[2]
        temp = [(wr[arc1], wr[arc2]), (wi[arc1], wi[arc2]), (wr[arc1], wi[arc2]), (wi[arc1], wr[arc2])]
        temp2 = [(wr_bar[arc1], wr_bar[arc2]), (wi_bar[arc1], wi_bar[arc2]), (wr_bar[arc1], wi_bar[arc2]), (wi_bar[arc1], wr_bar[arc2])]
        for j in 1:length(temp)
          populate_bieq(b_ieq, temp[j], temp2[j])
        end
      end
      for ba in busarcpairs
        bus = ba[1]
        arc = (ba[2], ba[3])
        temp = [(w[bus], wr[arc]), (w[bus], wi[arc])]
        temp2 = [(w_bar[bus], wr_bar[arc]), (w_bar[bus], wi_bar[arc])]
        for j in 1:length(temp)
          populate_bieq(b_ieq, temp[j], temp2[j])
        end
      end
    end
  end
  if params["cycle_relax"] == "mc" && length(cyc) == 4
    if cs_w == "cs"
      cs_bar = round_densearray(JuMP.value.(cs))
      si_bar = round_densearray(JuMP.value.(si))
      arcpairs = []
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[3], cyc[4])))
      push!(arcpairs, ((cyc[1], cyc[4]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[4]), (cyc[3], cyc[4])))
      push!(arcpairs, ((cyc[2], cyc[3]), (cyc[3], cyc[4])))
      push!(arcpairs, ((cyc[1], cyc[4]), (cyc[1], cyc[2])))
      for i in 1:(length(arcpairs)/2)
        ind = Int64(2 * i - 1)
        arc1 = arcpairs[ind][1]
        arc2 = arcpairs[ind][2]
        ind = Int64(2 * i)
        arc3 = arcpairs[ind][1]
        arc4 = arcpairs[ind][2]
        temp = [(cs[arc1], cs[arc2]), (si[arc1], si[arc2]), (cs[arc3], cs[arc4]), (si[arc3], si[arc4]), (cs[arc1], si[arc2]), (si[arc1], cs[arc2]), (cs[arc3], si[arc4]), (si[arc3], cs[arc4])]
        temp2 = [(cs_bar[arc1], cs_bar[arc2]), (si_bar[arc1], si_bar[arc2]), (cs_bar[arc3], cs_bar[arc4]), (si_bar[arc3], si_bar[arc4]), (cs_bar[arc1], si_bar[arc2]), (si_bar[arc1], cs_bar[arc2]), (cs_bar[arc3], si_bar[arc4]), (si_bar[arc3], cs_bar[arc4])]
        for j in 1:length(temp)
          populate_bieq(b_ieq, temp[j], temp2[j])
        end
      end
    end
  end
  return b_ieq
end

function create_cutvareq(params, cyc, cs_w)
  cut_var_eq = []
  if params["cycle_relax"] == "mc" && length(cyc) == 3
    if cs_w == "cs"
      # McCormick constraints
      append!(cut_var_eq, [cs[(cyc[1], cyc[3])], si[(cyc[1], cyc[3])], cs[(cyc[1], cyc[2])], si[(cyc[1], cyc[2])], cs[(cyc[2], cyc[3])], si[(cyc[2], cyc[3])]])
    elseif cs_w == "w"
      append!(cut_var_eq, [0 for i in 1:6])
    end
  end
  if params["cycle_relax"] == "mc" && length(cyc) == 4
    if cs_w == "cs"
      append!(cut_var_eq, [0 for i in 1:6])
    end
  end
  if params["cycle_relax"] == "epr" && length(cyc) == 3
    if cs_w == "cs"
      expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
      master_vars = [cs[(cyc[1], cyc[3])], si[(cyc[1], cyc[3])], cs[(cyc[1], cyc[2])], si[(cyc[1], cyc[2])], cs[(cyc[2], cyc[3])], si[(cyc[2], cyc[3])]]
      master_vars2 = [cs[(cyc[1], cyc[2])], cs[(cyc[2], cyc[3])], cs[(cyc[1], cyc[3])], si[(cyc[1], cyc[2])], si[(cyc[2], cyc[3])], si[(cyc[1], cyc[3])]]
      append!(cut_var_eq, master_vars)
      append!(cut_var_eq, 1)
      append!(cut_var_eq, master_vars2)
      append!(cut_var_eq, [0 for i in 1:length(expairs_3)])
    elseif cs_w == "w"
      expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
      master_vars2 = [wr[(cyc[1], cyc[2])], wr[(cyc[2], cyc[3])], wr[(cyc[1], cyc[3])], wi[(cyc[1], cyc[2])], wi[(cyc[2], cyc[3])], wi[(cyc[1], cyc[3])], w[cyc[1]], w[cyc[2]], w[cyc[3]]]
      append!(cut_var_eq, [0 for i in 1:6])
      append!(cut_var_eq, 1)
      append!(cut_var_eq, master_vars2)
      append!(cut_var_eq, [0 for i in 1:length(expairs_w_3)])
    end
  end
  if params["cycle_relax"] == "epr" && length(cyc) == 4
    if cs_w == "cs"
      expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
      master_vars2 = [cs[(cyc[1], cyc[2])], cs[(cyc[2], cyc[3])], cs[(cyc[3], cyc[4])], cs[(cyc[1], cyc[4])], si[(cyc[1], cyc[2])], si[(cyc[2], cyc[3])], si[(cyc[3], cyc[4])], si[(cyc[1], cyc[4])]]
      append!(cut_var_eq, [0 for i in 1:6])
      append!(cut_var_eq, 1)
      append!(cut_var_eq, master_vars2)
      append!(cut_var_eq, [0 for i in 1:length(expairs_4)])
    elseif cs_w == "w"
      expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
      num_cyc_con = 2
      master_vars2 = [wr[(cyc[1], cyc[2])], wr[(cyc[2], cyc[3])], wr[(cyc[3], cyc[4])], wr[(cyc[1], cyc[4])], wi[(cyc[1], cyc[2])], wi[(cyc[2], cyc[3])], wi[(cyc[3], cyc[4])], wi[(cyc[1], cyc[4])]]
      if params["rotate_4w"]
        append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
        num_cyc_con = 6
        append!(master_vars2, [w[cyc[1]], w[cyc[2]], w[cyc[3]], w[cyc[4]]])
      end
      append!(cut_var_eq, [0 for i in 1:num_cyc_con])
      append!(cut_var_eq, 1)
      append!(cut_var_eq, master_vars2)
      append!(cut_var_eq, [0 for i in 1:length(expairs_w_4)])
    end
  end
  return cut_var_eq
end

function create_cutvarieq(params, cyc, cs_w)
  cut_var_ieq = []
  if params["cycle_relax"] == "mc" && length(cyc) == 3
    if cs_w == "cs"
      arcpairs = []
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[2], cyc[3]), (cyc[1], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[3]), (cyc[1], cyc[2])))
      for ap in arcpairs
        arc1 = ap[1]
        arc2 = ap[2]
        temp = [(cs[arc1], cs[arc2]), (si[arc1], si[arc2]), (cs[arc1], si[arc2]), (si[arc1], cs[arc2])]
        for j in 1:4
          y1 = temp[j][1]
          y2 = temp[j][2]
          lb = [lower_bound(y1), lower_bound(y2)]
          ub = [upper_bound(y1), upper_bound(y2)]
          push!(cut_var_ieq, lb[1]*y2 + lb[2]*y1 - lb[1]*lb[2])
          push!(cut_var_ieq, ub[1]*y2 + ub[2]*y1 - ub[1]*ub[2])
          push!(cut_var_ieq, - lb[1]*y2 - ub[2]*y1 + lb[1]*ub[2])
          push!(cut_var_ieq, - ub[1]*y2 - lb[2]*y1 + ub[1]*lb[2])
        end
      end
    elseif cs_w == "w"
      arcpairs = []
      busarcpairs = []
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[2], cyc[3]), (cyc[1], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[3]), (cyc[1], cyc[2])))
      push!(busarcpairs, (cyc[2], cyc[1], cyc[3]))
      push!(busarcpairs, (cyc[3], cyc[1], cyc[2]))
      push!(busarcpairs, (cyc[1], cyc[2], cyc[3]))
      for ap in arcpairs
        arc1 = ap[1]
        arc2 = ap[2]
        temp = [(wr[arc1], wr[arc2]), (wi[arc1], wi[arc2]), (wr[arc1], wi[arc2]), (wi[arc1], wr[arc2])]
        for j in 1:length(temp)
          populate_cutvarieq(cut_var_ieq, temp[j])
        end
      end
      for ba in busarcpairs
        bus = ba[1]
        arc = (ba[2], ba[3])
        temp = [(w[bus], wr[arc]), (w[bus], wi[arc])]
        for j in 1:length(temp)
          populate_cutvarieq(cut_var_ieq, temp[j])
        end
      end
    end
  end
  if params["cycle_relax"] == "mc" && length(cyc) == 4
    if cs_w == "cs"
      cs_bar = round_densearray(JuMP.value.(cs))
      si_bar = round_densearray(JuMP.value.(si))
      arcpairs = []
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[3], cyc[4])))
      push!(arcpairs, ((cyc[1], cyc[4]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
      push!(arcpairs, ((cyc[1], cyc[4]), (cyc[3], cyc[4])))
      push!(arcpairs, ((cyc[2], cyc[3]), (cyc[3], cyc[4])))
      push!(arcpairs, ((cyc[1], cyc[4]), (cyc[1], cyc[2])))
      for i in 1:(length(arcpairs)/2)
        ind = Int64(2 * i - 1)
        arc1 = arcpairs[ind][1]
        arc2 = arcpairs[ind][2]
        ind = Int64(2 * i)
        arc3 = arcpairs[ind][1]
        arc4 = arcpairs[ind][2]
        temp = [(cs[arc1], cs[arc2]), (si[arc1], si[arc2]), (cs[arc3], cs[arc4]), (si[arc3], si[arc4]), (cs[arc1], si[arc2]), (si[arc1], cs[arc2]), (cs[arc3], si[arc4]), (si[arc3], cs[arc4])]
        for j in 1:length(temp)
          populate_cutvarieq(cut_var_ieq, temp[j])
        end
      end
    end
  end
  return cut_var_ieq
end

function populate_bieq(b_ieq, temp_j, temp2_j)
  y1 = temp_j[1]
  y2 = temp_j[2]
  lb = [lower_bound(y1), lower_bound(y2)]
  ub = [upper_bound(y1), upper_bound(y2)]
  y1_bar = temp2_j[1]
  y2_bar = temp2_j[2]
  push!(b_ieq, lb[1]*y2_bar + lb[2]*y1_bar - lb[1]*lb[2])
  push!(b_ieq, ub[1]*y2_bar + ub[2]*y1_bar - ub[1]*ub[2])
  push!(b_ieq, - lb[1]*y2_bar - ub[2]*y1_bar + lb[1]*ub[2])
  push!(b_ieq, - ub[1]*y2_bar - lb[2]*y1_bar + ub[1]*lb[2])
end

function populate_cutvarieq(cut_var_ieq, temp_j)
  y1 = temp_j[1]
  y2 = temp_j[2]
  lb = [lower_bound(y1), lower_bound(y2)]
  ub = [upper_bound(y1), upper_bound(y2)]
  push!(cut_var_ieq, lb[1]*y2 + lb[2]*y1 - lb[1]*lb[2])
  push!(cut_var_ieq, ub[1]*y2 + ub[2]*y1 - ub[1]*ub[2])
  push!(cut_var_ieq, - lb[1]*y2 - ub[2]*y1 + lb[1]*ub[2])
  push!(cut_var_ieq, - ub[1]*y2 - lb[2]*y1 + ub[1]*lb[2])
end

# In subproblem modifying the epr A_eq matrix for the constraints defining the polytope.
function modify_Aeq_epr_polytope(cyc, gamma_eq, x, cs_w)
  num_col = length(x)
  if cs_w == "cs" && length(cyc) == 3
    expairs_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6)]
    master_vars2 = [cs[(cyc[1], cyc[2])], cs[(cyc[2], cyc[3])], cs[(cyc[1], cyc[3])], si[(cyc[1], cyc[2])], si[(cyc[2], cyc[3])], si[(cyc[1], cyc[3])]] # The list of master problem variabels in the order of x^c_1 to x^c_6.
    modify_Aeq_util(cyc, gamma_eq, x, expairs_3, master_vars2, cs_w)
  elseif cs_w == "w" && length(cyc) == 3
    expairs_w_3 = [(1,2), (4,5), (1,5), (2,4), (2,3), (5,6), (2,6), (3,5), (1,3), (4,6), (3,4), (1,6), (8,3), (8,6), (9,1), (9,4), (7,2), (7,5)]
    master_vars2 = [wr[(cyc[1], cyc[2])], wr[(cyc[2], cyc[3])], wr[(cyc[1], cyc[3])], wi[(cyc[1], cyc[2])], wi[(cyc[2], cyc[3])], wi[(cyc[1], cyc[3])], w[cyc[1]], w[cyc[2]], w[cyc[3]]]
    modify_Aeq_util(cyc, gamma_eq, x, expairs_w_3, master_vars2, cs_w)
  elseif cs_w == "cs" && length(cyc) == 4
    expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
    master_vars2 = [cs[(cyc[1], cyc[2])], cs[(cyc[2], cyc[3])], cs[(cyc[3], cyc[4])], cs[(cyc[1], cyc[4])], si[(cyc[1], cyc[2])], si[(cyc[2], cyc[3])], si[(cyc[3], cyc[4])], si[(cyc[1], cyc[4])]]
    modify_Aeq_util(cyc, gamma_eq, x, expairs_4, master_vars2, cs_w)
  elseif cs_w == "w" && length(cyc) == 4
    expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
    master_vars2 = [wr[(cyc[1], cyc[2])], wr[(cyc[2], cyc[3])], wr[(cyc[3], cyc[4])], wr[(cyc[1], cyc[4])], wi[(cyc[1], cyc[2])], wi[(cyc[2], cyc[3])], wi[(cyc[3], cyc[4])], wi[(cyc[1], cyc[4])]]
    if params["rotate_4w"]
      append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
      append!(master_vars2, [w[cyc[1]], w[cyc[2]], w[cyc[3]], w[cyc[4]]])
    end
    modify_Aeq_util(cyc, gamma_eq, x, expairs_w_4, master_vars2, cs_w)
  end
end

# Used when modifying the A_eq in epr.
function modify_Aeq_util(cyc, gamma_eq, x, expairs, master_vars2, cs_w)
  start_λ = length(expairs)
  start_con = 6 + 1
  if cs_w == "w" && length(cyc) == 4 && !params["rotate_4w"]
    start_con = 2 + 1
  end
  # if cs_w == "w"
  #   start_con += specs["num_cons_li"][1]
  # end
  start_x = 0
  # if cs_w == "w"
  #   start_x += specs["num_vars_li"][1]
  # end
  lb, ub = mastervar2_bds(cyc, cs_w)
  # lb = [lower_bound(var) for var in master_vars2]
  # ub = [upper_bound(var) for var in master_vars2]
  X = cartesian_product(lb,ub)
  dim = length(master_vars2)
  for j in 1:dim
    cnt_1 = 0
    for i in 1:(2^dim)
      if X[i][j] == lb[j]
        set_normalized_coefficient(gamma_eq[start_con + j], x[start_x + start_λ + i], lb[j])
        cnt_1 += 1
      else
        set_normalized_coefficient(gamma_eq[start_con + j], x[start_x + start_λ + i], ub[j])
      end
    end
    # @assert cnt_1 == 2^dim / 2
  end
  num_lb_ub_con = 0
  if cs_w == "cs" && length(cyc) == 3
    num_lb_ub_con = 6
  elseif cs_w == "w" && length(cyc) == 3
    num_lb_ub_con = 9
  elseif cs_w == "cs" && length(cyc) == 4
    num_lb_ub_con = 8
  elseif cs_w == "w" && length(cyc) == 4
    if params["rotate_4w"]
      num_lb_ub_con = 12
    else
      num_lb_ub_con = 8
    end
  end
  for j in 1:length(expairs)
    ep = expairs[j]
    set_normalized_coefficient(gamma_eq[start_con + num_lb_ub_con + j], x[start_x + j], 1)
    for i in 1:(2^dim)
      set_normalized_coefficient(gamma_eq[start_con + num_lb_ub_con + j], x[start_x + start_λ + i], - X[i][ep[1]] * X[i][ep[2]])
    end
  end
  # println("cs_w: ", cs_w)
  # println("start_con: ", start_con)
  # println("start_x: ", start_x)
  # println("num_lb_ub_con: ", num_lb_ub_con)
  # println("gamma_eq: ", gamma_eq)
end

function sp_per_cyc(cyc_size, m_sp, params, cs_w)
  x = create_x(params, m_sp, cyc_size, cs_w)
  A_eq = create_Aeq(params, cyc_size, x, cs_w) # conefficient matrix for equality constraints
  num_constr_eq = size(A_eq)[1]
  b_eq = zeros(num_constr_eq)
  JuMP.@constraint(m_sp, gamma_eq[i = 1:num_constr_eq], sum(A_eq[i][j] * x[j] for j in 1:length(x)) == b_eq[i])

  A_ieq = nothing
  b_ieq = nothing
  num_constr_ieq = 0
  if params["cycle_relax"] == "mc"
    A_ieq = create_Aieq(params, cyc_size, x, cs_w) # coefficient matrix for inequality constraints.
    num_constr_ieq = size(A_ieq)[1]
    b_ieq = zeros(num_constr_ieq)
    JuMP.@constraint(m_sp, gamma_ieq[i = 1:num_constr_ieq], sum(A_ieq[i][j] * x[j] for j in 1:length(x)) >= b_ieq[i])
  end
  JuMP.@objective(m_sp, Min, 0)
  for cyc in cycle[string("cyc_", string(cyc_size))]
    b_eq = create_beq(params, cyc, cs_w)
    # println("cs_w: ", cs_w)
    # println("b_eq: ", b_eq)
    for i in 1:num_constr_eq
      set_normalized_rhs(gamma_eq[i], b_eq[i])
    end
    if params["cycle_relax"] == "mc"
      b_ieq = create_bieq(params, cyc, cs_w)
      for i in 1:num_constr_ieq
        set_normalized_rhs(gamma_ieq[i], b_ieq[i])
      end
    elseif params["cycle_relax"] == "epr"
      modify_Aeq_epr_polytope(cyc, gamma_eq, x, cs_w)
    end
    if params["print_model"] #&& cyc_size == 4
      Write_model(m_sp, "results/separation_model.txt")
    end
    specs["sp_time"] += @elapsed JuMP.optimize!(m_sp)
    sp_status = termination_status(m_sp)
    if sp_status != OPTIMAL
      println("Type of problem: ", cs_w)
      println("Termination status of sub: ", sp_status)
      # Write_model(m_sp, "results/separation_model.txt")
      specs["terminate"] = false
      gamma_eq_bar = round_vector(JuMP.dual.(gamma_eq))
      cut_var_eq = create_cutvareq(params, cyc, cs_w)
      gamma_ieq_bar = 0
      cut_var_ieq = 0
      if params["cycle_relax"] == "mc"
        gamma_ieq_bar = JuMP.dual.(gamma_ieq)
        cut_var_ieq = create_cutvarieq(params, cyc, cs_w)
      end
      # println("cyc_size", string(cyc_size))
      # println("cut_var_ieq: ", cut_var_ieq)
      # println("gamma_ieq_bar: ", gamma_ieq_bar)
      # println("cut LHS: ", sum(gamma_eq_bar[i] * cut_var_eq[i] for i in 1 : num_constr_eq) + sum(gamma_ieq_bar[i] * cut_var_ieq[i] for i in 1 : num_constr_ieq))
      # # if length(cyc) == 4
      # println("cut_var_eq: ", cut_var_eq)
      # println("length cut_var_eq: ", length(cut_var_eq))
      # println("gamma_eq_bar: ", gamma_eq_bar)
      # println("length gamma_eq_bar: ", length(gamma_eq_bar))
      # println("cut LHS value: ", sum(gamma_eq_bar[i] * b_eq[i] for i in 1 : num_constr_eq))
      # println("cut LHS eq: ", sum(gamma_eq_bar[i] * cut_var_eq[i] for i in 1 : num_constr_eq))
      # end
      if params["cycle_relax"] == "mc"
        JuMP.@constraint(m, sum(gamma_eq_bar[i] * cut_var_eq[i] for i in 1 : num_constr_eq) + sum(gamma_ieq_bar[i] * cut_var_ieq[i] for i in 1 : num_constr_ieq) <= - params["sep_tol"])
      elseif params["cycle_relax"] == "epr"
        JuMP.@constraint(m, sum(gamma_eq_bar[i] * cut_var_eq[i] for i in 1 : num_constr_eq) <= - params["sep_tol"])
      end
      specs["num_cuts"] += 1
      if cs_w == "cs"
        specs["num_cs_cuts"] += 1
      elseif cs_w == "w"
        specs["num_w_cuts"] += 1
      end
    end
  end
end

function solve_sp(cyc_size)
  if params["cycle_c_s_cuts"]
    m_sp = Get_sp_solver(params)
    sp_per_cyc(cyc_size, m_sp, params, "cs")
  end
  if params["cycle_wr_wi_cuts"]
    m_sp = Get_sp_solver(params)
    sp_per_cyc(cyc_size, m_sp, params, "w")
  end
end
