function create_constraints(_m, ref, params, bd_data)
    # relaxation of vm square. ots: (5)
    for (i, bus) in ref[:bus]
        quadratic_relax(_m, _m[:w][i], _m[:vm][i])
    end

    # JuMP.@constraint(_m, _m[:z][(1,2)] == 0)
    # JuMP.@constraint(_m, _m[:zc][(1,2,3)] == 0)

    if params["model"] == "ots" || params["model"] == "ots_relax"
        if params["ots_card_lb"] >= 0.5
            JuMP.@constraint(_m, sum((1 - _m[:z][bp]) for bp in keys(ref[:buspairs])) >= params["ots_card_lb"])
            # println("card constr: ", @build_constraint(sum((1 - _m[:z][bp]) for bp in keys(ref[:buspairs])) >= params["ots_card_lb"]))
        end
        if params["ots_card_ub"] != Inf
            JuMP.@constraint(_m, sum((1 - _m[:z][bp]) for bp in keys(ref[:buspairs])) <= params["ots_card_ub"])
        end
    end

    for (bp, buspair) in ref[:buspairs]
        i, j = bp
        td_max = buspair["angmax"]
        td_min = buspair["angmin"]
        td_u = bd_data["td_u"][bp]
        td_M = bd_data["td_M"]
        cos_min = bd_data["cos_min"][bp]
        cos_max = bd_data["cos_max"][bp]
        sin_min = bd_data["sin_min"][bp]
        sin_max = bd_data["sin_max"][bp]

       # Phase-angle diff constraint
       if params["ext_cons"]
           JuMP.@constraint(_m, _m[:va][i] - _m[:va][j] == _m[:td][bp])
       end

        # if params["model"] == "opf"
        #     JuMP.@constraint(_m, _m[:va][i] - _m[:va][j] == _m[:td][bp])
        # elseif params["model"] == "ots" || params["model"] == "ots_relax"
        #     JuMP.@constraint(_m, _m[:va][i] - _m[:va][j] - 2 * td_M * (1 - _m[:z][bp]) <= _m[:td][bp])
        #     JuMP.@constraint(_m, _m[:td][bp] <= _m[:va][i] - _m[:va][j] + 2 * td_M * (1 - _m[:z][bp]))
        # end

       # Cosine and sine relaxations.
        if params["model"] == "opf"
            relaxation_cos(_m, _m[:cs][bp], _m[:td][bp], td_min, td_max, td_u)

            relaxation_sin(_m, _m[:si][bp], _m[:td][bp], td_min, td_max, td_u)

            if params["ext_cons"] == false
                JuMP.@constraint(_m, _m[:va][i] - _m[:va][j] == _m[:td][bp])
            end
        elseif params["model"] == "ots" || params["model"] == "ots_relax"
            relaxation_cos_on_off(_m, bp, td_min, td_max, td_u, td_M, cos_min, cos_max)

            relaxation_sin_on_off(_m, bp, td_min, td_max, td_u, td_M, sin_min, sin_max)

            if params["ext_cons"]
                JuMP.@constraint(_m, _m[:td][bp] <= td_max * _m[:z][bp] + td_M * (1 - _m[:z][bp]))
                JuMP.@constraint(_m, _m[:td][bp] >= td_min * _m[:z][bp] - td_M * (1 - _m[:z][bp]))
            else
                # p.27 Hijazi:
                JuMP.@constraint(_m, _m[:va][i] - _m[:va][j] <= td_max * _m[:z][bp] + (1 - _m[:z][bp]) * td_M)
                JuMP.@constraint(_m, td_min * _m[:z][bp] - td_M * (1 - _m[:z][bp]) <= _m[:va][i] - _m[:va][j])
            end

            JuMP.@constraint(_m, _m[:wr][bp] <= bd_data["wr_max"][bp] * _m[:z][bp])
            JuMP.@constraint(_m, _m[:wr][bp] >= bd_data["wr_min"][bp] * _m[:z][bp])
            JuMP.@constraint(_m, _m[:wi][bp] <= bd_data["wi_max"][bp] * _m[:z][bp])
            JuMP.@constraint(_m, _m[:wi][bp] >= bd_data["wi_min"][bp] * _m[:z][bp])

            # JuMP.@constraint(_m, _m[:cs][bp] <= bd_data["cos_max"][bp] * _m[:z][bp])
            # JuMP.@constraint(_m, _m[:cs][bp] >= bd_data["cos_min"][bp] * _m[:z][bp])
            # JuMP.@constraint(_m, _m[:si][bp] <= bd_data["sin_max"][bp] * _m[:z][bp])
            # JuMP.@constraint(_m, _m[:si][bp] >= bd_data["sin_min"][bp] * _m[:z][bp])
        end

        if params["circ_cons"]
            if params["model"] == "opf"
                JuMP.@constraint(_m, _m[:cs][bp]^2 + _m[:si][bp]^2 <= 1)
            elseif params["model"] == "ots" || params["model"] == "ots_relax"
                JuMP.@constraint(_m, _m[:cs][bp]^2 + _m[:si][bp]^2 <= _m[:z][bp])
            end
        end

        # vals = fix_vals(_m)

        if params["relax"] == "tri"
          # Trilinear convex-hull relaxation
            if params["model"] == "opf"
                conv_tri_sum(_m, _m[:wr][bp], _m[:wi][bp], _m[:vm][i], _m[:vm][j], _m[:cs][bp], _m[:si][bp], _m[:λ_wr][bp,:], _m[:λ_wi][bp,:])
            elseif params["model"] == "ots" || params["model"] == "ots_relax"
                x3_y3_bd = [[bd_data["cos_min"][bp], bd_data["cos_max"][bp]], [bd_data["sin_min"][bp], bd_data["sin_max"][bp]]]
                conv_tri_sum_on_off(_m, _m[:wr][bp], _m[:wi][bp], _m[:vm][i], _m[:vm][j], _m[:cs][bp], _m[:si][bp], _m[:λ_wr][bp,:], _m[:λ_wi][bp,:], _m[:z][bp], x3_y3_bd)
                # JuMP.@constraint(_m, [l in 1:2], sum(_m[:λ_hat][bp, l, k] for k = 1:2) == 1 - z[bp])
            end
            if params["rlt_cuts"]
                conv_tri_sum_rlt_1(_m, _m[:wr][bp], _m[:wi][bp], _m[:vm][i], _m[:vm][j], _m[:cs][bp], _m[:si][bp], _m[:λ_wr][bp,:], _m[:λ_wi][bp,:], _m[:vc_i][bp], _m[:vc_j][bp], _m[:vs_i][bp], _m[:vs_j][bp], _m[:c_λ_wr][bp,:], _m[:s_λ_wi][bp,:])
                conv_tri_sum_rlt_2(_m, _m[:wr][bp], _m[:wi][bp], _m[:vm][i], _m[:vm][j], _m[:cs][bp], _m[:si][bp], _m[:λ_wr][bp,:], _m[:λ_wi][bp,:], _m[:wcc][bp], _m[:wss][bp], _m[:c_λ_wr][bp,:], _m[:s_λ_wi][bp,:])
            end
        elseif params["relax"] == "rmc"
            if params["model"] == "opf"
                mcc(_m, _m[:vv][bp], _m[:vm][i], _m[:vm][j])
                mcc(_m, _m[:wr][bp], _m[:vv][bp], _m[:cs][bp])
                mcc(_m, _m[:wi][bp], _m[:vv][bp], _m[:si][bp])
            elseif params["model"] == "ots" || params["model"] == "ots_relax"
                mcc_ots_wij(_m, bp)
            end
        end

       # Lifted nonlinear cut (k), (l)/ots: (10b), (10c)
        if params["lnc"]
            if params["model"] == "opf"
                relaxation_voltage_lnc(_m, _m[:wr][bp], _m[:wi][bp], _m[:w][i], _m[:w][j], buspair)
            elseif params["model"] == "ots" || params["model"] == "ots_relax"
                if params["ext_cons"]
                    relaxation_voltage_lnc_on_off(_m, _m[:wr][bp], _m[:wi][bp], _m[:w][i], _m[:w][j], buspair, _m[:z][bp], _m[:wz][bp], _m[:wz][(bp[2], bp[1])])
                end
            end
        end
    end

    # Current magnitude R-SOC constraint
    for (i, branch) in ref[:branch]
        bp = (branch["f_bus"], branch["t_bus"])
        buspair = ref[:buspairs][bp]
        tm = branch["tap"]

       # to prevent this constraint from being posted on multiple parallel branches
        if buspair["branch"] == i
          # extract quantities
            g, b = PowerModels.calc_branch_y(branch) # mutual admittance
            tr, ti = PowerModels.calc_branch_t(branch)
            g_fr = branch["g_fr"] # charging admittance
            b_fr = branch["b_fr"]

            if params["extra_rlt"]
                t_idx = (i, branch["t_bus"], branch["f_bus"])
                p_to = _m[:p][t_idx]
                q_to = _m[:q][t_idx]
            end

          # extract variables
            p_fr = _m[:p][(i, branch["f_bus"], branch["t_bus"])]
            q_fr = _m[:q][(i, branch["f_bus"], branch["t_bus"])]
            w_fr = _m[:w][branch["f_bus"]]
            w_to = _m[:w][branch["t_bus"]]

          # R-SOC constraint. Constraint (7d, opf)/(10d, ots)
            JuMP.@constraint(_m, p_fr^2 + q_fr^2 <= w_fr / tm^2 * _m[:cm][bp])
            if params["extra_rlt"]
                JuMP.@constraint(_m, p_to^2 + q_to^2 <= w_to * _m[:cm][(bp[2], bp[1])])
            end
            if params["model"] == "ots" || params["model"] == "ots_relax"
                cm_ub_br = bd_data["cm_ub"][(branch["f_bus"], branch["t_bus"])]
                # JuMP.@constraint(_m, p_fr^2 + q_fr^2 <= ref[:bus][branch["f_bus"]]["vmax"]^2 / tm^2 * _m[:cm][bp] * _m[:z][bp])
                # JuMP.@constraint(_m, p_fr^2 + q_fr^2 <= cm_ub_br / tm^2 * w_fr * _m[:z][bp])
                # Bounds for cm, constraint (t)
                if params["ext_cons"]
                    JuMP.@constraint(_m, _m[:cm][bp] <= cm_ub_br * _m[:z][bp])
                else
                    JuMP.@constraint(_m, _m[:cm][bp] <= cm_ub_br)
                end
                if params["extra_rlt"]
                    cm_ub_br = bd_data["cm_ub"][(branch["t_bus"], branch["f_bus"])]
                    if params["ext_cons"]
                        JuMP.@constraint(_m, _m[:cm][(bp[2], bp[1])] <= cm_ub_br * _m[:z][bp])
                    else
                        JuMP.@constraint(_m, _m[:cm][(bp[2], bp[1])] <= cm_ub_br)
                    end
                end
            end

            ym_sh_sqr = g_fr^2 + b_fr^2

            # constraint (7e) for opf/(n) for ots
            if params["model"] == "opf"
                JuMP.@constraint(_m, _m[:cm][bp] == (g^2 + b^2) * (w_fr / tm^2 + w_to - 2 * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) / tm^2) - ym_sh_sqr * (w_fr / tm^2) + 2 * (g_fr * p_fr - b_fr * q_fr))
                if params["extra_rlt"]
                    JuMP.@constraint(_m, _m[:cm][bp] == (g^2 + b^2) * (w_fr / tm^2 + w_to - 2 * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) / tm^2) + ym_sh_sqr * (w_fr / tm^2) + 2 * (g * g_fr + b * b_fr) * w_fr / tm^2 - 2 * ((g * g_fr + b * b_fr) * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) - (- b * g_fr + g * b_fr) * (- ti * _m[:wr][bp] + tr * _m[:wi][bp])) / tm^2)

                    JuMP.@constraint(_m, _m[:cm][(bp[2], bp[1])] == (g^2 + b^2) * (w_fr / tm^2 + w_to - 2 * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) / tm^2) + ym_sh_sqr * w_to + 2 * (g * g_fr + b * b_fr) * w_to - 2 * ((g * g_fr + b * b_fr) * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) + (- b * g_fr + g * b_fr) * (- ti * _m[:wr][bp] + tr * _m[:wi][bp])) / tm^2)
                end
            elseif params["model"] == "ots" || params["model"] == "ots_relax"
                # if params["ext_cons"]
                JuMP.@constraint(_m, _m[:cm][bp] == (g^2 + b^2) * (_m[:wz][bp] / tm^2 + _m[:wz][(bp[2], bp[1])] - 2 * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) / tm^2) - ym_sh_sqr * (_m[:wz][bp] / tm^2) + 2 * (g_fr * p_fr - b_fr * q_fr))
                # else
                #     JuMP.@constraint(_m, _m[:cm][bp] == (g^2 + b^2) * (_m[:wz][bp] + _m[:wz][(bp[2], bp[1])] - 2 * _m[:wr][bp]))
                # end
                if params["extra_rlt"]
                    JuMP.@constraint(_m, _m[:cm][bp] == (g^2 + b^2) * (_m[:wz][bp] / tm^2 + _m[:wz][(bp[2], bp[1])] - 2 * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) / tm^2) + ym_sh_sqr * (_m[:wz][bp] / tm^2) + 2 * (g * g_fr + b * b_fr) * _m[:wz][bp] / tm^2 - 2 * ((g * g_fr + b * b_fr) * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) - (- b * g_fr + g * b_fr) * (- ti * _m[:wr][bp] + tr * _m[:wi][bp])) / tm^2)

                    JuMP.@constraint(_m, _m[:cm][(bp[2], bp[1])] == (g^2 + b^2) * (_m[:wz][bp] / tm^2 + _m[:wz][(bp[2], bp[1])] - 2 * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) / tm^2) + ym_sh_sqr * _m[:wz][(bp[2], bp[1])] + 2 * (g * g_fr + b * b_fr) * _m[:wz][(bp[2], bp[1])] - 2 * ((g * g_fr + b * b_fr) * (tr * _m[:wr][bp] + ti * _m[:wi][bp]) + (- b * g_fr + g * b_fr) * (- ti * _m[:wr][bp] + tr * _m[:wi][bp])) / tm^2)
                end
                # Bounds for p and q, constraints (r) & (s) -> Possibly redundant. Removed.
                # JuMP.@constraint(_m, p_fr <= branch["rate_a"] * _m[:z][bp])
                # JuMP.@constraint(_m, q_fr <= branch["rate_a"] * _m[:z][bp])
                # JuMP.@constraint(_m, p_fr >= - branch["rate_a"] * _m[:z][bp])
                # JuMP.@constraint(_m, q_fr >= - branch["rate_a"] * _m[:z][bp])
            end
        end
        
    end

    # Reference bus - theta constraint
    for i in keys(ref[:ref_buses])
        JuMP.@constraint(_m, _m[:va][i] == 0)
    end

    # Kirchoffs current law constraints
    for (i, bus) in ref[:bus]
          # Bus KCL
        gs = 0.0
        bs = 0.0
        if length(ref[:bus_shunts][i]) > 0
            shunt_num = ref[:bus_shunts][i][1]
            gs = ref[:shunt][shunt_num]["gs"]
            bs = ref[:shunt][shunt_num]["bs"]
        end
        # ots (2b)
        JuMP.@constraint(_m, sum(_m[:p][a] for a in ref[:bus_arcs][i])  == sum(_m[:pg][g] for g in ref[:bus_gens][i]) - sum(load["pd"] * params["load_multiplier"] for (l, load) in ref[:load] if load["load_bus"] == i)  - gs * _m[:w][i])
        JuMP.@constraint(_m, sum(_m[:q][a] for a in ref[:bus_arcs][i]) == sum(_m[:qg][g] for g in ref[:bus_gens][i]) - sum(load["qd"] * params["load_multiplier"] for (l, load) in ref[:load] if load["load_bus"] == i) + bs * _m[:w][i])
    end

    # Valid inequality: if net load at a node > 0, then at least one line connected to it needs to be turned on.
    if params["valid_ieq"]
        for (l, load) in ref[:load]
            load_bus = load["load_bus"]
            load_diff_pd = 0.0
            if length(ref[:bus_gens][load_bus]) == 0
                load_diff_pd = load["pd"] * params["load_multiplier"]
            else
                load_diff_pd = load["pd"] * params["load_multiplier"] - sum(ref[:gen][g]["pmax"] for g in ref[:bus_gens][load_bus])
            end
            load_diff_qd = 0.0
            if length(ref[:bus_gens][load_bus]) == 0
                load_diff_qd = load["qd"] * params["load_multiplier"]
            else
                load_diff_qd = load["qd"] * params["load_multiplier"] - sum(ref[:gen][g]["qmax"] for g in ref[:bus_gens][load_bus])
            end
            load_diff = max(load_diff_pd, load_diff_qd)
            if load_diff > 1e-3
                expr = 0
                for bp in keys(ref[:buspairs])
                    if bp[1] == l || bp[2] == l
                        expr += _m[:z][bp]
                    end
                end
                JuMP.@constraint(_m, expr >= 1)
            end
        end
    end

    # Ohms laws, voltage-angle difference and thermal limits
    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])
        bp_idx2 = (branch["t_bus"], branch["f_bus"])

        p_fr = _m[:p][f_idx]
        q_fr = _m[:q][f_idx]
        p_to = _m[:p][t_idx]
        q_to = _m[:q][t_idx]

        w_fr = _m[:w][branch["f_bus"]]
        w_to = _m[:w][branch["t_bus"]]
        wr_br = _m[:wr][bp_idx]
        wi_br = _m[:wi][bp_idx]

       # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        g_to = branch["g_to"]
        b_fr = branch["b_fr"]
        b_to = branch["b_to"]
        tm = branch["tap"]

       # AC Line Flow Constraints
       if params["model"] == "opf"
            JuMP.@constraint(_m, p_fr ==  (g + g_fr) / tm^2 * w_fr + (-g * tr + b * ti) / tm^2 * wr_br + (-b * tr - g * ti) / tm^2 * wi_br )
            JuMP.@constraint(_m, q_fr == -(b + b_fr) / tm^2 * w_fr - (-b * tr - g * ti) / tm^2 * wr_br + (-g * tr + b * ti) / tm^2 * wi_br )

            JuMP.@constraint(_m, p_to ==  (g + g_to) * w_to + (-g * tr - b * ti) / tm^2 * wr_br + (-b * tr + g * ti) / tm^2 * -wi_br )
            JuMP.@constraint(_m, q_to == -(b + b_to) * w_to - (-b * tr + g * ti) / tm^2 * wr_br + (-g * tr - b * ti) / tm^2 * -wi_br )

            # Apparent PoLimit, From and To
            JuMP.@constraint(_m, p_fr^2 + q_fr^2 <= branch["rate_a"]^2)
            JuMP.@constraint(_m, p_to^2 + q_to^2 <= branch["rate_a"]^2)
        elseif params["model"] == "ots" || params["model"] == "ots_relax"
            wz_ij = _m[:wz][bp_idx]
            wz_ji = _m[:wz][bp_idx2]
            z_ij = _m[:z][bp_idx]
            # constraint (3a) (3b)
            JuMP.@constraint(_m, p_fr ==  (g + g_fr) / tm^2 * wz_ij + (-g * tr + b * ti) / tm^2 * wr_br + (-b * tr - g * ti) / tm^2 * wi_br )
            JuMP.@constraint(_m, q_fr == -(b + b_fr) / tm^2 * wz_ij - (-b * tr - g * ti) / tm^2 * wr_br + (-g * tr + b * ti) / tm^2 * wi_br )
            # constraint (3c)
            JuMP.@constraint(_m, wz_ij <= w_fr - (1 - z_ij) * ref[:bus][branch["f_bus"]]["vmin"]^2)
            JuMP.@constraint(_m, wz_ij >= w_fr - (1 - z_ij) * ref[:bus][branch["f_bus"]]["vmax"]^2)
            # constraint (3d)
            JuMP.@constraint(_m, wz_ji <= w_to - (1 - z_ij) * ref[:bus][branch["t_bus"]]["vmin"]^2)
            JuMP.@constraint(_m, wz_ji >= w_to - (1 - z_ij) * ref[:bus][branch["t_bus"]]["vmax"]^2)
            # if params["ext_cons"]
                # constraint (3e)
            JuMP.@constraint(_m, wz_ij >= ref[:bus][branch["f_bus"]]["vmin"]^2 * z_ij)
            JuMP.@constraint(_m, wz_ij <= ref[:bus][branch["f_bus"]]["vmax"]^2 * z_ij)
            # constraint (3f)
            JuMP.@constraint(_m, wz_ji >= ref[:bus][branch["t_bus"]]["vmin"]^2 * z_ij)
            JuMP.@constraint(_m, wz_ji <= ref[:bus][branch["t_bus"]]["vmax"]^2 * z_ij)
            # end
            # constraint (d)
            JuMP.@constraint(_m, p_to ==  (g + g_to) * wz_ji + (-g * tr - b * ti) / tm^2 * wr_br + (-b * tr + g * ti) / tm^2 * -wi_br )
            JuMP.@constraint(_m, q_to == -(b + b_to) * wz_ji - (-b * tr + g * ti) / tm^2 * wr_br + (-g * tr - b * ti) / tm^2 * -wi_br )

            # constraints (2j)
            JuMP.@constraint(_m, p_fr^2 + q_fr^2 <= branch["rate_a"]^2 * _m[:z][bp_idx]^2)
            JuMP.@constraint(_m, p_to^2 + q_to^2 <= branch["rate_a"]^2 * _m[:z][bp_idx]^2)
        end

        # if params["ext_cons"]
           # Phase Angle Difference Limit
        JuMP.@constraint(_m, wi_br <= tan(branch["angmax"]) * wr_br)
        JuMP.@constraint(_m, wi_br >= tan(branch["angmin"]) * wr_br)
        # end
    end

    # Relationship between cs/w variables and their reversed counterparts.
    for bp in keys(ref[:buspairs])
        JuMP.@constraint(_m, _m[:si][bp] == - _m[:si][(bp[2], bp[1])])
        JuMP.@constraint(_m, _m[:wi][bp] == - _m[:wi][(bp[2], bp[1])])
        JuMP.@constraint(_m, _m[:cs][bp] == _m[:cs][(bp[2], bp[1])])
        JuMP.@constraint(_m, _m[:wr][bp] == _m[:wr][(bp[2], bp[1])])
    end

    # heuristic: all branches in the spanning tree are turned on.
    if params["model"] == "ots" && params["span_tree"]
        span_tree_edges = Get_span_tree(ref)
        for bp in span_tree_edges 
            # The tree edges may have diffrent order than original edges.
            if bp in keys(ref[:buspairs])
                JuMP.@constraint(_m, _m[:z][bp] == 1)
            else
                JuMP.@constraint(_m, _m[:z][(bp[2], bp[1])] == 1)
            end
        end
    end

    if params["off_insights"]
        on_lines_li, off_lines_li = on_off_lines(instance, data, ref)
        for bp in on_lines_li
            JuMP.@constraint(_m, _m[:z][bp] == 1)
        end
        for bp in off_lines_li
            JuMP.@constraint(_m, _m[:z][bp] == 0)
        end
    end

    if params["debug"]
        # JuMP.@constraint(_m, _m[:vm][2] == 1.0100704619038228)
        # JuMP.@constraint(_m, _m[:vm][3] == 0.9928105071732412)
        # JuMP.@constraint(_m, _m[:vm][1] == 1.0912817221698126)
        JuMP.@constraint(_m, _m[:z][(2, 25)] == 0)
        JuMP.@constraint(_m, _m[:z][(12, 11)] == 0)
        JuMP.@constraint(_m, _m[:z][(3, 18)] == 0)
        JuMP.@constraint(_m, _m[:z][(4, 5)] == 0)
    end

    # JuMP.@constraint(_m, _m[:p][(33, 25, 24)] == 0)
    # JuMP.@constraint(_m, _m[:q][(33, 25, 24)] == 0)
end

# Debug PM.jl QC RM formulation.
function fix_vals(_m)
    vals = Dict()
    line_ind = Dict()
    vals[:va] = Dict()
    vals[:td] = Dict()
    vals[:cs] = Dict()
    vals[:si] = Dict()
    vals[:vm] = Dict()
    vals[:vm_fr] = Dict()
    vals[:vm_to] = Dict()
    vals[:vv] = Dict()
    vals[:wr] = Dict()
    vals[:wi] = Dict()
    vals[:z] = Dict()
    line_ind[1] = [1, 3]
    line_ind[2] = [3, 2]
    line_ind[3] = [1, 2]
    # vals[:ccm][1] = 0.29448786542959154
    # vals[:ccm][2] = 0.308641975308642
    # vals[:ccm][3] = 0.04507168667867171
    vals[:cs][1] = 0.9547660302269865
    vals[:cs][2] = 0.9469878431334768
    vals[:cs][3] = 0.9469878431334768
    # vals[:p][(1, 1, 3)] = 0.5458393612121765
    # vals[:p][(1, 3, 1)] = -0.527218013642359
    # vals[:p][(2, 2, 3)] = 0.4288568919368097
    # vals[:p][(2, 3, 2)] = -0.4227819863576409
    # vals[:p][(3, 1, 2)] = -0.19153614279130424
    # vals[:p][(3, 2, 1)] = 0.19307442092673677
    # vals[:pg][1] = 1.4543032184208722
    # vals[:pg][2] = 1.7219313128635465
    # vals[:pg][3] = 0.0
    # vals[:q][(1, 1, 3)] = -0.138989054545674
    # vals[:q][(1, 3, 1)] = -0.14842068968279487
    # vals[:q][(2, 2, 3)] = -0.257064023528016
    # vals[:q][(2, 3, 2)] = -0.26693714130794494
    # vals[:q][(3, 1, 2)] = -0.10895261806295053
    # vals[:q][(3, 2, 1)] = -0.17393988290488926
    # vals[:qg][1] = 0.15205832739137548
    # vals[:qg][2] = -0.031003906432905237
    # vals[:qg][3] = 0.08464216900926019
    vals[:si][1] = 0.22699315399039155
    vals[:si][2] = -0.32125191859422036
    vals[:si][3] = -0.0986205953093676
    vals[:td][1] = 0.2285890557274997
    vals[:td][2] = -0.32706963978247594
    vals[:td][3] = -0.09848058405497623
    vals[:va][1] = 0.0
    vals[:va][2] = 0.09848058405497623
    vals[:va][3] = -0.2285890557274997
    vals[:vm][1] = 1.0336615251440735
    vals[:vm][2] = 1.009962395540597
    vals[:vm][3] = 0.9929674075523719
    vals[:vm_fr][1] = 1.0336615251440735
    vals[:vm_fr][2] = 0.9929674075523719
    vals[:vm_fr][3] = 1.0336615251440735
    vals[:vm_to][1] = 0.9929674075523719
    vals[:vm_to][2] = 1.009962395540597
    vals[:vm_to][3] = 1.009962395540597
    vals[:vv][1] = 1.0324632435378231
    vals[:vv][2] = 0.9932227834022657
    vals[:vv][3] = 1.0379863127531375
    # vals[:w][1] = 1.077323050288147
    # vals[:w][2] = 1.0283809755117457
    # vals[:w][3] = 0.9894714022340392
    # vals[:w_fr][1] = 1.077323050288147
    # vals[:w_fr][2] = 0.9894714022340392
    # vals[:w_fr][3] = 1.077323050288147
    # vals[:w_to][1] = 0.9894714022340392
    # vals[:w_to][2] = 1.0283809755117457
    # vals[:w_to][3] = 1.0283809755117457
    vals[:wi][1] = 0.33169884288655405
    vals[:wi][2] = -0.3190709360050799
    vals[:wi][3] = -0.1745936537703452
    vals[:wr][1] = 0.9777301401124768
    vals[:wr][2] = 0.9405075647875042
    vals[:wr][3] = 1.0379863127531375
    vals[:z][1] = 1.0
    vals[:z][2] = 1.0
    vals[:z][3] = 1.0
    for i in 1:3
        fr_ind = line_ind[i][1]
        to_ind = line_ind[i][2]
        JuMP.@constraint(_m, _m[:va][i] == vals[:va][i])
        JuMP.@constraint(_m, _m[:td][(fr_ind, to_ind)] == vals[:td][i])

        JuMP.@constraint(_m, _m[:cs][(fr_ind, to_ind)] == vals[:cs][i])
        JuMP.@constraint(_m, _m[:si][(fr_ind, to_ind)] == vals[:si][i])
        JuMP.@constraint(_m, _m[:vv][(fr_ind, to_ind)] == vals[:vv][i])
        JuMP.@constraint(_m, _m[:wr][(fr_ind, to_ind)] == vals[:wr][i])
        # JuMP.@constraint(_m, _m[:wi][(fr_ind, to_ind)] == vals[:wi][i])

        JuMP.@constraint(_m, _m[:vm][i] == vals[:vm][i])
        JuMP.@constraint(_m, _m[:vmz][(fr_ind, to_ind)] == vals[:vm_fr][i])
        JuMP.@constraint(_m, _m[:vmz][(to_ind, fr_ind)] == vals[:vm_to][i])
        JuMP.@constraint(_m, _m[:z][(fr_ind, to_ind)] == vals[:z][i])
    end
    return vals
end

# function fix_vals_old(_m)
#     vals = Dict()
#     line_ind = Dict()
#     vals[:va] = Dict()
#     vals[:td] = Dict()
#     vals[:cs] = Dict()
#     vals[:si] = Dict()
#     vals[:vm] = Dict()
#     vals[:vm_fr] = Dict()
#     vals[:vm_to] = Dict()
#     vals[:vv] = Dict()
#     vals[:wr] = Dict()
#     vals[:wi] = Dict()
#     vals[:z] = Dict()
#     line_ind[1] = [1, 3]
#     line_ind[2] = [3, 2]
#     line_ind[3] = [1, 2]
#     # ccm[1] = 0.2774682418789966
#     # ccm[2] = 0.2673401396308076
#     # ccm[3] = 0.10498549961398765
#     vals[:cs][1] = 0.9469878431334768
#     vals[:cs][2] = 0.9678273132191462
#     vals[:cs][3] = 0.9634474671087025
#     # p[(1, 1, 3)] = 0.5616719939225907
#     # p[(1, 3, 1)] = -0.5425586483528695
#     # p[(2, 2, 3)] = 0.4121209265138007
#     # p[(2, 3, 2)] = -0.4074413516471305
#     # p[(3, 1, 2)] = -0.34907489060457825
#     # p[(3, 2, 1)] = 0.35421007040703434
#     # pg[1] = 1.3125971033180128
#     # pg[2] = 1.866330996920835
#     # pg[3] = 0.0
#     # q[(1, 1, 3)] = -0.09469355913120112
#     # q[(1, 3, 1)] = -0.21147567296848502
#     # q[(2, 2, 3)] = -0.2831204945088235
#     # q[(2, 3, 2)] = -0.2898141533008884
#     # q[(3, 1, 2)] = -0.03009449803881359
#     # q[(3, 2, 1)] = -0.19070836297411964
#     # qg[1] = 0.2752119428299853
#     # qg[2] = -0.0738288574829431
#     # qg[3] = -0.0012898262693734064
#     vals[:si][1] = 0.24756055534051882
#     vals[:si][2] = -0.30816419233084996
#     vals[:si][3] = -0.23479220437383674
#     vals[:td][1] = 0.32706963978247594
#     vals[:td][2] = -0.23649356847004327
#     vals[:td][3] = 0.09057607131243266
#     vals[:va][1] = 0.0
#     vals[:va][2] = -0.09057607131243266
#     vals[:va][3] = -0.32706963978247594
#     vals[:vm][1] = 1.0796452348982217
#     vals[:vm][2] = 1.018001142204249
#     vals[:vm][3] = 0.995868417741324
#     vals[:vm_fr][1] = 1.0796452348982217
#     vals[:vm_fr][2] = 0.995868417741324
#     vals[:vm_fr][3] = 1.0796452348982217
#     vals[:vm_to][1] = 0.995868417741324
#     vals[:vm_to][2] = 1.018001142204249
#     vals[:vm_to][3] = 1.0180011422042492
#     vals[:vv][1] = 1.0925204861247086
#     vals[:vv][2] = 1.0052565159401303
#     vals[:vv][3] = 1.0974110148127179
#     # w[1] = 1.1692904697964432
#     # w[2] = 1.0363257182596943
#     # w[3] = 1.0017368354826477
#     # w_fr[1] = 1.1692904697964432
#     # w_fr[2] = 1.0017368354826477
#     # w_fr[3] = 1.1692904697964432
#     # w_to[1] = 1.0017368354826477
#     # w_to[2] = 1.0363257182596943
#     # w_to[3] = 1.0363257182596943
#     vals[:wi][1] = 0.3372908444547613
#     vals[:wi][2] = -0.3071008572132988
#     vals[:wi][3] = -0.3202699625862079
#     vals[:wr][1] = 1.0283757763162156
#     vals[:wr][2] = 0.9663275649352973
#     vals[:wr][3] = 1.053182450014248
#     vals[:z][1] = 1.0
#     vals[:z][2] = 1.0
#     vals[:z][3] = 1.0
#     for i in 1:3
#         fr_ind = line_ind[i][1]
#         to_ind = line_ind[i][2]
#         JuMP.@constraint(_m, _m[:va][i] == vals[:va][i])
#         JuMP.@constraint(_m, _m[:td][(fr_ind, to_ind)] == vals[:td][i])

#         JuMP.@constraint(_m, _m[:cs][(fr_ind, to_ind)] == vals[:cs][i])
#         JuMP.@constraint(_m, _m[:si][(fr_ind, to_ind)] == vals[:si][i])
#         JuMP.@constraint(_m, _m[:vv][(fr_ind, to_ind)] == vals[:vv][i])
#         JuMP.@constraint(_m, _m[:wr][(fr_ind, to_ind)] == vals[:wr][i])
#         # JuMP.@constraint(_m, _m[:wi][(fr_ind, to_ind)] == vals[:wi][i])

#         JuMP.@constraint(_m, _m[:vm][i] == vals[:vm][i])
#         JuMP.@constraint(_m, _m[:vmz][(fr_ind, to_ind)] == vals[:vm_fr][i])
#         JuMP.@constraint(_m, _m[:vmz][(to_ind, fr_ind)] == vals[:vm_to][i])
#         JuMP.@constraint(_m, _m[:z][(fr_ind, to_ind)] == vals[:z][i])
#     end
# end
