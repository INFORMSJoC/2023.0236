# Get all constraints corresponding to the cycle constraints

function create_cyc_constr(_m, ref, params)
    if params["cycle_relax"] == "mc"
        if params["cycle_max_bnd"] >= 3
            arcpairs = []
            busarcpairs = []
            for cyc in cycle["cyc_3"]
                push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
                push!(busarcpairs, (cyc[2], cyc[1], cyc[3]))
                if params["rotate"]
                    push!(arcpairs, ((cyc[2], cyc[3]), (cyc[1], cyc[3])))
                    push!(arcpairs, ((cyc[1], cyc[3]), (cyc[1], cyc[2])))
                    push!(busarcpairs, (cyc[3], cyc[1], cyc[2]))
                    push!(busarcpairs, (cyc[1], cyc[2], cyc[3]))
                end
            end
            if params["cycle_c_s_cuts"]
            # McCormick constrs
                @constraints(_m, begin
                    [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[3])] == _m[:hcc][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] - _m[:hss][((cyc[1], cyc[2]), (cyc[2], cyc[3]))]
                    [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[3])] == _m[:hcs][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] + _m[:hsc][((cyc[1], cyc[2]), (cyc[2], cyc[3]))]
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[2])] == _m[:hcc][((cyc[2], cyc[3]), (cyc[1], cyc[3]))] + _m[:hss][((cyc[2], cyc[3]), (cyc[1], cyc[3]))]
                        [cyc in cycle["cyc_3"]], - _m[:si][(cyc[1], cyc[2])] == - _m[:hcs][((cyc[2], cyc[3]), (cyc[1], cyc[3]))] + _m[:hsc][((cyc[2], cyc[3]), (cyc[1], cyc[3]))]
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[2], cyc[3])] == _m[:hcc][((cyc[1], cyc[3]), (cyc[1], cyc[2]))] + _m[:hss][((cyc[1], cyc[3]), (cyc[1], cyc[2]))]
                        [cyc in cycle["cyc_3"]], - _m[:si][(cyc[2], cyc[3])] == _m[:hcs][((cyc[1], cyc[3]), (cyc[1], cyc[2]))] - _m[:hsc][((cyc[1], cyc[3]), (cyc[1], cyc[2]))]
                    end)
                end
                for ap in arcpairs
                    mcc(_m, _m[:hcc][ap], _m[:cs][ap[1]], _m[:cs][ap[2]])
                    mcc(_m, _m[:hss][ap], _m[:si][ap[1]], _m[:si][ap[2]])
                    mcc(_m, _m[:hcs][ap], _m[:cs][ap[1]], _m[:si][ap[2]])
                    mcc(_m, _m[:hsc][ap], _m[:si][ap[1]], _m[:cs][ap[2]])
                end
            end

            if params["cycle_wr_wi_cuts"]
                # McCormick constrs w version
                @constraints(_m, begin
                    [cyc in cycle["cyc_3"]], _m[:wwr][(cyc[2], cyc[1], cyc[3])] == _m[:wrr][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] - _m[:wii][((cyc[1], cyc[2]), (cyc[2], cyc[3]))]
                    [cyc in cycle["cyc_3"]], _m[:wwi][(cyc[2], cyc[1], cyc[3])] == _m[:wri][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] + _m[:wir][((cyc[1], cyc[2]), (cyc[2], cyc[3]))]
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_3"]], _m[:wwr][(cyc[3], cyc[1], cyc[2])] == _m[:wrr][((cyc[2], cyc[3]), (cyc[1], cyc[3]))] + _m[:wii][((cyc[2], cyc[3]), (cyc[1], cyc[3]))]
                        [cyc in cycle["cyc_3"]], - _m[:wwi][(cyc[3], cyc[1], cyc[2])] == - _m[:wri][((cyc[2], cyc[3]), (cyc[1], cyc[3]))] + _m[:wir][((cyc[2], cyc[3]), (cyc[1], cyc[3]))]
                        [cyc in cycle["cyc_3"]], _m[:wwr][(cyc[1], cyc[2], cyc[3])] == _m[:wrr][((cyc[1], cyc[3]), (cyc[1], cyc[2]))] + _m[:wii][((cyc[1], cyc[3]), (cyc[1], cyc[2]))]
                        [cyc in cycle["cyc_3"]], - _m[:wwi][(cyc[1], cyc[2], cyc[3])] == _m[:wri][((cyc[1], cyc[3]), (cyc[1], cyc[2]))] - _m[:wir][((cyc[1], cyc[3]), (cyc[1], cyc[2]))]
                    end)
                end
                for ba in busarcpairs
                    mcc(_m, _m[:wwr][ba], _m[:w][ba[1]], _m[:wr][(ba[2], ba[3])])
                    mcc(_m, _m[:wwi][ba], _m[:w][ba[1]], _m[:wi][(ba[2], ba[3])])
                end
                for ap in arcpairs
                    mcc(_m, _m[:wrr][ap], _m[:wr][ap[1]], _m[:wr][ap[2]])
                    mcc(_m, _m[:wii][ap], _m[:wi][ap[1]], _m[:wi][ap[2]])
                    mcc(_m, _m[:wri][ap], _m[:wr][ap[1]], _m[:wi][ap[2]])
                    mcc(_m, _m[:wir][ap], _m[:wi][ap[1]], _m[:wr][ap[2]])
                end
            end
            if params["cycle_sw_cuts"]
                # JuMP.@constraint(_m, )
            end
        end
        if params["cycle_max_bnd"] >= 4
            arcpairs = []
            busarctuples = []
            for cyc in cycle["cyc_4"]
                push!(arcpairs, ((cyc[1], cyc[2]), (cyc[3], cyc[4])))
                push!(arcpairs, ((cyc[1], cyc[4]), (cyc[2], cyc[3])))
                push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
                push!(arcpairs, ((cyc[1], cyc[4]), (cyc[3], cyc[4])))
                push!(arcpairs, ((cyc[2], cyc[3]), (cyc[3], cyc[4])))
                push!(arcpairs, ((cyc[1], cyc[4]), (cyc[1], cyc[2])))
                push!(busarctuples, (cyc[4], (cyc[1], cyc[2]), (cyc[2], cyc[3])))
                push!(busarctuples, (cyc[2], (cyc[1], cyc[4]), (cyc[3], cyc[4])))
                push!(busarctuples, (cyc[1], (cyc[2], cyc[3]), (cyc[3], cyc[4])))
                push!(busarctuples, (cyc[3], (cyc[1], cyc[4]), (cyc[1], cyc[2])))
            end
            if params["cycle_c_s_cuts"]
            # McCormick constrs
                @constraints(_m, begin
                    [cyc in cycle["cyc_4"]], _m[:hcc][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] - _m[:hss][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] == _m[:hcc][((cyc[1], cyc[4]), (cyc[2], cyc[3]))] + _m[:hss][((cyc[1], cyc[4]), (cyc[2], cyc[3]))]
                    [cyc in cycle["cyc_4"]], _m[:hcs][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] + _m[:hsc][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] == - _m[:hcs][((cyc[1], cyc[4]), (cyc[2], cyc[3]))] + _m[:hsc][((cyc[1], cyc[4]), (cyc[2], cyc[3]))]
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_4"]], _m[:hcc][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] - _m[:hss][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] == _m[:hcc][((cyc[1], cyc[4]), (cyc[3], cyc[4]))] + _m[:hss][((cyc[1], cyc[4]), (cyc[3], cyc[4]))]
                        [cyc in cycle["cyc_4"]], _m[:hcs][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] + _m[:hsc][((cyc[1], cyc[2]), (cyc[2], cyc[3]))] == - _m[:hcs][((cyc[1], cyc[4]), (cyc[3], cyc[4]))] + _m[:hsc][((cyc[1], cyc[4]), (cyc[3], cyc[4]))]
                        [cyc in cycle["cyc_4"]], _m[:hcc][((cyc[2], cyc[3]), (cyc[3], cyc[4]))] - _m[:hss][((cyc[2], cyc[3]), (cyc[3], cyc[4]))] == _m[:hcc][((cyc[1], cyc[4]), (cyc[1], cyc[2]))] + _m[:hss][((cyc[1], cyc[4]), (cyc[1], cyc[2]))]
                        [cyc in cycle["cyc_4"]], _m[:hcs][((cyc[2], cyc[3]), (cyc[3], cyc[4]))] + _m[:hsc][((cyc[2], cyc[3]), (cyc[3], cyc[4]))] == - _m[:hcs][((cyc[1], cyc[4]), (cyc[1], cyc[2]))] + _m[:hsc][((cyc[1], cyc[4]), (cyc[1], cyc[2]))]
                    end)
                end
                for ap in arcpairs
                    mcc(_m, _m[:hcc][ap], _m[:cs][ap[1]], _m[:cs][ap[2]])
                    mcc(_m, _m[:hss][ap], _m[:si][ap[1]], _m[:si][ap[2]])
                    mcc(_m, _m[:hcs][ap], _m[:cs][ap[1]], _m[:si][ap[2]])
                    mcc(_m, _m[:hsc][ap], _m[:si][ap[1]], _m[:cs][ap[2]])
                end
            end
            if params["cycle_wr_wi_cuts"]
                # McCormick constrs w version
                @constraints(_m, begin
                    [cyc in cycle["cyc_4"]], _m[:wrr][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] - _m[:wii][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] == _m[:wrr][((cyc[1], cyc[4]), (cyc[2], cyc[3]))] + _m[:wii][((cyc[1], cyc[4]), (cyc[2], cyc[3]))]
                    [cyc in cycle["cyc_4"]], _m[:wri][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] + _m[:wir][((cyc[1], cyc[2]), (cyc[3], cyc[4]))] == - _m[:wri][((cyc[1], cyc[4]), (cyc[2], cyc[3]))] + _m[:wir][((cyc[1], cyc[4]), (cyc[2], cyc[3]))]
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_4"]], _m[:wwrr][(cyc[4], (cyc[1], cyc[2]), (cyc[2], cyc[3]))] - _m[:wwii][(cyc[4], (cyc[1], cyc[2]), (cyc[2], cyc[3]))] == _m[:wwrr][(cyc[2], (cyc[1], cyc[4]), (cyc[3], cyc[4]))] + _m[:wwii][(cyc[2], (cyc[1], cyc[4]), (cyc[3], cyc[4]))]
                        [cyc in cycle["cyc_4"]], _m[:wwri][(cyc[4], (cyc[1], cyc[2]), (cyc[2], cyc[3]))] + _m[:wwir][(cyc[4], (cyc[1], cyc[2]), (cyc[2], cyc[3]))] == - _m[:wwri][(cyc[2], (cyc[1], cyc[4]), (cyc[3], cyc[4]))] + _m[:wwir][(cyc[2], (cyc[1], cyc[4]), (cyc[3], cyc[4]))]
                        # [cyc in cycle["cyc_4"]], _m[:wwrr][(cyc[1], (cyc[2], cyc[3]), (cyc[3], cyc[4]))] - _m[:wwii][(cyc[1], (cyc[2], cyc[3]), (cyc[3], cyc[4]))] == _m[:wwrr][(cyc[3], (cyc[1], cyc[4]), (cyc[1], cyc[2]))] + _m[:wwii][(cyc[3], (cyc[1], cyc[4]), (cyc[1], cyc[2]))]
                        # [cyc in cycle["cyc_4"]], _m[:wwri][(cyc[1], (cyc[2], cyc[3]), (cyc[3], cyc[4]))] + _m[:wwir][(cyc[1], (cyc[2], cyc[3]), (cyc[3], cyc[4]))] == - _m[:wwri][(cyc[3], (cyc[1], cyc[4]), (cyc[1], cyc[2]))] + _m[:wwir][(cyc[3], (cyc[1], cyc[4]), (cyc[1], cyc[2]))]
                    end)
                end
                for ap in arcpairs
                    mcc(_m, _m[:wrr][ap], _m[:wr][ap[1]], _m[:wr][ap[2]])
                    mcc(_m, _m[:wii][ap], _m[:wi][ap[1]], _m[:wi][ap[2]])
                    mcc(_m, _m[:wri][ap], _m[:wr][ap[1]], _m[:wi][ap[2]])
                    mcc(_m, _m[:wir][ap], _m[:wi][ap[1]], _m[:wr][ap[2]])
                end
                for ba in busarctuples
                    mcc(_m, _m[:wwrr][ba], _m[:w][ba[1]], _m[:wrr][(ba[2], ba[3])])
                    mcc(_m, _m[:wwii][ba], _m[:w][ba[1]], _m[:wii][(ba[2], ba[3])])
                    mcc(_m, _m[:wwri][ba], _m[:w][ba[1]], _m[:wri][(ba[2], ba[3])])
                    mcc(_m, _m[:wwir][ba], _m[:w][ba[1]], _m[:wir][(ba[2], ba[3])])
                end
            end
        end
    elseif params["cycle_relax"] == "epr"
        if params["model"] == "ots" || params["model"] == "ots_relax"
            link_zc_z(_m, params)
        end
        # Extrem point representation.
        if params["cycle_max_bnd"] >= 3
            expairs_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6)]
            expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
            if params["cycle_c_s_cuts"]
                if params["model"] == "opf"
                    @constraints(_m, begin
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[3])] == _m[:x_c][cyc, (1,2)] - _m[:x_c][cyc, (4,5)]
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[3])] == _m[:x_c][cyc, (1,5)] + _m[:x_c][cyc, (2,4)]
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[2])] == _m[:x_c][cyc, (2,3)] + _m[:x_c][cyc, (5,6)]
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[2])] == _m[:x_c][cyc, (2,6)] - _m[:x_c][cyc, (3,5)]
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[2], cyc[3])] == _m[:x_c][cyc, (1,3)] + _m[:x_c][cyc, (4,6)]
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[2], cyc[3])] == _m[:x_c][cyc, (1,6)] - _m[:x_c][cyc, (3,4)]
                        [cyc in cycle["cyc_3"]], sum(_m[:λ_c][cyc, i] for i in 1:(2^6)) == 1
                    end)
                elseif params["model"] == "ots" || params["model"] == "ots_relax"
                    @constraints(_m, begin
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[3])] - _m[:x_c][cyc, (1,2)] + _m[:x_c][cyc, (4,5)] <= (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[3])] - _m[:x_c][cyc, (1,2)] + _m[:x_c][cyc, (4,5)] >= - (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[3])] - _m[:x_c][cyc, (1,5)] - _m[:x_c][cyc, (2,4)] <= (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[3])] - _m[:x_c][cyc, (1,5)] - _m[:x_c][cyc, (2,4)] >= - (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[2])] - _m[:x_c][cyc, (2,3)] - _m[:x_c][cyc, (5,6)] <= (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[2])] - _m[:x_c][cyc, (2,3)] - _m[:x_c][cyc, (5,6)] >= - (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[2])] - _m[:x_c][cyc, (2,6)] + _m[:x_c][cyc, (3,5)] <= (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[2])] - _m[:x_c][cyc, (2,6)] + _m[:x_c][cyc, (3,5)] >= - (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[2], cyc[3])] - _m[:x_c][cyc, (1,3)] - _m[:x_c][cyc, (4,6)] <= (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[2], cyc[3])] - _m[:x_c][cyc, (1,3)] - _m[:x_c][cyc, (4,6)] >= - (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[2], cyc[3])] - _m[:x_c][cyc, (1,6)] + _m[:x_c][cyc, (3,4)] <= (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], _m[:si][(cyc[2], cyc[3])] - _m[:x_c][cyc, (1,6)] + _m[:x_c][cyc, (3,4)] >= - (1 - _m[:zc][Tuple(cyc)])
                        [cyc in cycle["cyc_3"]], sum(_m[:λ_c][cyc, i] for i in 1:(2^6)) == _m[:zc][Tuple(cyc)]
                    end)
                end
                for cyc in cycle["cyc_3"]
                    buspairs = [(cyc[1], cyc[2]), (cyc[2], cyc[3]), (cyc[1], cyc[3])]
                    x_li = [_m[:cs][buspairs[1]], _m[:cs][buspairs[2]], _m[:cs][buspairs[3]], _m[:si][buspairs[1]], _m[:si][buspairs[2]], _m[:si][buspairs[3]]]
                    if params["model"] == "opf"
                        conv_bi_x_li(_m, cyc, x_li, expairs_3, params, "cs")
                    elseif params["model"] == "ots" || params["model"] == "ots_relax"
                        ub = [bd_data["cos_max"][buspairs[1]], bd_data["cos_max"][buspairs[2]], bd_data["cos_max"][buspairs[3]], bd_data["sin_max"][buspairs[1]], bd_data["sin_max"][buspairs[2]], bd_data["sin_max"][buspairs[3]]]
                        lb = [bd_data["cos_min"][buspairs[1]], bd_data["cos_min"][buspairs[2]], bd_data["cos_min"][buspairs[3]], bd_data["sin_min"][buspairs[1]], bd_data["sin_min"][buspairs[2]], bd_data["sin_min"][buspairs[3]]]
                        ub2 = [max(ub[i], 0) for i in 1:6]
                        lb2 = [min(lb[i], 0) for i in 1:6]
                        conv_bi_x_li_on_off(_m, cyc, x_li, ub, lb, ub2, lb2, expairs_3, params, "cs")
                    end
                end
            end
            if params["cycle_wr_wi_cuts"]
                @constraints(_m, begin
                [cyc in cycle["cyc_3"]], _m[:x_w][cyc, (8,3)] == _m[:x_w][cyc, (1,2)] - _m[:x_w][cyc, (4,5)]
                [cyc in cycle["cyc_3"]], _m[:x_w][cyc, (8,6)] == _m[:x_w][cyc, (1,5)] + _m[:x_w][cyc, (2,4)]
                [cyc in cycle["cyc_3"]], _m[:x_w][cyc, (9,1)] == _m[:x_w][cyc, (2,3)] + _m[:x_w][cyc, (5,6)]
                [cyc in cycle["cyc_3"]], _m[:x_w][cyc, (9,4)] == _m[:x_w][cyc, (2,6)] - _m[:x_w][cyc, (3,5)]
                [cyc in cycle["cyc_3"]], _m[:x_w][cyc, (7,2)] == _m[:x_w][cyc, (1,3)] + _m[:x_w][cyc, (4,6)]
                [cyc in cycle["cyc_3"]], _m[:x_w][cyc, (7,5)] == _m[:x_w][cyc, (1,6)] - _m[:x_w][cyc, (3,4)]
                end)
                if params["model"] == "opf"
                    JuMP.@constraint(_m, [cyc in cycle["cyc_3"]], sum(_m[:λ_w][cyc, i] for i in 1:(2^9)) == 1)
                elseif params["model"] == "ots" || params["model"] == "ots_relax"
                    JuMP.@constraint(_m, [cyc in cycle["cyc_3"]], sum(_m[:λ_w][cyc, i] for i in 1:(2^9)) == _m[:zc][Tuple(cyc)])
                end
                for cyc in cycle["cyc_3"]
                    buspairs = [(cyc[1], cyc[2]), (cyc[2], cyc[3]), (cyc[1], cyc[3])]
                    x_li = [_m[:wr][buspairs[1]], _m[:wr][buspairs[2]], _m[:wr][buspairs[3]], _m[:wi][buspairs[1]], _m[:wi][buspairs[2]], _m[:wi][buspairs[3]], _m[:w][cyc[1]], _m[:w][cyc[2]], _m[:w][cyc[3]]]
                    if params["model"] == "opf"
                        conv_bi_x_li(_m, cyc, x_li, expairs_w_3, params, "w")
                    elseif params["model"] == "ots" || params["model"] == "ots_relax"
                        ub = [bd_data["wr_max"][buspairs[1]], bd_data["wr_max"][buspairs[2]], bd_data["wr_max"][buspairs[3]], bd_data["wi_max"][buspairs[1]], bd_data["wi_max"][buspairs[2]], bd_data["wi_max"][buspairs[3]], bd_data["w_max"][cyc[1]], bd_data["w_max"][cyc[2]], bd_data["w_max"][cyc[3]]]
                        lb = [bd_data["wr_min"][buspairs[1]], bd_data["wr_min"][buspairs[2]], bd_data["wr_min"][buspairs[3]], bd_data["wi_min"][buspairs[1]], bd_data["wi_min"][buspairs[2]], bd_data["wi_min"][buspairs[3]], bd_data["w_min"][cyc[1]], bd_data["w_min"][cyc[2]], bd_data["w_min"][cyc[3]]]
                        ub2 = [max(ub[i], 0) for i in 1:6]
                        lb2 = [min(lb[i], 0) for i in 1:6]
                        append!(ub2, ub[7:end])
                        append!(lb2, lb[7:end])
                        conv_bi_x_li_on_off(_m, cyc, x_li, ub, lb, ub2, lb2, expairs_w_3, params, "w")
                    end
                end
            end
        end
        if params["cycle_max_bnd"] >= 4
            expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
            expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
            if params["rotate_4w"]
                append!(expairs_w_4, [(1,2, 12), (5,6, 12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)])
            end
            if params["cycle_c_s_cuts"]
                @constraints(_m, begin
                    [cyc in cycle["cyc_4"]], _m[:x_c][cyc, (1,3)] - _m[:x_c][cyc, (5,7)] == _m[:x_c][cyc, (2,4)] + _m[:x_c][cyc, (6,8)]
                    [cyc in cycle["cyc_4"]], _m[:x_c][cyc, (1,7)] + _m[:x_c][cyc, (3,5)] == -_m[:x_c][cyc, (4,6)] + _m[:x_c][cyc, (2,8)]
                    [cyc in cycle["cyc_4"]], _m[:x_c][cyc, (1,2)] - _m[:x_c][cyc, (5,6)] == _m[:x_c][cyc, (3,4)] + _m[:x_c][cyc, (7,8)]
                    [cyc in cycle["cyc_4"]], _m[:x_c][cyc, (1,6)] + _m[:x_c][cyc, (2,5)] == -_m[:x_c][cyc, (4,7)] + _m[:x_c][cyc, (3,8)]
                    [cyc in cycle["cyc_4"]], _m[:x_c][cyc, (2,3)] - _m[:x_c][cyc, (6,7)] == _m[:x_c][cyc, (1,4)] + _m[:x_c][cyc, (5,8)]
                    [cyc in cycle["cyc_4"]], _m[:x_c][cyc, (2,7)] + _m[:x_c][cyc, (3,6)] == -_m[:x_c][cyc, (4,5)] + _m[:x_c][cyc, (1,8)]
                end)
                if params["model"] == "opf"
                    JuMP.@constraint(_m, [cyc in cycle["cyc_4"]], sum(_m[:λ_c][cyc, i] for i in 1:(2^8)) == 1)
                elseif params["model"] == "ots" || params["model"] == "ots_relax"
                    JuMP.@constraint(_m, [cyc in cycle["cyc_4"]], sum(_m[:λ_c][cyc, i] for i in 1:(2^8)) == _m[:zc][Tuple(cyc)])
                end
                for cyc in cycle["cyc_4"]
                    buspairs = [(cyc[1], cyc[2]), (cyc[2], cyc[3]), (cyc[3], cyc[4]), (cyc[1], cyc[4])]
                    x_li = [_m[:cs][buspairs[1]], _m[:cs][buspairs[2]], _m[:cs][buspairs[3]], _m[:cs][buspairs[4]], _m[:si][buspairs[1]], _m[:si][buspairs[2]], _m[:si][buspairs[3]], _m[:si][buspairs[4]]]
                    if params["model"] == "opf"
                        conv_bi_x_li(_m, cyc, x_li, expairs_4, params, "cs")
                    elseif params["model"] == "ots" || params["model"] == "ots_relax"
                        ub = [bd_data["cos_max"][buspairs[1]], bd_data["cos_max"][buspairs[2]], bd_data["cos_max"][buspairs[3]], bd_data["cos_max"][buspairs[4]], bd_data["sin_max"][buspairs[1]], bd_data["sin_max"][buspairs[2]], bd_data["sin_max"][buspairs[3]], bd_data["sin_max"][buspairs[4]]]
                        lb = [bd_data["cos_min"][buspairs[1]], bd_data["cos_min"][buspairs[2]], bd_data["cos_min"][buspairs[3]], bd_data["cos_min"][buspairs[4]], bd_data["sin_min"][buspairs[1]], bd_data["sin_min"][buspairs[2]], bd_data["sin_min"][buspairs[3]], bd_data["sin_min"][buspairs[4]]]
                        ub2 = [max(ub[i], 0) for i in 1:length(ub)]
                        lb2 = [min(lb[i], 0) for i in 1:length(lb)]
                        conv_bi_x_li_on_off(_m, cyc, x_li, ub, lb, ub2, lb2, expairs_4, params, "cs")
                    end
                end
            end
            if params["cycle_wr_wi_cuts"]
                @constraints(_m, begin
                    [cyc in cycle["cyc_4"]], _m[:x_w][cyc, (1,3)] - _m[:x_w][cyc, (5,7)] == _m[:x_w][cyc, (2,4)] + _m[:x_w][cyc, (6,8)]
                    [cyc in cycle["cyc_4"]], _m[:x_w][cyc, (1,7)] + _m[:x_w][cyc, (3,5)] == -_m[:x_w][cyc, (4,6)] + _m[:x_w][cyc, (2,8)]
                end)
                if params["rotate_4w"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_4"]], _m[:x_w][cyc, (1,2,12)] - _m[:x_w][cyc, (5,6,12)] == _m[:x_w][cyc, (3,4,10)] + _m[:x_w][cyc, (7,8,10)]
                        [cyc in cycle["cyc_4"]], _m[:x_w][cyc, (1,6,12)] + _m[:x_w][cyc, (2,5,12)] == -_m[:x_w][cyc, (4,7,10)] + _m[:x_w][cyc, (3,8,10)]
                        [cyc in cycle["cyc_4"]], _m[:x_w][cyc, (2,3,9)] - _m[:x_w][cyc, (6,7,9)] == _m[:x_w][cyc, (1,4,11)] + _m[:x_w][cyc, (5,8,11)]
                        [cyc in cycle["cyc_4"]], _m[:x_w][cyc, (2,7,9)] + _m[:x_w][cyc, (3,6,9)] == -_m[:x_w][cyc, (4,5,11)] + _m[:x_w][cyc, (1,8,11)]
                    end)
                    if params["model"] == "opf"
                        JuMP.@constraint(_m, [cyc in cycle["cyc_4"]], sum(_m[:λ_w][cyc, i] for i in 1:(2^12)) == 1)
                    elseif params["model"] == "ots" || params["model"] == "ots_relax"
                        JuMP.@constraint(_m, [cyc in cycle["cyc_4"]], sum(_m[:λ_w][cyc, i] for i in 1:(2^12)) == _m[:zc][Tuple(cyc)])
                    end
                else
                    if params["model"] == "opf"
                        JuMP.@constraint(_m, [cyc in cycle["cyc_4"]], sum(_m[:λ_w][cyc, i] for i in 1:(2^8)) == 1)
                    elseif params["model"] == "ots" || params["model"] == "ots_relax"
                        JuMP.@constraint(_m, [cyc in cycle["cyc_4"]], sum(_m[:λ_w][cyc, i] for i in 1:(2^8)) == _m[:zc][Tuple(cyc)])
                    end
                end
                for cyc in cycle["cyc_4"]
                    buspairs = [(cyc[1], cyc[2]), (cyc[2], cyc[3]), (cyc[3], cyc[4]), (cyc[1], cyc[4])]
                    x_li = [_m[:wr][buspairs[1]], _m[:wr][buspairs[2]], _m[:wr][buspairs[3]], _m[:wr][buspairs[4]], _m[:wi][buspairs[1]], _m[:wi][buspairs[2]], _m[:wi][buspairs[3]], _m[:wi][buspairs[4]]]
                    if params["rotate_4w"]
                        append!(x_li, [_m[:w][cyc[1]], _m[:w][cyc[2]], _m[:w][cyc[3]], _m[:w][cyc[4]]])
                    end
                    if params["model"] == "opf"
                        conv_bi_x_li(_m, cyc, x_li, expairs_w_4, params, "w")
                    elseif params["model"] == "ots" || params["model"] == "ots_relax"
                        ub = [bd_data["wr_max"][buspairs[1]], bd_data["wr_max"][buspairs[2]], bd_data["wr_max"][buspairs[3]], bd_data["wr_max"][buspairs[4]], bd_data["wi_max"][buspairs[1]], bd_data["wi_max"][buspairs[2]], bd_data["wi_max"][buspairs[3]], bd_data["wi_max"][buspairs[4]]]
                        lb = [bd_data["wr_min"][buspairs[1]], bd_data["wr_min"][buspairs[2]], bd_data["wr_min"][buspairs[3]], bd_data["wr_min"][buspairs[4]], bd_data["wi_min"][buspairs[1]], bd_data["wi_min"][buspairs[2]], bd_data["wi_min"][buspairs[3]], bd_data["wi_min"][buspairs[4]]]
                        ub2 = [max(ub[i], 0) for i in 1:length(ub)]
                        lb2 = [min(lb[i], 0) for i in 1:length(lb)]
                        if params["rotate_4w"]
                            temp_li1 = [bd_data["w_max"][cyc[1]], bd_data["w_max"][cyc[2]], bd_data["w_max"][cyc[3]], bd_data["w_max"][cyc[4]]]
                            temp_li2 = [bd_data["w_min"][cyc[1]], bd_data["w_min"][cyc[2]], bd_data["w_min"][cyc[3]], bd_data["w_min"][cyc[4]]]
                            append!(ub, temp_li1)
                            append!(lb, temp_li2)
                            append!(ub2, temp_li1)
                            append!(ub2, temp_li2)
                        end
                        conv_bi_x_li_on_off(_m, cyc, x_li, ub, lb, ub2, lb2, expairs_w_4, params, "w")
                    end
                end
            end
        end
    elseif params["cycle_relax"] == "none"
        if params["cycle_c_s_cuts"]
            # Original bilinear constraints
            if params["cycle_max_bnd"] >= 3
                # td_bar = [-0.03196329287072686,  0.32706963978247594, 0.2951063469117491]
                arcpairs = [(1,2), (2,3), (1,3)]
                @constraints(_m, begin
                    [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[3])] == _m[:cs][(cyc[1], cyc[2])] * _m[:cs][(cyc[2], cyc[3])] - _m[:si][(cyc[1], cyc[2])] * _m[:si][(cyc[2], cyc[3])]
                    [cyc in cycle["cyc_3"]], _m[:si][(cyc[1], cyc[3])] == _m[:cs][(cyc[1], cyc[2])] * _m[:si][(cyc[2], cyc[3])] + _m[:si][(cyc[1], cyc[2])] * _m[:cs][(cyc[2], cyc[3])]
                    # [ap in arcpairs], td[ap] == td_bar[ap]
                    # [ap in arcpairs], _m[:cs][ap] == cs_bar[ap]
                    # [ap in arcpairs], _m[:si][ap] == si_bar[ap]
                end)
                # td_bar2 = [-0.029050451331921766, 0.3261192251358423, 0.29706877380392055]
                # arcpairs = []
                # arcpairs2 = []
                # for cyc in cycle["cyc_3"]
                #     # For cos
                #     push!(arcpairs, (cyc[2], cyc[3], cyc[1], cyc[3]))
                #     push!(arcpairs, (cyc[1], cyc[2], cyc[2], cyc[3]))
                #     push!(arcpairs, (cyc[1], cyc[2], cyc[1], cyc[3]))
                #     push!(arcpairs2, (cyc[1], cyc[2], cyc[2], cyc[3]))
                #     # For sin
                #     # push!(arcpairs2, (cyc[2], cyc[3], cyc[1], cyc[3]))
                #     # push!(arcpairs2, (cyc[1], cyc[2], cyc[1], cyc[3]))
                # end
                # JuMP.@variable(_m, hcc[ap in arcpairs2])
                # JuMP.@variable(_m, hss[ap in arcpairs])
                # @constraints(_m, begin
                #     [cyc in cycle["cyc_3"], ap in arcpairs2], _m[:hcc][(ap[1], ap[2], ap[3], ap[4])] == _m[:cs][(ap[1], ap[2])] * _m[:cs][(ap[3], ap[4])]
                #     [cyc in cycle["cyc_3"], ap in arcpairs], _m[:hss][(ap[1], ap[2], ap[3], ap[4])] == _m[:si][(ap[1], ap[2])] * _m[:si][(ap[3], ap[4])]
                #     [cyc in cycle["cyc_3"]], _m[:hcc][(cyc[1], cyc[2], cyc[2], cyc[3])] * _m[:cs][(cyc[1], cyc[3])] + _m[:cs][(cyc[1], cyc[2])] * _m[:hss][(cyc[2], cyc[3], cyc[1], cyc[3])] - _m[:hss][(cyc[1], cyc[2], cyc[2], cyc[3])] * _m[:cs][(cyc[1], cyc[3])] + _m[:hss][(cyc[1], cyc[2], cyc[1], cyc[3])] * _m[:cs][(cyc[2], cyc[3])] == 1
                #     # td[(1,2)] == td_bar2[(1,2)]
                #     # td[(2,3)] == td_bar2[(2,3)]
                #     # td[(1,3)] == td_bar2[(1,3)]
                #     # [cyc in cycle["cyc_3"]], _m[:hcc][(cyc[2], cyc[3], cyc[1], cyc[3])] * _m[:si][(cyc[1], cyc[2])] + _m[:si][(cyc[1], cyc[3])] * _m[:hss][(cyc[1], cyc[2], cyc[2], cyc[3])] + _m[:hcc][(cyc[1], cyc[2], cyc[1], cyc[3])] * _m[:si][(cyc[2], cyc[3])] - _m[:hcc][(cyc[1], cyc[2], cyc[2], cyc[3])] * _m[:si][(cyc[1], cyc[3])] == 0
                # end)
                # @NLconstraint(_m, [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[2])] * _m[:cs][(cyc[2], cyc[3])] * _m[:cs][(cyc[1], cyc[3])] + _m[:cs][(cyc[1], cyc[2])] * _m[:si][(cyc[2], cyc[3])] * _m[:si][(cyc[1], cyc[3])] - _m[:si][(cyc[1], cyc[2])] * _m[:si][(cyc[2], cyc[3])] * _m[:cs][(cyc[1], cyc[3])] + _m[:si][(cyc[1], cyc[2])] * _m[:cs][(cyc[2], cyc[3])] * _m[:si][(cyc[1], cyc[3])] == 1)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[1], cyc[2])] == _m[:cs][(cyc[2], cyc[3])] * _m[:cs][(cyc[1], cyc[3])] + _m[:si][(cyc[2], cyc[3])] * _m[:si][(cyc[1], cyc[3])]
                        [cyc in cycle["cyc_3"]], - _m[:si][(cyc[1], cyc[2])] == - _m[:cs][(cyc[2], cyc[3])] * _m[:si][(cyc[1], cyc[3])] + _m[:si][(cyc[2], cyc[3])] * _m[:cs][(cyc[1], cyc[3])]
                        [cyc in cycle["cyc_3"]], _m[:cs][(cyc[2], cyc[3])] == _m[:cs][(cyc[1], cyc[3])] * _m[:cs][(cyc[1], cyc[2])] + _m[:si][(cyc[1], cyc[3])] * _m[:si][(cyc[1], cyc[2])]
                        [cyc in cycle["cyc_3"]], - _m[:si][(cyc[2], cyc[3])] == _m[:cs][(cyc[1], cyc[3])] * _m[:si][(cyc[1], cyc[2])] - _m[:si][(cyc[1], cyc[3])] * _m[:cs][(cyc[1], cyc[2])]
                    end)
                end
            end
            if params["cycle_max_bnd"] >= 4
                @constraints(_m, begin
                    [cyc in cycle["cyc_4"]], _m[:cs][(cyc[1], cyc[2])] * _m[:cs][(cyc[3], cyc[4])] - _m[:si][(cyc[1], cyc[2])] * _m[:si][(cyc[3], cyc[4])] == _m[:cs][(cyc[1], cyc[4])] * _m[:cs][(cyc[2], cyc[3])] + _m[:si][(cyc[1], cyc[4])] * _m[:si][(cyc[2], cyc[3])]
                    [cyc in cycle["cyc_4"]], _m[:cs][(cyc[1], cyc[2])] * _m[:si][(cyc[3], cyc[4])] + _m[:si][(cyc[1], cyc[2])] * _m[:cs][(cyc[3], cyc[4])] == - _m[:cs][(cyc[1], cyc[4])] * _m[:si][(cyc[2], cyc[3])] + _m[:si][(cyc[1], cyc[4])] * _m[:cs][(cyc[2], cyc[3])]
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_4"]], _m[:cs][(cyc[1], cyc[2])] * _m[:cs][(cyc[2], cyc[3])] - _m[:si][(cyc[1], cyc[2])] * _m[:si][(cyc[2], cyc[3])] == _m[:cs][(cyc[1], cyc[4])] * _m[:cs][(cyc[3], cyc[4])] + _m[:si][(cyc[1], cyc[4])] * _m[:si][(cyc[3], cyc[4])]
                        [cyc in cycle["cyc_4"]], _m[:cs][(cyc[1], cyc[2])] * _m[:si][(cyc[2], cyc[3])] + _m[:si][(cyc[1], cyc[2])] * _m[:cs][(cyc[2], cyc[3])] == - _m[:cs][(cyc[1], cyc[4])] * _m[:si][(cyc[3], cyc[4])] + _m[:si][(cyc[1], cyc[4])] * _m[:cs][(cyc[3], cyc[4])]
                        [cyc in cycle["cyc_4"]], _m[:cs][(cyc[2], cyc[3])] * _m[:cs][(cyc[3], cyc[4])] - _m[:si][(cyc[2], cyc[3])] * _m[:si][(cyc[3], cyc[4])] == _m[:cs][(cyc[1], cyc[4])] * _m[:cs][(cyc[1], cyc[2])] + _m[:si][(cyc[1], cyc[4])] * _m[:si][(cyc[1], cyc[2])]
                        [cyc in cycle["cyc_4"]], _m[:cs][(cyc[2], cyc[3])] * _m[:si][(cyc[3], cyc[4])] + _m[:si][(cyc[2], cyc[3])] * _m[:cs][(cyc[3], cyc[4])] == - _m[:cs][(cyc[1], cyc[4])] * _m[:si][(cyc[1], cyc[2])] + _m[:si][(cyc[1], cyc[4])] * _m[:cs][(cyc[1], cyc[2])]
                    end)
                end
            end
        end

        if params["cycle_wr_wi_cuts"]
            # Original bilinear constraints w version
            if params["cycle_max_bnd"] >= 3
                # @constraints(_m, begin
                #     [cyc in cycle["cyc_3"]], _m[:w][cyc[2]] * _m[:wr][(cyc[1], cyc[3])] == _m[:wr][(cyc[1], cyc[2])] * _m[:wr][(cyc[2], cyc[3])] - _m[:wi][(cyc[1], cyc[2])] * _m[:wi][(cyc[2], cyc[3])]
                #     [cyc in cycle["cyc_3"]], _m[:w][cyc[2]] * _m[:wi][(cyc[1], cyc[3])] == _m[:wr][(cyc[1], cyc[2])] * _m[:wi][(cyc[2], cyc[3])] + _m[:wi][(cyc[1], cyc[2])] * _m[:wr][(cyc[2], cyc[3])]
                # end)
                arcpairs = []
                arcpairs2 = []
                arcpairs3 = []
                for cyc in cycle["cyc_3"]
                    # For cos
                    push!(arcpairs, (cyc[2], cyc[3], cyc[1], cyc[3]))
                    push!(arcpairs, (cyc[1], cyc[2], cyc[2], cyc[3]))
                    push!(arcpairs, (cyc[1], cyc[2], cyc[1], cyc[3]))
                    push!(arcpairs2, (cyc[1], cyc[2], cyc[2], cyc[3]))
                    push!(arcpairs3, (cyc[1], cyc[2]))
                    # For sin
                    # push!(arcpairs2, (cyc[2], cyc[3], cyc[1], cyc[3]))
                    # push!(arcpairs2, (cyc[1], cyc[2], cyc[1], cyc[3]))
                end
                JuMP.@variable(_m, wwr[ap in arcpairs2])
                JuMP.@variable(_m, wwi[ap in arcpairs])
                JuMP.@variable(_m, ww[ap in arcpairs3])
                @constraints(_m, begin
                    [cyc in cycle["cyc_3"], ap in arcpairs2], _m[:wwr][(ap[1], ap[2], ap[3], ap[4])] == _m[:wr][(ap[1], ap[2])] * _m[:wr][(ap[3], ap[4])]
                    [cyc in cycle["cyc_3"], ap in arcpairs], _m[:wwi][(ap[1], ap[2], ap[3], ap[4])] == _m[:wi][(ap[1], ap[2])] * _m[:wi][(ap[3], ap[4])]
                    [cyc in cycle["cyc_3"], ap in arcpairs3], _m[:ww][(ap[1], ap[2])] == _m[:w][ap[1]] * _m[:w][ap[2]]
                    [cyc in cycle["cyc_3"]], _m[:wwr][(cyc[1], cyc[2], cyc[2], cyc[3])] * _m[:wr][(cyc[1], cyc[3])] + _m[:wr][(cyc[1], cyc[2])] * _m[:wwi][(cyc[2], cyc[3], cyc[1], cyc[3])] - _m[:wwi][(cyc[1], cyc[2], cyc[2], cyc[3])] * _m[:wr][(cyc[1], cyc[3])] + _m[:wwi][(cyc[1], cyc[2], cyc[1], cyc[3])] * _m[:wr][(cyc[2], cyc[3])] == _m[:ww][(cyc[1], cyc[2])] * _m[:w][cyc[3]]
                    # [cyc in cycle["cyc_3"]], _m[:wwr][(cyc[2], cyc[3], cyc[1], cyc[3])] * _m[:wi][(cyc[1], cyc[2])] + _m[:wi][(cyc[1], cyc[3])] * _m[:wwi][(cyc[1], cyc[2], cyc[2], cyc[3])] + _m[:wwr][(cyc[1], cyc[2], cyc[1], cyc[3])] * _m[:wi][(cyc[2], cyc[3])] - _m[:wwr][(cyc[1], cyc[2], cyc[2], cyc[3])] * _m[:wi][(cyc[1], cyc[3])] == 0 # sin
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        [cyc in cycle["cyc_3"]], _m[:w][cyc[3]] * _m[:wr][(cyc[1], cyc[2])] == _m[:wr][(cyc[2], cyc[3])] * _m[:wr][(cyc[1], cyc[3])] + _m[:wi][(cyc[2], cyc[3])] * _m[:wi][(cyc[1], cyc[3])]
                        [cyc in cycle["cyc_3"]], - _m[:w][cyc[3]] * _m[:wi][(cyc[1], cyc[2])] == - _m[:wr][(cyc[2], cyc[3])] * _m[:wi][(cyc[1], cyc[3])] + _m[:wi][(cyc[2], cyc[3])] * _m[:wr][(cyc[1], cyc[3])]
                        [cyc in cycle["cyc_3"]], _m[:w][cyc[1]] * _m[:wr][(cyc[2], cyc[3])] == _m[:wr][(cyc[1], cyc[3])] * _m[:wr][(cyc[1], cyc[2])] + _m[:wi][(cyc[1], cyc[3])] * _m[:wi][(cyc[1], cyc[2])]
                        [cyc in cycle["cyc_3"]], - _m[:w][cyc[1]] * _m[:wi][(cyc[2], cyc[3])] == _m[:wr][(cyc[1], cyc[3])] * _m[:wi][(cyc[1], cyc[2])] - _m[:wi][(cyc[1], cyc[3])] * _m[:wr][(cyc[1], cyc[2])]
                    end)
                end
            end
            if params["cycle_max_bnd"] >= 4
                @constraints(_m, begin
                    [cyc in cycle["cyc_4"]], _m[:wr][(cyc[1], cyc[2])] * _m[:wr][(cyc[3], cyc[4])] - _m[:wi][(cyc[1], cyc[2])] * _m[:wi][(cyc[3], cyc[4])] == _m[:wr][(cyc[1], cyc[4])] * _m[:wr][(cyc[2], cyc[3])] + _m[:wi][(cyc[1], cyc[4])] * _m[:wi][(cyc[2], cyc[3])]
                    [cyc in cycle["cyc_4"]], _m[:wr][(cyc[1], cyc[2])] * _m[:wi][(cyc[3], cyc[4])] + _m[:wi][(cyc[1], cyc[2])] * _m[:wr][(cyc[3], cyc[4])] == - _m[:wr][(cyc[1], cyc[4])] * _m[:wi][(cyc[2], cyc[3])] + _m[:wi][(cyc[1], cyc[4])] * _m[:wr][(cyc[2], cyc[3])]
                end)
                if params["rotate"]
                    @constraints(_m, begin
                        # [cyc in cycle["cyc_4"]], _m[:wr][(cyc[1], cyc[2])] * _m[:wr][(cyc[2], cyc[3])] - _m[:wi][(cyc[1], cyc[2])] * _m[:wi][(cyc[2], cyc[3])] == _m[:wr][(cyc[1], cyc[4])] * _m[:wr][(cyc[3], cyc[4])] + _m[:wi][(cyc[1], cyc[4])] * _m[:wi][(cyc[3], cyc[4])]
                        # [cyc in cycle["cyc_4"]], _m[:wr][(cyc[1], cyc[2])] * _m[:wi][(cyc[2], cyc[3])] + _m[:wi][(cyc[1], cyc[2])] * _m[:wr][(cyc[2], cyc[3])] == - _m[:wr][(cyc[1], cyc[4])] * _m[:wi][(cyc[3], cyc[4])] + _m[:wi][(cyc[1], cyc[4])] * _m[:wr][(cyc[3], cyc[4])]
                        # [cyc in cycle["cyc_4"]], _m[:wr][(cyc[2], cyc[3])] * _m[:wr][(cyc[3], cyc[4])] - _m[:wi][(cyc[2], cyc[3])] * _m[:wi][(cyc[3], cyc[4])] == _m[:wr][(cyc[1], cyc[4])] * _m[:wr][(cyc[1], cyc[2])] + _m[:wi][(cyc[1], cyc[4])] * _m[:wi][(cyc[1], cyc[2])]
                        # [cyc in cycle["cyc_4"]], _m[:wr][(cyc[2], cyc[3])] * _m[:wi][(cyc[3], cyc[4])] + _m[:wi][(cyc[2], cyc[3])] * _m[:wr][(cyc[3], cyc[4])] == - _m[:wr][(cyc[1], cyc[4])] * _m[:wi][(cyc[1], cyc[2])] + _m[:wi][(cyc[1], cyc[4])] * _m[:wr][(cyc[1], cyc[2])]
                    end)
                end
            end
        end
    end
end

function link_zc_z(_m, params)
    cycle_3_4 = []
    if params["cycle_max_bnd"] >= 3
        append!(cycle_3_4, cycle["cyc_3"])
    end
    if params["cycle_max_bnd"] >= 4
        append!(cycle_3_4, cycle["cyc_4"])
    end
    for cyc in cycle_3_4
        buspairs = cyc_to_buspair(cyc)
        cyc_tp = Tuple(cyc)
        cyc_size = length(cyc)
        @constraints(_m, begin
            _m[:zc][cyc_tp] >= 1 - sum((1 - _m[:z][bp]) for bp in buspairs) - 1e-3
            _m[:zc][cyc_tp] <= 1 / cyc_size * sum(_m[:z][bp] for bp in buspairs) + 1e-3
        end)
    end
end
