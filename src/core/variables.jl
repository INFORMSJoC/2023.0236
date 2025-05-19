function variables(_m, ref, params)

   wr_min, wr_max, wi_min, wi_max = PowerModels.ref_calc_voltage_product_bounds(ref[:buspairs])

   # computing bounds for cosine variables
   cos_min = Dict([(bp, -Inf) for bp in keys(ref[:buspairs])])
   cos_max = Dict([(bp, Inf) for bp in keys(ref[:buspairs])])

   for (bp, buspair) in ref[:buspairs]
      angmin = buspair["angmin"]
      angmax = buspair["angmax"]
      if angmin >= 0
         cos_max[bp] = cos(angmin)
         cos_min[bp] = cos(angmax)
      end
      if angmax <= 0
         cos_max[bp] = cos(angmax)
         cos_min[bp] = cos(angmin)
      end
      if angmin < 0 && angmax > 0
         cos_max[bp] = 1.0
         cos_min[bp] = min(cos(angmin), cos(angmax))
      end
   end

   # Computing bounds for sine variables.
   sin_min = Dict([(bp, -Inf) for bp in keys(ref[:buspairs])])
   sin_max = Dict([(bp, Inf) for bp in keys(ref[:buspairs])])
   for bp in keys(ref[:buspairs])
      sin_min[bp] = sin(ref[:buspairs][bp]["angmin"])
      sin_max[bp] = sin(ref[:buspairs][bp]["angmax"])
   end

   # Add bounds for variables corresponding to reversed buspairs.
   bp_li = [bp for bp in keys(ref[:buspairs])]
   new_li = [(bp[2], bp[1]) for bp in bp_li] # list of reversed buspairs.
   append!(bp_li, new_li)
   for bp in new_li
      sin_max[bp] = - sin_min[(bp[2], bp[1])]
      sin_min[bp] = - sin_max[(bp[2], bp[1])]
      wi_max[bp] = - wi_min[(bp[2], bp[1])]
      wi_min[bp] = - wi_max[(bp[2], bp[1])]
      cos_max[bp] = cos_max[(bp[2], bp[1])]
      cos_min[bp] = cos_min[(bp[2], bp[1])]
      wr_max[bp] = wr_max[(bp[2], bp[1])]
      wr_min[bp] = wr_min[(bp[2], bp[1])]
   end

   # current magnitude squared
   # compute upper bound
   cm_ub = Dict()
   for (bp, buspair) in ref[:buspairs]
      cm_ub[bp] = ((buspair["rate_a"] * buspair["tap"])/buspair["vm_fr_min"])^2
      if params["extra_rlt"]
         cm_ub[(bp[2], bp[1])] = (buspair["rate_a"] / buspair["vm_to_min"])^2
      end
   end
   # println("cm_ub: ", cm_ub)

   # Bounds for Variables
   bd_data["cm_ub"] = Dict()
   bd_data["td_u"] = Dict()
   bd_data["w_max"] = Dict()
   bd_data["w_min"] = Dict()
   bd_data["cos_min"] = cos_min
   bd_data["cos_max"] = cos_max
   bd_data["sin_min"] = sin_min
   bd_data["sin_max"] = sin_max
   bd_data["wr_min"] = wr_min
   bd_data["wr_max"] = wr_max
   bd_data["wi_min"] = wi_min
   bd_data["wi_max"] = wi_max
   for (bp, buspair) in ref[:buspairs]
      bd_data["cm_ub"][bp] = ((buspair["rate_a"] * buspair["tap"]) / buspair["vm_fr_min"])^2
      if params["extra_rlt"]
         bd_data["cm_ub"][(bp[2], bp[1])] = (buspair["rate_a"] /buspair["vm_to_min"])^2
      end

      bd_data["td_u"][bp] = max(abs(ref[:buspairs][bp]["angmin"]), abs(ref[:buspairs][bp]["angmax"]))
   end
   for i in keys(ref[:bus])
      bd_data["w_max"][i] = ref[:bus][i]["vmax"]^2
      bd_data["w_min"][i] = ref[:bus][i]["vmin"]^2
   end
   # bd_data["td_M"] = sum(bd_data["td_u"][bp] for bp in keys(ref[:buspairs]))
   td_u_sorted = []
   for bp in keys(ref[:buspairs])
      append!(td_u_sorted, bd_data["td_u"][bp])
   end
   sort!(td_u_sorted, rev = true)
   # bd_data["td_M"] = sum(td_u_sorted[i] for i in 1:(length(ref[:bus]) - 1))
   bd_data["td_M"] = sum(td_u_sorted[i] for i in 1:length(td_u_sorted))

   JuMP.@variables(_m, begin
      # voltage angle and magnitude
      va[i in keys(ref[:bus])]

      # generation pg and qg
      ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"]
      ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"]
   end)
   # println("bp_li: ", bp_li)
   if params["model"] == "opf"
      JuMP.@variables(_m, begin
         # voltage angle and magnitude
         ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]

         # voltage squared
         ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2

         # voltage product - lifted vars, wr = vi*vj*cs, wi = vi*vj*si
         wr_min[bp] <= wr[bp in bp_li] <= wr_max[bp]
         wi_min[bp] <= wi[bp in bp_li] <= wi_max[bp]
         # Cosine variables
         cos_min[bp] <= cs[bp in bp_li] <= cos_max[bp]
         # cos_min[bp] <= cs[bp in keys(ref[:buspairs])] <= cos_max[bp]
         # Sine variables
         sin_min[bp] <= si[bp in bp_li] <=  sin_max[bp]
         # sin(ref[:buspairs][bp]["angmin"]) <= si[bp in keys(ref[:buspairs])] <=  sin(ref[:buspairs][bp]["angmax"])

         # line flow variables
         -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"]
         -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"]

         # voltage-angle differences
         ref[:buspairs][bp]["angmin"] <= td[bp in keys(ref[:buspairs])]  <= ref[:buspairs][bp]["angmax"]

         # Current magnitude product variable
         # 0 <= cm[bp in keys(ref[:buspairs])] <= cm_ub[bp]
         # 0 <= cm[bp in bp_li] <= cm_ub[bp]
      end)
      if params["extra_rlt"]
         JuMP.@variable(_m, 0 <= cm[bp in bp_li] <= cm_ub[bp])
      else
         JuMP.@variable(_m, 0 <= cm[bp in keys(ref[:buspairs])] <= cm_ub[bp])
      end
   elseif params["model"] == "ots" || params["model"] == "ots_relax"
      JuMP.@variables(_m, begin
         ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]
         # vm[i in keys(ref[:bus])]

         # voltage squared
         ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2
         # w[i in keys(ref[:bus])]

         wr[bp in bp_li]
         wi[bp in bp_li]
         cs[bp in bp_li]
         si[bp in bp_li]

         # line flow variables
         p[(l,i,j) in ref[:arcs]]
         q[(l,i,j) in ref[:arcs]]

         # voltage-angle differences
         td[bp in keys(ref[:buspairs])]

         # Current magnitude product variable
         # cm[bp in keys(ref[:buspairs])]
         # cm[bp in bp_li]

         wz[bp in bp_li]
      end)
      if params["extra_rlt"]
         JuMP.@variable(_m, cm[bp in bp_li])
      else
         JuMP.@variable(_m, cm[bp in keys(ref[:buspairs])])
      end
      # Line On-Off variables
      if params["model"] == "ots" && params["turn_on_lines"] == false
         JuMP.@variable(_m, z[bp in keys(ref[:buspairs])], Bin)
      elseif params["model"] == "ots" && params["turn_on_lines"] == true
         JuMP.@variable(_m, z[bp in keys(ref[:buspairs])] == 1)
      elseif params["model"] == "ots_relax"
         JuMP.@variable(_m, 0 <= z[bp in keys(ref[:buspairs])] <= 1)
      end
   end

   if params["relax"] == "rmc"
      # VV variables
      if params["model"] == "opf"
         JuMP.@variable(_m, ref[:buspairs][bp]["vm_fr_min"] * ref[:buspairs][bp]["vm_to_min"] <= vv[bp in keys(ref[:buspairs])] <= ref[:buspairs][bp]["vm_fr_max"] * ref[:buspairs][bp]["vm_to_max"])
      elseif params["model"] == "ots" || params["model"] == "ots_relax"
         JuMP.@variable(_m, vv[bp in keys(ref[:buspairs])])
         JuMP.@variable(_m, vmz[bp in bp_li])
         # JuMP.@variable(_m, vvz[bp in keys(ref[:buspairs])])
         # JuMP.@variable(_m, csz[bp in bp_li])
         # JuMP.@variable(_m, siz[bp in bp_li])
      end
   elseif params["relax"] == "tri"
      # Non-neg lambda multipliers
      JuMP.@variable(_m, 0 <= λ_wr[bp in keys(ref[:buspairs]), 1:8] <= 1)
      JuMP.@variable(_m, 0 <= λ_wi[bp in keys(ref[:buspairs]), 1:8] <= 1)
      # if params["model"] == "ots" || params["model"] == "ots_relax"
      #    JuMP.@variable(_m, 0 <= λ_hat[(i,j) in keys(ref[:buspairs]), 1:2, 1:2] <= 1)
      # end
      if params["rlt_cuts"]

         JuMP.@variable(_m, wcc[bp in keys(ref[:buspairs])])
         JuMP.@variable(_m, wss[bp in keys(ref[:buspairs])])
         #JuMP.@variable(_m, c_λ_wr[bp in keys(ref[:buspairs]), 1:8])
         #JuMP.@variable(_m, s_λ_wi[bp in keys(ref[:buspairs]), 1:8])

         JuMP.@variable(_m, vc_i[bp in keys(ref[:buspairs])])
         JuMP.@variable(_m, vc_j[bp in keys(ref[:buspairs])])
         JuMP.@variable(_m, vs_i[bp in keys(ref[:buspairs])])
         JuMP.@variable(_m, vs_j[bp in keys(ref[:buspairs])])
         JuMP.@variable(_m, c_λ_wr[bp in keys(ref[:buspairs]), 1:8])
         JuMP.@variable(_m, s_λ_wi[bp in keys(ref[:buspairs]), 1:8])
      end
   end

   # # Lifted variables for OTS formulation
   # if (params["model"] == "ots") || (params["model"] == "ots_relax")
   #    JuMP.@variable(_m, 0 <= zw_fr[bp in keys(ref[:buspairs])] <= ref[:buspairs][bp]["vm_fr_max"]^2)
   #    JuMP.@variable(_m, 0 <= zw_to[bp in keys(ref[:buspairs])] <= ref[:buspairs][bp]["vm_to_max"]^2)
   # end
   cycle_3_4 = []
   if ((params["model"] == "ots") || (params["model"] == "ots_relax")) && params["cycle_cuts"]
      # cycle_3_4 = []
      append!(cycle_3_4, cycle["cyc_3"])
      if params["cycle_max_bnd"] >= 4
         append!(cycle_3_4, cycle["cyc_4"])
      end
      for i in 1:length(cycle_3_4)
         cycle_3_4[i] = Tuple(cycle_3_4[i])
      end

      if params["model"] == "ots"
         JuMP.@variable(_m, zc[cyc in cycle_3_4], Bin)
      elseif params["model"] == "ots_relax"
         JuMP.@variable(_m, 0 <= zc[cyc in cycle_3_4] <= 1)
      end
   end

   vars = Dict{String, Any}(
                            "va" => va,
                            "vm" => vm,
                            "w" => w,
                            "wr" => wr,
                            "wi" => wi,
                            "pg" => pg,
                            "qg" => qg,
                            "p" => p,
                            "q" => q,
                            "td" => td,
                            "cs" => cs,
                            "si" => si,
                            "cm" => cm
                           )

   if params["relax"] == "rmc"
      vars["vv"] = vv
      if params["model"] == "ots" || params["model"] == "ots_relax"
         vars["vmz"] = vmz
         # vars["vvz"] = vvz
         # vars["csz"] = csz
         # vars["siz"] = siz
      end
   elseif params["relax"] == "tri"
      vars["λ_wr"] = λ_wr
      vars["λ_wi"] = λ_wi
      # if (params["model"] == "ots") || (params["model"] == "ots_relax")
      #    # vars["λ_hat"] = λ_hat
      # end
      if params["rlt_cuts"]
         vars["wcc"] = wcc
         vars["wss"] = wss
         vars["vc_i"] = vc_i
         vars["vc_j"] = vc_j
         vars["vs_i"] = vs_i
         vars["vs_j"] = vs_j
         vars["c_λ_wr"] = c_λ_wr
         vars["s_λ_wi"] = s_λ_wi
      end
   end
   if (params["model"] == "ots") || (params["model"] == "ots_relax")
   #    vars["zw_fr"] = zw_fr
   #    vars["zw_to"] = zw_to
      vars["z"] = z
      vars["wz"] = wz
      if params["cycle_cuts"]
         vars["zc"] = zc
      end
   end

   if params["cycle_cuts"] && params["separation"] == false
      arcpairs = []
      busarcpairs = []
      busarctuples = [] # For 4-cycle constraints with trilinear terms.
      expairs_3 = []
      expairs_w_3 = []
      expairs_4 = []
      expairs_w_4 = []
      if params["cycle_max_bnd"] >= 3
         for cyc in cycle["cyc_3"]
            push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
            push!(arcpairs, ((cyc[2], cyc[3]), (cyc[1], cyc[3])))
            push!(arcpairs, ((cyc[1], cyc[3]), (cyc[1], cyc[2])))
            push!(busarcpairs, (cyc[2], cyc[1], cyc[3]))
            push!(busarcpairs, (cyc[3], cyc[1], cyc[2]))
            push!(busarcpairs, (cyc[1], cyc[2], cyc[3]))
         end
         expairs_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6)] # Pairs of extreme points for lambda formulation, 3-cycle, cs.
         expairs_w_3 = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6), (7,2), (7,5), (8,3), (8,6), (9,1), (9,4)]
      end
      if params["cycle_max_bnd"] >= 4
         for cyc in cycle["cyc_4"]
            if (((cyc[1], cyc[2]), (cyc[3], cyc[4])) in arcpairs) == false
               push!(arcpairs, ((cyc[1], cyc[2]), (cyc[3], cyc[4])))
            end
            if (((cyc[1], cyc[4]), (cyc[2], cyc[3])) in arcpairs) == false
               push!(arcpairs, ((cyc[1], cyc[4]), (cyc[2], cyc[3])))
            end
            if (((cyc[1], cyc[2]), (cyc[2], cyc[3])) in arcpairs) == false
               push!(arcpairs, ((cyc[1], cyc[2]), (cyc[2], cyc[3])))
            end
            if (((cyc[1], cyc[4]), (cyc[3], cyc[4])) in arcpairs) == false
               push!(arcpairs, ((cyc[1], cyc[4]), (cyc[3], cyc[4])))
            end
            if (((cyc[2], cyc[3]), (cyc[3], cyc[4])) in arcpairs) == false
               push!(arcpairs, ((cyc[2], cyc[3]), (cyc[3], cyc[4])))
            end
            if (((cyc[1], cyc[4]), (cyc[1], cyc[2])) in arcpairs) == false
               push!(arcpairs, ((cyc[1], cyc[4]), (cyc[1], cyc[2])))
            end
            push!(busarctuples, (cyc[4], (cyc[1], cyc[2]), (cyc[2], cyc[3])))
            push!(busarctuples, (cyc[2], (cyc[1], cyc[4]), (cyc[3], cyc[4])))
            push!(busarctuples, (cyc[1], (cyc[2], cyc[3]), (cyc[3], cyc[4])))
            push!(busarctuples, (cyc[3], (cyc[1], cyc[4]), (cyc[1], cyc[2])))
         end
         expairs_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2), (5,6), (3,4), (7,8), (1,6), (2,5), (4,7), (3,8), (2,3), (6,7), (1,4), (5,8), (2,7), (3,6), (4,5), (1,8)]
         if params["rotate_4w"]
            expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8), (1,2, 12), (5,6,12), (3,4,10), (7,8,10), (1,6,12), (2,5,12), (4,7,10), (3,8,10), (2,3,9), (6,7,9), (1,4,11), (5,8,11), (2,7,9), (3,6,9), (4,5,11), (1,8,11)]
         else
            expairs_w_4 = [(1,3), (5,7), (2,4), (6,8), (1,7), (3,5), (4,6), (2,8)]
         end
      end
      if params["cycle_relax"] == "mc"
         if params["cycle_c_s_cuts"]
            JuMP.@variable(_m, hcc[ap in arcpairs]) # lifted variable for cs * cs, for McCormick
            JuMP.@variable(_m, hss[ap in arcpairs])
            JuMP.@variable(_m, hcs[ap in arcpairs])
            JuMP.@variable(_m, hsc[ap in arcpairs])
            vars["hcc"] = hcc
            vars["hss"] = hss
            vars["hcs"] = hcs
            vars["hsc"] = hsc
         end

         if params["cycle_wr_wi_cuts"]
            JuMP.@variable(_m, wwr[ba in busarcpairs]) # lifted variable for w * wr, for McCormick
            JuMP.@variable(_m, wwi[ba in busarcpairs])
            vars["wwr"] = wwr
            vars["wwi"] = wwi

            wrr_bounds = Dict()
            wii_bounds = Dict()
            wri_bounds = Dict()
            wir_bounds = Dict()
            for ap in arcpairs
               wrr_bounds[ap] = Get_bounds([wr[ap[1]], wr[ap[2]]])
               wii_bounds[ap] = Get_bounds([wi[ap[1]], wi[ap[2]]])
               wri_bounds[ap] = Get_bounds([wr[ap[1]], wi[ap[2]]])
               wir_bounds[ap] = Get_bounds([wi[ap[1]], wr[ap[2]]])
            end
            JuMP.@variables(_m, begin
               wrr_bounds[ap][1] <= wrr[ap in arcpairs] <= wrr_bounds[ap][2] # lifted variable for wr * wr, for McCormick
               wii_bounds[ap][1] <= wii[ap in arcpairs] <= wii_bounds[ap][2]
               wri_bounds[ap][1] <= wri[ap in arcpairs] <= wri_bounds[ap][2]
               wir_bounds[ap][1] <= wir[ap in arcpairs] <= wir_bounds[ap][2]
            end)
            vars["wrr"] = wrr
            vars["wii"] = wii
            vars["wri"] = wri
            vars["wir"] = wir

            JuMP.@variable(_m, wwrr[ba in busarctuples])
            JuMP.@variable(_m, wwii[ba in busarctuples])
            JuMP.@variable(_m, wwri[ba in busarctuples])
            JuMP.@variable(_m, wwir[ba in busarctuples])
            vars["wwrr"] = wwrr
            vars["wwii"] = wwii
            vars["wwri"] = wwri
            vars["wwir"] = wwir
         end
      elseif params["cycle_relax"] == "epr"
         cyc_i = Dict()
         cyc_ep = Dict()
         cycle_3_4 = []
         append!(cycle_3_4, cycle["cyc_3"])
         if params["cycle_max_bnd"] >= 4
            append!(cycle_3_4, cycle["cyc_4"])
         end
         if params["cycle_c_s_cuts"]
            cyc_i[3] = 1:2^6
            cyc_i[4] = 1:2^8
            cyc_ep[3] = expairs_3
            cyc_ep[4] = expairs_4
            JuMP.@variable(_m, λ_c[cyc in cycle_3_4, cyc_i[length(cyc)]] >= 0)
            JuMP.@variable(_m, x_c[cyc in cycle_3_4, ep in cyc_ep[length(cyc)]])
            vars["λ_c"] = λ_c
            vars["x_c"] = x_c
         end
         if params["cycle_wr_wi_cuts"]
            cyc_i[3] = 1:(2^9)
            cyc_i[4] = 0
            if params["rotate_4w"]
               cyc_i[4] = 1:(2^12)
            else
               cyc_i[4] = 1:(2^8)
            end
            cyc_ep[3] = expairs_w_3
            cyc_ep[4] = expairs_w_4
            JuMP.@variable(_m, λ_w[cyc in cycle_3_4, cyc_i[length(cyc)]] >= 0)
            JuMP.@variable(_m, x_w[cyc in cycle_3_4, ep in cyc_ep[length(cyc)]])
            vars["λ_w"] = λ_w
            vars["x_w"] = x_w
         end
      end
   end

   if params["warm_start"]
      if params["model"] == "opf"
         for i in keys(ref[:bus])
            # JuMP.set_start_value(va[i], )
            JuMP.set_start_value(vm[i], (ref[:bus][i]["vmax"] + ref[:bus][i]["vmin"]) / 2)
         end
         for bp in bp_li
            JuMP.set_start_value(cs[bp], (cos_min[bp] + cos_max[bp]) / 2)
            JuMP.set_start_value(si[bp], (sin_min[bp] + sin_max[bp]) / 2)
         end
      elseif params["model"] == "ots" || params["model"] == "ots_relax"
         if params["cycle_cuts"] == false
            warm_start_ots(vars, bp_li)
         else
            warm_start_ots_cyc(vars, bp_li, cycle_3_4)
         end
      end
   elseif params["partial_start"]
      for bp in keys(ref[:buspairs])
         JuMP.set_start_value(z[bp], 1)
      end
   end

   if params["preprocess_tightening"]
      preprocess_pqz(_m)
   end

   return vars, bd_data
end

# This warm start function does not deal with model with cycle. See next function for warm start with cycle constraints.
function warm_start_ots(vars, bp_li)
   string_cycle = ""
   string_rmc = ""
   if params["relax"] == "rmc"
      string_rmc = "rmc_"
   end
   input_path = string("results/opf_sol/", params["instance"], "_", string_cycle, string_rmc, "sol.txt")
   str = read(input_path, String)
   sol = eval(Meta.parse(str))
   # Initialize solutions for ots-specific variables.
   sol["z"] = Dict()
   for bp in keys(ref[:buspairs])
      sol["z"][bp] = 1
   end
   sol["wz"] = Dict()
   for bp in bp_li
      sol["wz"][bp] = sol["w"][bp[1]]
   end
   if params["relax"] == "rmc"
      sol["vmz"] = Dict()
      for bp in keys(ref[:buspairs])
         sol["vmz"][bp] = sol["vm"][bp[1]]
         sol["vmz"][(bp[2], bp[1])] = sol["vm"][bp[2]]
      end
   end

   for item in keys(vars)
      println(item)
      the_var = vars[item]
      the_sol = sol[item]
      for ind in keys(the_var)
         if length(size(the_var)) == 1
            JuMP.set_start_value(the_var[ind[1]], the_sol[ind[1]])
         elseif length(size(the_var)) == 2
            JuMP.set_start_value(the_var[ind[1], ind[2]], the_sol[ind[1], ind[2]])
         end
      end
   end
end

function warm_start_ots_cyc(vars, bp_li, cycle_3_4)
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
   input_path = string("results/opf_sol/", params["instance"], "_", string_cycle, string_load, string_circ2, "sol.txt")
   str = read(input_path, String)
   sol = eval(Meta.parse(str))

   # Initialize solutions for ots-specific variables.
   sol["z"] = Dict()
   for bp in keys(ref[:buspairs])
      sol["z"][string(vars["z"][bp])] = 1
      # sol["z"][vars["z"][bp]] = 1
   end
   sol["wz"] = Dict()
   for bp in bp_li
      sol["wz"][string(vars["wz"][bp])] = sol["w"][string(vars["w"][bp[1]])]
   end
   sol["zc"] = Dict()
   for cyc in cycle_3_4
      sol["zc"][string(vars["zc"][Tuple(cyc)])] = 1 # Equals 1 iff all lines in a cycle are turned on.
   end

   for item in keys(vars)
      the_var = vars[item]
      the_sol = sol[item]
      for item2 in the_var
         # set_lower_bound(item2, the_sol[string(item2)])
         # set_upper_bound(item2, the_sol[string(item2)])
         JuMP.set_start_value(item2, the_sol[string(item2)])
      end
   end
end

function preprocess_pqz(_m)
   tightened = false
    for (l, load) in ref[:load]
        load_bus = load["load_bus"]
        net_pd = load["pd"]
        net_qd = load["qd"]
        if length(ref[:bus_gens][load_bus]) >= 1
           net_pd -= sum(ref[:gen][g]["pmax"] for g in ref[:bus_gens][load_bus])
           net_qd -= sum(ref[:gen][g]["qmax"] for g in ref[:bus_gens][load_bus])
        end
        sum_rate_a = sum(ref[:branch][a[1]]["rate_a"] for a in ref[:bus_arcs][load_bus])
        for a in ref[:bus_arcs][load_bus]
           branch_idx = a[1]
           bp_idx = (a[2], a[3])
           max_load = max(net_pd - (sum_rate_a - ref[:branch][branch_idx]["rate_a"]), net_qd - (sum_rate_a - ref[:branch][branch_idx]["rate_a"]))
           println("load_bus: ", load_bus, ", bp_idx: ", bp_idx, ", max_load: ", max_load)
           if max_load > 1e-3
             set_lower_bound(_m[:z][bp_idx], 1)
             tightened = true
           end
        end
    end
    if tightened == true
       println("Bound tightened!")
    else
      println("Bound not tightened..")
    end
end
