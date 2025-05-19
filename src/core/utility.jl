function Get_bounds(v)
   v_l = [lower_bound(v[1]), upper_bound(v[1])]
   v_u = [lower_bound(v[2]), upper_bound(v[2])]
   M = v_l*v_u'

   if length(v) == 2
      # Bilinear
      return minimum(M), maximum(M)

   elseif (length(v) == 3) || (length(v) == 4)
      M1 = zeros(Float64, (2,2,2))
      M1[:,:,1] = M*lower_bound(v[3])
      M1[:,:,2] = M*upper_bound(v[3])
      if length(v) == 3
         # Trilinear
         return minimum(M1), maximum(M1)
      elseif length(v) == 4
         # Quadrilinear
         M2 = zeros(Float64, (2,2,4))
         M2[:,:,1:2] = M1*lower_bound(v[4])
         M2[:,:,3:4] = M1*upper_bound(v[4])
         return minimum(M2), maximum(M2)
      end
   end
end

function Get_cycle_adjacency(c)
   A = Set()
   for i=1:(length(c)-1)
      I = min(c[i],c[i+1])
      J = max(c[i],c[i+1])
      # @printf("I: %s, J: %s\n", I, J)
      # @assert (I,J) in keys(ref[:buspairs])
      push!(A, (I,J))
   end
   I = min(c[1],c[end])
   J = max(c[1],c[end])
   # @assert (I,J) in keys(ref[:buspairs])
   push!(A, (I,J))
   return A
end

function Get_unique_cycles(cyc)
   for i=1:(length(cyc)-1)
      for j=(i+1):length(cyc)
         A_i = Get_cycle_adjacency(cyc[i])
         A_j = Get_cycle_adjacency(cyc[j])
         if A_i == A_j
            if cyc[j] != cyc[i]
               cyc[j] = cyc[i]
            end
         end
      end
   end
   return unique(cyc)
end

function Get_cycles(ref, params)
   #num_bus = length(ref[:bus])
   num_bus = length(data["bus"])
   num_edges = length(ref[:branch])

   if params["LightGraphs"] == false
      println(">>>>>>> Invoking custom labeling code for cycle enumeration <<<<<<<<")
      A = zeros(Int, num_bus, num_bus)
      #for i=1:length(ref[:branch])
      for i in keys(ref[:branch])
         A[ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"]] = 1
         A[ref[:branch][i]["t_bus"], ref[:branch][i]["f_bus"]] = 1
      end

      # Get all 3-cycles
      if params["cycle_max_bnd"] <= 4
         cyc_3 = Any[]
         #for i=1:length(ref[:branch])
         for i in keys(ref[:branch])
            e = (ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"])
            buses = setdiff(keys(ref[:bus]), e)
            for j in buses
               v3 = [e[1],e[2],j]
               if A[j, e[1]] == 1 && A[j, e[2]] == 1 && (length(unique(v3)) == 3)
                  push!(cyc_3, v3)
               end
            end
         end
         #@show cyc_3
         cyc_3 = Get_unique_cycles(unique(cyc_3))

         println("Total number of 3 cycles: ", length(cyc_3))
         if params["cycle_max_bnd"] == 3
            return cyc_3
         end
      end
   elseif params["LightGraphs"] == true
      println(">>>>>>> Invoking LightGraphs for cycle enumeration <<<<<<<<")
      A = zeros(Int, num_bus, num_bus)
      for i in keys(ref[:branch])
         A[ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"]] = 1
         A[ref[:branch][i]["t_bus"], ref[:branch][i]["f_bus"]] = 1
      end
      D = SimpleDiGraph(A)
      cyc_time = @elapsed cyc_2_6 = simplecycles_limited_length(D, params["cycle_max_bnd"])
      @printf("# Time for finding cycles: %s secs\n", cyc_time)

      cyc_3 = Any[]
      cyc_4 = Any[]
      cyc_5 = Any[]
      cyc_6 = Any[]
      for li in cyc_2_6
         if length(li) == 3
            push!(cyc_3, li)
         elseif length(li) == 4
            push!(cyc_4, li)
         elseif length(li) == 5
            push!(cyc_5, li)
         elseif length(li) == 6
            push!(cyc_6, li)
         end
      end
      cyc = [Any[] for i in 3:6]
      @elapsed cyc[1] = Get_unique_cycles(unique(cyc_3))
      @elapsed cyc[2] = Get_unique_cycles(unique(cyc_4))
      @elapsed cyc[3] = Get_unique_cycles(unique(cyc_5))
      @elapsed cyc[4] = Get_unique_cycles(unique(cyc_6))

      println("Total number of 3 cycles: ", length(cyc[1]))
      println("Total number of 4 cycles: ", length(cyc[2]))
      println("Total number of 5 cycles: ", length(cyc[3]))
      println("Total number of 6 cycles: ", length(cyc[4]))
      cycle = Dict()
      for i in 3:params["cycle_max_bnd"]
         str = string("cyc_", string(i))
         cycle[str] = cyc[i - 2]
      end
   end
   return cycle
end

function Get_cycle_basis(ref, params)
   #num_bus = length(ref[:bus])
   num_bus = length(data["bus"])
   num_edges = length(ref[:branch])

   if params["LightGraphs"] == true
      println(">>>>>>> Invoking LightGraphs for cycle basis enumeration <<<<<<<<")
      elist = []
      for i in keys(ref[:branch])
         push!(elist, (ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"]))
      end
      G = SimpleGraph(Edge.(elist))
      cyc_time = @elapsed cyc_basis = cycle_basis(G)
      @printf("# Time for finding cycles basis: %s secs\n", cyc_time)
      cyc_3 = Any[]
      cyc_4 = Any[]
      cyc_5 = Any[]
      cyc_6 = Any[]
      for li in cyc_basis #cyc_2_6
         if length(li) == 3
            push!(cyc_3, li)
         elseif length(li) == 4
            push!(cyc_4, li)
         elseif length(li) == 5
            push!(cyc_5, li)
         elseif length(li) == 6
            push!(cyc_6, li)
         end
      end
      cyc = [Any[] for i in 3:6]
      # The following makes sure each edge first index is smaller than the second.
      cyc[1] = [sort!(cyc) for cyc in cyc_3]
      cyc[2] = [sort!(cyc) for cyc in cyc_4]
      cyc[3] = [sort!(cyc) for cyc in cyc_5]
      cyc[4] = [sort!(cyc) for cyc in cyc_6]

      println("Total number of 3 cycles: ", length(cyc[1]))
      println("Total number of 4 cycles: ", length(cyc[2]))
      println("Total number of 5 cycles: ", length(cyc[3]))
      println("Total number of 6 cycles: ", length(cyc[4]))
      cycle = Dict()
      for i in 3:params["cycle_max_bnd"]
         str = string("cyc_", string(i))
         cycle[str] = cyc[i - 2]
      end
   end
   return cycle
end

function Get_span_tree(ref)
   println(">>>>>>> Invoking LightGraphs for getting spanning tree <<<<<<<<")
      elist = []
      for i in keys(ref[:branch])
         if params["cal_weight"]
            weight = Get_span_tree_weight(ref)
            push!(elist, (ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"], weight[i]))
         else
            push!(elist, (ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"], 1.0))
         end
      end
      G = SimpleGraph(Edge.(elist))
      comps = connected_components(G)

      span_tree_edges = []
      tree_time = 0.0
      for comp in comps
         # Create a subgraph for the current component
         subgraph, vmap = induced_subgraph(G, comp)
         
         # Compute the minimum spanning tree using Kruskal's algorithm
         if params["cal_weight"]
            tree_time = @elapsed tree = kruskal_mst(subgraph; minimize = false)
         else
            tree_time = @elapsed tree = kruskal_mst(subgraph)
         end
     
         for e in tree
            push!(span_tree_edges, (vmap[src(e)], vmap[dst(e)]))
         end
      end
      @printf("# Time for finding spanning tree: %s secs\n", tree_time)

      return span_tree_edges
end

function Get_span_tree_weight(ref)
   nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes")
   result = solve_ac_opf(data, nlp_solver)

   # minlp_solver = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => params["time_limit"])
   # result = solve_opf(data, QCLSPowerModel, minlp_solver)


   weight = Dict()
   for (i, branch) in ref[:branch]
      qf = result["solution"]["branch"][string(i)]["qf"]
      qt = result["solution"]["branch"][string(i)]["qt"]
      pf = result["solution"]["branch"][string(i)]["pf"]
      pt = result["solution"]["branch"][string(i)]["pt"]
      weight[i] = max((qf^2 + qt^2)/branch["rate_a"]^2, (pf^2 + pt^2)/branch["rate_a"]^2)
   end

   return weight
end

function Write_model(_m, _file_path)
   open(_file_path, "w") do io
      println(io, _m)
   end
end

# function Write_LP(_m, _lp_path)
#    # lp_file = MathOptFormat.LP.Model()
#    lp_file = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_LP)
#    MOI.copy_to(lp_file, backend(_m))
#    MOI.write_to_file(lp_file, _lp_path)
# end

function cartesian_product(lb, ub)
   v = []
   n = length(lb)
   for i in 1:n
      push!(v, (lb[i], ub[i]))
   end
   return collect(Base.product(v...)) # The splat operator, ..., takes out the entries in v.
end

# function cartesian_product(lb, ub)
#     # Input: lb and ub are the vectors of lower and upper bounds of the variables
#     # Output: Generate cartesian product of sets of 2 elements, i.e., lower and upper bounds of variables
#
#     @assert length(lb) == length(ub)
#     n = length(lb)
#     # @assert n <= 8
#
#     (typeof(lb[1]) == String) && (cart_prod = Matrix{String}(undef, 2^n, n))
#     (typeof(lb[1]) != String) && (cart_prod = Matrix{Float64}(undef, 2^n, n))
#     ii = 1
#
#     if n == 2
#        v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]];
#        for i=1:2
#           for j=1:2
#             cart_prod[ii,:] =  [v1[j], v2[i]]
#             ii += 1
#           end
#        end
#
#     elseif n == 3
#        v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]]; v3 = [lb[3], ub[3]];
#        for i=1:2
#           for j=1:2
#              for k=1:2
#                 cart_prod[ii,:] =  [v1[k], v2[j], v3[i]]
#                 ii += 1
#              end
#           end
#        end
#
#     elseif n == 4
#        v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]]; v3 = [lb[3], ub[3]]; v4 = [lb[4], ub[4]];
#        for i=1:2
#           for j=1:2
#              for k=1:2
#                 for l=1:2
#                    cart_prod[ii,:] =  [v1[l], v2[k], v3[j], v4[i]]
#                    ii += 1
#                 end
#              end
#           end
#        end
#
#     elseif n == 5
#        v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]]; v3 = [lb[3], ub[3]]; v4 = [lb[4], ub[4]]; v5 = [lb[5], ub[5]];
#        for i=1:2
#           for j=1:2
#              for k=1:2
#                 for l=1:2
#                    for m=1:2
#                       cart_prod[ii,:] =  [v1[m], v2[l], v3[k], v4[j], v5[i]]
#                       ii += 1
#                    end
#                 end
#              end
#           end
#        end
#
#     elseif n == 6
#        v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]]; v3 = [lb[3], ub[3]]; v4 = [lb[4], ub[4]]; v5 = [lb[5], ub[5]]; v6 = [lb[6], ub[6]];
#        for i=1:2
#           for j=1:2
#              for k=1:2
#                 for l=1:2
#                    for m=1:2
#                       for o=1:2
#                          cart_prod[ii,:] =  [v1[o], v2[m], v3[l], v4[k], v5[j], v6[i]]
#                          ii += 1
#                       end
#                    end
#                 end
#              end
#           end
#        end
#
#     elseif n == 7
#        v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]]; v3 = [lb[3], ub[3]]; v4 = [lb[4], ub[4]]; v5 = [lb[5], ub[5]]; v6 = [lb[6], ub[6]]; v7 = [lb[7], ub[7]];
#        for i=1:2
#           for j=1:2
#              for k=1:2
#                 for l=1:2
#                    for m=1:2
#                       for o=1:2
#                          for p=1:2
#                             cart_prod[ii,:] =  [v1[p], v2[o], v3[m], v4[l], v5[k], v6[j], v7[i]]
#                             ii += 1
#                          end
#                       end
#                    end
#                 end
#              end
#           end
#        end
#
#     elseif n == 8
#         v1 = [lb[1], ub[1]]; v2 = [lb[2], ub[2]]; v3 = [lb[3], ub[3]]; v4 = [lb[4], ub[4]]; v5 = [lb[5], ub[5]]; v6 = [lb[6], ub[6]]; v7 = [lb[7], ub[7]]; v8 = [lb[8], ub[8]];
#         for i=1:2
#            for j=1:2
#               for k=1:2
#                  for l=1:2
#                     for m=1:2
#                        for o=1:2
#                           for p=1:2
#                             for q=1:2
#                                 cart_prod[ii,:] = [v1[q], v2[p], v3[o], v4[m], v5[l], v6[k], v7[j], v8[i]]
#                                 ii += 1
#                             end
#                           end
#                        end
#                     end
#                  end
#               end
#            end
#         end
#
#      else
#         v = []
#         for i in 1:n
#            push!(v, (lb[i], ub[i]))
#         end
#         cart_prod = collect(Base.product(v...))
#     end
#
#     return cart_prod
# end

# No more need this function after defining cs and w for reversed arcs.
# For cycle constraints, sort the index of an arc then return the corresponding  trigonometric variable. Input: _arc: tuple, an arc; "trigo": string, "cs" or "si".
function trigo_correct_ap(_arc, _trigo)
   if _arc[1] < _arc[2]
      if _trigo == "cs"
         return cs[_arc]
      elseif _trigo == "si"
         return si[_arc]
      end
   else
      if _trigo == "cs"
         return cs[(_arc[2], _arc[1])]
      elseif _trigo == "si"
         return -si[(_arc[2], _arc[1])]
      end
   end
end

# Getting the bounds of cs, si, wr, and wi variables.
function get_bounds_csw(arg_min, arg_max, vm_min, vm_max)
   wr_min = Dict((bp, -Inf) for bp in keys(ref[:buspairs]))
   wr_max = Dict((bp,  Inf) for bp in keys(ref[:buspairs]))
   wi_min = Dict((bp, -Inf) for bp in keys(ref[:buspairs]))
   wi_max = Dict((bp,  Inf) for bp in keys(ref[:buspairs]))
   cos_min = Dict([(bp, -Inf) for bp in keys(ref[:buspairs])])
   cos_max = Dict([(bp, Inf) for bp in keys(ref[:buspairs])])
   sin_min = Dict([(bp, -Inf) for bp in keys(ref[:buspairs])])
   sin_max = Dict([(bp, Inf) for bp in keys(ref[:buspairs])])

   for bp in keys(ref[:buspairs])
     i,j = bp
     sin_min[bp] = sin(arg_min[bp])
     sin_max[bp] = sin(arg_max[bp])

     if arg_min[bp] >= 0
         cos_max[bp] = cos(arg_min[bp])
         cos_min[bp] = cos(arg_max[bp])
         wr_max[bp] = vm_max[i] * vm_max[j] * cos_min[bp]
         wr_min[bp] = vm_min[i] * vm_min[j] * cos_max[bp]
         wi_max[bp] = vm_max[i] * vm_max[j] * sin_max[bp]
         wi_min[bp] = vm_min[i] * vm_min[j] * sin_min[bp]
     end
     if arg_max[bp] <= 0
         cos_max[bp] = cos(arg_max[bp])
         cos_min[bp] = cos(arg_min[bp])
         wr_max[bp] = vm_max[i] * vm_max[j] * cos_max[bp]
         wr_min[bp] = vm_min[i] * vm_min[j] * cos_min[bp]
         wi_max[bp] = vm_min[i] * vm_min[j] * sin_max[bp]
         wi_min[bp] = vm_max[i] * vm_max[j] * sin_min[bp]
     end
     if arg_min[bp] < 0 && arg_max[bp] > 0
         cos_max[bp] = 1.0
         cos_min[bp] = min(cos(arg_min[bp]), cos(arg_max[bp]))
         wr_max[bp] = vm_max[i] * vm_max[j] * 1.0
         wr_min[bp] = vm_min[i] * vm_min[j] * min(cos_min[bp], cos_max[bp])
         wi_max[bp] = vm_max[i] * vm_max[j] * sin_max[bp]
         wi_min[bp] = vm_max[i] * vm_max[j] * sin_min[bp]
     end
   end

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

   return wr_min, wr_max, wi_min, wi_max, cos_min, cos_max, sin_min, sin_max
end

# Round down small entries of a vector to zero.
function round_vector(vec)
   for i in 1:length(vec)
      if abs(vec[i]) <= 1e-4
         vec[i] = 0
      end
   end
   return vec
end

function round_densearray(vec)
   # for ind in keys(vec)
   #    if abs(vec[ind]) <= 1e-4
   #       vec[ind.I[1]] = 0
   #    end
   # end
   return vec
end

function leaf_gens(ref)
   ref[:gens_on_leaf] = Dict()
   ref[:gens_not_on_leaf] = []
   temp = Dict()
   for (l, load) in ref[:load]
      temp[load["load_bus"]] = Dict()
      temp[load["load_bus"]]["pd"] = load["pd"]
      temp[load["load_bus"]]["qd"] = load["qd"]
   end
   bus_load = Dict()
   for bus in keys(ref[:bus])
      bus_load[bus] = Dict()
      if bus in keys(temp)
         bus_load[bus]["pd"] = temp[bus]["pd"]
         bus_load[bus]["qd"] = temp[bus]["qd"]
      else
         bus_load[bus]["pd"] = 0
         bus_load[bus]["qd"] = 0
      end
   end
   for bus in keys(ref[:bus])
      # println("bus: ", bus)
      # println("ref[:bus_arc]: ", ref[:bus_arcs])
      arcs = ref[:bus_arcs][bus]
      if length(arcs) == 1 && bus_load[bus]["pd"] == 0 && bus_load[bus]["qd"] == 0
         for gen in ref[:bus_gens][bus]
            ind, i, j = ref[:bus_arcs][bus][1]
            # Making sure the order of arc is the same as the order of buspair, as z is indexed by buspair.
            if (i, j) in keys(ref[:buspairs])
               ref[:gens_on_leaf][gen] = (i, j)
            else
               ref[:gens_on_leaf][gen] = (j, i)
            end
         end
      else
         append!(ref[:gens_not_on_leaf], ref[:bus_gens][bus])
      end
   end
end

function mcc_ots_wij(m, bp)
   i, j = bp
   lb = Dict()
   ub = Dict()
   lb["vmz"] = Dict()
   lb["vmz"][i] = ref[:bus][i]["vmin"]
   lb["vmz"][j] = ref[:bus][j]["vmin"]
   ub["vmz"] = Dict()
   ub["vmz"][i] = ref[:bus][i]["vmax"]
   ub["vmz"][j] = ref[:bus][j]["vmax"]

   lb["vvz"] = ref[:buspairs][bp]["vm_fr_min"] * ref[:buspairs][bp]["vm_to_min"]
   ub["vvz"] = ref[:buspairs][bp]["vm_fr_max"] * ref[:buspairs][bp]["vm_to_max"]
   #
   lb["csz"] = bd_data["cos_min"][bp]
   ub["csz"] = bd_data["cos_max"][bp]
   #
   lb["siz"] = bd_data["sin_min"][bp]
   ub["siz"] = bd_data["sin_max"][bp]

   if params["relax"] == "rmc"
      lb["vvz"] = 0
      lb["csz"] = 0
   end

   function bds_on_off(x, xz, lb, ub)
      # Hijazi et. al. do not have the following two constraints.
      @constraint(m, x - (1 - m[:z][bp]) * ub <= xz)
      @constraint(m, xz <= x - (1 - m[:z][bp]) * lb)
      if params["ext_cons"]
         @constraint(m, lb * m[:z][bp] <= xz)
         @constraint(m, xz <= ub * m[:z][bp])
      end
   end

   function mcc_on_off(xy, xz, yz, xlb, xub, ylb, yub)
      @constraint(m, xy >= xlb * yz + ylb * xz - xlb * ylb * m[:z][bp])
      @constraint(m, xy >= xub * yz + yub * xz - xub * yub * m[:z][bp])
      @constraint(m, xy <= xlb * yz + yub * xz - xlb * yub * m[:z][bp])
      @constraint(m, xy <= xub * yz + ylb * xz - xub * ylb * m[:z][bp])
   end

   bds_on_off(m[:vm][i], m[:vmz][bp], lb["vmz"][i], ub["vmz"][i])
   bds_on_off(m[:vm][j], m[:vmz][(j, i)], lb["vmz"][j], ub["vmz"][j])

   @constraint(m, m[:vv][bp] <= ub["vmz"][i] * ub["vmz"][j] * m[:z][bp])
   @constraint(m, lb["vmz"][i] * lb["vmz"][j] * m[:z][bp] <= m[:vv][bp])

   mcc_on_off(m[:vv][bp], m[:vmz][bp], m[:vmz][(j,i)], lb["vmz"][i], ub["vmz"][i], lb["vmz"][j], ub["vmz"][j])

   mcc_on_off(m[:wr][bp], m[:vv][bp], m[:cs][bp], lb["vvz"], ub["vvz"], lb["csz"], ub["csz"])
   mcc_on_off(m[:wi][bp], m[:vv][bp], m[:si][bp], lb["vvz"], ub["vvz"], lb["siz"], ub["siz"])
end

# Turn a tuple for cycle to the set of buspairs in the cycle.
function cyc_to_buspair(cyc)
  buspairs = []
  for i in 1:length(cyc)
    ind1 = i
    ind2 = i + 1
    if i == length(cyc)
      ind2 = 1
    end
    bp = (cyc[ind1], cyc[ind2])
    if bp in keys(ref[:buspairs])
      push!(buspairs, bp)
    else
      push!(buspairs, (cyc[ind2], cyc[ind1]))
    end
  end
  return buspairs
end

# Checking the feasibility of McCormick relaxation with PM results.
function check_pm_mcc(result)
   # Need variables: wr -> , wi, vv, cs, si, vmz, vm
   for (bp, buspair) in ref[:buspairs]

   end
end
