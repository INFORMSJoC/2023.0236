
function network_data(instance, params)

   if instance[end-2:end] == "api"
      if instance[1:5] == "pglib"
      # if params["instance_type"] == "pglib"
         path = "data/pglib/api/"
      elseif instance[1:5] == "nesta"
      # elseif params["instance_type"] == "nesta"
         path = "data/nesta/api/"
      end
   elseif instance[end-2:end] == "sad"
      if instance[1:5] == "pglib"
      # if params["instance_type"] == "pglib"
         path = "data/pglib/sad/"
      elseif instance[1:5] == "nesta"
      # elseif params["instance_type"] == "nesta"
         path = "data/nesta/sad/"
      end
   elseif instance[end-2:end] == "nco"
         path = "data/nesta/nco/"
   elseif instance[end-2:end] == "utl"
         path = "data/nesta/utl/"
   elseif instance[end-2:end] == "rad"
         path = "data/nesta/rad/"
   else
      if instance[1:5] == "pglib"
      # if params["instance_type"] == "pglib"
         path = "data/pglib/"
      elseif instance[1:5] == "nesta"
      # elseif params["instance_type"] == "nesta"
         path = "data/nesta/"
      end
   end

   println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
   println("Running $instance")
   #println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
   data = PowerModels.parse_file(string(path,instance,".m"))

   #------------------------------------------;
   # Set bus numbers to linear ordering
   #------------------------------------------;
   if params["correct_bus_numbering"]
      bus_ids = Any[]
      for i in keys(data["bus"])
         push!(bus_ids, data["bus"]["$i"]["bus_i"])
      end
      if (minimum(bus_ids) != 1) || (maximum(bus_ids) != length(keys(data["bus"])))
         println(">>>>>>> Creating bus map to have a linear ordering <<<<<<<")

         data_new = deepcopy(data)
         # Update data["bus"]
         for i in keys(data_new["bus"])
            delete!(data_new["bus"], "$i")
         end

         k=1
         for i in keys(sort(data["bus"]))
            data_new["bus"]["$k"] = Dict{String, Any}("zone"  => data["bus"]["$i"]["zone"],
                                                      "bus_i" => k,
                                                      "bus_type" => data["bus"]["$i"]["bus_type"],
                                                      "vmax" => data["bus"]["$i"]["vmax"],
                                                      "source_id" => Any["bus", k],
                                                      "area" => data["bus"]["$i"]["area"],
                                                      "vmin" => data["bus"]["$i"]["vmin"],
                                                      "index" => k,
                                                      "va" => data["bus"]["$i"]["va"],
                                                      "vm" => data["bus"]["$i"]["vm"],
                                                      "base_kv" => data["bus"]["$i"]["base_kv"]
                                                     )
            data["bus"]["$i"]["map"] = k  #Create a map of bus ids
            k+=1
         end

         # Update data["gen"]
         for i in keys(data_new["gen"])
            bus_id = data["gen"]["$i"]["gen_bus"]
            data_new["gen"]["$i"]["gen_bus"] = data["bus"]["$bus_id"]["map"]
         end

         # Update data["branch"]
         for i in keys(data_new["branch"])
            f_bus_id = data["branch"]["$i"]["f_bus"]
            t_bus_id = data["branch"]["$i"]["t_bus"]
            data_new["branch"]["$i"]["f_bus"] = data["bus"]["$f_bus_id"]["map"]
            data_new["branch"]["$i"]["t_bus"] = data["bus"]["$t_bus_id"]["map"]
         end

         # Update data["load"]
         for i in keys(data_new["load"])
            load_bus_id = data["load"]["$i"]["load_bus"]
            source_id_ = data["load"]["$i"]["source_id"][2]
            data_new["load"]["$i"]["load_bus"] = data["bus"]["$load_bus_id"]["map"]
            data_new["load"]["$i"]["source_id"][2] = data["bus"]["$source_id_"]["map"]
         end

         # Update data["shunt"]
         for i in keys(data_new["shunt"])
            shunt_bus_id = data["shunt"]["$i"]["shunt_bus"]
            source_id_ = data["shunt"]["$i"]["source_id"][2]
            data_new["shunt"]["$i"]["shunt_bus"] = data["bus"]["$shunt_bus_id"]["map"]
            data_new["shunt"]["$i"]["source_id"][2] = data["bus"]["$source_id_"]["map"]
         end

         data = deepcopy(data_new)
      end
   end

   #----------------------------------------------;
   # Set from and to bus numbering for branches
   #----------------------------------------------;
   # if params["cycle_cuts"]
   #    println("Activating correct branch directions since cycle_cuts are active")
   #    params["correct_branch_directions"] = true
   # end

   if params["correct_branch_directions"]
      for i=1:length(data["branch"])
         branch_orginal = copy(data["branch"]["$i"])
         branch = data["branch"]["$i"]

         f = branch_orginal["f_bus"]
         t = branch_orginal["t_bus"]
         if f > t
            # println("Reversing the orientation of branch ($f,$t) to be consistent with cycle cuts")

            branch["f_bus"] = branch_orginal["t_bus"]
            branch["t_bus"] = branch_orginal["f_bus"]
            branch["g_to"] = branch_orginal["g_fr"] .* branch_orginal["tap"]'.^2
            branch["b_to"] = branch_orginal["b_fr"] .* branch_orginal["tap"]'.^2
            branch["g_fr"] = branch_orginal["g_to"] ./ branch_orginal["tap"]'.^2
            branch["b_fr"] = branch_orginal["b_to"] ./ branch_orginal["tap"]'.^2
            branch["tap"] = 1 ./ branch_orginal["tap"]
            branch["br_r"] = branch_orginal["br_r"] .* branch_orginal["tap"]'.^2
            branch["br_x"] = branch_orginal["br_x"] .* branch_orginal["tap"]'.^2
            branch["shift"] = -branch_orginal["shift"]
            branch["angmin"] = -branch_orginal["angmax"]
            branch["angmax"] = -branch_orginal["angmin"]

         end
      end
   end

   # data["branch"]["1"]["br_status"] = 0

   if params["turn_off_lines"] # || params["off_insights"]
      turn_off_lines(instance, data)
   end

   return data
end

function turn_off_lines(instance, data)
   # result_path = "results/fix_dec/rmc.csv"
   result_path = "results/fix_dec/ext.csv"
   # result_path = "results/fix_dec/cyc.csv"

   # if params["cycle_cuts"] == false
   #    result_path = "results/noCycle_ots_results_true_false-grb.csv"
   # else
   #    # result_path = "results/ots_qc_pm_results_true_false.csv"
   #    result_path = "results/ots_cb_4cycle_results_epr_both_noRotation4wepr-grb.csv"
   # end
   df = CSV.read(result_path, DataFrame)
   df[!, :off_lines] = convert(Array{Any,1}, df[!, :off_lines])

   ind = findall(df[!,:Instance][1:72] .== instance)[1]
   # if params["cycle_cuts"] == false
   #    ind = findall(df[!,:Instance][1:51] .== instance)[1]
   # else
   #    ind = findall(df[!,:Instance][1:10] .== instance)[1]
   # end
   z_bar = eval(Meta.parse(df[!, :off_lines][ind]))

   br_off_ind = bp_to_branch(z_bar, data)

   for br in br_off_ind
      data["branch"][br]["br_status"] = 0
   end
end

# For getting on and off lines of an instance.
function on_off_lines(instance, data, ref)
   # result_path = "results/fix_dec/rmc.csv"
   result_path = "results/fix_dec/ext.csv"
   # result_path = "results/fix_dec/cyc.csv"

   # if params["cycle_cuts"] == false
   #    result_path = "results/noCycle_ots_results_true_false-grb.csv"
   # else
   #    # result_path = "results/ots_qc_pm_results_true_false.csv"
   #    result_path = "results/ots_cb_4cycle_results_epr_both_noRotation4wepr-grb.csv"
   # end
   df = CSV.read(result_path, DataFrame)
   df[!, :off_lines] = convert(Array{Any,1}, df[!, :off_lines])

   ind = findall(df[!,:Instance][1:72] .== instance)[1]
   # if params["cycle_cuts"] == false
   #    ind = findall(df[!,:Instance][1:51] .== instance)[1]
   # else
   #    ind = findall(df[!,:Instance][1:10] .== instance)[1]
   # end
   z_bar = eval(Meta.parse(df[!, :off_lines][ind]))

   on_lines_li = []
   off_lines_li = []
   for i in keys(ref[:branch])
      br = (ref[:branch][i]["f_bus"], ref[:branch][i]["t_bus"]) 
      if br in z_bar
         push!(off_lines_li, br)
      else
         push!(on_lines_li, br)
      end
   end

   println("On lines: ", on_lines_li)
   println("Off lines: ", off_lines_li)

   # br_off_ind = bp_to_branch(z_bar, data)

   # for br in br_off_ind
   #    data["branch"][br]["br_status"] = 0
   # end

   return on_lines_li, off_lines_li
end

# Change an array of buspairs to branch indices
function bp_to_branch(bps, data)
   br_ind = []
   for i=1:length(data["branch"])
      branch = data["branch"]["$i"]

      f = branch["f_bus"]
      t = branch["t_bus"]

      if ((f, t) in bps) || ((t, f) in bps)
         push!(br_ind, "$i")
      end
   end
   return br_ind
end
