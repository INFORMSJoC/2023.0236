va = vars["va"]
vm = vars["vm"]
w  = vars["w"]
wr = vars["wr"]
wi = vars["wi"]
pg = vars["pg"]
qg = vars["qg"]
p  = vars["p"]
q  = vars["q"]
td = vars["td"]
cs = vars["cs"]
si = vars["si"]
cm = vars["cm"]

if params["relax"] == "tri"
   λ_wr = vars["λ_wr"]
   λ_wi = vars["λ_wi"]
   # if (params["model"] == "ots") || (params["model"] == "ots_relax")
   #    λ_hat = vars["λ_hat"]
   # end
   if params["rlt_cuts"]
      wcc = vars["wcc"]
      wss = vars["wss"]
      vc_i = vars["vc_i"]
      vc_j = vars["vc_j"]
      vs_i = vars["vs_i"]
      vs_j = vars["vs_j"]
      c_λ_wr = vars["c_λ_wr"]
      s_λ_wi = vars["s_λ_wi"]
   end
elseif params["relax"] == "rmc"
   vv = vars["vv"]
   if params["model"] == "ots" || params["model"] == "ots_relax"
      vmz = vars["vmz"]
      # vvz = vars["vvz"]
      # csz = vars["csz"]
      # siz = vars["siz"]
   end
end

if (params["model"] == "ots") || (params["model"] == "ots_relax")
#    zw_fr = vars["zw_fr"]
#    zw_to = vars["zw_to"]
   z = vars["z"]
   wz = vars["wz"]
   if params["cycle_cuts"]
      zc = vars["zc"]
   end
end

if params["cycle_cuts"] && params["separation"] == false
   if params["cycle_relax"] == "mc"
      if params["cycle_c_s_cuts"]
         hcc = vars["hcc"]
         hss = vars["hss"]
         hcs = vars["hcs"]
         hsc = vars["hsc"]
      end
      if params["cycle_wr_wi_cuts"]
         wwr = vars["wwr"]
         wwi = vars["wwi"]
         wrr = vars["wrr"]
         wii = vars["wii"]
         wri = vars["wri"]
         wir = vars["wir"]
         wwrr = vars["wwrr"]
         wwii = vars["wwii"]
         wwri = vars["wwri"]
         wwir = vars["wwir"]
      end
   elseif params["cycle_relax"] == "epr"
      if params["cycle_c_s_cuts"]
         λ_c = vars["λ_c"]
         x_c = vars["x_c"]
      end
      if params["cycle_wr_wi_cuts"]
         λ_w = vars["λ_w"]
         x_w = vars["x_w"]
      end
   end
end
