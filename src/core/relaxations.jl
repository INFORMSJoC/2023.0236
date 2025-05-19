function mcc(m,xy,x,y)
  lb = [lower_bound(x), lower_bound(y)]
  ub = [upper_bound(x), upper_bound(y)]
  JuMP.@constraint(m, xy >= lb[1]*y + lb[2]*x - lb[1]*lb[2])
  JuMP.@constraint(m, xy >= ub[1]*y + ub[2]*x - ub[1]*ub[2])
  JuMP.@constraint(m, xy <= lb[1]*y + ub[2]*x - lb[1]*ub[2])
  JuMP.@constraint(m, xy <= ub[1]*y + lb[2]*x - ub[1]*lb[2])
  return m
end

# No more need this function after defining cs and w for reversed arcs.
# McCormick relaxation for term xy at the arcpair _ap, _trigo_i (string) term indicate whether the ith term is a cosine function or a sine function.
# Used in the 4-cycle c-s constraints.
function mcc_correct_ap(m, xy, _ap, _trigo1, _trigo2)
   x = trigo_correct_ap(_ap[1], _trigo1)
   y = trigo_correct_ap(_ap[2], _trigo2)
   lb_x = 0
   lb_y = 0
   ub_x = 0
   ub_y = 0
   if _trigo1 == "cs" ||(_ap[1][1] < _ap[1][2])
      lb_x = lower_bound(x)
      ub_x = upper_bound(x)
   else
      lb_x = - upper_bound(si[(_ap[1][2], _ap[1][1])])
      ub_x = - lower_bound(si[(_ap[1][2], _ap[1][1])])
   end
   if _trigo2 == "cs" ||(_ap[2][1] < _ap[2][2])
      lb_y = lower_bound(y)
      ub_y = upper_bound(y)
   else
      lb_y = - upper_bound(si[(_ap[2][2], _ap[2][1])])
      ub_y = - lower_bound(si[(_ap[2][2], _ap[2][1])])
   end
   lb = [lb_x, lb_y]
   ub = [ub_x, ub_y]
   JuMP.@constraint(m, xy >= lb[1]*y + lb[2]*x - lb[1]*lb[2])
   JuMP.@constraint(m, xy >= ub[1]*y + ub[2]*x - ub[1]*ub[2])
   JuMP.@constraint(m, xy <= lb[1]*y + ub[2]*x - lb[1]*ub[2])
   JuMP.@constraint(m, xy <= ub[1]*y + lb[2]*x - ub[1]*lb[2])
   return m
end

# function mcc_on_off(m,xy,x,y,xz,yz,z)
#   lb = [lower_bound(x), lower_bound(y)]
#   ub = [upper_bound(x), upper_bound(y)]
#   JuMP.@variable(m, w)
#   JuMP.@constraint(m, xy >= lb[1]*yz + lb[2]*xz - lb[1]*lb[2]*z)
#   JuMP.@constraint(m, xy >= ub[1]*yz + ub[2]*xz - ub[1]*ub[2]*z)
#   JuMP.@constraint(m, xy <= lb[1]*yz + ub[2]*xz - lb[1]*ub[2]*z)
#   JuMP.@constraint(m, xy <= ub[1]*yz + lb[2]*xz - ub[1]*lb[2]*z)
# end

function quadratic_relax(m,xy,x)
  lb = lower_bound(x)
  ub = upper_bound(x)
  JuMP.@constraint(m, xy >= x^2)
  JuMP.@constraint(m, xy <= (lb + ub)*x - lb*ub)
  return m
end

function conv_tri(m,x_lift,x1,x2,x3,λ)
   lb = [lower_bound(x1), lower_bound(x2), lower_bound(x3)]
   ub = [upper_bound(x1), upper_bound(x2), upper_bound(x3)]
   ex_pt = [lb[1]*lb[2]*lb[3],
            lb[1]*ub[2]*lb[3],
            ub[1]*lb[2]*lb[3],
            ub[1]*ub[2]*lb[3],
            lb[1]*lb[2]*ub[3],
            lb[1]*ub[2]*ub[3],
            ub[1]*lb[2]*ub[3],
            ub[1]*ub[2]*ub[3]]
   # println("lb: ", lb)
   # println("ub: ", ub)
   # println("ex_pt: ", ex_pt)
   JuMP.@constraint(m, x_lift == sum(ex_pt[i]*λ[i] for i=1:8))
   JuMP.@constraint(m, x1 == (λ[1] + λ[2] + λ[5] + λ[6])*lb[1] + (λ[3] + λ[4] + λ[7] + λ[8])*ub[1])
   JuMP.@constraint(m, x2 == (λ[1] + λ[3] + λ[5] + λ[7])*lb[2] + (λ[2] + λ[4] + λ[6] + λ[8])*ub[2])
   JuMP.@constraint(m, x3 == (λ[1] + λ[2] + λ[3] + λ[4])*lb[3] + (λ[5] + λ[6] + λ[7] + λ[8])*ub[3])
   JuMP.@constraint(m, sum(λ[i] for i=1:8) == 1)

   return m
end

function conv_tri_on_off(m, x_lift, x1, x2, x3, λ, z, x3_bd)
   lb = [lower_bound(x1), lower_bound(x2), x3_bd[1]]
   ub = [upper_bound(x1), upper_bound(x2), x3_bd[2]]
   ex_pt = [lb[1]*lb[2]*lb[3],
            lb[1]*ub[2]*lb[3],
            ub[1]*lb[2]*lb[3],
            ub[1]*ub[2]*lb[3],
            lb[1]*lb[2]*ub[3],
            lb[1]*ub[2]*ub[3],
            ub[1]*lb[2]*ub[3],
            ub[1]*ub[2]*ub[3]]
   # println("lb: ", lb)
   # println("ub: ", ub)
   # println("ex_pt: ", ex_pt)
   JuMP.@constraint(m, x_lift == sum(ex_pt[i]*λ[i] for i=1:8))

   JuMP.@constraint(m, x1 <= (λ[1] + λ[2] + λ[5] + λ[6])*lb[1] + (λ[3] + λ[4] + λ[7] + λ[8])*ub[1] + ub[1] * (1 - z))
   JuMP.@constraint(m, x1 >= (λ[1] + λ[2] + λ[5] + λ[6])*lb[1] + (λ[3] + λ[4] + λ[7] + λ[8])*ub[1] + lb[1] * (1 - z))
   JuMP.@constraint(m, x2 <= (λ[1] + λ[3] + λ[5] + λ[7])*lb[2] + (λ[2] + λ[4] + λ[6] + λ[8])*ub[2] + ub[2] * (1 - z))
   JuMP.@constraint(m, x2 >= (λ[1] + λ[3] + λ[5] + λ[7])*lb[2] + (λ[2] + λ[4] + λ[6] + λ[8])*ub[2] + lb[2] * (1 - z))

   # JuMP.@constraint(m, x1 == λ_hat[1, 1] * lb[1] + λ_hat[1, 2] * ub[1] + (λ[1] + λ[2] + λ[5] + λ[6])*lb[1] + (λ[3] + λ[4] + λ[7] + λ[8])*ub[1])
   # JuMP.@constraint(m, x2 == λ_hat[2, 1] * lb[2] + λ_hat[2, 2] * ub[2] + (λ[1] + λ[3] + λ[5] + λ[7])*lb[2] + (λ[2] + λ[4] + λ[6] + λ[8])*ub[2])

   JuMP.@constraint(m, x3 == (λ[1] + λ[2] + λ[3] + λ[4])*lb[3] + (λ[5] + λ[6] + λ[7] + λ[8])*ub[3])
   JuMP.@constraint(m, sum(λ[i] for i = 1:8) == z)

   return m
end

function conv_tri_sum(m,x_lift,y_lift,x1,x2,x3,y3,λ_x,λ_y)
   # trilinear term x1*x2*x3 relaxation
   lb_x = [lower_bound(x1), lower_bound(x2), lower_bound(x3)]
   ub_x = [upper_bound(x1), upper_bound(x2), upper_bound(x3)]
   conv_tri(m, x_lift, x1, x2, x3, λ_x)

   # trilinear term x1*x2*y3 relaxation
   lb_y = [lower_bound(x1), lower_bound(x2), lower_bound(y3)]
   ub_y = [upper_bound(x1), upper_bound(x2), upper_bound(y3)]
   conv_tri(m, y_lift, x1, x2, y3, λ_y)

   # tying constraint
   val = [lb_x[1]*lb_x[2],
          lb_x[1]*ub_x[2],
          ub_x[1]*lb_x[2],
          ub_x[1]*ub_x[2],
          lb_x[1]*lb_x[2],
          lb_x[1]*ub_x[2],
          ub_x[1]*lb_x[2],
          ub_x[1]*ub_x[2]]
   JuMP.@constraint(m, sum(λ_x[i]*val[i] for i=1:8) == sum(λ_y[i]*val[i] for i=1:8))

   return m
end

# For ACOTS: x_lift -> w^R, y_lift -> w^I, x1 -> v_i, x2 -> v_j, x3 -> cs, y_3 -> si, λ_hat -> λ_hat_ij
function conv_tri_sum_on_off(m, x_lift, y_lift, x1, x2, x3, y3, λ_x, λ_y, z, x3_y3_bd)
   # trilinear term x1*x2*x3 relaxation
   lb_x = [lower_bound(x1), lower_bound(x2), x3_y3_bd[1][1]]
   ub_x = [upper_bound(x1), upper_bound(x2), x3_y3_bd[1][2]]

   # trilinear term x1*x2*y3 relaxation
   lb_y = [lower_bound(x1), lower_bound(x2), x3_y3_bd[2][1]]
   ub_y = [upper_bound(x1), upper_bound(x2), x3_y3_bd[2][2]]

   conv_tri_on_off(m, x_lift, x1, x2, x3, λ_x, z, x3_y3_bd[1])
   conv_tri_on_off(m, y_lift, x1, x2, y3, λ_y, z, x3_y3_bd[2])

   # tying constraint
   val = [lb_x[1]*lb_x[2],
          lb_x[1]*ub_x[2],
          ub_x[1]*lb_x[2],
          ub_x[1]*ub_x[2],
          lb_x[1]*lb_x[2],
          lb_x[1]*ub_x[2],
          ub_x[1]*lb_x[2],
          ub_x[1]*ub_x[2]]
   JuMP.@constraint(m, sum(λ_x[i]*val[i] for i=1:8) == sum(λ_y[i]*val[i] for i=1:8))

   return m
end


 # RLT cuts to capture (v_i*c_ij)*c_ij + (v_i*s_ij)*s_ij = v_i and
                     # (v_j*c_ij)*c_ij + (v_j*s_ij)*s_ij = v_j

function conv_tri_sum_rlt_1(m,x_lift,y_lift,x1,x2,x3,y3,λ_x,λ_y,x3_lift_i,x3_lift_j,y3_lift_i,y3_lift_j,λ_x3,λ_y3)
   # Additional RLT constraints to capture (v_i*c_ij)*c_ij + (v_i*s_ij)*s_ij = v_i and (v_j*c_ij)*c_ij + (v_j*s_ij)*s_ij = v_j
   for i=1:length(λ_x3)
      mcc(m, λ_x3[i], λ_x[i], x3)
      mcc(m, λ_y3[i], λ_y[i], y3)
   end

   lb = [lower_bound(x1), lower_bound(x2), lower_bound(x3)]
   ub = [upper_bound(x1), upper_bound(x2), upper_bound(x3)]
   ex_pt = [lb[1]*lb[3],
            lb[1]*lb[3],
            ub[1]*lb[3],
            ub[1]*lb[3],
            lb[1]*ub[3],
            lb[1]*ub[3],
            ub[1]*ub[3],
            ub[1]*ub[3]]
   JuMP.@constraint(m, x3_lift_i == sum(ex_pt[i]*λ_x3[i] for i=1:8))

   ex_pt = [lb[2]*lb[3],
            ub[2]*lb[3],
            lb[2]*lb[3],
            ub[2]*lb[3],
            lb[2]*ub[3],
            ub[2]*ub[3],
            lb[2]*ub[3],
            ub[2]*ub[3]]
   JuMP.@constraint(m, x3_lift_j == sum(ex_pt[i]*λ_x3[i] for i=1:8))

   lb = [lower_bound(x1), lower_bound(x2), lower_bound(y3)]
   ub = [upper_bound(x1), upper_bound(x2), upper_bound(y3)]
   ex_pt = [lb[1]*lb[3],
            lb[1]*lb[3],
            ub[1]*lb[3],
            ub[1]*lb[3],
            lb[1]*ub[3],
            lb[1]*ub[3],
            ub[1]*ub[3],
            ub[1]*ub[3]]
   JuMP.@constraint(m, y3_lift_i == sum(ex_pt[i]*λ_y3[i] for i=1:8))

   ex_pt = [lb[2]*lb[3],
            ub[2]*lb[3],
            lb[2]*lb[3],
            ub[2]*lb[3],
            lb[2]*ub[3],
            ub[2]*ub[3],
            lb[2]*ub[3],
            ub[2]*ub[3]]
   JuMP.@constraint(m, y3_lift_j == sum(ex_pt[i]*λ_y3[i] for i=1:8))

   JuMP.@constraint(m, x1  == x3_lift_i + y3_lift_i)
   JuMP.@constraint(m, x2  == x3_lift_j + y3_lift_j)

   return m
end


 # RLT cuts to capture wc_ij*c_ij + ws_ij*s_ij = vi*vj

function conv_tri_sum_rlt_2(m,x_lift,y_lift,x1,x2,x3,y3,λ_x,λ_y,x3_lift,y3_lift,λ_x3,λ_y3)

   lb_x = [lower_bound(x1), lower_bound(x2), lower_bound(x3)]
   ub_x = [upper_bound(x1), upper_bound(x2), upper_bound(x3)]
   lb_y = [lower_bound(x1), lower_bound(x2), lower_bound(y3)]
   ub_y = [upper_bound(x1), upper_bound(x2), upper_bound(y3)]
   val = [lb_x[1]*lb_x[2],
          lb_x[1]*ub_x[2],
          ub_x[1]*lb_x[2],
          ub_x[1]*ub_x[2],
          lb_x[1]*lb_x[2],
          lb_x[1]*ub_x[2],
          ub_x[1]*lb_x[2],
          ub_x[1]*ub_x[2]]

   # Additional RLT constraints to capture wc_ij*c_ij + ws_ij*s_ij = vi*vj
   for i=1:length(λ_x3)
      mcc(m, λ_x3[i], λ_x[i], x3)
      mcc(m, λ_y3[i], λ_y[i], y3)
   end

   lb = [lower_bound(x1), lower_bound(x2), lower_bound(x3)]
   ub = [upper_bound(x1), upper_bound(x2), upper_bound(x3)]
   ex_pt = [lb[1]*lb[2]*lb[3],
            lb[1]*ub[2]*lb[3],
            ub[1]*lb[2]*lb[3],
            ub[1]*ub[2]*lb[3],
            lb[1]*lb[2]*ub[3],
            lb[1]*ub[2]*ub[3],
            ub[1]*lb[2]*ub[3],
            ub[1]*ub[2]*ub[3]]
   JuMP.@constraint(m, x3_lift == sum(ex_pt[i]*λ_x3[i] for i=1:8))

   lb = [lower_bound(x1), lower_bound(x2), lower_bound(y3)]
   ub = [upper_bound(x1), upper_bound(x2), upper_bound(y3)]
   ex_pt = [lb[1]*lb[2]*lb[3],
            lb[1]*ub[2]*lb[3],
            ub[1]*lb[2]*lb[3],
            ub[1]*ub[2]*lb[3],
            lb[1]*lb[2]*ub[3],
            lb[1]*ub[2]*ub[3],
            ub[1]*lb[2]*ub[3],
            ub[1]*ub[2]*ub[3]]
   JuMP.@constraint(m, y3_lift == sum(ex_pt[i]*λ_y3[i] for i=1:8))

   # RLT constraint
   #mcc(m, x3_lift, x3, x_lift)
   #mcc(m, y3_lift, y3, y_lift)
   JuMP.@constraint(m, sum(λ_x[i]*val[i] for i=1:8) == x3_lift + y3_lift)

   return m
end

function relaxation_cos(m, cs, td, td_min, td_max, td_u)
   JuMP.@constraint(m, cs <= 1 - (1-cos(td_u))/(td_u*td_u)*(td^2))
   JuMP.@constraint(m, cs >= (cos(td_min) - cos(td_max))/(td_min-td_max)*(td - td_min) + cos(td_min))
   return m
end

function relaxation_sin(m, si, td, td_min, td_max, td_u)
   if td_min < 0 && td_max > 0
      JuMP.@constraint(m, si <= cos(td_u/2)*(td - td_u/2) + sin(td_u/2))
      JuMP.@constraint(m, si >= cos(td_u/2)*(td + td_u/2) - sin(td_u/2))
   end
   if td_max <= 0
      JuMP.@constraint(m, si <= (sin(td_min) - sin(td_max))/(td_min-td_max)*(td - td_min) + sin(td_min))
      JuMP.@constraint(m, si >= cos(td_u/2)*(td + td_u/2) - sin(td_u/2))
   end
   if td_min >= 0
      JuMP.@constraint(m, si <= cos(td_u/2)*(td - td_u/2) + sin(td_u/2))
      JuMP.@constraint(m, si >= (sin(td_min) - sin(td_max))/(td_min-td_max)*(td - td_min) + sin(td_min))
   end
   return m
end

function relaxation_cos_on_off(m, bp, td_min, td_max, td_u, td_M, cos_min, cos_max)
   i = bp[1]
   j = bp[2]
   cs = m[:cs][bp]
   td = m[:td][bp]
   z = m[:z][bp]
   va_i = m[:va][i]
   va_j = m[:va][j]
   @assert td_min >= -pi/2 && td_max <= pi/2
   A = (cos(td_max) - cos(td_min)) / (td_max - td_min)
   B = (1 - cos(td_u)) / (td_u)^2
   if params["ext_cons"]
      JuMP.@constraint(m, cs <= cos_max * z)
   end
   JuMP.@constraint(m, cs >= cos_min * z)
   # can this be integrated?
   #JuMP.JuMP.@constraint(m, y >= (cos(lb) - cos(ub))/(lb-ub)*(x - lb) + cos(lb))
   if params["ext_cons"]
      # (8a)
      JuMP.@constraint(m, - cs + A * td <= (A * td_min - cos(td_min)) * z + abs(A) * td_M * (1 - z))
   end

   # (8b)
   if params["ext_cons"]
      JuMP.@constraint(m, cs <= z - B * (td^2) + B * td_M^2 * (1-z))

      # S = {1} or {2}
      # if A >= 0
      #    JuMP.@constraint(m, - cs <= (A * td_max - cos(td_min) - A * td_min) * z)
      # else
      #    JuMP.@constraint(m, - cs <= (A * td_max - cos(td_min) - A * td_max) * z)
      # end
      #
      # JuMP.@constraint(m, A * td <= (A * td_max - cos(td_min) + cos_max) * z + abs(A) * td_M * (1 - z))
   else
      JuMP.@constraint(m, cs <= z - B * ((va_i - va_j)^2) + B * td_M^2 * (1-z))
   end

   return m
end

function relaxation_sin_on_off(m, bp, td_min, td_max, td_u, td_M, sin_min, sin_max)
    i = bp[1]
    j = bp[2]
    si = m[:si][bp]
    td = m[:td][bp]
    z = m[:z][bp]
    va_i = m[:va][i]
    va_j = m[:va][j]
    @assert td_min >= -pi/2 && td_max <= pi/2
    C = cos(td_u / 2)
    D = (sin(td_min) - sin(td_max)) / (td_min - td_max)
    if params["ext_cons"]
       JuMP.@constraint(m, si <= sin_max * z)
       JuMP.@constraint(m, si >= sin_min * z)
    end

    if td_min < 0 && td_max > 0
      if params["ext_cons"]
         # JuMP.@constraint(m, si <= cos(td_u/2)*(td - td_u/2) + sin(td_u/2))
         # JuMP.@constraint(m, si >= cos(td_u/2)*(td + td_u/2) - sin(td_u/2))
          JuMP.@constraint(m,  si - C * td <= (- C * td_u / 2 + sin(td_u / 2)) * z + C * td_M * (1 - z))
          JuMP.@constraint(m, -si + C * td <= (- C * td_u / 2 + sin(td_u / 2)) * z + C * td_M * (1 - z))
          # println("(- C * td_u/2 + sin(td_u / 2)): ", (- C * td_u/2 + sin(td_u / 2)))
          # println("C * td_M: ", C * td_M)
          # S = {1} or {2}:
          # JuMP.@constraint(m, si <= (- C * td_u / 2 + sin(td_u / 2) + C * td_max) * z)
          # JuMP.@constraint(m, - si <= (- C * td_u / 2 + sin(td_u / 2) - C * td_min) * z)
          # JuMP.@constraint(m, -C * td <= (- C * td_u / 2 + sin(td_u / 2) - sin_min) * z + (C * td_M) * (1 - z))
          # JuMP.@constraint(m, C * td <= (- C * td_u / 2 + sin(td_u / 2) + sin_max) * z + (C * td_M) * (1 - z))
       else
          # JuMP.@constraint(m, si <= cos(td_u/2) * (va_i - va_j - td_u/2) + sin(td_u/2))
          # JuMP.@constraint(m, si >= cos(td_u/2) * (va_i - va_j + td_u/2) - sin(td_u/2))
          JuMP.@constraint(m,  si - C * (va_i - va_j) <= (- C * td_u / 2 + sin(td_u / 2)) * z + C * td_M * (1 - z))
          JuMP.@constraint(m, -si + C * (va_i - va_j) <= (- C * td_u / 2 + sin(td_u / 2)) * z + C * td_M * (1 - z))
       end
    end
    if params["ext_cons"]
       if td_max <= 0
         JuMP.@constraint(m, si - D * td <= (-D * td_min + sin(td_min)) * z + D * td_M * (1 - z))
         JuMP.@constraint(m, -si + C * td <= (- C * td_u/2 + sin(td_u / 2)) * z + C * td_M * (1-z))
         JuMP.@constraint(m, si <= (- D * td_min + sin(td_min) + D * td_max) * z)
         JuMP.@constraint(m, -D * td <= (- D * td_min + sin(td_min) - sin_min) * z + (D * td_M) * (1 - z))
       end
       if td_min >= 0
         JuMP.@constraint(m,  si - C * td <= (- C * td_u / 2 + sin(td_u / 2)) * z + C * td_M * (1 - z))
         JuMP.@constraint(m, -si + D * td <= (D * td_min - sin(td_min)) * z + D * td_M * (1 - z))
         JuMP.@constraint(m, - si <= (D * td_min - sin(td_min) - D * td_min) * z)
         JuMP.@constraint(m, -D * td <= (D * td_min - sin(td_min) + sin_max) * z + (D * td_M) * (1 - z))
       end
    end

    # JuMP.JuMP.@constraint(m, y <= z*(sin(max_ad/2) + cos(max_ad/2)*max_ad/2))
    # JuMP.JuMP.@constraint(m, -y <= z*(sin(max_ad/2) + cos(max_ad/2)*max_ad/2))

    # JuMP.JuMP.@constraint(m, cos(max_ad/2)*x <= z*(sin(max_ad/2) - cos(max_ad/2)*max_ad/2 + sin(max_ad)) + (1-z)*(cos(max_ad/2)*M_x))
    # JuMP.JuMP.@constraint(m, -cos(max_ad/2)*x <= z*(sin(max_ad/2) - cos(max_ad/2)*max_ad/2 + sin(max_ad)) + (1-z)*(cos(max_ad/2)*M_x))
    return m
end

function relaxation_voltage_lnc(m, wr, wi, w_i, w_j, buspair)
   vfub = buspair["vm_fr_max"]
   vflb = buspair["vm_fr_min"]
   vtub = buspair["vm_to_max"]
   vtlb = buspair["vm_to_min"]
   tdub = buspair["angmax"]
   tdlb = buspair["angmin"]

   phi = (tdub + tdlb)/2
   d   = (tdub - tdlb)/2

   sf = vflb + vfub
   st = vtlb + vtub

   # Voltage Product Relaxation Upper-bound
   JuMP.@constraint(m, sf*st*(cos(phi)*wr + sin(phi)*wi) - vtub*cos(d)*st*w_i - vfub*cos(d)*sf*w_j >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub))
   JuMP.@constraint(m, sf*st*(cos(phi)*wr + sin(phi)*wi) - vtlb*cos(d)*st*w_i - vflb*cos(d)*sf*w_j >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub))
end

function relaxation_voltage_lnc_on_off(m, wr, wi, w_i, w_j, buspair, z, wz_ij, wz_ji)
   vfub = buspair["vm_fr_max"]
   vflb = buspair["vm_fr_min"]
   vtub = buspair["vm_to_max"]
   vtlb = buspair["vm_to_min"]
   tdub = buspair["angmax"]
   tdlb = buspair["angmin"]

   phi = (tdub + tdlb)/2
   d   = (tdub - tdlb)/2

   sf = vflb + vfub
   st = vtlb + vtub

   # Voltage Product Relaxation Upper-bound
   JuMP.@constraint(m, sf*st*(cos(phi)*wr + sin(phi)*wi) - vtub*cos(d)*st*wz_ij - vfub*cos(d)*sf*wz_ji >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub) * z)
   JuMP.@constraint(m, sf*st*(cos(phi)*wr + sin(phi)*wi) - vtlb*cos(d)*st*wz_ij - vflb*cos(d)*sf*wz_ji >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub) * z)
end

# Extreme point constraints for a list of variables x_li and their bilinear combination as listed in expairs.
function conv_bi_x_li(m, cyc, x_li, expairs, params, x_type, x_li_bar=:x_li_bar)
   λ = nothing
   x = nothing
   if x_type == "cs"
      λ = m[:λ_c]
      x = m[:x_c]
   elseif x_type == "w"
      λ = m[:λ_w]
      x = m[:x_w]
   end
   lb = [lower_bound(var) for var in x_li]
   ub = [upper_bound(var) for var in x_li]
   X = cartesian_product(lb,ub)
   dim = length(x_li)
   # print(X)
   beta = []
   for j in 1:dim
     expr_1 = 0
     expr_2 = 0
     cnt_1 = 0
     for i in 1:(2^dim)
      # println("lb: $(lb[j]), ub: $(ub[j])")
      if X[i][j] == lb[j]
         # if X[i, j] == lb[j]
         expr_1 += λ[cyc, i]
         cnt_1 += 1
      else
         expr_2 += λ[cyc, i]
      end
     end
     # println("cnt_1: $(cnt_1), total: $(2^dim)")
     # @assert cnt_1 == 2^dim / 2
     if params["separation"]
        push!(beta, JuMP.@constraint(m, lb[j] * expr_1 + ub[j] * expr_2 == x_li_bar[j]))
     else
        JuMP.@constraint(m, x_li[j] == lb[j] * expr_1 + ub[j] * expr_2)
     end
   end
   if x_type == "w" && length(cyc) == 4 && params["rotate_4w"]
      expairs_1 = expairs[1:8]
      expairs_2 = expairs[9:end]
      JuMP.@constraint(m, [ep in expairs_1], x[cyc, ep] == sum(λ[cyc, i] * X[i][ep[1]] * X[i][ep[2]] for i in 1:(2^dim)))
      JuMP.@constraint(m, [ep in expairs_2], x[cyc, ep] == sum(λ[cyc, i] * X[i][ep[1]] * X[i][ep[2]] * X[i][ep[3]] for i in 1:(2^dim)))
   else
      JuMP.@constraint(m, [ep in expairs], x[cyc, ep] == sum(λ[cyc, i] * X[i][ep[1]] * X[i][ep[2]] for i in 1:(2^dim)))
   end
   if params["separation"]
      return beta
   end
end

function conv_bi_x_li_on_off(m, cyc, x_li, ub, lb, ub2, lb2, expairs, params, x_type, x_li_bar=:x_li_bar)
   λ = nothing
   x = nothing
   if x_type == "cs"
      λ = m[:λ_c]
      x = m[:x_c]
   elseif x_type == "w"
      λ = m[:λ_w]
      x = m[:x_w]
   end
   X = cartesian_product(lb,ub)
   dim = length(x_li)
   # print(X)
   beta = []
   for j in 1:dim
     expr_1 = 0
     expr_2 = 0
     cnt_1 = 0
     for i in 1:(2^dim)
      # println("lb: $(lb[j]), ub: $(ub[j])")
      if X[i][j] == lb[j]
         # if X[i, j] == lb[j]
         expr_1 += λ[cyc, i]
         cnt_1 += 1
      else
         expr_2 += λ[cyc, i]
      end
     end
     # println("cnt_1: $(cnt_1), total: $(2^dim)")
     # @assert cnt_1 == 2^dim / 2
     JuMP.@constraint(m, x_li[j] >= lb[j] * expr_1 + ub[j] * expr_2 + lb2[j] * (1 - m[:zc][Tuple(cyc)]))
     JuMP.@constraint(m, x_li[j] <= lb[j] * expr_1 + ub[j] * expr_2 + ub2[j] * (1 - m[:zc][Tuple(cyc)]))
   end
   if x_type == "w" && length(cyc) == 4 && params["rotate_4w"]
      expairs_1 = expairs[1:8]
      expairs_2 = expairs[9:end]
      JuMP.@constraint(m, [ep in expairs_1], x[cyc, ep] == sum(λ[cyc, i] * X[i][ep[1]] * X[i][ep[2]] for i in 1:(2^dim)))
      JuMP.@constraint(m, [ep in expairs_2], x[cyc, ep] == sum(λ[cyc, i] * X[i][ep[1]] * X[i][ep[2]] * X[i][ep[3]] for i in 1:(2^dim)))
   else
      JuMP.@constraint(m, [ep in expairs], x[cyc, ep] == sum(λ[cyc, i] * X[i][ep[1]] * X[i][ep[2]] for i in 1:(2^dim)))
   end
end
