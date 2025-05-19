# Comparing the ext and mc relaxations for cycle constraints.
using JuMP
using Gurobi
using CSV, DataFrames
# using LightGraphs
# import LightGraphs: cartesian_product
# using IterTools
using Printf
using Random, Distributions
using Plots
const GRB_ENV = Gurobi.Env()
path1 = "results/ext_bds.txt"
path2 = "results/rm_bds.txt"

nRounds = 5000
Random.seed!(123)
lb = -1.0
ub = 1.0
N = 6
pts = rand(Uniform(lb, ub), nRounds, N)
P = [(1,2), (1,3), (1,5), (1,6), (2,3),(2,4), (2,6),(3,4), (3,5),(4,5),(4,6),(5,6)]
function cartesian_product(lb, ub)
   v = []
   n = length(lb)
   for i in 1:n
      push!(v, (lb[i], ub[i]))
   end
   return collect(Base.product(v...)) # The splat operator, ..., takes out the entries in v.
end
X = cartesian_product(repeat([lb], N), repeat([ub], N))

function ext_gap(i)
   m = Model(() -> Gurobi.Optimizer(GRB_ENV))
   set_optimizer_attributes(m, "OutputFlag" => 0)

   # JuMP.@variable(m, lb <= x[k in 1:N] <= ub)
   JuMP.@variable(m, x_lift[p in P])
   JuMP.@variable(m, λ[l in 1:(2^N)] >= 0)

   JuMP.@objective(m, Max, sum(x_lift[p] for p in P))

   # @constraints(m, begin
   #      pts[i, 3] == x_lift[(1,2)] - x_lift[(4,5)]
   #      pts[i, 6] == x_lift[(1,5)] + x_lift[(2,4)]
   #      pts[i, 1] == x_lift[(2,3)] + x_lift[(5,6)]
   #      pts[i, 4] == x_lift[(2,6)] - x_lift[(3,5)]
   #      pts[i, 2] == x_lift[(1,3)] + x_lift[(4,6)]
   #      pts[i, 5] == x_lift[(1,6)] - x_lift[(3,4)]
   #      sum(λ[l] for l in 1:(2^N)) == 1
   # end)
   JuMP.@constraint(m, sum(λ[l] for l in 1:(2^N)) == 1)
   for k in 1:N
    expr_1 = 0
    expr_2 = 0
    cnt_1 = 0
    for l in 1:(2^N)
     if X[l][k] == lb
        expr_1 += λ[l]
        cnt_1 += 1
     else
        expr_2 += λ[l]
     end
    end
    JuMP.@constraint(m, pts[i, k] == lb * expr_1 + ub * expr_2)
   end
   JuMP.@constraint(m, [p in P], x_lift[p] == sum(λ[l] * X[l][p[1]] * X[l][p[2]] for l in 1:(2^N)))
   # print(m)
   JuMP.optimize!(m)
   max_obj = JuMP.objective_value(m)

   JuMP.@objective(m, Min, sum(x_lift[p] for p in P))
   JuMP.optimize!(m)
   min_obj = JuMP.objective_value(m)

   gap = max_obj - min_obj

   return gap
end

function rm_gap(i)
   m = Model(() -> Gurobi.Optimizer(GRB_ENV))
   set_optimizer_attributes(m, "OutputFlag" => 0)

   # JuMP.@variable(m, lb <= x[k in 1:N] <= ub)
   JuMP.@variable(m, x_lift[p in P])
   JuMP.@variable(m, λ[l in 1:(2^N)] >= 0)

   JuMP.@objective(m, Max, sum(x_lift[p] for p in P))

   function mcc(m,xy,x,y)
     JuMP.@constraint(m, xy >= lb*y + lb*x - lb*lb)
     JuMP.@constraint(m, xy >= ub*y + ub*x - ub*ub)
     JuMP.@constraint(m, xy <= lb*y + ub*x - lb*ub)
     JuMP.@constraint(m, xy <= ub*y + lb*x - ub*lb)
     return m
   end

   for p in P
      mcc(m, x_lift[p], pts[i, p[1]], pts[i, p[2]])
   end

   JuMP.optimize!(m)
   max_obj = JuMP.objective_value(m)

   JuMP.@objective(m, Min, sum(x_lift[p] for p in P))
   JuMP.optimize!(m)
   min_obj = JuMP.objective_value(m)

   gap = max_obj - min_obj

   return gap
end

ext_gap_li = []
rm_gap_li = []
for r in 1:nRounds
   push!(ext_gap_li, ext_gap(r))
   push!(rm_gap_li, rm_gap(r))
end

open(path1, "w") do f
   write(f, string(ext_gap_li))
end
open(path2, "w") do f
   write(f, string(rm_gap_li))
end

str1 = read(path1, String)
ext_gap_li = eval(Meta.parse(str1))
str2 = read(path2, String)
rm_gap_li = eval(Meta.parse(str2))

p = plot(ext_gap_li, rm_gap_li, seriestype = :scatter, legend = false, grid = false, xlim = (0, 20), ylim = (0, 20), xlabel = "Gap of extreme-point representation", guide = "Gap of McCormick relaxation", markersize = 2, markerstrokewidth = 0, markershape = :none, msc=:auto, tickfontsize = 10, guidefontsize = 15)#,markercolor = :blue , markerstrockcolor = :blue; msc automatically matching marker strock color.
plot!([0,20], [0,20], color = :lightgray)
savefig(p, "results/cyc_ext_rm_largefont.pdf")
