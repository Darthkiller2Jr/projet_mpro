using JuMP, CPLEX, Plots, BenchmarkTools
include("optim.jl")

function visualize_solution(x, obj, coordinates, n, C, T)
    # Create a new plot for each solution
    filename = "fig/Sol_n$(n)_C$(C)_T$(T).png"
    plt = scatter(coordinates[:, 1], coordinates[:, 2], legend=false, title="VRP Solution: C=$C, T=$T, obj=$obj", xlabel="X", ylabel="Y", markersize=8)

    # Draw the arcs (paths) that are selected in the solution
    for i in 1:n
        for j in 1:n
            if i != j && x[(i, j)] == 1
                plot!([coordinates[i, 1], coordinates[j, 1]], [coordinates[i, 2], coordinates[j, 2]], color=:blue, linewidth=2, label="")
            end
        end
    end

    # Display the plot
    savefig(plt, filename)
end

"""
# Define sets and parameters
n = 10  # Number of vertices
V = 1:n  # Set of vertices (adjust n as needed)
A = [(i, j) for i in V, j in V if i != j]  # Set of arcs
t_hat = rand(n)*10  # Example nominal travel times (replace with actual data)

coordinates = rand(n, 2) * 10

# Initialize travel time matrix t
t = fill(0.0, n, n)

# Compute t[i, j] as the Euclidean distance between point i and point j
for (i, j) in A
    t[i, j] = sqrt((coordinates[i, 1] - coordinates[j, 1])^2 + (coordinates[i, 2] - coordinates[j, 2])^2)
end

d = rand(n)*10  # Example demands (replace with actual data)

x1,obj1 = robust_opt(V,A,t_hat,t,d,20,0)
visualize_solution(x1,obj1, coordinates, n, 20, 0)

x2,obj2 = robust_opt(V,A,t_hat,t,d,40,0)
visualize_solution(x2,obj2, coordinates, n, 40, 0)

x3,obj3 = robust_opt(V,A,t_hat,t,d,20,1)
visualize_solution(x3,obj3, coordinates, n, 20, 1)

x4,obj4 = robust_opt(V,A,t_hat,t,d,20,3)
visualize_solution(x4,obj4, coordinates, n, 20, 3)"""

"""include("data/n_10-euclidean_true")
x2,obj2 = robust_opt(V,A,th,t,d,C,0)
println(obj2)"""
i=1
obj_false=0
obj_true=0
"""for i in 5:20
    println("n = $(i)")
    include("data/n_$(i)-euclidean_false")
    x,obj_false = robust_opt(n,th,t,d,C,T)
    println("objectif (non euclidien): $obj_false")
    include("data/n_$(i)-euclidean_true")
    x,obj_true = robust_opt(n,th,t,d,C,T)
    println("objectif (euclidien): $obj_true")
end"""

#include("data/n_10-euclidean_true")
#print(robust_clark_wright(n,t,th,d,C))

max_cost = true
"""include("data/n_10-euclidean_true")
routes_CW = robust_clark_wright(n, t, th, d, C, max_cost=max_cost)
println("Routes CW: ",routes_CW)
#println(routes_to_x(routes_CW,n))
println("Borne sup cout total CW : ", total_cost(routes_CW,t,th,max_cost=max_cost))
println("Cout réel de la sol CW: ", real_cost(routes_CW,n,th,t,T))
routes_LK = lin_kernighan_VRP(n, t, th, d, C, max_cost=max_cost, two_opt=false)
println("Routes LK: ",routes_LK)
#println(routes_to_x(routes_LK,n))
println("Borne sup cout total LK : ", total_cost(routes_LK,t,th,max_cost=max_cost))
println("Cout réel de la sol LK: ", real_cost(routes_LK,n,th,t,T))
routes_LK_2 = hybrid_heuristic(n, t, th, d, C, max_cost; two_opt=false)
println("Routes LK^2: ",routes_LK_2)
println("Borne sup cout total LK entre routes : ", total_cost(routes_LK_2,t,th,max_cost=max_cost))
println("Cout réel de la sol LK entre routes : ", real_cost(routes_LK_2,n,th,t,T))
sol_opt = robust_opt(n,th,t,d,C,T;verbose=false)
println("Optimum réel : ",sol_opt[2])
println("Routes optimale : ", x_to_routes(sol_opt[1],n))"""


function comparaison_with_without_heuristique(n::Int, euclidien::Bool;two_opt::Bool=true)
    file = "data/n_$n-euclidean_$euclidien"
    println(file)
    include(file)
    max_cost=true

    time_heuristic = @elapsed begin
        routes_LK_2 = hybrid_heuristic(n, t, th, d, C, euclidien;max_cost=max_cost, two_opt=two_opt)
    end
    #println("Routes LK^2: ",routes_LK_2)
    println("Cout réel de la sol heuristique : ", real_cost(routes_LK_2,n,th,t,T))
    println("Time of heuristic: ", time_heuristic, " seconds")

    heuristic_x = routes_to_x(routes_LK_2, n)
    # With heuristic
    #println("Solving with heuristic:")
    sol_opt_with_heuristic = robust_opt(n, th, t, d, C, T; heuristic_solution=heuristic_x, verbose=false)
    Gap = (sol_opt_with_heuristic[2]-sol_opt_with_heuristic[3])/sol_opt_with_heuristic[2]
    println(sol_opt_with_heuristic[3])
    println("Meilleure solution du solver : ", sol_opt_with_heuristic[2], " (Gap : $Gap)")
    println("Temps avec heuristique: ", sol_opt_with_heuristic[4], " seconds")
    #println("Routes optimale (with heuristic): ", x_to_routes(sol_opt_with_heuristic[1], n))

    # Without heuristic
    #println("\nSolving without heuristic:")
    sol_opt_without_heuristic = robust_opt(n, th, t, d, C, T; verbose=false)
    Gap_without_heuristique = (sol_opt_without_heuristic[2]-sol_opt_without_heuristic[3])/sol_opt_without_heuristic[2]
    println(sol_opt_without_heuristic[3])
    println("Meilleure solution du solver : ", sol_opt_without_heuristic[2], " (Gap : $Gap_without_heuristique)")
    println("Temps avec heuristique: ", sol_opt_without_heuristic[4], " seconds")
    #println("Routes optimale (without heuristic): ", x_to_routes(sol_opt_without_heuristic[1], n))

    # Calculate speedup
    speedup = sol_opt_without_heuristic[4] / sol_opt_with_heuristic[4]
    println("Speedup : ", speedup, "x")
    Gap_increase = (Gap_without_heuristique - Gap)/Gap_without_heuristique
    println("Gap increase : ", Gap_increase,"\n")
end

function comparaison_with_without_heuristique_all_i()
    for i in 5:20
        euclidien = true
        comparaison_with_without_heuristique(i,euclidien)
        euclidien = false
        comparaison_with_without_heuristique(i,euclidien)
    end
end

function test_heuristiques(n::Int, euclidien::Bool, max_cost::Bool; time_limit::Float64=10.0)
    filename = "data/n_$n-euclidean_$euclidien"
    if isfile(filename)
        include(filename)
        println("File: ", filename)

        # Timing robust_clark_wright
        cw_time = @elapsed begin 
            routes_CW = robust_clark_wright(n, t, th, d, C, max_cost=max_cost)
        end
        println("robust_clark_wright : $(real_cost_smart(routes_CW, n, th, t, T)) ($(cw_time)s).")

        # Timing sous_tours_heuristic with LK (LKH)
        lk_time = @elapsed begin
            routes_LK = sous_tours_heuristic(n, t, th, d, C, max_cost, euclidien; two_opt=false, LK=true)
        end
        println("sous_tours_heuristic (LK) : $(real_cost_smart(routes_LK, n, th, t, T)) ($(lk_time)s).")

        # Timing sous_tours_heuristic with 2-opt
        twoopt_time = @elapsed begin
            route_2opt = sous_tours_heuristic(n, t, th, d, C, max_cost, euclidien; two_opt=true, LK=false)
        end
        println("sous_tours_heuristic (2-opt) : $(real_cost_smart(route_2opt, n, th, t, T)) ($(twoopt_time)s).")
        
        # Timing sous_tours_heuristic with 3-opt
        threeopt_time = @elapsed begin
            route_3opt = sous_tours_heuristic(n, t, th, d, C, max_cost, euclidien; two_opt=false, LK=false)
        end
        println("sous_tours_heuristic (3-opt) : $(real_cost_smart(route_3opt, n, th, t, T)) ($(threeopt_time)s).")

        # Timing hybrid_heuristic: 3-opt on sub-tours + swap 2-opt between tours (max cost)
        hybrid_swap2_max_time = @elapsed begin
            routes_3opt_swap_max_cost = hybrid_heuristic(n, t, th, d, C, euclidien; max_cost=max_cost, LK=false, real_cost=false, two_opt=false)
        end
        println("hybrid_heuristic (3-opt + swap 2-opt, max cost) : $(real_cost_smart(routes_3opt_swap_max_cost, n, th, t, T)) ($(hybrid_swap2_max_time)s).")

        # Timing hybrid_heuristic: 3-opt on sub-tours + swap 2-opt between tours (real cost)
        hybrid_swap2_real_time = @elapsed begin
            routes_3opt_swap_real_cost = hybrid_heuristic(n, t, th, d, C, euclidien; max_cost=max_cost, LK=false, real_cost=true, two_opt=false)
        end
        println("hybrid_heuristic (3-opt + swap 2-opt, real cost) : $(real_cost_smart(routes_3opt_swap_real_cost, n, th, t, T)) ($(hybrid_swap2_real_time)s).")

        # Timing hybrid_heuristic: 3-opt on sub-tours + swap 3-opt between tours (real cost)
        #hybrid_swap3_real_time = @elapsed begin
        #    routes_3opt_swap3_real_cost = hybrid_heuristic(n, t, th, d, C, euclidien; max_cost=max_cost, LK=false, real_cost=true, two_opt=false, three_opt_swap=true)
        #end
        #println("hybrid_heuristic (3-opt + swap 3-opt, real cost) : $(real_cost_smart(routes_3opt_swap3_real_cost, n, th, t, T)) ($(hybrid_swap3_real_time)s).")

        # lk_swap_time = @elapsed begin
        #     routes_LK_swap_max_cost = hybrid_heuristic(n, t, th, d, C, euclidien; max_cost=max_cost, LK=true, real_cost=false, two_opt=false)
        # end
        # println("hybrid_heuristic (LK + swap 2-opt, max cost) took $(lk_swap_time) seconds.")
        # println("Objectif la sol explo LK sur les sous-tours et swap 2opt entre tours: ", real_cost_smart(routes_LK_swap_max_cost, n, th, t, T))

        # Timing the MILP optimization
        milp_time = @elapsed begin
            sol_opt = robust_opt(n, th, t, d, C, T, time_limit=time_limit)
        end
        println("Meilleure sol par MILP : ", sol_opt[2], ", Meilleure borne inf : ", sol_opt[3], " ($(milp_time)s)")
        
        # MILP avec warm start
        milp_WS_time = @elapsed begin
            heuristic_x = routes_to_x(routes_3opt_swap_real_cost, n)
            sol_opt_WS = robust_opt(n, th, t, d, C, T, heuristic_solution=heuristic_x, time_limit=time_limit)
        end
        println("Meilleure sol par MILP (warm start): ", sol_opt_WS[2], ", Meilleure borne inf : ", sol_opt_WS[3], " ($(milp_WS_time)s)\n")
    end
end

function test_heuristiques_all(euclidien::Bool)
    
    for i in 5:100
        filename = "data/n_$i-euclidean_$euclidien"
        if !isfile(filename)
            continue
        end
        println(i)
        include(filename)
        
        routes_LK_2_real_cost = hybrid_heuristic(n, t, th, d, C, euclidien; real_cost=true, two_opt=false)
        println("Cout réel de la sol LK entre routes (cout réel): ", real_cost(routes_LK_2_real_cost,n,th,t,T))
        
        routes_LK_2 = hybrid_heuristic(n, t, th, d, C, euclidien; real_cost=false, two_opt=false)
        println("Cout réel de la sol LK entre routes (cout max_cost): ", real_cost(routes_LK_2,n,th,t,T),"\n")
    
    end
end

for n in 40:40
    test_heuristiques(n,true,true,time_limit=2.0)
end

n = 5
euclidien = true
filename = "data/n_$n-euclidean_$euclidien"
include(filename)
relache = relache_continu(n, th, t, d, C, T, verbose=false)
println("Meilleure solution du solver : ", relache[2])
println("Temps : ", relache[4], " seconds")