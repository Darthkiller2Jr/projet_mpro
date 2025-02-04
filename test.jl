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

x1,obj1 = simple_opt(V,A,t_hat,t,d,20,0)
visualize_solution(x1,obj1, coordinates, n, 20, 0)

x2,obj2 = simple_opt(V,A,t_hat,t,d,40,0)
visualize_solution(x2,obj2, coordinates, n, 40, 0)

x3,obj3 = simple_opt(V,A,t_hat,t,d,20,1)
visualize_solution(x3,obj3, coordinates, n, 20, 1)

x4,obj4 = simple_opt(V,A,t_hat,t,d,20,3)
visualize_solution(x4,obj4, coordinates, n, 20, 3)"""

"""include("data/n_10-euclidean_true")
x2,obj2 = simple_opt(V,A,th,t,d,C,0)
println(obj2)"""
i=1
obj_false=0
obj_true=0
"""for i in 5:20
    println("n = $(i)")
    include("data/n_$(i)-euclidean_false")
    x,obj_false = simple_opt(n,th,t,d,C,T)
    println("objectif (non euclidien): $obj_false")
    include("data/n_$(i)-euclidean_true")
    x,obj_true = simple_opt(n,th,t,d,C,T)
    println("objectif (euclidien): $obj_true")
end"""

#include("data/n_10-euclidean_true")
#print(robust_clark_wright(n,t,th,d,C))

"""max = true
include("data/n_20-euclidean_true")
routes_CW = robust_clark_wright(n, t, th, d, C, max)
println("Routes CW: ",routes_CW)
#println(routes_to_x(routes_CW,n))
println("Borne sup cout total CW : ", total_cost(routes_CW,t,th,max))
println("Cout réel de la sol CW: ", real_cost(routes_CW,n,th,t,T))
routes_LK = lin_kernighan_VRP(n, t, th, d, C, max; two_opt=false)
println("Routes LK: ",routes_LK)
#println(routes_to_x(routes_LK,n))
println("Borne sup cout total LK : ", total_cost(routes_LK,t,th,max))
println("Cout réel de la sol LK: ", real_cost(routes_LK,n,th,t,T))
routes_LK_2 = hybrid_heuristic(n, t, th, d, C, max; two_opt=false)
println("Routes LK^2: ",routes_LK_2)
println("Borne sup cout total LK entre routes : ", total_cost(routes_LK_2,t,th,max))
println("Cout réel de la sol LK entre routes : ", real_cost(routes_LK_2,n,th,t,T))
sol_opt = simple_opt(n,th,t,d,C,T;verbose=false)
println("Optimum réel : ",sol_opt[2])
println("Routes optimale : ", x_to_routes(sol_opt[1],n))"""


function comparaison_with_without_heuristique(file::String,euclidien::Bool;two_opt::Bool=true)
    println(file)
    include(file)
    max=true

    time_heuristic = @elapsed begin
        routes_LK_2 = hybrid_heuristic(n, t, th, d, C, max, euclidien; two_opt=two_opt)
    end
    #println("Routes LK^2: ",routes_LK_2)
    println("Cout réel de la sol heuristique : ", real_cost(routes_LK_2,n,th,t,T))
    println("Time of heuristic: ", time_heuristic, " seconds")

    heuristic_x = routes_to_x(routes_LK_2, n)
    # With heuristic
    #println("Solving with heuristic:")
    time_with_heuristic = @elapsed begin
        sol_opt_with_heuristic = simple_opt(n, th, t, d, C, T; heuristic_solution=heuristic_x, verbose=false)
        #lb = JuMP.objective_bound(m)
    end
    println("Optimum réel : ", sol_opt_with_heuristic[2])
    println("Time with heuristic: ", time_with_heuristic, " seconds")
    #println("Routes optimale (with heuristic): ", x_to_routes(sol_opt_with_heuristic[1], n))

    # Without heuristic
    #println("\nSolving without heuristic:")
    time_without_heuristic = @elapsed begin
        sol_opt_without_heuristic = simple_opt(n, th, t, d, C, T; verbose=false)
    end
    println("Time without heuristic: ", time_without_heuristic, " seconds")
    #println("Routes optimale (without heuristic): ", x_to_routes(sol_opt_without_heuristic[1], n))

    # Calculate speedup
    speedup = time_without_heuristic / time_with_heuristic
    println("Speedup: ", speedup, "x")
    println()
end

function comparaison_with_without_heuristique_all_i()
    for i in 5:15
        euclidien = true
        file = "data/n_$i-euclidean_$euclidien"
        comparaison_with_without_heuristique(file,euclidien)
        euclidien = false
        file = "data/n_$i-euclidean_$euclidien"
        comparaison_with_without_heuristique(file,euclidien)
    end
end
