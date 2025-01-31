using JuMP, CPLEX, Plots
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
visualize_solution(x4,obj4, coordinates, n, 20, 3)