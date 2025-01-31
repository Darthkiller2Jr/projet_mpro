using JuMP, CPLEX, Plots

function simple_opt(V,A,t_hat,t,d,C,T)

    m = Model(CPLEX.Optimizer) # Define the model

    @variable(m, theta)
    @variable(m, x[(i, j) in A], Bin)
    @variable(m, u[i in V])

    # Dualisation pour l'ensemble d'incertitude polyhedrique
    @variable(m,eta_1)
    @variable(m,eta_2)
    @variable(m, eta_3[(i,j) in A])
    @variable(m, eta_4[(i,j) in A])

    # Contraintes du dual de l'incertitude sur R
    @constraint(m, [(i,j) in A], eta_1 + eta_3[(i,j)] >= (t_hat[i] + t_hat[j]) * x[(i,j)])
    @constraint(m, [(i,j) in A], eta_2 + eta_4[(i,j)] >= (t_hat[i] * t_hat[j]) * x[(i,j)])
    @constraint(m,[(i,j) in A], eta_3[(i,j)]>=0)
    @constraint(m,[(i,j) in A], eta_4[(i,j)]>=0)

    # sup en dessous de theta - sum(t[i,j]*t[i,j])
    @constraint(m, eta_1 * T + eta_2 * T^2 + sum(eta_3[(i,j)] + 2 * eta_4[(i,j)] for (i,j) in A) <= theta - sum(t[i,j]*x[(i,j)] for (i,j) in A))    
    
    # Contraintes du problème de base
    @constraint(m, [j in V; j != 1], sum(x[(i, j)] for i in V if i != j) == 1)
    @constraint(m, [i in V; i != 1], sum(x[(i, j)] for j in V if i != j) == 1)
    @constraint(m,[(i,j) in A; i != 1 && d[i] + d[j] <= C], u[j]-u[i]>=d[j]-C*(1-x[(i,j)]))
    @constraint(m,[i in V; i!=1], u[i]<=C)
    @constraint(m,[i in V; i!=1], u[i]>=d[i])
    
    @objective(m,Min,theta)

    optimize!(m)

    return value.(x), JuMP.objective_value(m)
end

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