using JuMP, CPLEX, LinearAlgebra

function simple_opt(n::Int,t_hat::Vector{Int},t::Matrix{Int},d::Vector{Int},C::Int,T::Int;verbose=false)

    m = Model(CPLEX.Optimizer) # Define the model
    if !verbose
        set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
    end
    set_optimizer_attribute(m, "CPX_PARAM_TILIM", 30)
    
    # param
    V = 1:n
    A = [(i,j) for i in V for j in V if i!=j]

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

function savings_temps_max(n::Int, t::Matrix{Int},t_hat::Vector{Int})
    savings = Vector{Tuple{Int, Int, Int}}()
    for i in 1:n-1
        for j in i+1:n
            # Utilisation des temps max
            savings_ij = t[i,j] + (t_hat[i] + t_hat[j]) + 2 * t_hat[j] * t_hat[i]
            push!(savings, (savings_ij, i, j))
        end
    end
    return savings
end

function savings_temps_nominal(n::Int, t::Matrix{Int})
    savings = Vector{Tuple{Int, Int, Int}}()
    for i in 1:n-1
        for j in i+1:n
            # Utilisation des temps nominaux
            savings_ij = t[i,j]
            push!(savings, (savings_ij, i, j))
        end
    end
    return savings
end

function robust_clark_wright(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int)
    # Initialisation : une route/client
    routes = [[i] for i in 1:n]
    
    # Calculs des savings
    savings = savings_temps_nominal(n, t)
    
    # Tri des savings
    sort!(savings, by=x -> x[1], rev=true)
    
    # Merge des routes
    for (saving, i, j) in savings
        println(i)
        println(j)
        # Routes terminant/commencant par i/j
        route_i = find_route(routes, i, false)
        route_j = find_route(routes, j, true)
        println(routes)
        #println(route_i)
        #println(route_j)
        
        # Check if merging is feasible
        if route_i !== nothing && route_j !== nothing && route_i != route_j && is_feasible(route_i, route_j, d, C)
            # Merge routes
            merged_route = merge_routes(route_i, route_j)
            filter!(x -> x != route_i, routes)
            filter!(x -> x != route_j, routes)
            filter!(x -> x != reverse(route_i), routes)
            filter!(x -> x != reverse(route_j), routes)
            push!(routes, merged_route)
        end
    end
    
    return routes
end

function find_route(routes, customer::Int, deb::Bool)
    for route in routes
        if (customer == route[1] && deb) || (customer == route[end] && !deb)
            return route
        end
        if (customer == route[1] && !deb) || (customer == route[end] && deb)
            return reverse(route)
        end
    end
    println("Pas de route trouvée")
    return nothing
end

function is_feasible(route_i, route_j, d, C)
    # Check if the merged route satisfies capacity constraints
    merged_demand = sum(d[i] for i in route_i) + sum(d[j] for j in route_j)
    return merged_demand <= C
end

# Merge de 2 routes
function merge_routes(route_i, route_j)
    return vcat(route_i, route_j)
end