using LinearAlgebra

function cost(t::Matrix{Int},t_hat::Vector{Int},i::Int,j::Int,max::Bool)
    if max
        return t[i,j] + (t_hat[i] + t_hat[j]) + 2 * t_hat[j] * t_hat[i]
    else
        return t[i,j]
    end
end

function calc_savings(n::Int, t::Matrix{Int},t_hat::Vector{Int},max::Bool)
    savings = Vector{Tuple{Int, Int, Int}}()
    for i in 1:n-1
        for j in i+1:n
            # Utilisation des temps max
            savings_ij = cost(t,t_hat,i,j,max)
            push!(savings, (savings_ij, i, j))
        end
    end
    return savings
end

function find_route(routes::Vector{Vector{Int}}, customer::Int, deb::Bool)
    for route in routes
        if (customer == route[1] && deb) || (customer == route[end] && !deb)
            return route
        end
        if (customer == route[1] && !deb) || (customer == route[end] && deb)
            return reverse(route)
        end
    end
    #println("Pas de route trouvée")
    return nothing
end

function is_feasible(route_i::Vector{Int}, route_j::Vector{Int}, d::Vector{Int}, C::Int)
    # Check if the merged route satisfies capacity constraints
    merged_demand = sum(d[i] for i in route_i) + sum(d[j] for j in route_j)
    return merged_demand <= C
end

# Merge de 2 routes
function merge_routes(route_i::Vector{Int}, route_j::Vector{Int})
    return vcat(route_i, route_j)
end

function robust_clark_wright(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, max::Bool)
    # Initialisation : une route/client
    routes = [[i] for i in 1:n]

    # Calculs des savings
    savings = calc_savings(n, t, t_hat, max)

    # Tri des savings
    sort!(savings, by=x -> x[1], rev=true)

    # Merge des routes
    for (saving, i, j) in savings
        #println(i)
        #println(j)
        # Routes terminant/commencant par i/j
        route_i = find_route(routes, i, false)
        route_j = find_route(routes, j, true)
        #println(routes)
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

# Lin-Kernighan Procedure for TSP (applied to each route)
function lin_kernighan(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    n = length(route)
    improvement = true
    while improvement
        improvement = false
        for i in 1:n
            for j in i+1:n
                # Perform a 2-opt move
                new_route = two_opt_swap(route, i, j)
                # Calculate the cost of the new route
                # A améliorer
                current_cost = calculate_route_cost(route, t, t_hat, max)
                new_cost = calculate_route_cost(new_route, t, t_hat, max)
                # If the new route is better, accept it
                if new_cost < current_cost
                    route = new_route
                    improvement = true
                end
            end
        end
    end
    return route
end

# Helper function to perform a 2-opt swap
function two_opt_swap(route::Vector{Int}, i::Int, j::Int)
    # Reverse the segment of the route between i and j
    return vcat(route[1:i-1], reverse(route[i:j]), route[j+1:end])
end

# Calcul du coup d'un swap 2-opt
function swap_cost(route::Vector{Int},i::Int,j::Int,max::Bool)
    
end

# Helper function to calculate the cost of a route
function calculate_route_cost(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    route_cost = 0.0
    for k in 1:length(route)-1
        route_cost += cost(t, t_hat, route[k], route[k+1], max)
    end
    # Add the cost of returning to the depot
    route_cost += cost(t, t_hat, route[end], route[1], max)
    return route_cost
end

function total_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    c = 0.0
    for route in routes
        c+=calculate_route_cost(route, t, t_hat, max)
    end
    return c
end

# Main function to solve VRP using Clark and Wright + Lin-Kernighan
function solve_vrp(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, max::Int)
    # Clark and Wright pour sol initiale
    routes = robust_clark_wright(n, t, t_hat, d, C, max)
    println("Coût du meilleur routing après Clark Wright : ", total_cost(routes, t, t_hat, max))
    # Lin-Kernighen sur chaque tour
    optimized_routes = []
    for route in routes
        optimized_route = lin_kernighan(route, t, t_hat, true)
        push!(optimized_routes, optimized_route)
        print(optimized_route)
    end
    println("Coût d'une meilleure route après Lin-Kernighan : ", total_cost(optimized_routes, t, t_hat, max))
    return optimized_routes
end

include("data/n_10-euclidean_true")
optimized_routes = solve_vrp(n, t, th, d, C,true)
println("Optimized Routes: ", optimized_routes)