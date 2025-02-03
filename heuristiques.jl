using JuMP, CPLEX

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
    routes = [[i] for i in 2:n]

    # Calculs des savings
    savings = calc_savings(n, t, t_hat, max)

    # Tri des savings
    sort!(savings, by=x -> x[1], rev=true)

    # Merge des routes
    for (saving, i, j) in savings

        # Routes terminant/commencant par i/j
        route_i = find_route(routes, i, false)
        route_j = find_route(routes, j, true)
        
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
function lin_kernighan_one_route(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool; two_opt::Bool=false)
    n = length(route)
    improvement = true
    #println(route,calculate_route_cost(route,t,t_hat,max))
    while improvement
        improvement = false
        for i in 1:n-2
            for j in i+1:n-1
                if two_opt
                    swap_cost = swap_2opt_cost(route,i,j,t,t_hat,max)
                    if swap_cost<0
                        route = two_opt_swap(route, i, j)
                        improvement = true
                        #println(route,calculate_route_cost(route,t,t_hat,max))
                    end
                else
                    for k in j+1:n
                        new_route, improvement = three_opt_swap_best(route, i, j, k, t, t_hat, max)
                        if improvement
                            route = new_route
                            #println(route, calculate_route_cost(route,t,t_hat,max))
                        end
                    end
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

function three_opt_swap_best(route::Vector{Int}, i::Int, j::Int, k::Int, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    
    if i == 1
        end_segment_1 = 1
    else
        end_segment_1 = route[i-1]
    end
    if k == length(route)
        deb_segment_4 = 1
    else
        deb_segment_4 = route[k+1]
    end

    # Define the segments
    segment1 = route[1:i-1]
    segment2 = route[i:j]
    segment3 = route[j+1:k]
    segment4 = route[k+1:end]
    
    #println("Segments:", segment1, segment2, segment3, segment4)

    # Original route cost
    original_cost = cost(t, t_hat, end_segment_1, route[i], max) +
                cost(t, t_hat, route[j], route[j+1], max) +
                cost(t, t_hat, route[k], deb_segment_4, max)

    # All possible reconnections and their costs
    reconnections = [
        (route, original_cost),  # Original route
        (vcat(segment1, reverse(segment2), segment3, segment4), 
        cost(t, t_hat, end_segment_1, segment2[end], max) + cost(t, t_hat, segment2[1], segment3[1], max) + cost(t, t_hat, segment3[end], deb_segment_4, max)),
        (vcat(segment1, segment2, reverse(segment3), segment4), 
        cost(t, t_hat, end_segment_1, segment2[1], max) + cost(t, t_hat, segment2[end], segment3[end], max) + cost(t, t_hat, segment3[1], deb_segment_4, max)),
        (vcat(segment1, reverse(segment2), reverse(segment3), segment4), 
        cost(t, t_hat, end_segment_1, segment2[end], max) + cost(t, t_hat, segment2[1], segment3[end], max) + cost(t, t_hat, segment3[1], deb_segment_4, max)),
        (vcat(segment1, segment3, segment2, segment4), 
        cost(t, t_hat, end_segment_1, segment3[1], max) + cost(t, t_hat, segment3[end], segment2[1], max) + cost(t, t_hat, segment2[end], deb_segment_4, max)),
        (vcat(segment1, segment3, reverse(segment2), segment4), 
        cost(t, t_hat, end_segment_1, segment3[1], max) + cost(t, t_hat, segment3[end], segment2[end], max) + cost(t, t_hat, segment2[1], deb_segment_4, max)),
        (vcat(segment1, reverse(segment3), segment2, segment4), 
        cost(t, t_hat, end_segment_1, segment3[end], max) + cost(t, t_hat, segment3[1], segment2[1], max) + cost(t, t_hat, segment2[end], deb_segment_4, max)),
        (vcat(segment1, reverse(segment3), reverse(segment2), segment4), 
        cost(t, t_hat, end_segment_1, segment3[end], max) + cost(t, t_hat, segment3[1], segment2[end], max) + cost(t, t_hat, segment2[1], deb_segment_4, max))
    ]

    # Find the best reconnection
    best_index = argmin([rc[2] for rc in reconnections])
    best_route, best_cost = reconnections[best_index]

    # Check if the best cost is an improvement
    if best_index == 1
        improvement = false
    else
        improvement = true
    end

    #println(best_index)
    #println([rc[2] for rc in reconnections])
    #println(improvement)
    return best_route, improvement
end

# Calcul du coup d'un swap 2-opt
function swap_2opt_cost(route::Vector{Int}, i::Int, j::Int, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    if i == 1
        end_segment_1 = 1
    else
        end_segment_1 = route[i-1]
    end
    if j == length(route)
        deb_segment_3 = 1
    else
        deb_segment_3 = route[j+1]
    end
    old_cost = cost(t, t_hat, end_segment_1, route[i], max) + cost(t, t_hat, route[j], deb_segment_3, max)
    new_cost = cost(t, t_hat, end_segment_1, route[j], max) + cost(t, t_hat, route[i], deb_segment_3, max)
    return new_cost - old_cost
end

function merge_routes_2opt(route_i::Vector{Int},route_j::Vector{Int},k::Int,l::Int)
    segment1 = route_i[1:k-1]
    segment2 = route_i[k:end]
    segment3 = route_j[1:l-1]
    segment4 = route_j[l:end]
    route_1 = vcat(segment1,segment4)
    route_2 = vcat(segment3,segment2)
    return route_1,route_2
end

function cost_merge_routes_2opt(route_i::Vector{Int},route_j::Vector{Int},k::Int,l::Int,  t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    if k == 1
        end_segment_1 = 1
    else
        end_segment_1 = route_i[k-1]
    end
    if l == length(route_j)
        deb_segment_3 = 1
    else
        deb_segment_3 = route[l+1]
    end
    old_cost = cost(t, t_hat, end_segment_1, route_i[k], max) + cost(t, t_hat, route_j[l], deb_segment_3, max)
    new_cost = cost(t, t_hat, end_segment_1, route_j[l], max) + cost(t, t_hat, route_i[k], deb_segment_3, max)
    return new_cost - old_cost

end

# Helper function to calculate the cost of a route
function calculate_route_cost(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    # Dépot vers ville 1
    route_cost = cost(t, t_hat, route[1], 1, max)
    for k in 1:length(route)-1
        route_cost += cost(t, t_hat, route[k], route[k+1], max)
    end
    # Dernière ville vers dépot
    route_cost += cost(t, t_hat, route[end], 1, max)
    return route_cost
end

function total_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
    c = 0.0
    for route in routes
        c+=calculate_route_cost(route, t, t_hat, max)
    end
    return c
end

function routes_to_x(routes::Vector{Vector{Int}},n::Int)
    x = Dict((i, j) => 0 for i in 1:n for j in 1:n if i != j)
    for route in routes
        for k in 1:length(route)-1
            i = route[k]
            j = route[k+1]
            x[(i, j)] = 1
        end
        # Dépot vers ville 1
        x[(1, route[1])] = 1
        # Dernière ville vers dépot
        x[(route[end],1)] = 1
    end
    
    return x
end

function real_cost(routes::Vector{Vector{Int}}, n::Int, t_hat::Vector{Int}, t::Matrix{Int}, T::Int)
    # Define the model
    m0 = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m0, "CPX_PARAM_SCRIND", 0)  # Suppress CPLEX output

    # Define sets
    V = 1:n
    A = [(i, j) for i in V for j in V if i != j]

    # Convert routes to x[(i, j)]
    x = routes_to_x(routes,n)

    # Define variables
    @variable(m0, delta_1[(i, j) in A] >= 0)
    @variable(m0, delta_2[(i, j) in A] >= 0)

    # Add constraints
    @constraint(m0, [(i, j) in A], delta_1[(i, j)] <= 1)
    @constraint(m0, [(i, j) in A], delta_2[(i, j)] <= 2)

    @constraint(m0, sum(delta_1[(i, j)] for (i, j) in A) <= T)
    @constraint(m0, sum(delta_2[(i, j)] for (i, j) in A) <= T^2)

    @objective(m0, Max, sum((delta_1[(i, j)] * (t_hat[i] + t_hat[j]) + delta_2[(i, j)] * t_hat[i] * t_hat[j]) * x[(i, j)] for (i, j) in A))

    # Solve the model
    optimize!(m0)

    # Return the objective value plus the fixed cost term
    return objective_value(m0) + sum(t[i, j] * x[(i, j)] for (i, j) in A)
end

# 2 opt sur les routes de CW
function lin_kernighan_VRP(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, max::Bool; two_opt::Bool=false)
    # Clark and Wright pour sol initiale
    routes = robust_clark_wright(n, t, t_hat, d, C, max)
    #println("Coût du meilleur routing après Clark Wright : ", total_cost(routes, t, t_hat, max))
    # Lin-Kernighen sur chaque tour
    optimized_routes = Vector{Vector{Int}}()
    for route in routes
        optimized_route = lin_kernighan(route, t, t_hat, true; two_opt=two_opt)
        push!(optimized_routes, optimized_route)
        #print(optimized_route)
    end
    #println("Coût d'une meilleure route après Lin-Kernighan : ", total_cost(optimized_routes, t, t_hat, max))
    return optimized_routes
end

function hybrid_heuristic(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, max::Bool; two_opt::Bool=false)
    
    # Clark and Wright pour sol initiale
    CW_routes = robust_clark_wright(n, t, t_hat, d, C, max)

    # Lin-Kernighen sur chaque tour
    LK_routes = Vector{Vector{Int}}()
    for route in CW_routes
        optimized_route = lin_kernighan_one_route(route, t, t_hat, true; two_opt=two_opt)
        push!(optimized_routes, optimized_route)
    end

    # Merge entre routes


end