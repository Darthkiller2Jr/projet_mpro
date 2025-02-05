using JuMP, CPLEX

function cost(t::Matrix{Int},t_hat::Vector{Int},i::Int,j::Int;max::Bool=true)
    if max
        return t[i,j] + (t_hat[i] + t_hat[j]) + 2 * t_hat[j] * t_hat[i]
    else
        return t[i,j]
    end
end

function calc_savings(n::Int, t::Matrix{Int},t_hat::Vector{Int};max::Bool=true)
    savings = Vector{Tuple{Int, Int, Int}}()
    for i in 1:n-1
        for j in i+1:n
            # Utilisation des temps max
            savings_ij = cost(t,t_hat,i,j,max=max)
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
    d_i = 0
    d_j = 0
    if length(route_i) != 0
        d_i = sum(d[i] for i in route_i)
    end
    if length(route_j) != 0
        d_j = sum(d[j] for j in route_j)
    end
    return d_i + d_j <= C
end

# Merge de 2 routes
function merge_routes(route_i::Vector{Int}, route_j::Vector{Int})
    return vcat(route_i, route_j)
end

function robust_clark_wright(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int; max::Bool=true)
    # Initialisation : une route/client
    routes = [[i] for i in 2:n]

    # Calculs des savings
    savings = calc_savings(n, t, t_hat, max=max)

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
function lin_kernighan_one_route(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, euclidien::Bool; max::Bool=true, two_opt::Bool=false)
    n = length(route)
    improvement = true
    #println(route,compute_route_cost(route,t,t_hat,max))
    old_cost = compute_route_cost(route,t,t_hat,max=max)
    while improvement
        improvement = false
        for i in 1:n-2
            for j in i+1:n-1
                #println("Route : $route, test : $i, $j")
                if two_opt
                    if euclidien
                        swap_cost = swap_2opt_cost(route,i,j,t,t_hat,max)
                        if swap_cost<0
                            #println("swap $i, $j ($swap_cost)")
                            route = two_opt_swap(route, i, j)
                            improvement = true
                            #println(route,compute_route_cost(route,t,t_hat,max))
                        end
                    else
                        new_route = two_opt_swap(route, i, j)
                        new_cost = compute_route_cost(new_route,t,t_hat,max=max)
                        if new_cost<old_cost
                            route = new_route
                            old_cost = new_cost
                        end
                    end   
                else
                    for k in j+1:n
                        if euclidien
                            new_route, improvement = three_opt_swap_best(route, i, j, k, t, t_hat, max)
                        else
                            new_route, improvement = three_opt_swap_best_non_eucl(route, i, j, k, t, t_hat, max)
                        end
                        if improvement
                            route = new_route
                            #println(route, compute_route_cost(route,t,t_hat,max))
                        end
                    end
                end
            end
        end
    end
    return route
end

function lin_kernighan_all_route(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C, euclidien::Bool; max::Bool=true, two_opt::Bool=false)
    improvement = true
    while improvement
        improvement = false
        num_routes = length(routes)

        for i in 1:num_routes
            for j in i+1:num_routes
                route_i = routes[i]
                route_j = routes[j]

                old_cost_ri = compute_route_cost(route_i, t, t_hat, max=max)
                old_cost_rj = compute_route_cost(route_j, t, t_hat, max=max)
                best_tot = old_cost_ri + old_cost_rj
                best_route_i = copy(route_i)
                best_route_j = copy(route_j)

                for k in 1:length(route_i)
                    for l in 1:length(route_j)
                        res_merge = merge_routes_2opt(route_i, route_j, k, l, d, C)
                        #println("$i, $j, $k, $l")
                        if res_merge !== nothing
                            #println("merge possible")
                            optimized_route_1 = lin_kernighan_one_route(res_merge[1], t, t_hat, euclidien; max=max, two_opt=two_opt)
                            #println("Descente route 1 terminée")
                            optimized_route_2 = lin_kernighan_one_route(res_merge[2], t, t_hat, euclidien; max=max,  two_opt=two_opt)
                            #println("Descente route 2 terminée")
                            new_cost_r1 = compute_route_cost(optimized_route_1, t, t_hat, max=max)
                            new_cost_r2 = compute_route_cost(optimized_route_2, t, t_hat, max=max)
                            if new_cost_r1 + new_cost_r2 < best_tot
                                improvement = true
                                #println("improvement")
                                best_tot = new_cost_r1 + new_cost_r2
                                #println(i," ",j," ",best_tot)
                                best_route_i = optimized_route_1
                                best_route_j = optimized_route_2
                            end
                        end
                    end
                end
                routes[i] = best_route_i
                routes[j] = best_route_j
            end
        end
    end

    return routes
end

function lin_kernighan_all_route_real_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C, euclidien::Bool=true; max::Bool=true, two_opt::Bool=false)
    improvement = true
    while improvement
        improvement = false
        num_routes = length(routes)
        old_cost = real_cost_smart(routes,n,t_hat,t,T)
        for i in 1:num_routes
            for j in i+1:num_routes
                route_i = routes[i]
                route_j = routes[j]

                best_route_i = copy(route_i)
                best_route_j = copy(route_j)

                for k in 1:length(route_i)
                    for l in 1:length(route_j)
                        res_merge = merge_routes_2opt(route_i, route_j, k, l, d, C)
                        #println("$i, $j, $k, $l")
                        if res_merge !== nothing
                            #println("merge possible")
                            optimized_route_1 = lin_kernighan_one_route(res_merge[1], t, t_hat, euclidien; max=max, two_opt=two_opt)
                            #println("Descente route 1 terminée")
                            optimized_route_2 = lin_kernighan_one_route(res_merge[2], t, t_hat, euclidien; max=max, two_opt=two_opt)
                            #println("Descente route 2 terminée")
                            new_routes = copy(routes)
                            new_routes[i] = optimized_route_1
                            new_routes[j] = optimized_route_2
                            merged_real_cost = real_cost_smart(new_routes,n,t_hat,t,T)
                            if merged_real_cost < old_cost
                                improvement = true
                                #println("improvement")
                                old_cost = merged_real_cost
                                #println(i," ",j," ",best_tot)
                                best_route_i = optimized_route_1
                                best_route_j = optimized_route_2
                            end
                        end
                    end
                end
                routes[i] = best_route_i
                routes[j] = best_route_j
            end
        end
    end

    return routes
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
    original_cost = cost(t, t_hat, end_segment_1, route[i], max=max) +
                cost(t, t_hat, route[j], route[j+1], max=max) +
                cost(t, t_hat, route[k], deb_segment_4, max=max)

    # All possible reconnections and their costs
    reconnections = [
        (route, original_cost),  # Original route
        (vcat(segment1, reverse(segment2), segment3, segment4), 
        cost(t, t_hat, end_segment_1, segment2[end], max=max) + cost(t, t_hat, segment2[1], segment3[1], max=max) + cost(t, t_hat, segment3[end], deb_segment_4, max=max)),
        (vcat(segment1, segment2, reverse(segment3), segment4), 
        cost(t, t_hat, end_segment_1, segment2[1], max=max) + cost(t, t_hat, segment2[end], segment3[end], max=max) + cost(t, t_hat, segment3[1], deb_segment_4, max=max)),
        (vcat(segment1, reverse(segment2), reverse(segment3), segment4), 
        cost(t, t_hat, end_segment_1, segment2[end], max=max) + cost(t, t_hat, segment2[1], segment3[end], max=max) + cost(t, t_hat, segment3[1], deb_segment_4, max=max)),
        (vcat(segment1, segment3, segment2, segment4), 
        cost(t, t_hat, end_segment_1, segment3[1], max=max) + cost(t, t_hat, segment3[end], segment2[1], max=max) + cost(t, t_hat, segment2[end], deb_segment_4, max=max)),
        (vcat(segment1, segment3, reverse(segment2), segment4), 
        cost(t, t_hat, end_segment_1, segment3[1], max=max) + cost(t, t_hat, segment3[end], segment2[end], max=max) + cost(t, t_hat, segment2[1], deb_segment_4, max=max)),
        (vcat(segment1, reverse(segment3), segment2, segment4), 
        cost(t, t_hat, end_segment_1, segment3[end], max=max) + cost(t, t_hat, segment3[1], segment2[1], max=max) + cost(t, t_hat, segment2[end], deb_segment_4, max=max)),
        (vcat(segment1, reverse(segment3), reverse(segment2), segment4), 
        cost(t, t_hat, end_segment_1, segment3[end], max=max) + cost(t, t_hat, segment3[1], segment2[end], max=max) + cost(t, t_hat, segment2[1], deb_segment_4, max=max))
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

function three_opt_swap_best_non_eucl(route::Vector{Int}, i::Int, j::Int, k::Int, t::Matrix{Int}, t_hat::Vector{Int}, max::Bool)
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
    original_cost = compute_route_cost(route,t,t_hat,max=max)

    # All possible reconnections and their costs
    reconnections = [
        route,  # Original route
        vcat(segment1, reverse(segment2), segment3, segment4), 
        vcat(segment1, segment2, reverse(segment3), segment4), 
        vcat(segment1, reverse(segment2), reverse(segment3), segment4), 
        vcat(segment1, segment3, segment2, segment4), 
        vcat(segment1, segment3, reverse(segment2), segment4), 
        vcat(segment1, reverse(segment3), segment2, segment4), 
        vcat(segment1, reverse(segment3), reverse(segment2), segment4), 
    ]


    new_costs = zeros(7)  # Preallocate the array
    for i in 1:7
        new_costs[i] = compute_route_cost(reconnections[i+1], t, t_hat, max=max)
    end

    # Find the best reconnection
    best_index = argmin(new_costs)

    # Check if the best cost is an improvement
    if new_costs[best_index]<original_cost
        improvement = true
        best_route = reconnections[best_index+1]
    else
        improvement = false
        best_route = route
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
    old_cost = cost(t, t_hat, end_segment_1, route[i], max=max) + cost(t, t_hat, route[j], deb_segment_3, max=max)
    new_cost = cost(t, t_hat, end_segment_1, route[j], max=max) + cost(t, t_hat, route[i], deb_segment_3, max=max)
    return new_cost - old_cost
end

function merge_routes_2opt(route_i::Vector{Int},route_j::Vector{Int},k::Int,l::Int,d::Vector{Int},C::Int)
    segment1 = route_i[1:k-1]
    segment2 = route_i[k:end]
    segment3 = route_j[1:l-1]
    segment4 = route_j[l:end]
    if is_feasible(segment1,segment4,d,C) && is_feasible(segment3,segment2,d,C)
        route_1 = vcat(segment1,segment4)
        route_2 = vcat(segment3,segment2)
        return route_1,route_2
    else
        return nothing
    end
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
    old_cost = cost(t, t_hat, end_segment_1, route_i[k], max=max) + cost(t, t_hat, route_j[l], deb_segment_3, max=max)
    new_cost = cost(t, t_hat, end_segment_1, route_j[l], max=max) + cost(t, t_hat, route_i[k], deb_segment_3, max=max)
    return new_cost - old_cost

end

# Helper function to calculate the cost of a route
function compute_route_cost(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}; max::Bool=true)
    # Dépot vers ville 1
    route_cost = cost(t, t_hat, route[1], 1, max=max)
    for k in 1:length(route)-1
        route_cost += cost(t, t_hat, route[k], route[k+1], max=max)
    end
    # Dernière ville vers dépot
    route_cost += cost(t, t_hat, route[end], 1, max=max)
    return route_cost
end

function total_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}; max::Bool=true)
    c = 0.0
    for route in routes
        c+=compute_route_cost(route, t, t_hat, max=max)
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

function x_to_routes(x::Dict{Tuple{Int, Int}, Float64}, n::Int)
    V = 1:n
    routes = Vector{Vector{Int}}()
    next_node = Dict{Int, Int}()
    
    for ((i, j), val) in x
        if val > 0.5 && i != 1 && j != 1
            next_node[i] = j
        end
    end
    
    # Collect all edges starting from the depot (node 1)
    depot_edges = Int[]
    for j in 1:n
        if j != 1 && get(x, (1, j), 0.0) > 0.5
            push!(depot_edges, j)
        end
    end
    
    # Build each route starting from the depot
    for j in depot_edges
        route = [j]
        current = j
        
        while true
            if haskey(next_node, current)
                nxt = next_node[current]
                push!(route, nxt)
                current = nxt
            else
                if get(x, (current, 1), 0.0) > 0.5
                    break
                else
                    error("Invalid x: node $current does not return to depot")
                end
            end
        end
        push!(routes, route)
    end
    return routes
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

function real_cost_smart(routes::Vector{Vector{Int}}, n::Int, t_hat::Vector{Int}, t::Matrix{Int}, T::Int)

    V = 1:n
    A = [(i, j) for i in V for j in V if i != j]
    x = routes_to_x(routes, n)

    # Create coefficient dictionaries for delta_1 and delta_2 over all arcs
    coeffs_delta_1 = Dict{Tuple{Int,Int}, Int}()
    coeffs_delta_2 = Dict{Tuple{Int,Int}, Int}()
    for (i, j) in A
        coeffs_delta_1[(i, j)] = (t_hat[i] + t_hat[j]) * x[(i, j)]
        coeffs_delta_2[(i, j)] = t_hat[i] * t_hat[j] * x[(i, j)]
    end

    sorted_arcs_delta_1 = sort(collect(keys(coeffs_delta_1)), by = k -> coeffs_delta_1[k], rev = true)
    sorted_arcs_delta_2 = sort(collect(keys(coeffs_delta_2)), by = k -> coeffs_delta_2[k], rev = true)

    # Initialize delta dictionaries with zero values for every arc
    delta_1 = Dict{Tuple{Int,Int}, Int}((i, j) => 0 for (i, j) in A)
    delta_2 = Dict{Tuple{Int,Int}, Int}((i, j) => 0 for (i, j) in A)

    # Allocate delta_1 values until a total of T is reached.
    sum_delta_1 = 0
    idx = 1  # Using 1-indexing for arrays in Julia.
    while sum_delta_1 < T && idx <= length(sorted_arcs_delta_1)
        arc = sorted_arcs_delta_1[idx]
        # Here we assign 1 unit of delta_1 for each arc in order.
        delta_1[arc] = 1
        sum_delta_1 += 1
        idx += 1
    end
    # Pas de reste possible

    # Allocate delta_2 values until a total of T^2 is reached.
    sum_delta_2 = 0
    idx = 1
    while sum_delta_2 < T^2 && idx <= length(sorted_arcs_delta_2)
        arc = sorted_arcs_delta_2[idx]
        delta_2[arc] = 2
        sum_delta_2 += 2
        idx += 1
    end
    # If there is a remainder, allocate it to the next arc (if one exists)
    if sum_delta_2 < T^2 && idx <= length(sorted_arcs_delta_2)
        arc = sorted_arcs_delta_2[idx]
        delta_2[arc] = T^2 - sum_delta_2
        sum_delta_2 = T^2
    end

    # Compute and return the overall cost.
    cost = 0
    for (i, j) in A
        cost += (t[i, j] + delta_1[(i, j)] * (t_hat[i] + t_hat[j]) +
                 delta_2[(i, j)] * t_hat[i] * t_hat[j]) * x[(i, j)]
    end

    return cost
end


# 2 opt sur les routes de CW
function lin_kernighan_VRP(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, max::Bool, euclidien; two_opt::Bool=false)
    # Clark and Wright pour sol initiale
    routes = robust_clark_wright(n, t, t_hat, d, C, max=max)
    #println("Coût du meilleur routing après Clark Wright : ", total_cost(routes, t, t_hat, max))
    # Lin-Kernighen sur chaque tour
    optimized_routes = Vector{Vector{Int}}()
    for route in routes
        optimized_route = lin_kernighan_one_route(route, t, t_hat,euclidien; max=max, two_opt=two_opt)
        push!(optimized_routes, optimized_route)
        #print(optimized_route)
    end
    #println("Coût d'une meilleure route après Lin-Kernighan : ", total_cost(optimized_routes, t, t_hat, max))
    return optimized_routes
end

function hybrid_heuristic(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, euclidien::Bool=true; max::Bool=true, real_cost=false, two_opt::Bool=false)
    
    # Clark and Wright pour sol initiale
    CW_routes = robust_clark_wright(n, t, t_hat, d, C, max=max)
    #println(CW_routes,real_cost(CW_routes,n,th,t,T))

    # Lin-Kernighen sur chaque tour
    LK_routes = Vector{Vector{Int}}()
    for route in CW_routes
        LK_route = lin_kernighan_one_route(route, t, t_hat, euclidien; max=max, two_opt=two_opt)
        push!(LK_routes, LK_route)
    end
    #println(LK_routes,real_cost(LK_routes,n,th,t,T))
    # Merge entre routes
    if real_cost
        LK_all_routes = lin_kernighan_all_route_real_cost(LK_routes, t, t_hat, d, C, euclidien;max=max, two_opt=two_opt)
    else
        LK_all_routes = lin_kernighan_all_route(LK_routes, t, t_hat, d, C, euclidien;max=max, two_opt=two_opt)
    end
    #println(LK_all_routes,real_cost(LK_all_routes,n,th,t,T))

    return LK_all_routes
end