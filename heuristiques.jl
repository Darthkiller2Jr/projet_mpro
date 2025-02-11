using JuMP, CPLEX

function cost(t::Matrix{Int},t_hat::Vector{Int},i::Int,j::Int;max_cost::Bool=true)
    if max_cost
        return t[i,j] + (t_hat[i] + t_hat[j]) + 2 * t_hat[j] * t_hat[i]
    else
        return t[i,j]
    end
end

function calc_savings(n::Int, t::Matrix{Int},t_hat::Vector{Int};euclidien::Bool=false, max_cost::Bool=true)
    savings = Vector{Tuple{Int, Int, Int}}()
    
    function interval(i::Int)
        if euclidien
            return i+1:n
        else
            return 1:n
        end
    end

    for i in 1:n
        for j in interval(i)
            if i != j
                # Utilisation des temps max_cost
                savings_ij = cost(t,t_hat,1,j,max_cost=max_cost) + cost(t,t_hat,i,1,max_cost=max_cost) - cost(t,t_hat,i,j,max_cost=max_cost)
                push!(savings, (savings_ij, i, j))
            end
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

function is_feasible(route::Vector{Int}, d::Vector{Int}, C::Int)
    d_i = 0
    if length(route) != 0
        d_i = sum(d[i] for i in route)
    end
    return d_i <= C
end

# Merge de 2 routes
function merge_routes(route_i::Vector{Int}, route_j::Vector{Int})
    return vcat(route_i, route_j)
end

function robust_clark_wright(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int; max_cost::Bool=true)
    # Initialisation : une route/client
    routes = [[i] for i in 2:n]

    # Calculs des savings
    savings = calc_savings(n, t, t_hat, max_cost=max_cost)

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

function explo_2_3opt(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, euclidien::Bool; max_cost::Bool=true, two_opt::Bool=false)
    n = length(route)
    improvement = true
    #println(route,compute_route_cost(route,t,t_hat,max_cost))
    old_cost = compute_route_cost(route,t,t_hat, d, C, max_cost=max_cost)
    while improvement
        improvement = false
        for i in 1:n-2
            for j in i+1:n-1
                #println("Route : $route, test : $i, $j")
                if two_opt
                    if euclidien
                        swap_cost = swap_2opt_cost(route,i,j,t,t_hat,max_cost)
                        if swap_cost<0
                            #println("swap $i, $j ($swap_cost)")
                            route = two_opt_swap(route, i, j)
                            improvement = true
                            #println(route,compute_route_cost(route,t,t_hat,max_cost))
                        end
                    else
                        new_route = two_opt_swap(route, i, j)
                        new_cost = compute_route_cost(new_route,t,t_hat, d, C, max_cost=max_cost)
                        if new_cost<old_cost
                            route = new_route
                            old_cost = new_cost
                        end
                    end   
                else
                    for k in j+1:n
                        if euclidien
                            new_route, improvement = three_opt_swap_best(route, i, j, k, t, t_hat, max_cost)
                        else
                            new_route, improvement = three_opt_swap_best_non_eucl(route, i, j, k, t, t_hat, max_cost)
                        end
                        if improvement
                            route = new_route
                            #println(route, compute_route_cost(route,t,t_hat,max_cost))
                        end
                    end
                end
            end
        end
    end
    return route
end

function lin_kernighan_one_route(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, euclidien::Bool; max_cost::Bool=true)
    best_route = copy(route)
    best_cost = compute_route_cost(best_route, t, t_hat, d, C, max_cost=max_cost)

    # We will try improvements until no further gain is possible.
    improvement = true
    while improvement
        improvement = false
        # Try each edge in the current best tour as a candidate for removal.
        n = length(best_route)
        for i in 1:n
            # t1 is the node where we start a sequence of moves.
            t1 = best_route[i]
            # t2 is the neighbor following t1 (wrap-around at the end)
            t2 = (i < n) ? best_route[i+1] : best_route[1]
            # The gain from removing edge (t1,t2) is our starting point.
            G0 = t[t1, t2]
            # used will keep track of nodes already involved in the move sequence
            used = falses(length(t))
            used[t1] = true
            used[t2] = true

            # Start a recursive search for a sequence of moves
            candidate_route, candidate_cost, candidate_gain = lk_move!(best_route, i, t1, t2, G0, used, t, t_hat; max_cost=max_cost)
            if candidate_gain < 0 && candidate_cost < best_cost
                best_route = candidate_route
                best_cost = candidate_cost
                improvement = true
                break  # Accept the move and restart search over all edges
            end
        end
    end

    return best_route
end

# Recursive function for variable-depth move search.
# Parameters:
# - current_route: the current tour (vector of node indices)
# - i: index in the route where t1 resides (for reference)
# - t1: starting node of the move sequence
# - t_prev: the node last removed from the tour (the "old" neighbor)
# - G: accumulated gain so far (we try to make G negative)
# - used: Boolean vector marking nodes already used in the sequence
#
# Returns: (new_route, new_cost, gain) 
function lk_move!(current_route::Vector{Int}, i::Int, t1::Int, t_prev::Int, G::Int, used::AbstractVector{Bool},
                    t::Matrix{Int}, t_hat::Vector{Int}; max_cost::Bool=true, depth::Int=1, max_costDepth::Int=5)

    best_route = current_route
    best_cost = compute_route_cost(current_route, t, t_hat, d, C, max_cost=max_cost)
    best_gain = G

    n = length(current_route)
    # Loop over candidate nodes t_new that are not used yet
    for j in 1:n
        t_new = current_route[j]
        # Do not allow t_new to be t1 or already used in this sequence.
        if used[t_new] || t_new == t1
            continue
        end

        # Consider adding edge (t_prev, t_new) as a new connection.
        # Compute the gain from adding (t_prev, t_new)
        delta = t[t_prev, t_new]  # gain: cost of edge we add
        newG = G - delta
        # If newG is not promising, skip.
        if newG >= 0
            continue
        end

        # Now, try to complete a move sequence by reconnecting t_new to t1.
        # The candidate move is: remove (t1,?) and add (t_prev, t_new), then add (t_new,t1).
        reconnectionCost = t[t_new, t1]
        totalGain = newG + reconnectionCost

        # If this move results in an overall improvement, try to update.
        if totalGain < 0
            # Construct a new route by performing a 2-opt move.
            # We assume that the removal of edge (t1,t_prev) and addition of (t_prev,t_new) plus closing edge (t_new,t1)
            # is equivalent to reversing an appropriate segment. (The precise details depend on the structure of the tour.)
            # Here, we pick indices so that performing a 2-opt move between positions corresponding to t_prev and t_new yields the new tour.
            pos_t_prev = findfirst(==(t_prev), current_route)
            pos_t_new = findfirst(==(t_new), current_route)
            if pos_t_prev !== nothing && pos_t_new !== nothing && pos_t_prev < pos_t_new
                new_route = two_opt_swap(current_route, pos_t_prev, pos_t_new)
                new_cost = compute_route_cost(new_route, t, t_hat, d, C, max_cost=max_cost)
                if new_cost < best_cost
                    best_route = new_route
                    best_cost = new_cost
                    best_gain = totalGain
                end
            end
        end

        # If we haven't reached our recursion depth limit, try to extend the sequence.
        if depth < max_costDepth
            # Mark t_new as used and search deeper.
            used_copy = copy(used)
            used_copy[t_new] = true
            rec_route, rec_cost, rec_gain = lk_move!(current_route, i, t1, t_new, newG, used_copy, t, t_hat; max_cost=max_cost, depth=depth+1, max_costDepth=max_costDepth)
            if rec_gain < best_gain && rec_cost < best_cost
                best_route = rec_route
                best_cost = rec_cost
                best_gain = rec_gain
            end
        end
    end

    return best_route, best_cost, best_gain
end


function swap_between_routes(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C, euclidien::Bool; max_cost::Bool=true, LK::Bool=false, two_opt::Bool=false)
    improvement = true
    while improvement
        improvement = false
        num_routes = length(routes)

        for i in 1:num_routes
            for j in i+1:num_routes
                route_i = routes[i]
                route_j = routes[j]

                old_cost_ri = compute_route_cost(route_i, t, t_hat, d, C, max_cost=max_cost)
                old_cost_rj = compute_route_cost(route_j, t, t_hat, d, C, max_cost=max_cost)
                best_tot = old_cost_ri + old_cost_rj
                best_route_i = copy(route_i)
                best_route_j = copy(route_j)

                for k in 1:length(route_i)
                    for l in 1:length(route_j)
                        res_merge = merge_routes_2opt(route_i, route_j, k, l, d, C)
                        #println("$i, $j, $k, $l")
                        
                        if LK
                            optimized_route_1 = lin_kernighan_one_route(res_merge[1], t, t_hat, d, C, euclidien; max_cost=max_cost)
                            optimized_route_2 = lin_kernighan_one_route(res_merge[2], t, t_hat, d, C, euclidien; max_cost=max_cost)
                        else
                            optimized_route_1 = explo_2_3opt(res_merge[1], t, t_hat, d, C, euclidien ; two_opt=two_opt, max_cost=max_cost)
                            optimized_route_2 = explo_2_3opt(res_merge[2], t, t_hat, d, C, euclidien; two_opt=two_opt, max_cost=max_cost)
                        end
                        new_cost_r1 = compute_route_cost(optimized_route_1, t, t_hat, d, C, max_cost=max_cost)
                        new_cost_r2 = compute_route_cost(optimized_route_2, t, t_hat, d, C, max_cost=max_cost)
                        if new_cost_r1 + new_cost_r2 < best_tot && is_feasible(optimized_route_1,d,C) && is_feasible(optimized_route_2,d,C)
                            improvement = true
                            #println("improvement")
                            best_tot = new_cost_r1 + new_cost_r2
                            #println(i," ",j," ",best_tot)
                            best_route_i = optimized_route_1
                            best_route_j = optimized_route_2
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

function three_opt_permutations(routes::Vector{Vector{Int}}, route_i::Int, route_j::Int, route_k::Int, i::Int, j::Int, k::Int)

    # Extract route segments
    si_1 = routes[route_i][1:i-1]
    si_2 = routes[route_i][i:end]
    sj_1 = routes[route_j][1:j-1]
    sj_2 = routes[route_j][j:end]
    sk_1 = routes[route_k][1:k-1]
    sk_2 = routes[route_k][k:end]

    # Generate all possible permutations
    permutations = Vector{Vector{Vector{Int}}}()

    # Permutation 1
    new_routes_p1 = copy(routes)
    new_routes_p1[route_i] = vcat(si_1, sj_2)
    new_routes_p1[route_j] = vcat(sj_1, si_2)
    new_routes_p1[route_k] = vcat(sk_1, sk_2)
    push!(permutations, new_routes_p1)

    # Permutation 2
    new_routes_p2 = copy(routes)
    new_routes_p2[route_i] = vcat(si_1, sj_2)
    new_routes_p2[route_j] = vcat(sj_1, sk_2)
    new_routes_p2[route_k] = vcat(sk_1, si_2)
    push!(permutations, new_routes_p2)

    # Permutation 3
    new_routes_p3 = copy(routes)
    new_routes_p3[route_i] = vcat(si_1, sk_2)
    new_routes_p3[route_j] = vcat(sj_1, si_2)
    new_routes_p3[route_k] = vcat(sk_1, sj_2)
    push!(permutations, new_routes_p3)

    # Permutation 4
    new_routes_p4 = copy(routes)
    new_routes_p4[route_i] = vcat(si_1, sk_2)
    new_routes_p4[route_j] = vcat(sj_1, sj_2)
    new_routes_p4[route_k] = vcat(sk_1, si_2)
    push!(permutations, new_routes_p4)

    # Permutation 5
    new_routes_p5 = copy(routes)
    new_routes_p5[route_i] = vcat(sj_1, si_2)
    new_routes_p5[route_j] = vcat(si_1, sj_2)
    new_routes_p5[route_k] = vcat(sk_1, sk_2)
    push!(permutations, new_routes_p5)

    # Permutation 6
    new_routes_p6 = copy(routes)
    new_routes_p6[route_i] = vcat(si_1, si_2)
    new_routes_p6[route_j] = vcat(sj_1, sk_2)
    new_routes_p6[route_k] = vcat(sk_1, sj_2)
    push!(permutations, new_routes_p6)

    # Permutation 7
    new_routes_p8 = copy(routes)
    new_routes_p8[route_i] = vcat(sk_1, sj_2)
    new_routes_p8[route_j] = vcat(si_1, sk_2)
    new_routes_p8[route_k] = vcat(sj_1, si_2)
    push!(permutations, new_routes_p8)

    return permutations
end

function swap_between_routes_3opt(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C, euclidien::Bool; max_cost::Bool=true, LK::Bool=false, two_opt::Bool=false)
    improvement = true
    while improvement
        improvement = false
        num_routes = length(routes)

        for i in 1:num_routes
            for j in i+1:num_routes
                for k in j+1:num_routes
                    route_i = routes[i]
                    route_j = routes[j]
                    route_k = routes[k]

                    old_cost_ri = compute_route_cost(route_i, t, t_hat, d, C, max_cost=max_cost)
                    old_cost_rj = compute_route_cost(route_j, t, t_hat, d, C, max_cost=max_cost)
                    old_cost_rk = compute_route_cost(route_k, t, t_hat, d, C, max_cost=max_cost)
                    best_tot = old_cost_ri + old_cost_rj + old_cost_rk

                    best_route_i = copy(route_i)
                    best_route_j = copy(route_j)
                    best_route_k = copy(route_k)

                    for l in 1:length(route_i)
                        for m in 1:length(route_j)
                            for n in 1:length(route_k)

                                permutations = three_opt_permutations(routes,i,j,k,l,m,n)
                                
                                for routes_p in permutations
                                    if LK
                                        optimized_route_1 = lin_kernighan_one_route(routes_p[i], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_2 = lin_kernighan_one_route(routes_p[j], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_3 = lin_kernighan_one_route(routes_p[k], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                    else
                                        optimized_route_1 = explo_2_3opt(routes_p[i], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_2 = explo_2_3opt(routes_p[j], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_3 = explo_2_3opt(routes_p[k], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                    end

                                    new_cost_r1 = compute_route_cost(optimized_route_1, t, t_hat, d, C, max_cost=max_cost)
                                    new_cost_r2 = compute_route_cost(optimized_route_2, t, t_hat, d, C, max_cost=max_cost)
                                    new_cost_r3 = compute_route_cost(optimized_route_3, t, t_hat, d, C, max_cost=max_cost)

                                    if new_cost_r1 + new_cost_r2 + new_cost_r3 < best_tot && is_feasible(optimized_route_1,d,C) && is_feasible(optimized_route_2,d,C) && is_feasible(optimized_route_3,d,C)
                                        improvement = true
                                        #println("improvement")
                                        best_tot = new_cost_r1 + new_cost_r2
                                        #println(i," ",j," ",best_tot)
                                        best_route_i = optimized_route_1
                                        best_route_j = optimized_route_2
                                        best_route_k = optimized_route_3
                                    end
                                end
                            end
                        end
                    end
                    routes[i] = best_route_i
                    routes[j] = best_route_j
                    routes[k] = best_route_k
                end
            end
        end
    end

    return routes
end

function swap_between_routes_3opt_real_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C, euclidien::Bool; max_cost::Bool=true, LK::Bool=false, two_opt::Bool=false)
    improvement = true
    while improvement
        improvement = false
        num_routes = length(routes)
        old_cost = real_cost_smart(routes,n,t_hat,t,T)
        for i in 1:num_routes
            for j in i+1:num_routes
                for k in j+1:num_routes
                    route_i = routes[i]
                    route_j = routes[j]
                    route_k = routes[k]

                    best_route_i = copy(route_i)
                    best_route_j = copy(route_j)
                    best_route_k = copy(route_k)

                    for l in 1:length(route_i)
                        for m in 1:length(route_j)
                            for n in 1:length(route_k)

                                permutations = three_opt_permutations(routes,i,j,k,l,m,n)
                                
                                for routes_p in permutations
                                    if LK
                                        optimized_route_1 = lin_kernighan_one_route(routes_p[i], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_2 = lin_kernighan_one_route(routes_p[j], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_3 = lin_kernighan_one_route(routes_p[k], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                    else
                                        optimized_route_1 = explo_2_3opt(routes_p[i], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_2 = explo_2_3opt(routes_p[j], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                        optimized_route_3 = explo_2_3opt(routes_p[k], t, t_hat, d, C, euclidien; max_cost=max_cost)
                                    end

                                    new_routes = copy(routes)
                                    new_routes[i] = optimized_route_1
                                    new_routes[j] = optimized_route_2
                                    new_routes[k] = optimized_route_3

                                    merged_real_cost = real_cost_smart(new_routes,n,t_hat,t,T)

                                    if merged_real_cost < old_cost && is_feasible(optimized_route_1,d,C) && is_feasible(optimized_route_2,d,C) && is_feasible(optimized_route_3,d,C)
                                        improvement = true
                                        old_cost = merged_real_cost
                                        best_route_i = optimized_route_1
                                        best_route_j = optimized_route_2
                                        best_route_k = optimized_route_3
                                    end
                                end
                            end
                        end
                    end
                end
                routes[i] = best_route_i
                routes[j] = best_route_j
                routes[k] = best_route_k
            end
        end
    end

    return routes
end

function swap_between_routes_real_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, T::Int,euclidien::Bool=true; max_cost::Bool=true, LK::Bool=false, two_opt::Bool=false)
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

                for k in 2:length(route_i)-1
                    for l in 2:length(route_j)-1
                        res_merge = merge_routes_2opt(route_i, route_j, k, l, d, C)
                        #println("$i, $j, $k, $l")
                        if res_merge !== nothing
                            if LK
                                optimized_route_1 = lin_kernighan_one_route(res_merge[1], t, t_hat, d, C,euclidien; max_cost=max_cost)
                                optimized_route_2 = lin_kernighan_one_route(res_merge[2], t, t_hat, d, C, euclidien; max_cost=max_cost)
                            else
                                route_1,route_2,route_1_bis,route_2bis = res_merge[1],res_merge[2], res_merge[3],res_merge[4]

                                optimized_route_1 = explo_2_3opt(route_1, t, t_hat, d, C, euclidien; two_opt=two_opt, max_cost=max_cost)
                                optimized_route_2 = explo_2_3opt(route_2, t, t_hat, d, C, euclidien; two_opt=two_opt, max_cost=max_cost)

                                optimized_route_1bis = explo_2_3opt(route_1_bis, t, t_hat, d, C, euclidien; two_opt=two_opt, max_cost=max_cost)
                                optimized_route_2bis = explo_2_3opt(route_2bis, t, t_hat, d, C, euclidien; two_opt=two_opt, max_cost=max_cost)
                            end

                            new_routes = copy(routes)
                            new_routes[i] = optimized_route_1
                            new_routes[j] = optimized_route_2
                            merged_real_cost = real_cost_smart(new_routes,n,t_hat,t,T)

                            new_routes_bis = copy(routes)
                            new_routes_bis[i] = optimized_route_1bis
                            new_routes_bis[j] = optimized_route_2bis
                            merged_real_cost_bis = real_cost_smart(new_routes_bis,n,t_hat,t,T)

                            if merged_real_cost_bis < merged_real_cost
                                new_routes = new_routes_bis
                                merged_real_cost = merged_real_cost_bis
                                #println(new_routes)
                                #println(merged_real_cost)
                            end

                            if merged_real_cost < old_cost && is_feasible(new_routes[i],d,C) && is_feasible(new_routes[j],d,C)
                                improvement = true
                                #println("improvement")
                                old_cost = merged_real_cost
                                #println(i," ",j," ",best_tot)
                                best_route_i = new_routes[i]
                                best_route_j = new_routes[j]
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

function three_opt_swap_best(route::Vector{Int}, i::Int, j::Int, k::Int, t::Matrix{Int}, t_hat::Vector{Int}, max_cost::Bool)
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
    original_cost = cost(t, t_hat, end_segment_1, route[i], max_cost=max_cost) +
                cost(t, t_hat, route[j], route[j+1], max_cost=max_cost) +
                cost(t, t_hat, route[k], deb_segment_4, max_cost=max_cost)

    # All possible reconnections and their costs
    reconnections = [
        (route, original_cost),  # Original route
        (vcat(segment1, reverse(segment2), segment3, segment4), 
        cost(t, t_hat, end_segment_1, segment2[end], max_cost=max_cost) + cost(t, t_hat, segment2[1], segment3[1], max_cost=max_cost) + cost(t, t_hat, segment3[end], deb_segment_4, max_cost=max_cost)),
        (vcat(segment1, segment2, reverse(segment3), segment4), 
        cost(t, t_hat, end_segment_1, segment2[1], max_cost=max_cost) + cost(t, t_hat, segment2[end], segment3[end], max_cost=max_cost) + cost(t, t_hat, segment3[1], deb_segment_4, max_cost=max_cost)),
        (vcat(segment1, reverse(segment2), reverse(segment3), segment4), 
        cost(t, t_hat, end_segment_1, segment2[end], max_cost=max_cost) + cost(t, t_hat, segment2[1], segment3[end], max_cost=max_cost) + cost(t, t_hat, segment3[1], deb_segment_4, max_cost=max_cost)),
        (vcat(segment1, segment3, segment2, segment4), 
        cost(t, t_hat, end_segment_1, segment3[1], max_cost=max_cost) + cost(t, t_hat, segment3[end], segment2[1], max_cost=max_cost) + cost(t, t_hat, segment2[end], deb_segment_4, max_cost=max_cost)),
        (vcat(segment1, segment3, reverse(segment2), segment4), 
        cost(t, t_hat, end_segment_1, segment3[1], max_cost=max_cost) + cost(t, t_hat, segment3[end], segment2[end], max_cost=max_cost) + cost(t, t_hat, segment2[1], deb_segment_4, max_cost=max_cost)),
        (vcat(segment1, reverse(segment3), segment2, segment4), 
        cost(t, t_hat, end_segment_1, segment3[end], max_cost=max_cost) + cost(t, t_hat, segment3[1], segment2[1], max_cost=max_cost) + cost(t, t_hat, segment2[end], deb_segment_4, max_cost=max_cost)),
        (vcat(segment1, reverse(segment3), reverse(segment2), segment4), 
        cost(t, t_hat, end_segment_1, segment3[end], max_cost=max_cost) + cost(t, t_hat, segment3[1], segment2[end], max_cost=max_cost) + cost(t, t_hat, segment2[1], deb_segment_4, max_cost=max_cost))
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

function three_opt_swap_best_non_eucl(route::Vector{Int}, i::Int, j::Int, k::Int, t::Matrix{Int}, t_hat::Vector{Int}, max_cost::Bool)
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
    original_cost = compute_route_cost(route,t,t_hat,d,C,max_cost=max_cost)

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
        new_costs[i] = compute_route_cost(reconnections[i+1], t, t_hat, d, C, max_cost=max_cost)
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
function swap_2opt_cost(route::Vector{Int}, i::Int, j::Int, t::Matrix{Int}, t_hat::Vector{Int}, max_cost::Bool)
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
    old_cost = cost(t, t_hat, end_segment_1, route[i], max_cost=max_cost) + cost(t, t_hat, route[j], deb_segment_3, max_cost=max_cost)
    new_cost = cost(t, t_hat, end_segment_1, route[j], max_cost=max_cost) + cost(t, t_hat, route[i], deb_segment_3, max_cost=max_cost)
    return new_cost - old_cost
end

function merge_routes_2opt(route_i::Vector{Int},route_j::Vector{Int},k::Int,l::Int,d::Vector{Int},C::Int)
    segment1 = route_i[1:k-1]
    segment2 = route_i[k:end]
    segment3 = route_j[1:l-1]
    segment4 = route_j[l:end]

    route_1 = vcat(segment1,segment4)
    route_2 = vcat(segment3,segment2)

    route_1bis = vcat(segment1,reverse(segment3))
    route_2bis = vcat(reverse(segment4),segment2)
    #println(k,l,route_i,route_j,route_1,route_2, route_1bis, route_2bis)
    return route_1,route_2, route_1bis, route_2bis
end

function cost_merge_routes_2opt(route_i::Vector{Int},route_j::Vector{Int},k::Int,l::Int,  t::Matrix{Int}, t_hat::Vector{Int}, max_cost::Bool)
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
    old_cost = cost(t, t_hat, end_segment_1, route_i[k], max_cost=max_cost) + cost(t, t_hat, route_j[l], deb_segment_3, max_cost=max_cost)
    new_cost = cost(t, t_hat, end_segment_1, route_j[l], max_cost=max_cost) + cost(t, t_hat, route_i[k], deb_segment_3, max_cost=max_cost)
    return new_cost - old_cost

end

# Helper function to calculate the cost of a route
function compute_route_cost(route::Vector{Int}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int; max_cost::Bool=true)
    # Dépot vers ville 1
    route_cost = cost(t, t_hat, route[1], 1, max_cost=max_cost)
    for k in 1:length(route)-1
        route_cost += cost(t, t_hat, route[k], route[k+1], max_cost=max_cost)
    end
    # Dernière ville vers dépot
    route_cost += cost(t, t_hat, route[end], 1, max_cost=max_cost)

    # Non réalisabilitée 
    non_realisabilite = 0
    if length(route) != 0
        non_realisabilite = max(0,sum(d[j] for j in route)-C)
    end

    return route_cost + non_realisabilite*200
end

function total_cost(routes::Vector{Vector{Int}}, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int ; max_cost::Bool=true)
    c = 0.0
    for route in routes
        c+=compute_route_cost(route, t, t_hat, d, C, max_cost=max_cost)
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

function collect_active_arcs(routes::Vector{Vector{Int}})
    # Collect active arcs directly from routes
    active_arcs = Tuple{Int, Int}[]
    for route in routes
        len = length(route)
        for i in 1:(len-1)
            push!(active_arcs, (route[i], route[i+1]))
        end
        push!(active_arcs, (1, route[1]))
        push!(active_arcs, (route[end], 1))
    end
    num_arcs = length(active_arcs)
    
    return active_arcs, num_arcs
end

function real_cost_smart(routes::Vector{Vector{Int}}, n::Int, t_hat::Vector{Int}, t::Matrix{Int}, T::Int)
    
    active_arcs, num_arcs = collect_active_arcs(routes)

    # Precompute necessary values for active arcs
    t_vals = Vector{Int}(undef, num_arcs)
    t_hat_sum = Vector{Int}(undef, num_arcs)
    t_hat_prod = Vector{Int}(undef, num_arcs)
    @inbounds for idx in 1:num_arcs
        i, j = active_arcs[idx]
        t_vals[idx] = t[i, j]
        t_hat_sum[idx] = t_hat[i] + t_hat[j]
        t_hat_prod[idx] = t_hat[i] * t_hat[j]
    end

    # Sort indices based on precomputed coefficients
    sorted_order_delta1 = sortperm(t_hat_sum, rev=true)
    sorted_order_delta2 = sortperm(t_hat_prod, rev=true)

    # Allocate delta_1 efficiently
    delta_1 = zeros(Int, num_arcs)
    T_alloc = min(T, num_arcs)
    if T_alloc > 0
        @views delta_1[sorted_order_delta1[1:T_alloc]] .= 1
    end

    # Allocate delta_2 with possible remainder
    delta_2 = zeros(Int, num_arcs)
    max_cost_assign = min(num_arcs, T^2 ÷ 2)
    sum_d2 = 2 * max_cost_assign
    if max_cost_assign > 0
        @views delta_2[sorted_order_delta2[1:max_cost_assign]] .= 2
    end

    # Handle remainder for delta_2
    remainder = T - sum_d2
    if remainder > 0 && max_cost_assign < num_arcs
        delta_2[sorted_order_delta2[max_cost_assign + 1]] += remainder
    end

    # Compute the total cost using SIMD and inbounds for max_costimum speed
    cost = 0
    @inbounds @simd for idx in 1:num_arcs
        cost += t_vals[idx] + delta_1[idx] * t_hat_sum[idx] + delta_2[idx] * t_hat_prod[idx]
    end

    return cost
end

# 2 opt sur les routes de CW
function sous_tours_heuristic(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, max_cost::Bool, euclidien; two_opt::Bool=false, LK::Bool=false)
    # Clark and Wright pour sol initiale
    routes = robust_clark_wright(n, t, t_hat, d, C, max_cost=max_cost)
    #println("Coût du meilleur routing après Clark Wright : ", total_cost(routes, t, t_hat, max_cost))
    # Lin-Kernighen sur chaque tour
    optimized_routes = Vector{Vector{Int}}()
    for route in routes
        if LK
            optimized_route = lin_kernighan_one_route(route, t, t_hat, d, C, euclidien; max_cost=max_cost)
        else
            optimized_route = explo_2_3opt(route, t, t_hat, d, C, euclidien; two_opt=two_opt,max_cost=max_cost)
        end
        push!(optimized_routes, optimized_route)
        #print(optimized_route)
    end
    #println("Coût d'une meilleure route après Lin-Kernighan : ", total_cost(optimized_routes, t, t_hat, max_cost))
    return optimized_routes
end

function hybrid_heuristic(n::Int, t::Matrix{Int}, t_hat::Vector{Int}, d::Vector{Int}, C::Int, T::Int, euclidien::Bool=true; max_cost::Bool=true, real_cost=false, LK::Bool=false, two_opt::Bool=false, three_opt_swap::Bool=false)
    
    # Clark and Wright pour sol initiale
    CW_routes = robust_clark_wright(n, t, t_hat, d, C, max_cost=max_cost)
    #println(CW_routes,real_cost(CW_routes,n,th,t,T))

    # explo sur chaque tour
    LK_routes = Vector{Vector{Int}}()
    for route in CW_routes
        if LK
            LK_route = lin_kernighan_one_route(route, t, t_hat, d, C, euclidien; max_cost=max_cost)
        else
            LK_route = explo_2_3opt(route, t, t_hat, d, C, euclidien; two_opt=two_opt, max_cost=max_cost)
        end
        push!(LK_routes, LK_route)
    end
    #println(LK_routes,real_cost(LK_routes,n,th,t,T))
    # Merge entre routes
    if real_cost
        if three_opt_swap
            LK_all_routes = swap_between_routes_3opt(LK_routes, t, t_hat, d, C, euclidien; max_cost=max_cost, LK=LK, two_opt=two_opt)
        else
            LK_all_routes = swap_between_routes(LK_routes, t, t_hat, d, C, euclidien; max_cost=max_cost, LK=LK, two_opt=two_opt)
        end
    else
        if three_opt_swap
            LK_all_routes = swap_between_routes_3opt_real_cost(LK_routes, t, t_hat, d, C, euclidien; max_cost=max_cost, LK=LK, two_opt=two_opt)
        else
            LK_all_routes = swap_between_routes_real_cost(LK_routes, t, t_hat, d, C, T, euclidien; max_cost=max_cost, LK=LK, two_opt=two_opt)
        end
    end
    #println(LK_all_routes,real_cost(LK_all_routes,n,th,t,T))

    return LK_all_routes
end