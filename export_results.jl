include("heuristiques.jl")
include("dual.jl")
include("plans_coupants.jl")
include("optim.jl")

function comparaison_heuristiques_n(n::Int, euclidien::Bool)
    
    filename = "data/n_$n-euclidean_$euclidien"
    if isfile(filename)
        include(filename)
        time_CW = @elapsed begin 
            sol_CW = real_cost_smart(robust_clark_wright(n, t, th, d, C, max_cost=true),n,th,t,T)
        end
        time_LK = @elapsed begin
            sol_LK = real_cost_smart(sous_tours_heuristic(n, t, th, d, C, true, euclidien; two_opt=false, LK=true),n,th,t,T)
        end
        time_2opt = @elapsed begin
            sol_2opt = real_cost_smart(sous_tours_heuristic(n, t, th, d, C, true, euclidien; two_opt=true, LK=false),n,th,t,T)
        end
        time_3opt = @elapsed begin
            sol_3opt = real_cost_smart(sous_tours_heuristic(n, t, th, d, C, true, euclidien; two_opt=false, LK=false),n,th,t,T)
        end
        time_hybrid_max_cost = @elapsed begin
            sol_hybrid_max_cost = real_cost_smart(hybrid_heuristic(n, t, th, d, C, euclidien; max_cost=true, LK=false, real_cost=false, two_opt=false),n,th,t,T)
        end
        time_hybrid_real_cost = @elapsed begin
            sol_hybrid_real_cost = real_cost_smart(hybrid_heuristic(n, t, th, d, C, euclidien; max_cost=true, LK=false, real_cost=true, two_opt=false),n,th,t,T)
        end

        relache = relache_continu(n, th, t, d, C, T, time_limit=20)
        if relache !== nothing
            borne_inf = relache[2]
        else
            borne_inf = robust_opt(n, th, t, d, C, T, time_limit=20)[2]
        end

        sol_values = [sol_CW, sol_LK, sol_2opt, sol_3opt, sol_hybrid_real_cost, sol_hybrid_max_cost]
        times = [time_CW, time_LK, time_2opt, time_3opt, time_hybrid_real_cost, time_hybrid_max_cost]
        return times, sol_values, borne_inf
    else
        return nothing
    end
    
end

function comparaison_heuristiques(euclidien::Bool)
    filename_time = "results_heuristiques/time_heuristic_$euclidien.csv"
    filename_value = "results_heuristiques/value_heuristic_$euclidien.csv"
    
    file_time = open(filename_time, "w")
    println(file_time, "n,CW,LK,2opt,3opt,hybrid_rc,hybrid_mc")
    
    file_value = open(filename_value, "w")
    println(file_value, "n,CW,LK,2opt,3opt,hybrid_rc,hybrid_mc,borne_inf")
    
    for n in 6:100        
        res = comparaison_heuristiques_n(n, euclidien)
        if res !== nothing
            times, sol_values, borne_inf = res

            println(file_time, "$n,$(times[1]),$(times[2]),$(times[3]),$(times[4]),$(times[5]),$(times[6])")
            flush(file_time)
            
            println(file_value, "$n,$(sol_values[1]),$(sol_values[2]),$(sol_values[3]),$(sol_values[4]),$(sol_values[5]),$(sol_values[6]),$borne_inf")
            flush(file_value)
        end
    end
    
    close(file_time)
    close(file_value)
end

function to_csv(time_limit::Float64)

    file_name_gap = "results_exact_methods/results_gap.csv"
    file_name_time = "results_exact_methods/results_time.csv"

    file_gap = open(file_name_gap, "w")
    println(file_gap, "n,euclidien,plans coupants,B&C,B&C with heuristic,dual,dual with warm start")

    file_time = open(file_name_time, "w")
    println(file_time, "n,euclidien,plans coupants,B&C,B&C with heuristic,dual,dual with warm start")

    for i in 5:100
        for euclidien in ["true","false"]
            instance = "data/n_$i-euclidean_$euclidien"
            if !isfile(instance)
                continue
            end
            
            # plans coupants classique
            sol_pc = plans_coupants(instance,slave_heur=false,time_limit=time_limit)
            gap_pc = 1 - sol_pc[2]/sol_pc[1]
            time_pc = sol_pc[3]

            # plans coupants heuristic
            #sol_pch = plans_coupants(instance,slave_heur=true,time_limit=time_limit)
            #gap_pch = 1 - sol_pch[2]/sol_pch[1]
            #time_pch = sol_pch[3]

            # B&C classique
            sol_bc = branch_and_cut(instance,slave_heur=false,time_limit=time_limit)
            gap_bc = 1 - sol_bc[2]/sol_bc[1]
            time_bc = sol_bc[3]

            # B&C heuristic
            sol_bch = branch_and_cut(instance,slave_heur=true,time_limit=time_limit)
            gap_bch = 1 - sol_bch[2]/sol_bch[1]
            time_bch = sol_bch[3]

            # dual
            sol_dual = dual(instance,time_limit=time_limit)
            gap_dual = 1 - sol_dual[2]/sol_dual[1]
            time_dual = sol_dual[3]

            # dual with warm start
            sol_dual_ws = dual(instance,warm_start=true,time_limit=time_limit)
            gap_dual_ws = 1 - sol_dual_ws[2]/sol_dual_ws[1]
            time_dual_ws = sol_dual_ws[3]

            println(file_gap,"$n,$euclidien,$gap_pc,$gap_bc,$gap_bch,$gap_dual,$gap_dual_ws")
            println(file_time,"$n,$euclidien,$time_pc,$time_bc,$time_bch,$time_dual,$time_dual_ws")

            flush(file_gap)
            flush(file_time)
        end
    end
    close(file_gap)
    close(file_time)
end