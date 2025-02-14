using JuMP
using CPLEX
include("heuristiques.jl")

function dual(file::String; warm_start::Bool=false, time_limit = 30.0)
    include(file)

    # 
    m = Model(CPLEX.Optimizer)

    # Limite de temps
    
    # MOI.set(m, MOI.NumberOfThreads(), 1)
    set_silent(m)
    # # Désactive le presolve (simplification automatique du modèle)
    # set_optimizer_attribute(m, "CPXPARAM_Preprocessing_Presolve", 0)
    # # Désactive la génération de coupes automatiques
    # set_optimizer_attribute(m, "CPXPARAM_MIP_Limits_CutsFactor", 0)
    # # Désactive la génération de solutions entières à partir de solutions
    # # fractionnaires
    # set_optimizer_attribute(m, "CPXPARAM_MIP_Strategy_FPHeur",-1)
    # # Désactive les sorties de CPLEX (optionnel)
    # set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

    # Definition des variables du problème dual
    @variable(m, x[i in 1:n, j in 1:n], Bin)  # x_ij ∈ {0,1}
    @variable(m, u[i in 2:n] >= 0)  # u_i ≥ 0
    @variable(m, α1 >= 0)
    @variable(m, α2 >= 0)
    @variable(m, β1[i in 1:n, j in 1:n] >= 0)
    @variable(m, β2[i in 1:n, j in 1:n] >= 0)

    time_used=0
    euclidean = occursin("true",file)
    if warm_start
        start_warm = time()
        # heuristic_routes = hybrid_heuristic(n,t,th,d,C,T,true)
        # heuristic_x = routes_to_x(heuristic_routes, n)
        sol_3opt = sous_tours_heuristic(n, t, th, d, C, true, euclidean; two_opt=false, LK=false)
        heuristic_x = routes_to_x(sol_3opt, n)
        time_used = time()-start_warm
        for i in 1:n
            for j in 1:n
                set_start_value(x[i,j], get(heuristic_x, (i,j), 0.0))
            end
        end
    end

    set_time_limit_sec(m, max(1,time_limit-time_used))
    # Fonction objectif
    @objective(m, Min, 
    sum(t[i,j] * x[i,j] + β1[i,j] + 2 * β2[i,j] for i in 1:n, j in 1:n if j!=i) + α1 * T + α2 * T^2
    )

    # Contraintes
    @constraint(m, [i in 2:n], sum(x[j,i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[i,1] for i in 2:n) == sum(x[1,j] for j in 2:n))

    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n,  j != i], u[j] - u[i] >= d[i] - C * (1 - x[i,j]))

    @constraint(m, [i in 1:n, j in 1:n,  j != i], α1 + β1[i,j] >= (th[i] + th[j]) * x[i,j])
    @constraint(m, [i in 1:n, j in 1:n,  j != i], α2 + β2[i,j] >= th[i]th[j] * x[i,j])

    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1,j]))

    start = time()
    optimize!(m)
    comput_time = time()- start
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
        # Récupération des valeurs d’une variable
        vX = JuMP.value.(x)
        bound = ceil(JuMP.objective_bound(m))
        println("File: ", file, "\t Valeur de l’objectif : ", JuMP.objective_value(m), "\t Meilleure borne : ", bound, "\t time : ", comput_time) 
        return round(JuMP.objective_value(m)), bound, comput_time+time_used, vX
    end
end

function main_dual(time_limit = 10, warm_start = false)
    if warm_start
        name_results = "resultats_dual_ws_"*string(time_limit)*"s.txt"
        name_solution = "solutions_1_duales_ws.txt"
    else
        name_results = "resultats_dual_"*string(time_limit)*"s.txt"
        name_solution = "solutions_1_duales.txt"
    end

    results_file = open("results/"*name_results, "w")
    sol_file = open("results/"*name_solution, "w")

    println(results_file, "file \t comput time \t limit time \t val \t gap")
    nb_resolue = 0
    for file in readdir("data")
        file_name = "data"*"/"*file
        val, bound, comput_time, x = dual(file_name; warm_start=warm_start, time_limit=time_limit)
        gap =(val-bound)/bound*100
        println(results_file, file_name,"\t", comput_time, "\t", time_limit,"\t", val, "\t", gap)

        if gap<1e-2
            println(sol_file, file_name, "*********************************************************************")
            for i in findall(x.> 1-1e-4)
                print(sol_file, "(",i[1],",", i[2],")"," ")
            end
            println(sol_file)
            nb_resolue +=1
        end
    end   
    close(results_file)
    println(sol_file, "nb instances résolues : ", nb_resolue)
    close(sol_file)
end
    