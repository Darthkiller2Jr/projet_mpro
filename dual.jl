using JuMP
using CPLEX

function dual(file::String, time_limit = 30.0)
    include(file)

    # 
    m = Model(CPLEX.Optimizer)

    # Limite de temps
    set_time_limit_sec(m, time_limit)
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
        println("File: ", file, "\t Valeur de l’objectif : ", JuMP.objective_value(m), "\t Meilleure borne : ", bound) 
        return round(JuMP.objective_value(m)), bound, comput_time
    end
end

function main_dual(time_limit = 10)
    name_results = "resultats_dual_"*string(time_limit)*"s.txt"
    results_file = open("results/"*name_results, "w")
    println(results_file, "file \t comput time \t limit time \t gap")
    for file in readdir("data")
        file_name = "data"*"/"*file
        val, bound, comput_time = dual(file_name, time_limit)
        gap =(val-bound)/bound*100
        println(results_file, file_name,"\t", comput_time, "\t", time_limit, "\t", gap)
    end   
    close(results_file)
end
