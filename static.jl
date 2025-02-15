using JuMP
using CPLEX

# resout le pb statique de l'instance file avec un temps limite
function static(file, time_limit =10)
    include(file)

    m = Model(CPLEX.Optimizer)

    set_silent(m)

    # Variables
    @variable(m, x[i in 1:n, j in 1:n], Bin)  # x_ij ∈ {0,1}
    @variable(m, u[i in 2:n] >= 0)  # u_i ≥ 0

    set_time_limit_sec(m, time_limit)

    # Contraintes 
    @constraint(m, [i in 2:n], sum(x[j,i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[i,1] for i in 2:n) == sum(x[1,j] for j in 2:n))

    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n,  j != i], u[j] - u[i] >= d[i] - C * (1 - x[i,j]))
    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1,j]))
    @objective(m, Min, 
    sum(t[i,j] * x[i,j] for i in 1:n, j in 1:n if j!=i))

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
        return round(JuMP.objective_value(m)), bound, comput_time, vX
    end
end

# la fonction résout le problème statique sur toutes les instances avec la limite de temps indiquée.
function main_static(time_limit = 10)


    name_results = "resultats_static_"*string(time_limit)*"s.txt"
    name_solution = "solutions_static.txt"

    results_file = open("results/"*name_results, "w")
    sol_file = open("results/"*name_solution, "w")

    println(results_file, "file \t comput time \t limit time \t val \t gap")
    nb_resolue = 0
    for file in readdir("data")
        file_name = "data"*"/"*file
        # résolution
        val, bound, comput_time, x = static(file_name, time_limit)
        gap =(val-bound)/bound*100
        # écriture dans le fichier
        println(results_file, file_name,"\t", comput_time, "\t", time_limit,"\t", val, "\t", gap)

        if gap<1e-2
            # écriture de la solution si résolue.
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
    