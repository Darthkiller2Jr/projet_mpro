using JuMP
using CPLEX

function master(n, C, d, U, time_limit)
    m = Model(CPLEX.Optimizer)

    # MOI.set(m, MOI.NumberOfThreads(), 1)
    set_silent(m)
    set_time_limit_sec(m, time_limit)

    # Variables
    @variable(m, x[i in 1:n, j in 1:n], Bin)  # x_ij ∈ {0,1}
    @variable(m, u[i in 2:n] >= 0)  # u_i ≥ 0
    @variable(m, z)

    # Fonction objectif
    @objective(m, Min, z)

    # Contraintes
    @constraint(m, [i in 2:n], sum(x[j,i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[i,1] for i in 2:n) == sum(x[1,j] for j in 2:n))

    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n,  j != i], u[j] - u[i] >= d[i] - C * (1 - x[i,j]))
    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1,j]))
    @constraint(m, [t_prim in U], sum(t_prim[i,j]*x[i,j] for i in 1:n, j in 1:n if j!=i) - z <= 0)

    optimize!(m)
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
        # Récupération des valeurs d’une variable
        vx = JuMP.value.(x)
        vz = JuMP.value.(z)
        return vz, vx
    end
end

function slave(n, x, t, th, T, time_limit = 30)
    m = Model(CPLEX.Optimizer)

    # MOI.set(m, MOI.NumberOfThreads(), 1)
    set_silent(m)
    set_time_limit_sec(m, time_limit)

    # Variables
    @variable(m, δ1[i in 1:n, j in 1:n] >=0)
    @variable(m, δ2[i in 1:n, j in 1:n] >=0)

    @objective(m, Max,
    sum((t[i,j] + δ1[i,j]*(th[i]+th[j])+ δ2[i,j]*th[i]*th[j])*x[i,j] for i in 1:n, j in 1:n if j!=i))

    @constraint(m, sum(δ1[i,j] for i in 1:n, j in 1:n if j!=i) <= T)
    @constraint(m, sum(δ2[i,j] for i in 1:n, j in 1:n if j!=i) <= T^2)
    @constraint(m, [i in 1:n, j in 1:n], δ1[i,j] <= 1)
    @constraint(m, [i in 1:n, j in 1:n], δ2[i,j] <= 2)

    optimize!(m)
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
        # Récupération des valeurs d’une variable
        vδ1 = JuMP.value.(δ1)
        vδ2 = JuMP.value.(δ2)
        t_star = Matrix(zeros(n,n))
        for i in 1:n
            for j in 1:n
                if j!=i
                    t_star[i,j] = t[i,j] + vδ1[i,j]*(th[i]+th[j])+ vδ2[i,j]*th[i]*th[j]
                end
            end
        end
        return JuMP.objective_value(m), t_star
    end
end

function slave_fast(n, x, t, th, T, time_limit = 1)

    indexes = [(i,j) for i in 1:n for j in 1:n]

    d1 = sort(indexes, by = idx -> (th[idx[1]]+th[idx[2]])*x[idx[1],idx[2]], rev = true)
    d2 = sort(indexes, by = idx -> th[idx[1]]*th[idx[2]]*x[idx[1],idx[2]], rev = true)
    # display([[(th[i]+th[j])*x[i,j] for i in 1:n] for j in 1:n])
    # display([[(th[i]*th[j])*x[i,j] for i in 1:n] for j in 1:n])
    # for i in 1:(n^2)
    #     println("d1 : ", d1[i], (th[d1[i][1]]+th[d1[i][2]])*x[d1[i][1],d1[i][2]] )
    #     println("d2 : ", d2[i], th[d2[i][1]]*th[d2[i][2]]*x[d2[i][1],d2[i][2]] )
    # end

    δ1 = Matrix(zeros(n,n))
    δ2 = Matrix(zeros(n,n))

    for k in 1:T
        if k <= length(d1)
            idx = d1[k]
            δ1[idx[1], idx[2]] = 1
        end
    end
    k = 1
    while k <= T*T/2 && k <= length(d2)
        idx = d2[k]
        δ2[idx[1], idx[2]] = 2
        k+=1
    end
    if mod(T*T,2)==1
        idx = d2[k]
        δ2[idx[1], idx[2]] = 1
    end

    t_star = Matrix(zeros(n,n))
    for i in 1:n
        for j in 1:n
            if j!=i
                t_star[i,j] = t[i,j] + δ1[i,j]*(th[i]+th[j])+ δ2[i,j]*th[i]*th[j]
            end
        end
    end
    z = sum(t_star[i,j]*x[i,j] for i in 1:n for j in 1:n)
    # display(δ1)
    # display(δ2)
    # println("__________________________________________________________________________")
    return z, t_star
end

function plans_coupants(file, slave_heur = true, time_limit = 30)
    include(file)
    # choix de la methode de résolution du sous pb
    if slave_heur
        slave_fn = slave_fast
    else
        slave_fn = slave
    end

    val_slave = 1
    val_z = 0

    U = Vector{Matrix{Float64}}(undef, 1)
    U[1] = t
    start = time()
    val_z, x_star = master(n,C,d,U, time_limit)
    s = time()
    val_slave, t_star = slave_fn(n,x_star,t,th,T, time_limit)

    # Condition d'arêt:
    #   optimalité si le sous problème donne une meilleure solution que le problème maître
    #   le temps limite est dépassé
    while val_slave > val_z && time()-start < time_limit
        
        # Résolution du problème maître
        time_left = time_limit - (time()-start)
        if time_left < 0.5
            break
        end
        val_z, x_star = master(n,C,d,U, time_left)
        # Résolution du sous-problème
        time_left = time_limit - (time()-start)
        if time_left < 0.5
            break
        end
            val_slave, t_star = slave_fn(n,x_star,t,th,T, time_left)
            
        # On ajoute la contrainte en augmentant U avec le pire scénario trouvé
        push!(U,t_star)
    end
    
    # la valeur du ss pb est la valeur d'une solution réalisable, la valeur du probleme maître une borne inf (relaxation de contraintes)
    println("File: ", file, "\t Valeur de l’objectif : ", val_slave, "\t Meilleure borne : ", val_z)
    return val_slave, val_z, time()-start
end

function branch_and_cut(file, time_limit = 1)
    print("File: ", file)
    include(file)

    # if slave_heur
    #     slave_fn = slave_fast
    # else
    #     slave_fn = slave
    # end

    m = Model(CPLEX.Optimizer)

    MOI.set(m, MOI.NumberOfThreads(), 1)
    # set_silent(m)
    set_time_limit_sec(m, time_limit)

    # Variables
    @variable(m, x[i in 1:n, j in 1:n], Bin)  # x_ij ∈ {0,1}
    @variable(m, u[i in 2:n] >= 0)  # u_i ≥ 0
    @variable(m, z)

    # Fonction objectif
    @objective(m, Min, z)

    # Contraintes
    @constraint(m, [i in 2:n], sum(x[j,i] for j in 1:n if j != i) == 1)
    @constraint(m, [i in 2:n], sum(x[i,j] for j in 1:n if j != i) == 1)
    @constraint(m, sum(x[i,1] for i in 2:n) == sum(x[1,j] for j in 2:n))

    @constraint(m, [i in 2:n], u[i] <= C - d[i])
    @constraint(m, [i in 2:n, j in 2:n,  j != i], u[j] - u[i] >= d[i] - C * (1 - x[i,j]))
    @constraint(m, [j in 2:n], u[j] <= C * (1 - x[1,j]))

    function callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
        # If the callback is called because a feasible integer solution is found
        if context_id == CPX_CALLBACKCONTEXT_CANDIDATE
            
            # This line must be called before getting the feasible integer solution
            CPLEX.load_callback_variable_primal(cb_data, context_id)

            x_star = callback_value.((cb_data,), x)
            val_z = callback_value(cb_data, z)

            val_slave, t_star = slave_fast(n, x_star, t, th, T)
            println("coupe ?")
            if val_z < val_slave
                cstr = @build_constraint(z >= sum(t_star[i,j]*x[i,j] for i in 1:n, j in 1:n if j!=i))
                MOI.submit(m, MOI.LazyConstraint(cb_data), cstr)
                println("coupe")
            end
        end
    end
    MOI.set(m, CPLEX.CallbackFunction(), callback)
    start = time()
    optimize!(m)
    comput_time = time()-start
    feasibleSolutionFound = primal_status(m) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(m) == MOI.OPTIMAL
    if feasibleSolutionFound
        bound = ceil(JuMP.objective_bound(m))
        println("\t Valeur de l’objectif : ", JuMP.objective_value(m), "\t Meilleure borne : ", bound) 
        return round(JuMP.objective_value(m)), bound, comput_time
    else
        println(file, "no feasible solutions found")
        return 0,0,comput_time
    end
end

function main_plans_coupants(time_limit = 10)
    name_results = "resultats_plans_coupants_"*string(time_limit)*"s.txt"
    results_file = open("results/"*name_results, "w")
    println(results_file, "file \t comput time \t limit time \t gap")
    for file in readdir("data")
        file_name = "data"*"/"*file
        val, bound, comput_time = plans_coupants(file_name, time_limit)
        gap =(val-bound)/bound*100
        println(results_file, file_name,"\t", comput_time, "\t", time_limit, "\t", gap)
    end   
    close(results_file)
end

function main_branch_and_cut(time_limit = 1)
    name_results = "resultats_bnc_"*string(time_limit)*"s.txt"
    results_file = open("results/"*name_results, "w")
    println(results_file, "file \t comput time \t limit time \t gap")
    for file in readdir("data")
        file_name = "data/"*file
        val, bound, comput_time = branch_and_cut(file_name, time_limit)
        gap =(val-bound)/bound*100
        println(results_file, file_name,"\t", comput_time, "\t", time_limit, "\t", gap)
    end   
    close(results_file)
end