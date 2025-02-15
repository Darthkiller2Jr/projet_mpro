include("dual.jl")
include("plans_coupants.jl")
include("static.jl")

function main(time_limit)

    println("Dualisation")
    main_dual(time_limit, false)

    println("Dualisation avec ws")
    main_dual(time_limit, true)

    println("Algorithme de plans coupants avec heuristique")
    main_plans_coupants(time_limit, true)

    println("Algorithme de branch and cut avec heuristique")
    main_branch_and_cut(time_limit, true,true)

    println("statique")
    main_static(time_limit)
end