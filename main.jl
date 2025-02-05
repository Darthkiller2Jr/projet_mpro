include("dual.jl")
include("plans_coupants.jl")

function main(time_limit)
    println("Dualisation")
    main_dual(time_limit)

    println("Algorithme de plans coupants")
    main_plans_coupants(time_limit)

    println("Algorithme de branch and cut sans heuristique")
    main_branch_and_cut(time_limit, false)

    println("Algorithme de branch and cut avec heuristique")
    main_branch_and_cut(time_limit, true)
end