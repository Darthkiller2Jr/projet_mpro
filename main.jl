include("dual.jl")
include("plans_coupants.jl")

function main(limite_dual, limite_pc)
    println("Dualisation")
    main_dual(limite_dual)

    println("Algorithme de plans coupants")
    main_plans_coupants(limite_pc)
end