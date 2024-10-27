using Revise
using BirthDeathShiftSim
using Distributions


model = BirthDeathShiftSim.BirthDeathSpeciationExponential(0.3, 0.20, 0.05)
mid = quantile(Exponential(0.3), 0.5)

n_trees = 1_000_000

survivals = Bool[]
for _ in 1:n_trees
    try 
        tree = simulate(model, mid, 5.0, 10_000)

        for branch in tree.children
            if number_of_taxa(branch.outbounds) > 0
                push!(survivals, true)
            else
                push!(survivals, false)
            end
        end
    catch
        #push!(survivals, true)
    end
end

#n_taxa = [number_of_taxa(tree) for tree in trees]

survival_frequency = sum(survivals) / (length(survivals))
extinction_frequency = 1.0 - survival_frequency


model1 = BDSExponential(0.3, 0.20, 0.1)

E = extinct_probability(mid, 5.0, model1)[end][1]


using CairoMakie

fig = Figure();
ax = Axis(
    fig[1,1],
    xlabel = L"\text{x rate }(\lambda)",
    ylabel = L"\text{extinction, E}(\lambda, t)",
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
)

violin!(ax, [1.0], [extinction_frequency], label = L"\eta = 0.01", color = :gray)
lines!(ax, [0.5, 1.5], [E, E], linestyle = :dash, color = :red)

fig



