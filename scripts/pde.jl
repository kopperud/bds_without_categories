using ModelingToolkit
using MethodOfLines
using OrdinaryDiffEq
using DomainSets
using CairoMakie

@parameters λ t
@variables E(..) D(..)

Dt = Differential(t)
Dx = Differential(λ)
Dxx = Differential(λ)^2

# domains for speciation rate and time
domains = [
           λ ∈ Interval(0.0, 1.5),
           t ∈ Interval(0.0, 15.0)
          ]

bcs = [
       Dx(E(0.0, t)) ~ 0,
       Dx(E(1.5, t)) ~ 0,
       Dx(D(0.0, t)) ~ 0,
       Dx(D(1.5, t)) ~ 0,
       E(λ, 0) ~ 0,
       D(λ, 0) ~ 1,
      ]

μ = 0.4
σ = 0.12

ηs = [0.0001, 0.001, 0.01, 0.05]

sols = []

plts = []

for η in ηs
    eq = [
          ## extinction probability
        Dt(E(λ,t)) ~ μ + λ * E(λ,t)^2 -
        (η + λ + μ) * E(λ,t) +
        η * E(λ,t) +
        Dx(E(λ,t)) * (1/σ - λ) * η +
        η * Dxx(E(λ,t)) * (1/(σ^2) - λ/σ + λ^2),
          ## observation probability
        Dt(D(λ,t)) ~ -(λ+μ+η) * (D(λ,t)) +
        2 * λ * E(λ,t) * D(λ,t) +
        η * D(λ,t) +
        η * Dx(D(λ,t)) * (1/σ - λ) +
        η * Dxx(D(λ,t)) * (1/(σ^2) - λ/σ + λ^2)
       ]



    @named pdesys = PDESystem(eq, bcs, domains, [λ,t],[E(λ,t), D(λ,t)])

    dλ = 0.05
    order = 2

    discretization = MOLFiniteDifference([λ => dλ], t)

    # convert the PDE into an ODE
    prob = discretize(pdesys, discretization)

    ## solve the ODE problem
    #sol = solve(prob, Tsit5(), saveat =1)
    sol = solve(prob, Tsit5())
    push!(sols, sol)

end
#=
    using Plots
    plt = heatmap(
        discrete_t, discrete_λ, solu, 
        xlabel = "time before present (t)",
        ylabel = "speciation rate (λ)",
       xflip = true)
    push!(plts, plt)

p = plot(plts...,
    title = "extinction probability, E(λ,t)",
    size = (800, 400)
        )
p
=#

sols1 = deepcopy(sols)

fig = Figure()
axs = []
hms = []

for i in 1:2
    for j in 1:2
        eta = pop!(ηs)
        ax = Axis(fig[i,j],
                  xlabel = "time before present (t)",
                  ylabel = "speciation rate (λ)",
                  xreversed = true,
                  title = string("eta = ", eta),
                 )
        
        s = sols[q]

        # Plot results and compare with exact solution
        #sol = sols[q]
        sol = pop!(sols)
        discrete_λ = sol[λ]
        discrete_t = sol[t]
        solu = sol[E(λ, t)]

        #hm = heatmap!(ax, discrete_t, discrete_λ, log.(solu'))
        hm = heatmap!(ax, discrete_t, discrete_λ, solu')
        #contour!(ax, discrete_t, discrete_λ, solu')

        
        push!(axs, ax)
        push!(hms, hm)
    end
end

linkaxes!(axs...)

cb = Colorbar(fig[:, end+1], hms[1], label = "Extinction probability, E(λ,t)")
save("figures/extinction_prob.pdf", fig)
#cb = Colorbar(fig[:, end+1], hms[1], label = "Observation probability, log(D(λ,t))")
#save("figures/observation_prob.pdf", fig)


fig
