struct BirthDeathShiftDiscrete
    λ::Vector{Float64}
    μ::Vector{Float64}
    η::Float64
end

η = 0.05

function discrete_extinct_prob(dE, E, p, t)
    model, K = p
    λ = model.λ
    μ = model.μ
    η = model.η

    sumE = sum(E)

    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E .* E .+ (η/(K-1)) .* (sumE .- E) 
    #dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E .* E .+ (η/K) .* (sumE) 
end

function extinction_probability(
    model::BirthDeathShiftDiscrete,
    state::Int64,
    t::Float64,
    ) 
    alg = OrdinaryDiffEq.Tsit5()

    K = length(model.λ)
    pE = (model, K)

    tspan = (0.0, t)
    E0 = repeat([0.0], K)
    
    pr = OrdinaryDiffEq.ODEProblem(discrete_extinct_prob, E0, tspan, pE);
    
    E = OrdinaryDiffEq.solve(pr, alg)
    res = E[end][state]
    return(res)
end

function make_quantiles(d, k)
    quantiles = zeros(k)
    step = 0.5
    for i in 1:k
        p = (i-step)/k
        quantiles[i] = Distributions.quantile(d, p)
    end
    return(quantiles)
end


using Distributions

function discrete_ext(n)
    K = 2*n + 1
    μ = 0.2
    λmean = 0.3

    dλ = Exponential(λmean)

    λquantiles = make_quantiles(dλ, K)

    discrete_model = BirthDeathShiftDiscrete(λquantiles, [μ for _ in 1:K], η)

    E = extinction_probability(discrete_model, n, 5.0)
    return(E)
end 

mid = quantile(Exponential(0.3), 0.5)


continuous_model = BDSExponential(0.3, 0.2, η)
E_continuous = extinct_probability(mid, 5.0, continuous_model)

#ns = collect(range(1, 1000; step = 1))
ns = [1, 2, 4, 5, 8, 15, 50, 200, 1000]
ks = [(n*2)+1 for n in ns]
Es = [discrete_ext(n) for n in ns]

fig = Figure(size = (400, 400));
ax = Axis(fig[1,1], xscale = log10, 
    xlabel = L"\text{number of rate categories (K)}",
    ylabel = L"\text{extinction probability E}(\lambda=0.207,t=5.0)",
    title = L"\hat{\lambda} = 0.3, \mu = 0.2, \eta = %$η",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    )
scatter!(ax, ks, Es, color = :black, label = "discretized (Pesto/RB)")
pts = [(ki, Ei) for (ki, Ei) in zip(ks, Es)]
text_labels = [string(x) for x in ks]
text!(ax, pts, text = text_labels)
lines!(ax, ks, Es, color = :black)
lines!(ax, [extrema(ks)...], [E_continuous, E_continuous], label = "continous (1st order approx)")
lines!(ax, [extrema(ks)...], [extinction_frequency, extinction_frequency], label = "simulation frequency (2 mil repl)")
axislegend(ax)
#ylims!(ax, (0.0, 1.0))
xlims!(ax, (2, maximum(ks)*2.3)) 

fig


save("figures/validation-single.pdf", fig)



