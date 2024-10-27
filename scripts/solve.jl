using OrdinaryDiffEq

struct BDSExponential
    λmean::Float64 ## mean of base distribution for speciation rate
    μ::Float64 ## constant extinction rate
    η::Float64 ## shift rate
end

function ode(du, u, p, s)
    model, t1, λ1 = p
    μ = model.μ

    λ = find_lambda(model, s, t1, λ1)
    
    # du[1] is the extinction probability at E(r,s)
    du[1] = μ - (λ + μ) * u[1] + λ * u[1]^2
end

function ode2(du, u, p, t)
    σ, η = p
    # du[1] is d\lambda / dt
    du[1] = η * (u[1] - (1/σ))
end

function find_initial_lambda_ode(model, t::Float64, λ1::Float64)
    σ = 1.0 / model.λmean
    η = model.η

    tspan = (t, 0.0)
    p = (σ, η)

    u1 = [λ1]

    prob = ODEProblem(ode2, u1, tspan, p)
    sol = solve(prob)
    return(sol.u[end][1])
end


function find_initial_lambda(model, t::Float64, λ1)
    σ = 1 / model.λmean
    η = model.η

    #λ0 = model.λmean * (1 + exp((0.0 - t)/η)*(σ*λ1 - 1.0))
    λ0 = model.λmean * (1 + exp((0.0 - t)*η)*(σ*λ1 - 1.0))

    return(λ0)
end

function find_lambda(model, t::Float64, t1::Float64, λ1)
    σ = 1 / model.λmean
    η = model.η

    λ0 = model.λmean * (1 + exp((t - t1)*η)*(σ*λ1 - 1.0))

    return(λ0)
end

find_lambda(model, 0.0, 5.0, 0.207)

find_initial_lambda(models[1], 0.1, 0.5)

function extinct_probability(
    λ1::Float64, ## the speciation rate
    t1::Float64, ## the time span
    model::BDSExponential, ## model
    )
    λmean = model.λmean

    tspan = (0.0, t1)

    p = (model, t1, λ1)

    λ0 = find_initial_lambda_ode(models[1], t1, λ1)
    z0 = 0.0

    u0 = [0.0]

    prob = ODEProblem(ode, u0, tspan, p)

    sol = solve(prob, Tsit5())
    E = sol.u[end][1]
    return(E)
end


find_initial_lambda(models[2], 1.0, 0.5)
find_initial_lambda_ode(models[2], 1.0, 0.5)

sol = extinct_probability(0.0, 5.0, models[2])

find_initial_lambda(models[2], 1.0, 0.4)

using CairoMakie
using LaTeXStrings


extinct_probability(0.5, 10.0, models[2])

λs = collect(range(0.0, 2.0; length = 20))
ts = collect(range(0.0, 10.0; length = 20))

models = [
    BDSExponential(0.3, 0.2, 0.0),
    BDSExponential(0.3, 0.2, 1.0),
    BDSExponential(0.3, 0.2, 0.1),
]

fig = Figure(size = (650, 500));
## (left, right, bottom, top)
protrusions = [
    (30,50, 50, 10),
    (30,50, 50, 30), 
    (30,60, 50, 0),
]

axs = []

is = [
    1 3
    2 missing
]

surfaces = []

for (i, j) in Iterators.product(1:2, 1:2)
    idx = is[i,j]
    println(i, ", ", j, ", ", idx)
    if !ismissing(idx)
        model = models[idx]
        η = model.η
        Es = [extinct_probability(λi, ti, model) for (λi, ti) in Iterators.product(λs, ts)]

        ax = Axis3(fig[i,j],
        xlabel = L"\text{speciation rate }(λ)",
        ylabel = L"\text{time }(t)",
        zlabel = L"E(\lambda,t)",
        title = L"\hat{\lambda} = 0.3, \mu = 0.2, \eta = %$(η)",
        azimuth = -0.3*π, 
        viewmode = :stretch,
        protrusions = protrusions[idx],
        )
        zlims!(ax, 0.0, 1.0)

        s = surface!(ax, λs, ts, Es, colormap = :heat, colorrange = (0.0, 1.0))
        
        push!(surfaces, s)
        push!(axs, ax)
    end
end

rowgap!(fig.layout, 0.0)
colgap!(fig.layout, 0.0)

for i in 1:2
    colsize!(fig.layout, i, Relative(0.5))
    rowsize!(fig.layout, i, Relative(0.5))
end

fig

## draw a few lines
## along the curves
## z(r, t)
ls = []
for (i, model) in enumerate(models)
    for λ0 in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        model = models[i]

        #r = λ0
        t1 = maximum(ts)
        r0 = find_initial_lambda_ode(model, t1, λ0)

        #t1 = 5.0
        r1 = model.λmean * (1.0 + exp(-t1 * model.η) * ((1/model.λmean)*λ0 - 1.0))
        #r0 = model.λmean * (1.0 + exp(-t0 * model.η) * ((1/model.λmean)*λ0 - 1.0))
        #r0 = find_initial_lambda(model, t1, λ0)

        #rs = [r for _ in 1:length(ts)]
        rs = collect(range(r0, r1; length = length(ts)))

        λs_r = [model.λmean * (1 + ((1/model.λmean)*ri - 1.0)*exp(ti*model.η)) for (ri, ti) in zip(rs, ts)]
        Es = [extinct_probability(λi, ti, model) for (λi, ti) in zip(λs_r, ts)]

        #Es = [extinct_prob_r(ri, ti, model) for (ri, ti) in zip(rs, ts)]

        l = lines!(axs[i], λs_r, ts, Es, color = :black, linestyle = :dash)
        push!(ls, l)
    end
end

fig

ax4 = Axis(fig[2,2],
    xminorgridvisible = false,
    yminorgridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
    leftspinevisible = false,
    bottomspinevisible = false,
    )
ylims!(ax4, 0.0, 0.3)
hidedecorations!(ax4)


#elem_1 = PolyElement(color = "#3d2e72", strokewidth = 1)
elem_1 = PolyElement(color = :orange, strokewidth = 1)
elem_2 = LineElement(color = :black, strokewidth = 1, linestyle = :dash)

Legend(fig[2,2], 
    [elem_1, elem_2],
    [L"\text{E}(\lambda, t)", L"z(r, s)",], 
    position = :cc,
    labelsize = 30,
    patchsize = (40 ,20))

text!(
    ax4, 
    [0.0],
    [0.04],
    text = L"\text{E}(\lambda,t) = z(r(\lambda,t), s(t))",
    fontsize = 23,
    align = (:center, :center),
)

fig

save("figures/characteristic-lines.png", fig)

