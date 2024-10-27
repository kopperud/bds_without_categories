using CairoMakie
using Distributions
using LaTeXStrings


fig = Figure(size = (350, 250));

ax = Axis(
    fig[1,1],
    xlabel = L"\text{speciation rate }(\lambda)",
    ylabel = L"\text{probability density}",
    topspinevisible = false,
    rightspinevisible = false,
    xgridvisible = false,
    ygridvisible = false,
)

λmean = 0.3
σ = 1 / λmean

g(λ) = σ*exp(-σ*λ)

λs = collect(range(0.0, 1.0; length = 50))
y = [g(λi) for λi in λs]

lines!(ax, λs, y ,
    label = L"g(\lambda',t|\lambda, t + \Delta t)",
    color = :black)

lines!(ax, [λmean, λmean], [0.0, g(λmean)], linestyle = :dash, color = :gray, label = L"\text{mean }\lambda = 0.3")
axislegend(ax)

fig

save("figures/exponential.pdf", fig)




