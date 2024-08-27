using ModelingToolkit
using MethodOfLines
using OrdinaryDiffEq
using DomainSets
using CairoMakie

function ode_problem(μ, σ, η, λ_lower, λ_upper; n = 30)
    @parameters λ t
    @variables E(..) D(..)

    Dt = Differential(t)
    Dx = Differential(λ)
    Dxx = Differential(λ)^2
    #Dxxx = Differential(λ)^3

    t0 = 0.0
    t1 = 1.0

    # domains for speciation rate and time
    domains = [
               λ ∈ Interval(λ_lower, λ_upper),
               t ∈ Interval(t0, t1)
              ]

    # boundary conditions
    # we assume that the
    # derivative of E and D with respect to
    # the speciation rate is zero at the boundaries
    # (Neumann boundary condition)
    bcs = [
           Dx(E(λ_lower, t)) ~ 0,
           Dx(E(λ_upper, t)) ~ 0,
           Dx(D(λ_lower, t)) ~ 0,
           Dx(D(λ_upper, t)) ~ 0,
           E(λ, 0) ~ 0,
           D(λ, 0) ~ 1,
          ]

    # the partial differential equations E(λ,t) and D(λ,t)
    eq = [
          ## extinction probability
        Dt(E(λ,t)) ~ μ + λ * E(λ,t)^2 -
        (η + λ + μ) * E(λ,t) +
        η * E(λ,t) +
        Dx(E(λ,t)) * (1/σ - λ) * η +
        η * Dxx(E(λ,t)) * (1/(σ^2) - λ/σ + λ^2),
        #η * Dxxx(E(λ,t)) * (1/(σ^3) - λ^2/σ + λ/(σ^2) - λ^3),
          ## observation probability
        Dt(D(λ,t)) ~ -(λ+μ+η) * (D(λ,t)) +
        2 * λ * E(λ,t) * D(λ,t) +
        η * D(λ,t) +
        η * Dx(D(λ,t)) * (1/σ - λ) +
        η * Dxx(D(λ,t)) * (1/(σ^2) - λ/σ + λ^2)
        #η * Dxxx(D(λ,t)) * (1/(σ^3) - λ^2/σ + λ/(σ^2) - λ^3)
       ]

    @named pdesys = PDESystem(eq, bcs, domains, [λ,t],[E(λ,t), D(λ,t)])

    # convert the PDE into an ODE
    # by discretizing λ
    λ_span = λ_upper - λ_lower
    dλ = λ_span / (n-1)
    @assert dλ > 0.0
    λs = collect(range(λ_lower, λ_upper; step = dλ))

    discretization = MOLFiniteDifference([λ => dλ], t)
    #discretization = MOLFiniteDifference([λ => dλ], t; advection_scheme = WENOScheme())
    prob = discretize(pdesys, discretization)

    return(prob, λs)
end


function foo(prob, t0, t1, E_initial, D_initial)
    prob2 = remake(prob; tspan = (t0, t1))
    n_discretizations = Int64(length(prob.u0)/2)

    prob2.u0[1:n_discretizations] = E_initial
    prob2.u0[(n_discretizations+1):end] = D_initial

    # solve using Tsit5 ODE solver
    sol = solve(prob, Tsit5(),
                save_everystep = true
               )

    #D_sol = sol.u[D(λ, t)][:,end]
    #E_sol = sol.u[E(λ, t)][:,end]

    #res = hcat(E_sol, D_sol)
    
    return(sol)
end


μ = 0.4
σ = 0.12
η = 0.0001
λ_lower = 0.0001
λ_upper = 1.5
t0 = 0.0
t1 = 15.0

prob, λs = ode_problem(μ, σ, η, λ_lower, λ_upper; n = 50)
n_discretizations = Int64(length(prob.u0)/2)

E_initial = [0.0 for _ in 1:n_discretizations]
D_initial = [1.0 for _ in 1:n_discretizations]

sol = foo(prob, t0, t1, E_initial, D_initial)


using Pesto


notneg(u,p,t) = any(x -> x < 0, u)

function logl(data, σ, μ, η, K; λ_lower = 0.0, λ_upper = 2.5)
    prob, λs = ode_problem(μ, σ, η, λ_lower, λ_upper; n = K)

    nrows = size(data.edges, 1)
    Ntip = length(data.tiplab)
    descendants = make_descendants(data)

    root_index = Ntip+1
    left_edge, right_edge = descendants[root_index]

    local sf_left, sf_right
    local u_left, u_right


    u_left, sf_left = logl_po(left_edge, prob, data, descendants, Ntip, K, λs)
    u_right, sf_right = logl_po(right_edge, prob, data, descendants, Ntip, K, λs)

    E = u_left[1:K]
    D = u_left[(K+1):end] .* u_right[(K+1):end] .* λs

    c = sum(D)
    D = D ./ c
    sf = sf_left + sf_right

    if c > 0.0
        sf += log(c)
    else
        sf -= Inf
    end

    if !(sol.retcode == OrdinaryDiffEq.ReturnCode.Success)
        sf -= Inf
    end

    nonextinct = (1.0 .- E) .^2
    x = D ./ (nonextinct .* λs)

    p = sum(D) / K
    logl = log(p) + sum(sf)
    return(logl)
end

function logl_po(edge_index, prob, data::SSEdata, descendants, Ntip, K, λs)
    anc, dec = data.edges[edge_index,:]
    alg = OrdinaryDiffEq.Tsit5()
    node_age = data.node_depth[dec]
    parent_node_age = data.node_depth[anc]
    tspan = (node_age, parent_node_age)

    local sf_left, sf_right
    if dec <= Ntip
        E0 = zeros(K)
        D0 = ones(K)
        u0 = vcat(E0, D0)
        sf_left = 0.0
        sf_right = 0.0
    else
        left_edge, right_edge = descendants[dec]

        local u_left, u_right
        begin
            u_left, sf_left = logl_po(left_edge, prob, data, descendants, Ntip, K, λs)
            u_right, sf_right = logl_po(right_edge, prob, data, descendants, Ntip, K, λs)
        end

        E0 = u_left[1:K]
        D0 = u_left[(K+1):end] .* u_right[(K+1):end] .* λs
        u0 = vcat(E0, D0)
    end

    sf = sf_left + sf_right

    if isinf(sf) | isnan(sf) | in(u0, NaN)
        u = ones(K*2)
        sf -= Inf
    else
        prob = OrdinaryDiffEq.remake(prob, u0 = u0, tspan = tspan)
        #prob.u0 = u0

        sol = OrdinaryDiffEq.solve(prob, Tsit5(), save_everystep=true, 
                                   isoutofdomain = notneg, reltol = 1e-3)

        D, E = values(sol.u)
        n_knots = size(D)[2]
        D = D[:,end]
        E = E[:,end]

        ## small hack
        ## some times the boundary is negative, e.g. -1e-23, even if we have notneg?
        D[D .< 0] .= minimum(abs.(D))
        E[E .< 0] .= minimum(abs.(E))

        c = sum(D)
        D = D ./ c

        if c > 0.0
            sf += log(c)
        else
            sf -= Inf
        end

        if !(sol.retcode == OrdinaryDiffEq.ReturnCode.Success)
            sf -= Inf
        end

        u = vcat(E, D)
        Δt = tspan[2] - tspan[1]
        println("edge_index: $edge_index, \t scaling factor: $(round(sf, digits = 2)), n_knots: $n_knots", "  Δt: $(round(Δt, digits = 2)),   sum(D): $(round(c, digits = 5))")
    end


    return(u, sf)
end

## sanity check

K = 100
prob, λs = ode_problem(μ, σ, η, 0.00001, 30.0; n = K)

## one solve
u0 = vcat(zeros(K), ones(K))
prob1 = remake(prob, tspan = (0.0, 20.0), u0 = u0)
sol1 = OrdinaryDiffEq.solve(prob1, Tsit5(), save_everystep=false, 
                           isoutofdomain = notneg, reltol = 1e-8)
D, E = values(sol1.u)
D = D[:,end]
E = E[:,end]
res1 = vcat(E, D)

## two solves
## a) from 0.0 to 1.0
u0 = vcat(zeros(K), ones(K))
prob2 = remake(prob, tspan = (0.0, 10.0), u0 = u0)
sol2 = OrdinaryDiffEq.solve(prob2, Tsit5(), save_everystep=false, 
                           isoutofdomain = notneg, reltol = 1e-8)
D, E = values(sol2.u)
D = D[:,end]
E = E[:,end]
res2 = vcat(E, D)


## b) from 0.0 to 1.0
prob3 = remake(prob, tspan = (10.0, 20.0), u0 = res2)
sol3 = OrdinaryDiffEq.solve(prob3, Tsit5(), save_everystep=false, 
                           isoutofdomain = notneg, reltol = 1e-8)
D, E = values(sol3.u)
D = D[:,end]
E = E[:,end]
res3 = vcat(E, D)



tree = readtree(Pesto.path("primates.tre"))
data = SSEdata(tree, 1.0)

logl(data, 0.20, 0.13, 0.001, 100; λ_lower = 0.00001, λ_upper = 20.0)




