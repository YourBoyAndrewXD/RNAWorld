include("SolutionVol1.2.jl")
using InfiniteOpt, Ipopt

### ModelingToolkit
D = ModelingToolkit.D_nounits
eqs = [D(s) ~ nₜ[i] for (i, s) in enumerate(n)]
cost = -[RNA(200.)]
@named sys = ODESystem(eqs, t; costs = cost)
sys = mtkcompile(sys; inputs = [V])
n0 = [0.0, 0.0, 100.0, 0.0, 0.0]    # initial state
prob = JuMPDynamicOptProblem(sys, Dict([n .=> n0; V => 1.]), (0., 200.); dt = 1.)
sol = solve(prob, JuMPCollocation(Ipopt.Optimizer))

p1 = plot(sol.sol)
p2 = plot(sol.input_sol)
plot(p1, p2, layout = (1, 2))
