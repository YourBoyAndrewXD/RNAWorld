include("rna_model.jl")
using OptimizationBBO
using OrdinaryDiffEqVerner
using OrdinaryDiffEqRosenbrock

init_controller = Returns(1.)
const V_lo = 0.1
const V_hi = 10.
n0 = [10.0, 10.0, 10.0, 10.0, 0.0]    # initial state
λ0 = [0.0, 0.0, 0.0, 0.0, 1.0]   # final costate
ks = [4.55e-3, 0.183, 0.167, 0.1, 3.33e-3]
p = transpose(hcat(n0, λ0))
tspan = (0.0, 20.0)

controller, state_trajectory, costate_trajectory = ForwardBackwardSweep(init_controller, Rodas5P(), n0, λ0, ks, tspan; maxiters = 200);

function plot_rna_trajectories(ks, V_hi, V_lo, x0, λ0)
    ode_func_hi!(dy, y, p, t) = state_dynamics!(dy, y, ks, V_hi)
    prob_hi = ODEProblem(ode_func_hi!, x0, tspan)
    sol_hi = solve(prob_hi, Rodas5P())

    ode_func_low!(dy, y, p, t) = state_dynamics!(dy, y, ks, V_lo)
    prob_low = ODEProblem(ode_func_low!, x0, tspan)
    sol_low = solve(prob_low, Rodas5P())

    controller, sol_optim, sol_co = ForwardBackwardSweep(init_controller, Vern7(), x0, λ0, ks, tspan)

    p1 = plot(controller, label="True V(t)", lw = 2, size = (1200, 900))
    p2 = plot(sol_hi; idxs = 5, label = "RNA (optimal control)", lw=2, xlabel="Time", ylabel="RNA",
          title="RNA Production: Optimal V(t) vs Constant V", legend=:bottomright, size = (1200, 900))
    plot!(p2, sol_low; idxs = 5, label = "RNA (V = $V_lo)", lw = 2, linestyle = :dash)
    plot!(p2, sol_optim; idxs = 5, label = "RNA (V = $V_hi)", lw = 2, linestyle = :dash)
    plot(p1, p2)
end


####################################
###### Parameter optimization ######
####################################

ode_func_low!(du, u, p, t) = state_dynamics!(du, u, ks, V_lo)
ode_func_hi!(du, u, p, t) = state_dynamics!(du, u, ks, V_hi)

function run_trial(ks, ps)
    ks = exp.(ks)
    x0, λ0 = ps[1, :], ps[2, :]
    controller, sol_optim, sol_co = ForwardBackwardSweep(init_controller, Vern7(), x0, λ0, ks, tspan)
    prob_low = ODEProblem(ode_func_low!, x0, tspan)
    sol_lo = solve(prob_low, Vern7())

    return sol_lo.u[end][end] - sol_optim.u[end][end]
end

prob = OptimizationProblem(run_trial, log.(ks), p, lb = fill(-7, 5), ub = fill(7, 5))
sol = solve(prob, BBO_xnes())
