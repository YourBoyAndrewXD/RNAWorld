using ForwardDiff
using OrdinaryDiffEq
using DataInterpolations
using Symbolics
using ModelingToolkit
import CasADi
using LinearAlgebra
using Optim
using Plots
import Symbolics: derivative

# Constants
k1, k2, k3, k4, k5 = 4.55e-3, 0.183, 0.167, 0.1, 3.33e-3
const c6 = 1.06
const Kd = 24.7

t = ModelingToolkit.t_nounits
@variables α(t) Cp(t) αCp(t) CppC(t) RNA(..)
@variables V(t) [input = true, bounds = (0.1, 10.)]
@variables λ₁ λ₂ λ₃ λ₄ λ₅

const ∂ = Symbolics.derivative
const r1 = k1*αCp^2 / V^2
const r2 = k2*α*CppC / V^2
const r3 = k3*CppC / V
const r4 = k4*α*Cp / V^2
const r5 = k5*αCp / V
const r6 = c6 * (1 + (Cp + αCp) / (Kd * V)) * CppC / V

const H = V * (
        λ₁ * (r1 + r5 - r2 - r4) +
        λ₂ * (r3 + r5 - r4) +
        λ₃ * (2r2 - 2r1 + r3 + r4 - r5) +
        λ₄ * (r1 - r2 - r3 - r6) + 
        λ₅ * r6
    )
const λₜ = -[∂(H, α), ∂(H, Cp), ∂(H, αCp), ∂(H, CppC), ∂(H, RNA(t))]
const nₜ = [∂(H, λ₁), ∂(H, λ₂), ∂(H, λ₃), ∂(H, λ₄), ∂(H, λ₅)]

λ = [λ₁, λ₂, λ₃, λ₄, λ₅]
n = [α, Cp, αCp, CppC, RNA(t)]

Hamiltonian = Symbolics.build_function(H, n, λ, V; expression = Val{false})
adjoint_dynamics, adjoint_dynamics! = Symbolics.build_function(λₜ, n, λ, V; expression = Val{false})
state_dynamics, state_dynamics! = Symbolics.build_function(nₜ, n, V; expression = Val{false})

# Optimal V(t) maximizes H(n, λ, V)
function optimal_V(n, λ)
    obj(V) = -Hamiltonian(n, λ, V)  # negative for maximization
    result = optimize(obj, 1.0, 100.0)
    return Optim.minimizer(result)
end

##########################################
### Solve Pontryagin Maximum Principle ###
##########################################

function ForwardBackwardSweep(init_controller, solver, x0, λ0, tspan; tol = 1e-5, maxiters = 1000)
    # x_i, λ_i, u_i 
    # x_{i+1}, λ_{i+1}, u_{i+1}
    
    uᵢ = init_controller
    function state!(du, u, p, t)
        state_dynamics!(du, u, uᵢ(t))
    end
    prob_fwd = ODEProblem(state!, x0, tspan)
    xᵢ = solve(prob_fwd, solver; dtmax = 1.)

    function costate!(dλ, λ, p, t)
        adjoint_dynamics!(dλ, xᵢ(t), λ, uᵢ(t))
    end
    prob_bwd = ODEProblem(costate!, λ0, (tspan[2], tspan[1]))
    λᵢ = solve(prob_bwd, solver; dtmax = 1.)
    uᵢ = construct_controller(xᵢ, λᵢ)

    for sweep in 1:maxiters
        println("Sweep $sweep")
        x_prev = xᵢ
        λ_prev = λᵢ
        u_prev = uᵢ

        xᵢ = solve(prob_fwd, solver; dtmax = 1.)
        λᵢ = solve(prob_bwd, solver; dtmax = 1.)
        uᵢ = construct_controller(xᵢ, λᵢ)

        below_tol(xᵢ, x_prev; tol) && below_tol(λᵢ, λ_prev; tol) && below_tol(uᵢ, u_prev; tol) && break
    end
    return uᵢ, xᵢ, λᵢ
end

function plot_H(n, λ)
    vals = [-Hamiltonian(n, λ, V) for V in 0.1:0.1:10.]
    plot(0.1:0.1:10., vals)
end

function below_tol(sol, sol_prev; tol = 1e-4)
    tsteps = sol.t
    err = 0
    for (i, t) in enumerate(tsteps)
        err += norm(sol(t) - sol_prev(t))
    end
    @show err
    err < tol
end

function construct_controller(x, λ)
    u_vals = [optimal_V(x(t), λ(t)) for t in x.t]
    ConstantInterpolation(u_vals, x.t; extrapolation = ExtrapolationType.Extension)
end
