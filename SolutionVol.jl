using DifferentialEquations
using ForwardDiff
using Optim
using OrdinaryDiffEq
using Plots
using Interpolations
using Polynomials
using Catalyst

# Constants
const k1 = 1
const k2 = 1
const k3 = 1
const k4 = 1
const k5 = 1
const c6 = 1.06
const Kd = 24.7

# keff and its partial with respect to V
function keff(n, V)
    Cp, αCp = n[2], n[3]
    return c6 * (1 + (Cp + αCp) / (Kd * V))
end

function dkeff_dV(n, V)
    Cp, αCp = n[2], n[3]
    return -c6 * (Cp + αCp) / (Kd * V^2)
end

# State dynamics
function state_dynamics!(du, n, p, t, λ, V)
    α, Cp, αCp, CppC, RNA = n
    r1 = k1 * αCp^2
    r2 = k2 * CppC * α
    r3 = k3 * CppC
    r4 = k4 * α * Cp
    r5 = k5 * αCp
    r6 = keff(n, V) * CppC

    du[1] = r1 - r2 - r4
    du[2] = r3 + r5 - r4
    du[3] = -2r1 + 2r2 + r3 + r4 - r5
    du[4] = r1 - r2 - r3 - r6
    du[5] = r6
end

# Adjoint dynamics
function adjoint_dynamics!(dl, λ, n, V)
    α, Cp, αCp, CppC, _ = n
    L1, L2, L3, L4, L5 = λ
    dk_dV = dkeff_dV(n, V)

    dl[1] = -V * (L1*(-k2*CppC - k4*Cp) + L2*(-k4*Cp) + L3*(2k2*CppC + k4*Cp) - L4*(k2*CppC + k4*Cp))
    dl[2] = -V * (L1*(-k4*α) + L2*(-k4*α) + L3*(k4*α) - L4*(k4*α)) + L5 * CppC * dk_dV
    dl[3] = -V * (L1*(-2k1*αCp) + L2*(k5) + L3*(4k1*αCp - k5) - L4*(2k1*αCp)) + L5 * CppC * dk_dV
    dl[4] = -V * (L1*(-k2*α) + L2*(-k3) + L3*(2k2*α + k3) - L4*(k2*α + k3 + c6*(1 + (Cp + αCp)/(Kd * V)))) - L5 * c6*(1 + (Cp + αCp)/(Kd * V))
    dl[5] = 0.0
end

# Hamiltonian
function H(n, λ, V)
    α, Cp, αCp, CppC, _ = n
    L1, L2, L3, L4, L5 = λ
    term = V * (
        L1 * (k1*αCp^2 - k2*α*CppC - k4*α*Cp) +
        L2 * (k3*CppC + k5*αCp - k4*α*Cp) +
        L3 * (-2k1*αCp^2 + 2k2*α*CppC + k3*CppC + k4*α*Cp - k5*αCp) +
        L4 * (k1*αCp^2 - k2*α*CppC - k3*CppC - c6*CppC*(1 + (Cp + αCp)/(Kd*V)))
    )
    extra = L5 * c6 * CppC * (1 + (Cp + αCp)/(Kd*V))
    return term + extra
end

# Optimal V(t) solver: maximizes H(n, λ, V)
function optimal_V(n, λ)
    obj(V) = -H(n, λ, V)  # negative for maximization
    result = optimize(obj, 0.1, 10.0)
    return Optim.minimizer(result)
end

# Initial conditions
# n0 = [α, Cp, αCp, CppC, RNA] (initial concentrations)
n0 = [0.0, 0.0, 0.05, 0.0, 0.0]   # state
λ0 = [0.0, 0.0, 0.0, 0.0, 1.0]   # costate

# Time setup
tspan = (0.0, 1000.0)
tsteps = range(tspan[1], tspan[2], length=100)
V_vals = fill(1.0, length(tsteps))
rna_vals = zeros(length(tsteps))

# Forward-backward sweep
for sweep in 1:20
    println("Sweep $sweep")

    function forward_ode!(du, u, p, t)
        idx = findfirst(x -> x ≥ t, tsteps)
        V = V_vals[idx]
        λ_dummy = zeros(5)
        state_dynamics!(du, u, p, t, λ_dummy, V)
    end
    prob_fwd = ODEProblem(forward_ode!, n0, tspan)
    sol_n = solve(prob_fwd, Tsit5(), saveat=tsteps)

    function backward_ode!(dl, l, p, t)
        idx = findfirst(x -> x ≥ t, tsteps)
        n = sol_n(tsteps[idx])
        V = V_vals[idx]
        adjoint_dynamics!(dl, l, n, V)
    end
    prob_bwd = ODEProblem(backward_ode!, λ0, (tspan[2], tspan[1]))
    sol_λ = solve(prob_bwd, Tsit5(), saveat=reverse(tsteps))

    for (i, t) in enumerate(tsteps)
        n_t = sol_n(t)
        λ_t = sol_λ(t)
        V_vals[i] = optimal_V(n_t, λ_t)
        rna_vals[i] = n_t[5]
    end
end

# Plotting
plot(tsteps, V_vals, lw=2, xlabel="Time (t)", ylabel="V(t)",
     title="Optimal Control V(t)", legend=false)

#plot(tsteps, rna_vals,
 #    label = "RNA (n₅)",
  #   xlabel = "Time", ylabel = "Total RNA Produced",
   #  title = "RNA Production Over Time",
    # linewidth = 2, color = :green, legend = :topright, grid = true)

#display(plot!)
#println("Maximum RNA produced: ", round(maximum(rna_vals), digits=4))

# Interpolation and Polynomial Fit
V_interp = LinearInterpolation(tsteps, V_vals, extrapolation_bc=Line())
println("V at t=2.5 is: ", V_interp(2.5))

open("V_function.txt", "w") do io
    for (t, v) in zip(tsteps, V_vals)
        println(io, "V($t) = $v")
    end
end

pfit = fit(tsteps, V_vals, 5)
println("Approximate V(t) ≈ ", pfit)

#using V_spline 
V_spline = CubicSplineInterpolation(tsteps, V_vals, extrapolation_bc=Line())
plot(tsteps, V_vals, label="True V(t)", lw=2)
plot!(tsteps, V_spline.(tsteps), label="Cubic Spline V(t)", lw=2, linestyle=:dash)



# Define V(t) as the fitted polynomial function
V_fit(t) = pfit(t)





# Cubic spline interpolation of V(t)
V_spline = CubicSplineInterpolation(tsteps, V_vals, extrapolation_bc=Line())
V_fit(t) = V_spline(t)

# ODE system with spline-based V(t)
function ode_func_fit!(dy, y, p, t)
    α, Cp, αCp, CppC, RNA = y
    V = V_fit(t)
    r1 = k1 * αCp^2
    r2 = k2 * CppC * α
    r3 = k3 * CppC
    r4 = k4 * α * Cp
    r5 = k5 * αCp
    r6 = keff(y, V) * CppC

    dy[1] = r1 - r2 - r4
    dy[2] = r3 + r5 - r4
    dy[3] = -2r1 + 2r2 + r3 + r4 - r5
    dy[4] = r1 - r2 - r3 - r6
    dy[5] = r6
end

# Initial concentrations
n0 = [0.0, 0.0, .5, 0.0, 0.0]

# Time span and save points
tspan = (0.0, 1000.0)
saveat = tsteps

# Solve the ODE
prob_fit = ODEProblem(ode_func_fit!, n0, tspan)
sol_fit = solve(prob_fit, Tsit5(), saveat=saveat)

plot(sol_fit.t, sol_fit[1,:], label="α", lw=2)
plot!(sol_fit.t, sol_fit[2,:], label="Cp", lw=2)
plot!(sol_fit.t, sol_fit[3,:], label="αCp", lw=2)
plot!(sol_fit.t, sol_fit[4,:], label="CppC", lw=2)
plot!(sol_fit.t, sol_fit[5,:], label="RNA", lw=2, xlabel="Time", ylabel="Concentration",
      title="Species Concentrations Over Time", legend=:right)




# --- Cubic spline interpolation of V(t) ---
V_spline = CubicSplineInterpolation(tsteps, V_vals, extrapolation_bc=Line())
V_fit(t) = V_spline(t)

# --- Initial concentrations ---
n0 = [0.0, 0.0, 0.5, 0.0, 0.0]
tspan = (0.0, 1000.0)
saveat = tsteps

# --- ODE system with spline-based V(t) ---
function ode_func_spline!(dy, y, p, t)
    α, Cp, αCp, CppC, RNA = y
    V = V_fit(t)
    r1 = k1 * αCp^2
    r2 = k2 * CppC * α
    r3 = k3 * CppC
    r4 = k4 * α * Cp
    r5 = k5 * αCp
    r6 = keff(y, V) * CppC

    dy[1] = r1 - r2 - r4
    dy[2] = r3 + r5 - r4
    dy[3] = -2r1 + 2r2 + r3 + r4 - r5
    dy[4] = r1 - r2 - r3 - r6
    dy[5] = r6
end

prob_spline = ODEProblem(ode_func_spline!, n0, tspan)
sol_spline = solve(prob_spline, Tsit5(), saveat=saveat)

# --- ODE system with constant V=10 ---
V_constantW(t) = .1

function ode_func_constant!(dy, y, p, t)
    α, Cp, αCp, CppC, RNA = y
    V = V_constantW(t)
    r1 = k1 * αCp^2
    r2 = k2 * CppC * α
    r3 = k3 * CppC
    r4 = k4 * α * Cp
    r5 = k5 * αCp
    r6 = keff(y, V) * CppC

    dy[1] = r1 - r2 - r4
    dy[2] = r3 + r5 - r4
    dy[3] = -2r1 + 2r2 + r3 + r4 - r5
    dy[4] = r1 - r2 - r3 - r6
    dy[5] = r6
end

prob_constant = ODEProblem(ode_func_constant!, n0, tspan)
sol_constant = solve(prob_constant, Tsit5(), saveat=saveat)

# --- Plot RNA concentration for both cases ---
plot(sol_spline.t, sol_spline[5, :], label="RNA with spline V(t)", lw=2)
plot!(sol_constant.t, sol_constant[5, :], label="RNA with constant V=10", lw=2, linestyle=:dash)
xlabel!("Time")
ylabel!("RNA Concentration")
title!("RNA Concentration: Spline V(t) vs Constant V")
