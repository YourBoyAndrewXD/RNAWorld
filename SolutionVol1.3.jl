Pkg.activate("ControlOptimization")
using OrdinaryDiffEqFIRK
using ForwardDiff
using Optim
using OrdinaryDiffEq
using Plots
using Interpolations
using Polynomials
using Catalyst
# Constants
k1 = 1
k2 = 1
k3 = 1
k4 = 1 
k5 = 1
c6 = 1.06
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
    r1 = k1 * (αCp/V)^2
    r2 = k2 * CppC/V * α/V
    r3 = k3 * CppC/V
    r4 = k4 * α/V * Cp/V
    r5 = k5 * αCp/V
    r6 = keff(n, V) * CppC/V

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
    # ∂H/∂n = V * ( ... )

    dl[1] = - (L1*(-k2*CppC - k4*Cp) + L2*(-k4*Cp) + L3*(2k2*CppC + k4*Cp) - L4*(k2*CppC + k4*Cp))
    dl[2] = - (L1*(-k4*α) + L2*(-k4*α) + L3*(k4*α) - L4*(k4*α)) + L5 * CppC * dk_dV
    dl[3] = - (L1*(-2k1*αCp) + L2*(k5) + L3*(4k1*αCp - k5) - L4*(2k1*αCp)) + L5 * CppC * dk_dV
    dl[4] = - (L1*(-k2*α) + L2*(-k3) + L3*(2k2*α + k3) - L4*(k2*α + k3 + c6*(1 + (Cp + αCp)/(Kd * V)))) - L5 * c6*(1 + (Cp + αCp)/(Kd * V))
    dl[5] = 0.0
end

# Hamiltonian
function H(n, λ, V)
    α, Cp, αCp, CppC, _ = n
    L1, L2, L3, L4, L5 = λ
    term = (
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
n0 = [1e-4, 1e-4, 1 , 0.0, 0.0]   # state
λ0 = [0.0, 0.0, 0.0, 0.0, 1.0]   # costate

# Time setup
tspan = (0.0, 100.0)
tsteps = range(tspan[1], tspan[2], length=100)
V_vals = fill(1.0, length(tsteps))

rna_vals = zeros(length(tsteps))

# Forward-backward sweep
# Declare outside to persist after loop
sol_n = nothing
sol_λ = nothing

V0(t) = 1.0

function ForwardBackwardSweep(init_controller::Function; tolerance = 1e-5, maxiters = 1000)
    # x_i, λ_i, u_i 
    # x_{i+1}, λ_{i+1}, u_{i+1}
    uᵢ = init_controller
    function init_x0!(du, u, p, t)
        state_dynamics!(du, u, p, t, λ_dummy, V0(t))
    end
    prob_fwd = ODEProblem(init_x0!, n0, tspan)
    xᵢ = solve(prob_fwd, Tsit5(), saveat=tsteps)

    function init_λ0!(du, u, p, t)
        adjoint_dynamics!(du, u, x0(t), V0(t))
    end
    prob_bwd_local = ODEProblem(init_λ0!, λ0, (tspan[2], tspan[1]))
    λᵢ = solve(prob_bwd_local, Tsit5(), saveat=reverse(tsteps))

    λ_dummy = zeros(5)

    for sweep in 1:maxiters
        println("Sweep $sweep")
        function forward_ode!(du, u, p, t)
            state_dynamics!(du, u, p, t, λ_dummy, uᵢ(t))
        end

        prob_fwd = ODEProblem(forward_ode!, n0, tspan)
        sol_n = solve(prob_fwd, Tsit5(), saveat=tsteps)

        function backward_ode!(dl, l, p, t)
            idx = findfirst(x -> x ≥ t, tsteps)
            n = sol_n(tsteps[idx])
            V = V_vals[idx]
            adjoint_dynamics!(dl, l, n, V)
        end
        prob_bwd_local = ODEProblem(backward_ode!, λ0, (tspan[2], tspan[1]))
        sol_λ = solve(prob_bwd_local, Tsit5(), saveat=reverse(tsteps))
    end
end

for sweep in 1:100


    for (i, t) in enumerate(tsteps)
        n_t = sol_n(t)
        λ_t = sol_λ(tspan[2] - t)  # reverse time for correct 
        V_vals[i] = optimal_V(n_t, λ_t)
        rna_vals[i] = n_t[5]
    end
end

# Plot costate functions λ₁ through λ₅ over time
λ_matrix = [sol_λ(t) for t in tsteps] |> hcat  # shape: (5, length(tsteps))

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





# Cubic spline interpolation of V(t)
V_fit(t) = V_spline(t)

# Time span and save points
tspan = (0.0, 100.0)
saveat = tsteps

# --- ODE system with spline-based V(t) ---
function ode_func_spline!(dy, y, p, t)
    α, Cp, αCp, CppC, RNA = y
    V = V_spline(t)
    r1 = k1 * (αCp/V)^2
    r2 = k2 * CppC/V * α/V
    r3 = k3 * CppC/V
    r4 = k4 * α/V * Cp/V
    r5 = k5 * αCp/V
    r6 = keff(y, V) * CppC/V
    
    dy[1] = r1 - r2 - r4
    dy[2] = r3 + r5 - r4
    dy[3] = -2r1 + 2r2 + r3 + r4 - r5
    dy[4] = r1 - r2 - r3 - r6
    dy[5] = r6
end

prob_spline = ODEProblem(ode_func_spline!, n0, tspan)
sol_spline = solve(prob_spline, RadauIIA5())



#Test
# --- Safe wrapper around V_spline to clip values ---
function safe_V(t)
    V_raw = V_spline(t)
    return clamp(V_raw, 0.1, 10.0)  # V ∈ [0.1, 10]
end

# --- ODE system using stable spline-based V(t) ---
function ode_func_spline!(dy, y, p, t)
    α, Cp, αCp, CppC, RNA = y
    V = safe_V(t)
    r1 = k1 * (αCp / V)^2
    r2 = k2 * (CppC / V) * (α / V)
    r3 = k3 * (CppC / V)
    r4 = k4 * (α / V) * (Cp / V)
    r5 = k5 * (αCp / V)
    r6 = keff(y, V) * (CppC / V)

    dy[1] = r1 - r2 - r4
    dy[2] = r3 + r5 - r4
    dy[3] = -2r1 + 2r2 + r3 + r4 - r5
    dy[4] = r1 - r2 - r3 - r6
    dy[5] = r6
end

# --- Stable solve with tolerances ---
prob_spline = ODEProblem(ode_func_spline!, n0, tspan)
sol_spline = solve(prob_spline, Rodas5(), saveat=tsteps, reltol=1e-8, abstol=1e-10)

V_raw_vals = [V_spline(t) for t in tsteps]

plot(tsteps, V_raw_vals,
     xlabel = "Time (t)",
     ylabel = "V_raw(t)",
     title = "Raw Spline-Interpolated V(t)",
     label = "V_spline(t)",
     lw = 2,
     color = :blue,
     legend = :topright,
     grid = true)
#end test
# --- ODE system with constant V=10 ---
V_constant(t) = 10

function ode_func_constant!(dy, y, p, t)
    α, Cp, αCp, CppC, RNA = y
    V = V_constant(t)
    r1 = k1 * (αCp/V)^2
    r2 = k2 * CppC/V * α/V
    r3 = k3 * CppC/V
    r4 = k4 * α/V * Cp/V
    r5 = k5 * αCp/V
    r6 = keff(y, V) * CppC/V

    dy[1] = r1 - r2 - r4
    dy[2] = r3 + r5 - r4
    dy[3] = -2r1 + 2r2 + r3 + r4 - r5
    dy[4] = r1 - r2 - r3 - r6
    dy[5] = r6
end

# --- Solve constant V system ---
prob_constant = ODEProblem(ode_func_constant!, n0, tspan)
sol_constant = solve(prob_constant, Tsit5(), saveat=saveat)

# --- Extract RNA (n₅) values ---
rna_fit = [sol_spline(t)[5] for t in tsteps]
rna_constant = [sol_constant(t)[5] for t in tsteps]

# --- Plot comparison ---
plot(tsteps, rna_fit, label="RNA (V_fit(t))", lw=2, xlabel="Time", ylabel="RNA Amount",
     title="RNA Production: Fitted V(t) vs Constant V=0.1", legend=:bottomright)
plot!(tsteps, rna_constant, label="RNA (V = 0.1)", lw=2, linestyle=:dash)

# Optional: Print max RNA for each
println("Max RNA (V_fit): ", round(maximum(rna_fit), digits=4))
println("Max RNA (V = 0.1): ", round(maximum(rna_constant), digits=4))

# Extract concentrations from the spline solution
α_vals     = [sol_spline(t)[1] for t in tsteps]
Cp_vals    = [sol_spline(t)[2] for t in tsteps]
αCp_vals   = [sol_spline(t)[3] for t in tsteps]
RNA_vals   = [sol_spline(t)[5] for t in tsteps]

# Plot concentrations over time
plot(tsteps, Cp_vals, label="Cp", lw=2)
plot!(tsteps, αCp_vals, label="αCp", lw=2)
plot!(tsteps, α_vals, label="α", lw=2)
plot!(tsteps, RNA_vals, label="RNA", lw=2, xlabel="Time", ylabel="Concentration",
      title="System Species Over Time with V_fit", legend=:topright)



# Collect costate solutions λ₁ through λ₅ over time tsteps
λ_matrix = hcat([sol_λ(t) for t in tsteps]...)  # 5 rows × length(tsteps) columns



# Plot costates λ₁ to λ₅ over time
plot(
    tsteps, λ_matrix[1, :], label="λ₁", lw=2,
    xlabel="Time", ylabel="Costate Variables",
    title="Costate Variables λ₁ to λ₅ over Time",
    legend=:left, grid=true
)
plot!(tsteps, λ_matrix[2, :], label="λ₂", lw=2)
plot!(tsteps, λ_matrix[3, :], label="λ₃", lw=2)
plot!(tsteps, λ_matrix[4, :], label="λ₄", lw=2)
plot!(tsteps, λ_matrix[5, :], label="λ₅", lw=2)
