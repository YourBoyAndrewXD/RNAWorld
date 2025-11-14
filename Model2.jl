using Pkg
Pkg.activate("ControlOptimization")
using OrdinaryDiffEq, ForwardDiff, Optim, Plots

# Parameters
k1, k2, k3 = 1.8, 0.7 , 1.9
Ki = 2
Kd = 0.3
D = 2.5
Vmin, Vmax = 0.1, 100.0
ϵ = 1e-8  # small value to avoid division by zero

# Rates with stabilizers
rate1(A, C, V) = V * k1 * (A / (V + ϵ)) / (1 + (C / (V + ϵ)) / Ki)
rate2(B, V, D) = V * k2 * (B / (V + ϵ)) / (1 + (D / (V + ϵ)) / Kd)
rate3(C, V) = V * k3 * (C / (V + ϵ))

# System dynamics: n = [A, B, C]
function state_dynamics!(du, n, p, t, λ, V, D)
    A, B, C = n
    r1 = rate1(A, C, V)
    r2 = rate2(B, V, D)
    r3 = rate3(C, V)
    du[1] = -r1
    du[2] = r1 - r2
    du[3] = r2 - r3
end

# Hamiltonian
function H(n, λ, V, D)
    A, B, C = n
    L1, L2, L3 = λ
    r1 = rate1(A, C, V)
    r2 = rate2(B, V, D)
    r3 = rate3(C, V)
    return L1 * (-r1) + L2 * (r1 - r2) + L3 * (r2 - r3)
end

# Adjoint dynamics
function adjoint_dynamics!(dl, λ, n, V, D)
    # PMP: dl/dt = -∂H/∂n
    dl .= -ForwardDiff.gradient(n_ -> H(n_, λ, V, D), n)
end

# Pontryagin's Maximum Principle: maximize Hamiltonian w.r.t. V
function optimal_V(n, λ)
    obj(V) = -H(n, λ, V, D)  # minimize -H is same as maximize H
    res = optimize(obj, Vmin, Vmax, Brent())
    Vopt = Optim.minimizer(res)
    return clamp(Vopt, Vmin, Vmax)
end

# Forward-backward sweep
function forward_backward_sweep(n0, λT, tspan, tsteps; control_update=optimal_V)
    V_vals = fill((Vmin + Vmax) / 2, length(tsteps))
    sol_n, sol_λ = nothing, nothing
    for sweep in 1:30
        # Forward ODE
        function forward_ode!(du, u, p, t)
            idx = findfirst(x -> x ≥ t, tsteps)
            V = V_vals[clamp(idx, 1, length(V_vals))]
            λ_dummy = zeros(3)
            state_dynamics!(du, u, p, t, λ_dummy, V, D)
        end
        prob_fwd = ODEProblem(forward_ode!, n0, tspan)
        sol_n = solve(prob_fwd, Tsit5(), saveat=tsteps)
        # Backward ODE (adjoint)
        function backward_ode!(dl, l, p, t)
            idx = findfirst(x -> x ≥ t, reverse(tsteps))
            n = sol_n(tsteps[end - idx + 1])
            V = V_vals[end - idx + 1]
            adjoint_dynamics!(dl, l, n, V, D)
        end
        prob_bwd = ODEProblem(backward_ode!, λT, (tspan[2], tspan[1]))
        sol_λ = solve(prob_bwd, Tsit5(), saveat=reverse(tsteps))
        # Update control using optimal V(t)
        for (i, t) in enumerate(tsteps)
            n_t = sol_n(t)
            λ_t = sol_λ(tspan[2] - t)
            V_vals[i] = control_update(n_t, λ_t)
        end
    end
    return sol_n, sol_λ, V_vals
end

# Initial conditions
n0 = [200.0, 0.0, 0.0]  # [A, B, C]
λT = [0.0, 0.0, 1.0]   # final adjoint state: sensitivity to C

tspan = (0.0, 3.0)
tsteps = range(tspan[1], tspan[2], length=300)

# Run forward-backward sweep for optimal control
sol_n, sol_λ, V_vals = forward_backward_sweep(n0, λT, tspan, tsteps)

# Run for constant controls (all use tsteps of length 300)
function control_Vmax(n, λ)
    return Vmax
end
function control_Vmin(n, λ)
    return Vmin
end
function control_Vmean(n, λ)
    return (Vmin + Vmax) / 2
end

sol_n_Vmax, _, _ = forward_backward_sweep(n0, λT, tspan, tsteps; control_update=control_Vmax)
sol_n_Vmin, _, _ = forward_backward_sweep(n0, λT, tspan, tsteps; control_update=control_Vmin)
sol_n_Vmean, _, _ = forward_backward_sweep(n0, λT, tspan, tsteps; control_update=control_Vmean)

# Extract C concentrations
function safe_extract(sol, tsteps)
    # Use interpolation to get values at tsteps, regardless of sol.saveat
    return [sol(t)[3] for t in tsteps]
end
C_opt = safe_extract(sol_n, tsteps)
C_Vmax = safe_extract(sol_n_Vmax, tsteps)
C_Vmin = safe_extract(sol_n_Vmin, tsteps)
C_Vmean = safe_extract(sol_n_Vmean, tsteps)

# Plot all C curves together
pC = plot(tsteps, C_opt, label="Optimal V(t)", lw=2)
plot!(pC, tsteps, C_Vmax, label="Constant Vmax", lw=2, ls=:dash)
plot!(pC, tsteps, C_Vmin, label="Constant Vmin", lw=2, ls=:dot)
plot!(pC, tsteps, C_Vmean, label="Constant Vmean", lw=2, ls=:dashdot)
xlabel!("Time")
ylabel!("Concentration C")
title!("C(t) under different controls")

# Plot V(t) for optimal
pV = plot(tsteps, V_vals, label="Optimal V(t)", lw=2)
xlabel!("Time")
ylabel!("Volume V(t)")
title!("Optimal Control V(t)")

# Plot all together
plot(pC, pV, layout=(2,1), size=(800,600))

# Print max C for each control
println("Max C (Optimal): ", maximum(C_opt))
println("Max C (Vmax): ", maximum(C_Vmax))
println("Max C (Vmin): ", maximum(C_Vmin))
println("Max C (Vmean): ", maximum(C_Vmean))



# Detect times when Vmax or Vmin used in numerical optimal control

# If V_num, C_num, C_const are not defined, use V_vals, C_opt, C_Vmean as fallback
V_num = V_vals
C_num = C_opt
C_const = C_Vmean
times_Vmax = [tsteps[i] for i in 1:length(V_num) if isapprox(V_num[i], Vmax; atol=1e-4)]
times_Vmin = [tsteps[i] for i in 1:length(V_num) if isapprox(V_num[i], Vmin; atol=1e-4)]

println("Times when Vmax is used in numerical optimal control:")
println(times_Vmax)

println("Times when Vmin is used in numerical optimal control:")
println(times_Vmin)

# Find max C from all controls (numerical optimal and constant)

# Compute maxC_Vmax and maxC_Vmin for summary printout
maxC_num = maximum(C_num)
maxC_const = maximum(C_const)
maxC_Vmax = maximum(C_Vmax)
maxC_Vmin = maximum(C_Vmin)

println("Max concentration C (Numerical Optimal Control): ", maxC_num)
println("Max concentration C (Constant Control): ", maxC_const)

Pkg.add("RecipesBase")
using RecipesBase

function v_control_markers!(p, times_Vmax, times_Vmin)
    for t in times_Vmax
        vline!(p, [t], line=:dash, color=:red, label=false)
    end
    for t in times_Vmin
        vline!(p, [t], line=:dot, color=:blue, label=false)
    end
end

v_control_markers!(p2, times_Vmax, times_Vmin)

# Constant control at Vmax

# Already defined above, skip duplicate definitions

# Constant control at Vmin

# Already defined above, skip duplicate definitions

plot(tsteps, C_num, label="Numerical Optimal V(t)", lw=2)
plot!(tsteps, C_const, label="Constant V(t) = $(round((Vmin+Vmax)/2, digits=2))", lw=2, ls=:dot)
plot!(tsteps, C_Vmax, label="Constant V(t) = Vmax = $Vmax", lw=2, ls=:dashdot)
plot!(tsteps, C_Vmin, label="Constant V(t) = Vmin = $Vmin", lw=2, ls=:dashdot)
xlabel!("Time")
ylabel!("Concentration C")
title!("Comparison of C concentration under different controls")


println("\n======================")
println("Summary of Max C for Each Control Strategy")
println("======================")
println("Numerical Optimal V(t):      ", round(maxC_num, digits=4))
println("Constant V(t) = ", round((Vmin + Vmax)/2, digits=4), ":       ", round(maxC_const, digits=4))
println("Constant V(t) = Vmax = ", Vmax, ":    ", round(maxC_Vmax, digits=4))
println("Constant V(t) = Vmin = ", Vmin, ":    ", round(maxC_Vmin, digits=4))
println("======================\n")

# --- Parameter search for Model 2 ---
"""
Run the main simulation for a given parameter set and return max C for Vmax, Vmin, Vmean, and optimal.
"""
function run_main_simulation_for_params(params)
    k1, k2, k3, Ki, Kd, D = params
    Vmin = 0.1; Vmax = 100.0
    ϵ = 1e-8
    n0 = [200.0, 0.0, 0.0]
    λT = [0.0, 0.0, 1.0]
    tspan = (0.0, 3.0)
    tsteps = range(tspan[1], tspan[2], length=300)
    rate1(A, C, V) = V * k1 * (A / (V + ϵ)) / (1 + (C / (V + ϵ)) / Ki)
    rate2(B, V, D) = V * k2 * (B / (V + ϵ)) / (1 + (D / (V + ϵ)) / Kd)
    rate3(C, V) = V * k3 * (C / (V + ϵ))
    function state_dynamics!(du, n, p, t, λ, V, D)
        A, B, C = n
        r1 = rate1(A, C, V)
        r2 = rate2(B, V, D)
        r3 = rate3(C, V)
        du[1] = -r1
        du[2] = r1 - r2
        du[3] = r2 - r3
    end
    function H(n, λ, V, D)
        A, B, C = n
        L1, L2, L3 = λ
        r1 = rate1(A, C, V)
        r2 = rate2(B, V, D)
        r3 = rate3(C, V)
        return L1 * (-r1) + L2 * (r1 - r2) + L3 * (r2 - r3)
    end
    function adjoint_dynamics!(dl, λ, n, V, D)
        # PMP: dl/dt = -∂H/∂n
        dl .= -ForwardDiff.gradient(n_ -> H(n_, λ, V, D), n)
    end
    function optimal_V(n, λ)
        obj(V) = -H(n, λ, V, D)
        res = optimize(obj, Vmin, Vmax, Brent())
        Vopt = Optim.minimizer(res)
        return clamp(Vopt, Vmin, Vmax)
    end
    function control_Vmax(n, λ) Vmax end
    function control_Vmin(n, λ) Vmin end
    control_Vmean(n, λ) = (Vmin + Vmax) / 2
    function forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=optimal_V)
        V_vals = fill((Vmin + Vmax) / 2, length(tsteps))
        sol_n, sol_λ = nothing, nothing
        for sweep in 1:20
            function forward_ode!(du, u, p, t)
                idx = findfirst(x -> x ≥ t, tsteps)
                V = V_vals[clamp(idx, 1, length(V_vals))]
                λ_dummy = zeros(3)
                state_dynamics!(du, u, p, t, λ_dummy, V, D)
            end
            prob_fwd = ODEProblem(forward_ode!, n0, tspan)
            sol_n = solve(prob_fwd, Tsit5(), saveat=tsteps)
            function backward_ode!(dl, l, p, t)
                idx = findfirst(x -> x ≥ t, reverse(tsteps))
                n = sol_n(tsteps[end - idx + 1])
                V = V_vals[end - idx + 1]
                adjoint_dynamics!(dl, l, n, V, D)
            end
            prob_bwd = ODEProblem(backward_ode!, λT, (tspan[2], tspan[1]))
            sol_λ = solve(prob_bwd, Tsit5(), saveat=reverse(tsteps))
            for (i, t) in enumerate(tsteps)
                n_t = sol_n(t)
                λ_t = sol_λ(tspan[2] - t)
                V_vals[i] = control_update(n_t, λ_t)
            end
        end
        return sol_n, sol_λ, V_vals
    end
    sol_n_num, sol_λ_num, V_num = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=optimal_V)
    C_num = [sol_n_num(t)[3] for t in tsteps]
    maxC_num = maximum(C_num)
    sol_n_Vmax, _, _ = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=control_Vmax)
    sol_n_Vmin, _, _ = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=control_Vmin)
    sol_n_Vmean, _, _ = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=control_Vmean)
    C_Vmax = [sol_n_Vmax(t)[3] for t in tsteps]
    C_Vmin = [sol_n_Vmin(t)[3] for t in tsteps]
    C_Vmean = [sol_n_Vmean(t)[3] for t in tsteps]
    maxC_Vmax = maximum(C_Vmax)
    maxC_Vmin = maximum(C_Vmin)
    maxC_Vmean = maximum(C_Vmean)
    return maxC_num, maxC_Vmax, maxC_Vmin, maxC_Vmean
end
function simulate_and_compare_model2(params)
    k1, k2, k3, Ki, Kd, D = params
    Vmin = 0.1; Vmax = 100.0
    n0 = [200.0, 0.0, 0.0]  # match main script
    λT = [0.0, 0.0, 1.0]
    tspan = (0.0, 3.0)
    tsteps = range(tspan[1], tspan[2], length=300)
    # Use the same rate definitions as main script (with V factor)
    function rate1(A, C, V) V * k1 * (A / (V + 1e-8)) / (1 + (C / (V + 1e-8)) / Ki) end
    function rate2(B, V, D) V * k2 * (B / (V + 1e-8)) / (1 + (D / (V + 1e-8)) / Kd) end
    function rate3(C, V) V * k3 * (C / (V + 1e-8)) end
    function state_dynamics!(du, n, p, t, λ, V, D)
        A, B, C = n
        r1 = rate1(A, C, V)
        r2 = rate2(B, V, D)
        r3 = rate3(C, V)
        du[1] = -r1
        du[2] = r1 - r2
        du[3] = r2 - r3
    end
    function H(n, λ, V, D)
        A, B, C = n
        L1, L2, L3 = λ
        r1 = rate1(A, C, V)
        r2 = rate2(B, V, D)
        r3 = rate3(C, V)
        return L1 * (-r1) + L2 * (r1 - r2) + L3 * (r2 - r3)
    end
    function adjoint_dynamics!(dl, λ, n, V, D)
        # PMP: dl/dt = -∂H/∂n
        dl .= -ForwardDiff.gradient(n_ -> H(n_, λ, V, D), n)
    end
    function optimal_V(n, λ)
        obj(V) = -H(n, λ, V, D)
        res = optimize(obj, Vmin, Vmax, Brent())
        Vopt = Optim.minimizer(res)
        return clamp(Vopt, Vmin, Vmax)
    end
    function control_Vmax(n, λ) Vmax end
    function control_Vmin(n, λ) Vmin end
    control_Vmean(n, λ) = (Vmin + Vmax) / 2
    function forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=optimal_V)
        V_vals = fill((Vmin + Vmax) / 2, length(tsteps))
        sol_n, sol_λ = nothing, nothing
        for sweep in 1:20
            function forward_ode!(du, u, p, t)
                idx = findfirst(x -> x ≥ t, tsteps)
                V = V_vals[clamp(idx, 1, length(V_vals))]
                λ_dummy = zeros(3)
                state_dynamics!(du, u, p, t, λ_dummy, V, D)
            end
            prob_fwd = ODEProblem(forward_ode!, n0, tspan)
            sol_n = solve(prob_fwd, Tsit5(), saveat=tsteps)
            function backward_ode!(dl, l, p, t)
                idx = findfirst(x -> x ≥ t, reverse(tsteps))
                n = sol_n(tsteps[end - idx + 1])
                V = V_vals[end - idx + 1]
                adjoint_dynamics!(dl, l, n, V, D)
            end
            prob_bwd = ODEProblem(backward_ode!, λT, (tspan[2], tspan[1]))
            sol_λ = solve(prob_bwd, Tsit5(), saveat=reverse(tsteps))
            for (i, t) in enumerate(tsteps)
                n_t = sol_n(t)
                λ_t = sol_λ(tspan[2] - t)
                V_vals[i] = control_update(n_t, λ_t)
            end
        end
        return sol_n, sol_λ, V_vals
    end
    sol_n_num, sol_λ_num, V_num = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=optimal_V)
    C_num = [sol_n_num(t)[3] for t in tsteps]
    maxC_num = maximum(C_num)
    sol_n_Vmax, _, _ = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=control_Vmax)
    sol_n_Vmin, _, _ = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=control_Vmin)
    sol_n_Vmean, _, _ = forward_backward_sweep_local(n0, λT, tspan, tsteps; control_update=control_Vmean)
    C_Vmax = [sol_n_Vmax(t)[3] for t in tsteps]
    C_Vmin = [sol_n_Vmin(t)[3] for t in tsteps]
    C_Vmean = [sol_n_Vmean(t)[3] for t in tsteps]
    maxC_Vmax = maximum(C_Vmax)
    maxC_Vmin = maximum(C_Vmin)
    maxC_Vmean = maximum(C_Vmean)
    maxC_const_all = maximum([maxC_Vmax, maxC_Vmin, maxC_Vmean])
    # Ensure all returned vectors are of length 300
    if length(C_num) != 300 || length(C_Vmax) != 300 || length(C_Vmin) != 300 || length(C_Vmean) != 300
        error("ODE solution did not return 300 points. C_num: $(length(C_num)), C_Vmax: $(length(C_Vmax)), C_Vmin: $(length(C_Vmin)), C_Vmean: $(length(C_Vmean))")
    end
    return maxC_num - maxC_const_all, maxC_num, maxC_Vmax, maxC_Vmin, maxC_Vmean, params
end

println("\n--- Model 2 Parameter Search (diff > 0.5) ---")
best2 = (-Inf, 0, 0, 0, 0, zeros(6))
found = false
for trial in 1:500
    k1 = rand(0.1:0.1:2.0)
    k2 = rand(0.1:0.1:2.0)
    k3 = rand(0.1:0.1:2.0)
    Ki = rand(0.1:0.1:2.0)
    Kd = rand(0.1:0.1:2.0)
    D = rand(0.1:0.1:10.0)
    diff, maxC_num, maxC_Vmax, maxC_Vmin, maxC_Vmean, params = simulate_and_compare_model2([k1, k2, k3, Ki, Kd, D])
    if diff > 0.5 && isfinite(diff) && isfinite(maxC_num) && isfinite(maxC_Vmax) && isfinite(maxC_Vmin) && isfinite(maxC_Vmean)
        println("Params: k1=", params[1], ", k2=", params[2], ", k3=", params[3], ", Ki=", params[4], ", Kd=", params[5], ", D=", params[6])
        println("Numerical optimal max C: ", round(maxC_num, digits=4))
        println("Constant Vmax max C: ", round(maxC_Vmax, digits=4))
        println("Constant Vmin max C: ", round(maxC_Vmin, digits=4))
        println("Constant Vmean max C: ", round(maxC_Vmean, digits=4))
        println("Difference (numerical - best constant): ", round(diff, digits=4))
        if diff > best2[1]
            best2 = (diff, maxC_num, maxC_Vmax, maxC_Vmin, maxC_Vmean, params)
        end
        found = true
    end
end
if found
    println("\nBest found (diff > 0.5):")
    println("Difference: ", best2[1])
    println("Numerical optimal max C (parameter search): ", best2[2])
    println("Vmax (parameter search): ", best2[3], ", Vmin: ", best2[4], ", Vmean: ", best2[5])
    println("Params: k1=", best2[6][1], ", k2=", best2[6][2], ", k3=", best2[6][3], ", Ki=", best2[6][4], ", Kd=", best2[6][5], ", D=", best2[6][6])
    # Now run main simulation for these params and compare
    maxC_num_main, maxC_Vmax_main, maxC_Vmin_main, maxC_Vmean_main = run_main_simulation_for_params(best2[6])
    println("\n--- Main simulation with best parameter set ---")
    println("Numerical optimal max C (main): ", round(maxC_num_main, digits=4))
    println("Vmax (main): ", round(maxC_Vmax_main, digits=4), ", Vmin: ", round(maxC_Vmin_main, digits=4), ", Vmean: ", round(maxC_Vmean_main, digits=4))
    println("\nCompare: If these do not match, check for hidden differences in rates, ICs, or time grid.")
else
    println("\nNo parameter set found where numerical optimal C exceeds all constant controls by more than 0.5.")
end

# Plot A, B, and C across time for optimal control
A_opt = [sol_n(t)[1] for t in tsteps]
B_opt = [sol_n(t)[2] for t in tsteps]
C_opt = [sol_n(t)[3] for t in tsteps]
pABC = plot(tsteps, A_opt, label="A (Optimal)", lw=2)
plot!(pABC, tsteps, B_opt, label="B (Optimal)", lw=2)
plot!(pABC, tsteps, C_opt, label="C (Optimal)", lw=2)
xlabel!(pABC, "Time")
ylabel!(pABC, "Concentration")
title!(pABC, "A, B, C vs Time under Optimal Control")

# Plot numerical optimal V(t) below
pVnum = plot(tsteps, V_vals, label="Numerical Optimal V(t)", lw=2, color=:black)
xlabel!(pVnum, "Time")
ylabel!(pVnum, "V(t)")
title!(pVnum, "Numerical Optimal Control V(t)")

plot(pABC, pVnum, layout=(2,1), size=(800,600))
