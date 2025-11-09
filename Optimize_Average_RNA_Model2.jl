<<<<<<< HEAD
# ==========================================
# File: Optimize_Average_RNA_Model.jl
# Goal: Optimize control V(t) to maximize RNA Average (fast, stable version)
# ==========================================
using Pkg
Pkg.activate(".")
using DifferentialEquations, Optim, Trapz, Interpolations, Statistics, DataInterpolations
using Plots

# -------------------------
# Constants and Parameters
# -------------------------
const eps_div = 1e-9
Vmin, Vmax = 0.1, 100.0
tspan = (0.0, 100.0)
n0 = [50.0, 500.0, 50.0, 5.0, 0.0]          # [α, Cp, αCp, CppC, RNA]
params = (0.251188643150958, 1.0, 0.5011872336272722, 0.15848931924611134, 0.025118864315095794, 2.51188643150958, 0.7943282347242815, 0.025118864315095794)  # (k1,k2,k3,k4,k5,c6,Kd,D1)

# -------------------------
# Reaction Model
# -------------------------
function rates(n::AbstractVector{T}, V::T, p) where T
    α, Cp, αCp, CppC, RNA = n
    k1,k2,k3,k4,k5,c6,Kd,D1 = p
    Vsafe = max(V, T(eps_div))

    r1 = k1 * (αCp^2) / (Vsafe^2)
    r2 = k2 * (α * CppC) / (Vsafe^2)
    r3 = k3 * CppC / Vsafe
    r4 = k4 * (α * Cp) / (Vsafe^2)
    r5 = k5 * αCp / Vsafe
    r6 = c6 * (1 + (Cp + αCp) / (Kd * Vsafe)) * (CppC / Vsafe)
    r7 = D1 * RNA  # RNA degradation
    return r1,r2,r3,r4,r5,r6,r7
end

function model!(du, n, p, t, V)
    r1,r2,r3,r4,r5,r6,r7 = rates(n, V, p)
    du[1] = r1 + r5 - r2 - r4           # α̇
    du[2] = r3 + r5 - r4 + r7           # Cṗ
    du[3] = 2r2 - 2r1 + r3 + r4 - r5 + r6   # αCṗ
    du[4] = r1 - r2 - r3 - r6           # CppĊ
    du[5] = r6 - r7                     # RNȦ
end
# -------------------------



using DifferentialEquations, Optim, Trapz, Plots

# -------------------------
# Sinusoidal control
# -------------------------
function sinusoidal_V(t, A, ω, ϕ)
    Vmean = (Vmax + Vmin) / 2
    raw = Vmean + A * sin(ω * t + ϕ)
    # Keep it within [Vmin, Vmax]
    return clamp(raw, Vmin, Vmax)
end

# Simulation for given sinusoidal parameters
function simulate_sinusoidal(params, n0, tspan, A, ω, ϕ)
    f!(du, u, p, t) = model!(du, u, p, t, sinusoidal_V(t, A, ω, ϕ))
    prob = ODEProblem(f!, n0, tspan, params)
    sol = solve(prob, Tsit5(); abstol=1e-6, reltol=1e-6, saveat=1.0)
    return sol
end

# -------------------------
# Optimization for A, ω, ϕ
# -------------------------
function find_optimal_sinusoidal(params)
    println("🔹 Optimizing sinusoidal control V(t) = mean + A*sin(ωt + ϕ)...")

    obj(x) = begin
        A, ω, ϕ = x
        sol = simulate_sinusoidal(params, n0, tspan, A, ω, ϕ)
        if sol.retcode != :Success
            return 1e6
        end
        RNA_vals = getindex.(sol.u, 5)
        avg_RNA = trapz(sol.t, RNA_vals) / (sol.t[end] - sol.t[1])
        return -avg_RNA
    end

    # Reasonable parameter bounds
    lower = [0.0, 0.01, 0.0]     # A, ω, ϕ lower bounds
    upper = [(Vmax - Vmin)/2, 1.0, 2π]  # upper bounds
    x0 = [10.0, 0.1, 0.0]

    result = optimize(obj, lower, upper, x0,
                      Fminbox(NelderMead()),
                      Optim.Options(iterations=800, show_trace=true))

    println("✅ Optimization complete.")
    println("Best parameters: A = $(result.minimizer[1]), ω = $(result.minimizer[2]), ϕ = $(result.minimizer[3])")
    println("Best mean RNA = ", -result.minimum)

    return result.minimizer
end

# Run the optimization
opt_params = find_optimal_sinusoidal(params)
A_opt, ω_opt, ϕ_opt = opt_params

# Generate V(t) and RNA plots
tvals = range(tspan[1], tspan[2], length=300)
Vsig = [sinusoidal_V(t, A_opt, ω_opt, ϕ_opt) for t in tvals]
plot(tvals, Vsig, lw=2, color=:blue, xlabel="Time", ylabel="V(t)", 
     title="Optimized Sinusoidal Control V(t)", label="V(t)")
hline!([Vmin], linestyle=:dash, color=:red, label="Vmin")
hline!([Vmax], linestyle=:dash, color=:green, label="Vmax")

# Simulate with optimized parameters
sol_sin = simulate_sinusoidal(params, n0, tspan, A_opt, ω_opt, ϕ_opt)
RNA_sin = getindex.(sol_sin.u, 5)
plot(sol_sin.t, RNA_sin, lw=2, xlabel="Time", ylabel="[RNA]", 
     title="RNA Production (Optimized Sinusoidal Control)", color=:blue)



using DifferentialEquations, Optim, Trapz, Plots, Random

# -------------------------
# Constants for Parameter Search
# -------------------------
const eps_div = 1e-9
Npts = 100

# -------------------------
# Sinusoidal Control
# -------------------------
function sinusoidal_V(t, A, ω, ϕ)
    Vmean = (Vmax + Vmin) / 2
    raw = Vmean + A * sin(ω * t + ϕ)
    return clamp(raw, Vmin, Vmax)
end

function simulate_sinusoidal(params, n0, tspan, A, ω, ϕ)
    f!(du, u, p, t) = model!(du, u, p, t, sinusoidal_V(t, A, ω, ϕ))
    prob = ODEProblem(f!, n0, tspan, params)
    sol = solve(prob, Tsit5(); abstol=1e-6, reltol=1e-6, saveat=1.0)
    return sol
end

# -------------------------
# Optimization Function
# -------------------------
function find_optimal_sinusoidal(params)
    obj(x) = begin
        A, ω, ϕ = x
        sol = simulate_sinusoidal(params, n0, tspan, A, ω, ϕ)
        if sol.retcode != :Success
            return 1e6
        end
        RNA_vals = getindex.(sol.u, 5)
        avg_RNA = trapz(sol.t, RNA_vals) / (sol.t[end] - sol.t[1])
        return -avg_RNA
    end

    lower = [0.0, 0.01, 0.0]
    upper = [(Vmax - Vmin)/2, 1.0, 2π]
    x0 = [10.0, 0.1, 0.0]

    result = optimize(obj, lower, upper, x0,
                      Fminbox(NelderMead()),
                      Optim.Options(iterations=800, show_trace=false))

    return result.minimizer, -result.minimum
end

# -------------------------
# Parameter Search (Exponential Scaling)
# -------------------------
function random_parameters()
    return (
        10.0^(rand(-1.0:0.1:1.0)),  # k1
        10.0^(rand(-2.0:0.1:0.0)),  # k2
        10.0^(rand(-2.0:0.1:0.0)),  # k3
        10.0^(rand(-1.0:0.1:1.0)),  # k4
        10.0^(rand(-2.0:0.1:0.0)),  # k5
        10.0^(rand(-1.0:0.1:0.7)),  # c6
        10.0^(rand(-1.0:0.1:1.0)),  # Kd
        10.0^(rand(-3.0:0.1:-0.3))  # D1
    )
end

# --------------------------
# Simulate with constant volume
# --------------------------
# --- Constant-volume simulation function ---
function simulate_constant(params, n0, tspan, Vconst)
    f!(du, u, p, t) = model!(du, u, p, t, Vconst)
    prob = ODEProblem(f!, n0, tspan, params)
    sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=1.0)
    return sol
end

# --- Run for Vmin and Vmax ---
sol_min = simulate_constant(params, n0, tspan, Vmin)
sol_max = simulate_constant(params, n0, tspan, Vmax)

# --- Extract RNA (5th variable) ---
RNA_min = getindex.(sol_min.u, 5)
RNA_max = getindex.(sol_max.u, 5)

# --- Compute average RNA ---
avg_min = trapz(sol_min.t, RNA_min) / (sol_min.t[end] - sol_min.t[1])
avg_max = trapz(sol_max.t, RNA_max) / (sol_max.t[end] - sol_max.t[1])


# -------------------------
# Main Monte Carlo Loop 
# -------------------------
n_trials = 11025
successes = []

# -------------------------
# Main Monte Carlo Loop
# -------------------------
for i in 1:n_trials
    println("\n🧪 Trial $i / $n_trials")
    params_rand = random_parameters()

    try
        # Optimize sinusoidal V(t)
        Vparams, avg_RNA_sig = find_optimal_sinusoidal(params_rand)

        # Simulate constant Vmin and Vmax
        sol_min = simulate_constant(params_rand, n0, tspan, Vmin)
        sol_max = simulate_constant(params_rand, n0, tspan, Vmax)

        RNA_min = getindex.(sol_min.u, 5)
        RNA_max = getindex.(sol_max.u, 5)

        avg_RNA_min = trapz(sol_min.t, RNA_min) / (sol_min.t[end] - sol_min.t[1])
        avg_RNA_max = trapz(sol_max.t, RNA_max) / (sol_max.t[end] - sol_max.t[1])

        if avg_RNA_sig > avg_RNA_min + 10 && avg_RNA_sig > avg_RNA_max + 10
            push!(successes, (
                params = params_rand,
                sinusoidal_avg = avg_RNA_sig,
                Vmin_avg = avg_RNA_min,
                Vmax_avg = avg_RNA_max,
                Vopt_params = Vparams
            ))
            println("✅ Success! RNA gain ≥10 over both Vmin/Vmax.")
        else
            println("⚪ No significant gain.")
        end

    catch e
        println("❌ Skipped trial due to: ", e)
        continue
    end
end

println("\n--- SUMMARY ---")
println("Found $(length(successes)) parameter sets with ≥10 RNA advantage.")
for s in successes
    println(s)
end


# -------------------------
# Summary
# -------------------------
println("\n--- SUMMARY ---")
println("Found $(length(successes)) kinetic parameter sets with ≥10 RNA advantage.")
for s in successes
    println(s)
end


println("\n--- SUMMARY ---")
println("Found $(length(successes)) parameter sets with ≥10 RNA advantage.")
for s in successes
    println(s)
end


# --- Simulate optimized sinusoidal control ---
sol_sig = simulate_sinusoidal(params, n0, tspan, A_opt, ω_opt, ϕ_opt)
RNA_sig = getindex.(sol_sig.u, 5)
avg_sig = trapz(sol_sig.t, RNA_sig) / (sol_sig.t[end] - sol_sig.t[1])

# -------------------------
# Simulations for Comparison
# -------------------------
sol_sig = simulate_sinusoidal(params, n0, tspan, A_opt, ω_opt, ϕ_opt)
sol_min = simulate_constant(params, n0, tspan, Vmin)
sol_max = simulate_constant(params, n0, tspan, Vmax)

RNA_sig = getindex.(sol_sig.u, 5)
RNA_min = getindex.(sol_min.u, 5)
RNA_max = getindex.(sol_max.u, 5)

avg_sig = trapz(sol_sig.t, RNA_sig) / (sol_sig.t[end] - sol_sig.t[1])
avg_min = trapz(sol_min.t, RNA_min) / (sol_min.t[end] - sol_min.t[1])
avg_max = trapz(sol_max.t, RNA_max) / (sol_max.t[end] - sol_max.t[1])

# -------------------------
# Plot RNA Comparison
# -------------------------
plot(sol_sig.t, RNA_sig, lw=3, color=:blue, label="Optimized Sinusoidal")
plot!(sol_min.t, RNA_min, lw=2, color=:red, linestyle=:dash, label="Constant Vmin")
plot!(sol_max.t, RNA_max, lw=2, color=:green, linestyle=:dot, label="Constant Vmax")
xlabel!("Time")
ylabel!("[RNA]")
title!("RNA Production: Vsig vs Vmin vs Vmax")

annotate!(maximum(sol_sig.t)*0.6, maximum(RNA_sig)*0.8,
    text("⟨RNA⟩ Vsig = $(round(avg_sig, digits=2))", :blue, 10))
annotate!(maximum(sol_sig.t)*0.6, maximum(RNA_sig)*0.7,
    text("⟨RNA⟩ Vmin = $(round(avg_min, digits=2))", :red, 10))
annotate!(maximum(sol_sig.t)*0.6, maximum(RNA_sig)*0.6,
    text("⟨RNA⟩ Vmax = $(round(avg_max, digits=2))", :green, 10))

# --- Plot comparison ---
plot(
    sol_sig.t, RNA_sig, lw=3, color=:blue, label="Vsig (Sinusoidal)",
    xlabel="Time", ylabel="[RNA]",
    title="RNA Production: Vsig vs Vmin vs Vmax"
)
plot!(sol_min.t, RNA_min, lw=2, color=:red, linestyle=:dash, label="Vmin (Constant)")
plot!(sol_max.t, RNA_max, lw=2, color=:green, linestyle=:dot, label="Vmax (Constant)")

annotate!(maximum(sol_sig.t)*0.6, maximum(RNA_sig)*0.8,
    text("⟨RNA⟩: Vsig = $(round(avg_sig, digits=2))", :blue, 10))
annotate!(maximum(sol_sig.t)*0.6, maximum(RNA_sig)*0.7,
    text("⟨RNA⟩: Vmin = $(round(avg_min, digits=2))", :red, 10))
annotate!(maximum(sol_sig.t)*0.6, maximum(RNA_sig)*0.6,
    text("⟨RNA⟩: Vmax = $(round(avg_max, digits=2))", :green, 10))




    # --- Simulate with optimized Vsig control ---
sol_sig = simulate_sinusoidal(params, n0, tspan, A_opt, ω_opt, ϕ_opt)

# Extract each species (assuming order: α, Cp, αCp, CppC, RNA)
α_vals     = getindex.(sol_sig.u, 1)
Cp_vals    = getindex.(sol_sig.u, 2)
αCp_vals   = getindex.(sol_sig.u, 3)
CppC_vals  = getindex.(sol_sig.u, 4)
RNA_vals   = getindex.(sol_sig.u, 5)

# --- Plot all on one graph ---
plot(sol_sig.t, α_vals,    lw=2, label="α",     color=:red)
plot!(sol_sig.t, Cp_vals,  lw=2, label="Cp",    color=:orange)
plot!(sol_sig.t, αCp_vals, lw=2, label="αCp",   color=:green)
plot!(sol_sig.t, CppC_vals,lw=2, label="CppC",  color=:blue)
plot!(sol_sig.t, RNA_vals, lw=3, label="RNA",   color=:purple)

xlabel!("Time")
ylabel!("Concentration")
title!("Species Dynamics under Optimized Sinusoidal V(t)")

=======
# ==========================================
# File: Optimize_Average_RNA_Model.jl
# Goal: Optimize control V(t) to maximize RNA Average (fast, stable version)
# ==========================================
using Pkg
Pkg.activate(".")
using DifferentialEquations, Optim, Trapz, Interpolations, Statistics, DataInterpolations
using Plots

# -------------------------
# Constants and Parameters
# -------------------------
const eps_div = 1e-9
Vmin, Vmax = 0.1, 100.0
tspan = (0.0, 100.0)
n0 = [50.0, 500.0, 50.0, 5.0, 0.0]          # [α, Cp, αCp, CppC, RNA]
# include D1 as last parameter
params = (9.02, 0.1146, 0.2099, 9.21, 0.00206, 1.2, 3.0, 1)  # (k1,k2,k3,k4,k5,c6,Kd,D1)

# -------------------------
# Reaction Model
# -------------------------
function rates(n::AbstractVector{T}, V::T, p) where T
    α, Cp, αCp, CppC, RNA = n
    k1,k2,k3,k4,k5,c6,Kd,D1 = p
    Vsafe = max(V, T(eps_div))

    r1 = k1 * (αCp^2) / (Vsafe^2)
    r2 = k2 * (α * CppC) / (Vsafe^2)
    r3 = k3 * CppC / Vsafe
    r4 = k4 * (α * Cp) / (Vsafe^2)
    r5 = k5 * αCp / Vsafe
    r6 = c6 * (1 + (Cp + αCp) / (Kd * Vsafe)) * (CppC / Vsafe)
    r7 = D1 * RNA  # RNA degradation
    return r1,r2,r3,r4,r5,r6,r7
end

function model!(du, n, p, t, V)
    r1,r2,r3,r4,r5,r6,r7 = rates(n, V, p)
    du[1] = r1 + r5 - r2 - r4           # α̇
    du[2] = r3 + r5 - r4 + r7           # Cṗ
    du[3] = 2r2 - 2r1 + r3 + r4 - r5 + r6   # αCṗ
    du[4] = r1 - r2 - r3 - r6           # CppĊ
    du[5] = r6 - r7                     # RNȦ
end

# -------------------------
# Simulation Function
# -------------------------
function simulate(params, n0, tspan, Vvec)
    T = eltype(Vvec)
    tgrid = range(tspan[1], tspan[2], length=length(Vvec))

    # safeguard: avoid singular cubic spline if repeated values
    if length(unique(Vvec)) < 4
        Vvec .= Vvec .+ 1e-6 .* randn(length(Vvec))
    end

    # Correct argument order: CubicSpline(values, times)
    itp = CubicSpline(collect(Vvec), collect(tgrid))

    f!(du, u, p, t) = model!(du, u, p, t, itp(t))
    prob = ODEProblem(f!, T.(n0), tspan, params)
    sol = solve(prob, Tsit5(); abstol=1e-6, reltol=1e-6, saveat=1.0)
    return sol
end

# -------------------------
# Optimization Function
# -------------------------
function find_optimal_V(N::Int=25)
    V_init = fill((Vmin + Vmax)/2, N)
    lower = fill(Vmin, N)
    upper = fill(Vmax, N)

    println("🔹 Optimizing control function V(t)...")

    # Local objective using fixed params (D1 is constant in params)
    obj(Vvec) = begin
        sol = simulate(params, n0, tspan, Vvec)
        if sol.retcode != :Success
            return 1e6
        end
        RNA_vals = getindex.(sol.u, 5)
        avg_RNA = trapz(sol.t, RNA_vals) / (sol.t[end] - sol.t[1])
        return -avg_RNA
    end

    # Fminbox with Nelder-Mead
    result = optimize(obj, lower, upper, V_init,
                      Fminbox(NelderMead()),
                      Optim.Options(iterations=300, show_trace=true))

    println("✅ Optimization complete. Best mean RNA = ", -result.minimum)
    return result.minimizer
end

# -------------------------
# Run Optimization and Plot
# -------------------------
Vopt = find_optimal_V(100)  
sol_opt = simulate(params, n0, tspan, Vopt)
sol_min = simulate(params, n0, tspan, fill(Vmin, 20))
sol_max = simulate(params, n0, tspan, fill(Vmax, 20))

tvals = sol_opt.t
RNA_opt = getindex.(sol_opt.u, 5)
RNA_min = getindex.(sol_min.u, 5)
RNA_max = getindex.(sol_max.u, 5)

# Plotting
tgrid = range(tspan[1], tspan[2], length=length(Vopt))
p1 = plot(tgrid, Vopt, lw=2, color=:blue, xlabel="Time", ylabel="V(t)",
          title="Optimized Control Function V(t)", label="Vopt(t)")
hline!([Vmin], linestyle=:dash, color=:red, label="Vmin")
hline!([Vmax], linestyle=:dash, color=:green, label="Vmax")

p2 = plot(tvals, RNA_opt, lw=2, color=:blue, label="Vopt(t)",
          xlabel="Time", ylabel="[RNA]", title="RNA Production vs Time")
plot!(tvals, RNA_min, lw=2, color=:red, label="Vmin")
plot!(tvals, RNA_max, lw=2, color=:green, label="Vmax")

plot(p1, p2, layout=(2,1), size=(900,650))


# Compute average RNA for each case
avg_RNA_opt = trapz(sol_opt.t, RNA_opt) / (sol_opt.t[end] - sol_opt.t[1])
avg_RNA_min = trapz(sol_min.t, RNA_min) / (sol_min.t[end] - sol_min.t[1])
avg_RNA_max = trapz(sol_max.t, RNA_max) / (sol_max.t[end] - sol_max.t[1])

println("\n--- Average RNA Summary ---")
println("Vmin: $avg_RNA_min")
println("Vmax: $avg_RNA_max")
println("Vopt: $avg_RNA_opt")

# -------------------------
# Plotting
# -------------------------
# Plot V(t)
tgrid = range(tspan[1], tspan[2], length=length(Vopt))
p1 = plot(tgrid, Vopt, lw=2, color=:blue, xlabel="Time", ylabel="V(t)",
          title="Optimized Control Function V(t)", label="Vopt(t)")
hline!([Vmin], linestyle=:dash, color=:red, label="Vmin")
hline!([Vmax], linestyle=:dash, color=:green, label="Vmax")

# Plot RNA vs time for all three cases
p2 = plot(sol_opt.t, RNA_opt, lw=2, color=:blue, label="Vopt",
          xlabel="Time", ylabel="[RNA]", title="RNA Production vs Time")
plot!(sol_min.t, RNA_min, lw=2, color=:red, label="Vmin")
plot!(sol_max.t, RNA_max, lw=2, color=:green, label="Vmax")

# Combine plots vertically
plot(p1, p2, layout=(2,1), size=(900,650))


α_vals = getindex.(sol_opt.u, 1)
Cp_vals = getindex.(sol_opt.u, 2)
αCp_vals = getindex.(sol_opt.u, 3)
CppC_vals = getindex.(sol_opt.u, 4)

plot(sol_opt.t, α_vals, label="α")
plot!(sol_opt.t, Cp_vals, label="Cp")
plot!(sol_opt.t, αCp_vals, label="αCp")
plot!(sol_opt.t, CppC_vals, label="CppC")








using Random
using DifferentialEquations, Optim, Trapz, Interpolations, Statistics, DataInterpolations, Random, Plots

# -------------------------
# Constants
# -------------------------
const eps_div = 1e-9
Vmin, Vmax = 0.1, 100.0
tspan = (0.0, 100.0)
n0 = [50.0, 500.0, 5.0, 5.0, 0.0]  # [α, Cp, αCp, CppC, RNA]
Npts = 50  # number of points for V(t)

# -------------------------
# Reaction model
# -------------------------
function rates(n::AbstractVector{T}, V::T, p) where T
    α, Cp, αCp, CppC, RNA = n
    k1,k2,k3,k4,k5,c6,Kd,D1 = p
    Vsafe = max(V, T(eps_div))

    r1 = k1 * (αCp^2) / (Vsafe^2)
    r2 = k2 * (α * CppC) / (Vsafe^2)
    r3 = k3 * CppC / Vsafe
    r4 = k4 * (α * Cp) / (Vsafe^2)
    r5 = k5 * αCp / Vsafe
    r6 = c6 * (1 + (Cp + αCp) / (Kd * Vsafe)) * (CppC / Vsafe)
    r7 = D1 * RNA
    return r1,r2,r3,r4,r5,r6,r7
end

function model!(du, n, p, t, V)
    r1,r2,r3,r4,r5,r6,r7 = rates(n, V, p)
    du[1] = r1 + r5 - r2 - r4           # α̇
    du[2] = r3 + r5 - r4 + r7           # Cṗ
    du[3] = 2r2 - 2r1 + r3 + r4 - r5 + r6   # αCṗ
    du[4] = r1 - r2 - r3 - r6           # CppĊ
    du[5] = r6 - r7                     # RNȦ
end

# -------------------------
# Simulation function
# -------------------------
function simulate(params, n0, tspan, Vvec)
    T = eltype(Vvec)
    tgrid = range(tspan[1], tspan[2], length=length(Vvec))

    # safeguard to avoid singular cubic spline
    if length(unique(Vvec)) < 4
        Vvec .= Vvec .+ 1e-6 .* randn(length(Vvec))
    end

    # CubicSpline: values, times
    itp = CubicSpline(collect(Vvec), collect(tgrid))

    f!(du, u, p, t) = model!(du, u, p, t, itp(t))
    prob = ODEProblem(f!, T.(n0), tspan, params)
    sol = solve(prob, Tsit5(); abstol=1e-6, reltol=1e-6, saveat=1.0)
    return sol
end

# -------------------------
# Optimization function
# -------------------------
function find_optimal_V(params, N::Int=Npts)
    V_init = fill((Vmin + Vmax)/2, N)
    lower = fill(Vmin, N)
    upper = fill(Vmax, N)

    obj(Vvec) = begin
        sol = simulate(params, n0, tspan, Vvec)
        if sol.retcode != :Success
            return 1e6
        end
        RNA_vals = getindex.(sol.u, 5)
        avg_RNA = trapz(sol.t, RNA_vals) / (sol.t[end]-sol.t[1])
        return -avg_RNA
    end

    result = optimize(obj, lower, upper, V_init,
                  Fminbox(NelderMead()),
                  Optim.Options(iterations=10000, show_trace=true))

    return result.minimizer, -result.minimum
end

# -------------------------
# Random parameter search
# -------------------------
n_trials = 10  # number of random parameter sets
successes = []

for i in 1:n_trials
    # Random parameter ranges (adjust as needed)
    params_rand = (
        rand(0.1:0.1:20.0),  # k1
        rand(0.01:0.01:1.0), # k2
        rand(0.01:0.01:1.0), # k3
        rand(0.1:0.1:20.0),  # k4
        rand(0.01:0.01:1.0), # k5
        rand(0.1:0.1:5.0),   # c6
        rand(0.1:0.1:10.0),  # Kd
        rand(0.01:0.01:1.0)  # D1
    )

    # Optimize V(t)
    try
        Vopt, avg_RNA_opt = find_optimal_V(params_rand)
    catch
        continue
    end

    # Simulate fixed Vmin/Vmax
    sol_min = simulate(params_rand, n0, tspan, fill(Vmin, Npts))
    sol_max = simulate(params_rand, n0, tspan, fill(Vmax, Npts))
    avg_RNA_min = trapz(sol_min.t, getindex.(sol_min.u,5)) / (sol_min.t[end]-sol_min.t[1])
    avg_RNA_max = trapz(sol_max.t, getindex.(sol_max.u,5)) / (sol_max.t[end]-sol_max.t[1])

    # Check if Vopt > Vmin and Vmax
    if avg_RNA_opt > avg_RNA_min && avg_RNA_opt > avg_RNA_max
        push!(successes, (params=params_rand,
                          avg_RNA_opt=avg_RNA_opt,
                          avg_RNA_min=avg_RNA_min,
                          avg_RNA_max=avg_RNA_max))
    end
end

println("✅ Found ", length(successes), " parameter sets where Vopt > Vmin & Vmax")
for s in successes
    println(s)
end
# -------------------------
>>>>>>> ed8c9783785eaeb6248bc9506a046f52a4f90e9d
