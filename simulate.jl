include("rna_model.jl")

n0 = [0.0, 0.0, 100.0, 0.0, 0.0]    # initial state
λ0 = [0.0, 0.0, 0.0, 0.0, 1.0]   # final costate
tspan = (0.0, 200.0)
controller, state_trajectory, costate_trajectory = ForwardBackwardSweep(t -> 1., Rodas5P(), n0, λ0, tspan; maxiters = 200);

# plot controller
p1 = plot(controller, label="True V(t)", lw = 2, size = (1200, 900))

(Vl, Vh) = (1., 5.)

# --- ODE system with spline-based V(t) ---
ode_func_spline!(dy, y, p, t) = state_dynamics!(dy, y, controller(t))
prob_spline = ODEProblem(ode_func_spline!, n0, tspan)
sol_spline = solve(prob_spline, Rodas5P())

# --- ODE system with constant V=10 ---
ode_func_hi!(dy, y, p, t) = state_dynamics!(dy, y, Vh)
prob_hi = ODEProblem(ode_func_hi!, n0, tspan)
sol_hi = solve(prob_hi, Rodas5P())

ode_func_low!(dy, y, p, t) = state_dynamics!(dy, y, Vl)
prob_low = ODEProblem(ode_func_low!, n0, tspan)
sol_low = solve(prob_low, Rodas5P())

# --- Plot comparison ---
p2 = plot(sol_spline; idxs = 5, label="RNA (optimal control)", lw=2, xlabel="Time", ylabel="RNA",
          title="RNA Production: Optimal V(t) vs Constant V", legend=:bottomright, size = (1200, 900))
plot!(p2, sol_low; idxs = 5, label="RNA (V = $Vl)", lw=2, linestyle=:dash)
plot!(p2, sol_hi; idxs = 5, label="RNA (V = $Vh)", lw=2, linestyle=:dash)

# Plot concentrations over time
p3 = plot(sol_spline, label=["α" "Cp" "αCp" "CppC" "RNA"], lw=2, xlabel = "Time", ylabel = "Concentration", title = "System Species over time with optimal V", legend = :topright, size = (1200, 900))

# Plot costates λ₁ to λ₅ over time
p4 = plot(
          costate_trajectory, label=["λ₁" "λ₂" "λ₃" "λ₄" "λ₅"], lw=2,
    xlabel="Time", ylabel="Costate Variables",
    title="Costate Variables λ₁ to λ₅ over Time",
    legend=:right, grid=true, size = (1200, 900)
)
plot(p1, p2, p3, p4, layout = (2,2))
savefig("fig.png")
