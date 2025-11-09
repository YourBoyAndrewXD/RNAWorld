
###### Plot concentrations over time
# Plot utilities
p3 = plot(sol_spline, label=["α" "Cp" "αCp" "CppC" "RNA"], lw=2, xlabel = "Time", ylabel = "Concentration", title = "System Species over time with optimal V", legend = :topright, size = (1200, 900))

# Plot costates λ₁ to λ₅ over time
p4 = plot(
          costate_trajectory, label=["λ₁" "λ₂" "λ₃" "λ₄" "λ₅"], lw=2,
    xlabel="Time", ylabel="Costate Variables",
    title="Costate Variables λ₁ to λ₅ over Time",
    legend=:right, grid=true, size = (1200, 900)
)
