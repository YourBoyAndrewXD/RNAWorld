using Catalyst, Plots, OrdinaryDiffEqDefault

### Reactants: 
# AI: aminoimidazolium
# Cp_:
# Cp:
# Cp_pC:
#
# Primer extension only occurs when a dimer binds - both Cp and Cp_ competitively inhibit extension
# add flow reactions to model the pumping of Cp_pC

polym_model = @reaction_network begin
    @parameters begin
        K_M
        k_obsmax 
        K_i 
        K_eff ~ K_M*(1 + (Cp + Cp_)/K_i)
    end

    k1, Cp_ + Cp_ --> Cp_pC + AI
    k2, AI + Cp_pC --> Cp_ + Cp_ 
    k3, Cp_pC --> Cp_ + Cp
    k4, Cp_ --> Cp + AI
    mm(Cp_pC, k_obsmax, K_eff), Cp_pC + S --> S
    flow, ∅ --> Cp_
    flow, AI --> ∅
end

# Turn on flow events
t = default_t()
@parameters flow t_c
flowrate = 3.

onflow = (t == t_c) => [flow ~ 3.]
