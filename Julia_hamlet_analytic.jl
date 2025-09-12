using IdentityByDescentDispersal
using DataFrames, DelimitedFiles, CSV
using Turing, StatsBase, StatsPlots

#-------constant_plot--------------

L = 0.01         # Block length threshold (in Morgans)
#G = 1.0          # Genome length (in Morgans)
G = 13.372 * 4          # Hamlets - Genome length (in Morgans)
r_values = range(0.01, 200.0, length = 200);  # Distances
#
# p_constant = (plot(
#     r_values,
#     expected_ibd_blocks_constant_density.(r_values, 1.0, 0.5, L, G),
#     xlabel = "Geographic distance between individuals (r)",
#     ylabel = "E[#IBD blocks per unit of block length and per pair]",
#     title = "Constant effective density scenario",
#     label = "D=1.0, σ=0.5",
# ))
# plot!(
#     r_values,
#     expected_ibd_blocks_constant_density.(r_values, 2.0, 0.5, L, G),
#     label = "D=2.0, σ=0.5",
# )
# plot!(
#     r_values,
#     expected_ibd_blocks_constant_density.(r_values, 1.0, 0.8, L, G),
#     label = "D=1.0, σ=0.8",
# )
#
# savefig(p_constant, "output_analytical_constant.png")


#---------custom plot -------------------
D = 0.162672     # Effective population density
σ = 20.0         # Dispersal rate

function De(t, θ)
    D₀, α = θ
    max(D₀ * exp(-α * t), eps())
end
p_custom = (plot(
    r_values,
    [expected_ibd_blocks_custom(r, De, [D, 0.0], σ, L, G) for r in r_values],
    xlabel = "Geographic distance between individuals (r)",
    ylabel = "E[#IBD blocks per unit of block length and per pair]",
    title = "Constant effective density scenario",
    label = "growth rate α = 0.0",
))
plot!(
    r_values,
    [expected_ibd_blocks_custom(r, De, [D, 0.005], σ, L, G) for r in r_values],
    label = "growth rate α = 0.005",
)
plot!(
    r_values,
    [expected_ibd_blocks_custom(r, De, [D, -0.005], σ, L, G) for r in r_values],
    label = "growth rate α = -0.005",
)

savefig(p_custom, "output_analytical_custom.png")
