using IdentityByDescentDispersal
using DataFrames, DelimitedFiles, CSV, Plots, GLM, Statistics
using Turing, StatsBase, StatsPlots

#----inputs----
# reading input data

data_ibd_blocks = readdlm("input_ibd_blocks.txt", ',') #in Morgan
ibd_blocks = DataFrame(ID1=data_ibd_blocks[1, :], ID2=data_ibd_blocks[2, :], span=data_ibd_blocks[3, :])

data_distances = readdlm("input_distances.txt", ',')
individual_distances = DataFrame(ID1=data_distances[1, :], ID2=data_distances[2, :], distance=data_distances[3, :])

contig_lengths = vec(readdlm("input_contig_lengths.txt", ',', Float64))

bins = vec(readdlm("input_bins.txt", ',', Float64))

min_length = 0.002 #in Morgan




#-----


df = preprocess_dataset(ibd_blocks, individual_distances, bins, min_length);
open("long_df.txt", "w") do io
    show(io, MIME("text/plain"), df)
end


# ----- Plot empirical data -----

n_pairs = 291

IBD_array = zeros(n_pairs, 4) 
distance_list = zeros(n_pairs)

for bin in 1:4 # Julia loops are inclusive of both endpoints
  for i in 1:4:n_pairs # loops from 0 to n_pairs-1 with step size 4
    IBD_array[i,bin] = df[i+bin-1, 6]
    distance_list[i] = df[i+bin-1, 1]
  end
end

colors = [:red, :blue, :green, :orange]

plot(distance_list, IBD_array[:, 1], 
     xlabel="Distance", 
     ylabel="IBD", 
     title="IBD vs Distance", 
     marker=:circle, 
     seriestype=:scatter,
     label="0.2 - 0.8 cM",
     color=colors[1])

plot!(distance_list, IBD_array[:, 2],
     seriestype=:scatter,
     label="0.8 - 1.5 cM",
     color=colors[2])

plot!(distance_list, IBD_array[:, 3],
     seriestype=:scatter,
     label="1.5 - 2.9 cM",
     color=colors[3])

plot!(distance_list, IBD_array[:, 4],
     seriestype=:scatter,
     label="2.9 - 5.0 cM",
     color=colors[4])

L = [0.006, 0.007, 0.014, 0.021]   # length of the IBD block (in Morgans)
G = 13.372 * 4 * 10^-2             # Hamlets - Genome length (in Morgans)
D = 0.0264122                      # Effective Population Density
σ = 50                             # Root mean square dispersal distance per generation.
r_values = range(0.01, 200.0, length = 200)

for bin in 1:4
  plot!(
      r_values,
      expected_ibd_blocks_constant_density.(r_values, D, σ, L[bin], G),
      color=colors[bin],
      label="Expected IBD"
  )
end

savefig("output_plot_empirical.png")

println("saved fig")



#---- constant maxlike estimation
using Turing, StatsBase, StatsPlots
@model function constant_density(df, contig_lengths)
    D ~ Uniform(0, 100)
    σ ~ Uniform(0, 50)
    Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end

mle_estimate = maximum_likelihood(constant_density(df, contig_lengths))
mle_df = DataFrame(coeftable(mle_estimate)) # computed from the Fisher information matrix

# writing output file
open("output_MLE.txt", "w") do io
    show(io, MIME("text/plain"), mle_df)
end

# @model function power_density(df, contig_lengths)
#     D ~ Truncated(Normal(100, 20), 0, Inf)
#     σ ~ Truncated(Normal(1, 0.1), 0, Inf)
#     β ~ Normal(0, 0.5)
#     Turing.@addlogprob! composite_loglikelihood_power_density(D, β, σ, df, contig_lengths)
# end
# m = power_density(df, contig_lengths)
#
# chains = sample(m, NUTS(), MCMCSerial(), 1000, 4)
# baye = plot(chains)
# savefig(baye, "output_baye.png")

