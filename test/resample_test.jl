
using Distributions
using StatsBase
using Base.Test


# Weights

wk = [0.20; 0.50; 0.25; 0.05]

wk = wk / sum(wk)

# Perform resampling

res = Int64[]

for i = 1:100000

    append!(res, resample(wk))

end

# Check results

nr = counts(res, 1:maximum(res))

@test all(wk .== round(nr / sum(nr), 2))
