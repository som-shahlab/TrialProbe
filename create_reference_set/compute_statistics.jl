using Pkg
Pkg.activate(".")

using NPZ
using CSV
using StatsBase
using DataFrames
using Plots
using LaTeXStrings
using Empirikos
using HypothesisTests
using Hypatia
using ShapeConstrainedStats
using JSON

pgfplotsx()
empty!(PGFPlotsX.CUSTOM_PREAMBLE)

theme(
    :defaut;
    background_color_legend = :transparent,
    foreground_color_legend = :transparent,
    grid = nothing,
    frame = :box,
    legendfonthalign = :left,
    thickness_scaling = 2.0,
)

entries = readlines("unique_entries.txt")
entries = JSON.parse.(entries)

function to_int(val)
    Int(floor(val))
end

function get_sample(entry)
    table = entry["table"]
    Empirikos.NonCentralHypergeometricSample(
        to_int(table[1][1]),
        to_int(table[1][1] + table[1][2]),
        to_int(table[2][1] + table[2][2]),
        to_int(table[1][1] + table[2][1]),
    )
end

function get_p(entry)
    table = entry["table"]
    pvalue(FisherExactTest(
        to_int(table[1][1]),
        to_int(table[2][1]),
        to_int(table[1][2]),
        to_int(table[2][2]),
    ))
end

Zs_raw = get_sample.(entries)

Ps_raw = get_p.(entries)

# Require at least one event to be able to resolve hazard ratios
_idx = findall(getfield.(Zs_raw, :n) .> 0)
Zs_all = Zs_raw[_idx]
Ps_all = Ps_raw[_idx]

# Deconvolve
npmle_symm = NPMLE(
    convexclass = Empirikos.SymmetricDiscretePriorClass(; support = 0:0.01:8),
    solver = Hypatia.Optimizer,
)
npmle_symm_fit = fit(npmle_symm, Zs_all)

log_ORs = log.(Empirikos.odds_ratio.(Zs_all; offset = 0.000000001))
log_OR_cdf = ecdf([-log_ORs; log_ORs])

print(support(npmle_symm_fit.prior))
print('\n')
print(cdf(npmle_symm_fit.prior, -log(10)))
print("\nFOO\n")
print(log_OR_cdf(-log(10)))

x = support(npmle_symm_fit.prior)

npzwrite("../output/prior_support.npy", x)
npzwrite("../output/log_OR_cdf.npy", log_OR_cdf.(x))
npzwrite("../output/estimated_cdf.npy", cdf.(npmle_symm_fit.prior, x))

# Compute denoised log-odds ratios
postmeans = PosteriorMean.(Zs_all).(npmle_symm_fit.prior)
lower = exp.( quantile.(Empirikos.posterior.(Zs_all, Ref(npmle_symm_fit.prior)), 0.025))
upper = exp.( quantile.(Empirikos.posterior.(Zs_all, Ref(npmle_symm_fit.prior)), 0.975))

npzwrite("../output/log_ORs.npy", log_ORs)
npzwrite("../output/postmeans.npy", postmeans)

# Evaluation

for (i, postmean, lower, upper, p) in zip(_idx, postmeans, lower, upper, Ps_all)
    entries[i]["postmean"] = postmean
    entries[i]["lower"] = lower
    entries[i]["upper"] = upper
    entries[i]["p"] = p
end

first = _idx[1]
print(first, entries[first])
print(postmeans[1], Ps_all[1], Zs_all[1])

result_lines = JSON.json.(entries)

open("reference_set.txt", "w") do io
    for line in result_lines
        println(io, line)
    end
end
