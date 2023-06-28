using Pkg
Pkg.activate(".")

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

# Require at least one event per arm to avoid issues with infinite hazard ratios
_idx = findall(getfield.(Zs_raw, :Z) .> 0)
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

# Plot deconvolved distribution
plot(
    support(npmle_symm_fit.prior),
    x -> cdf(npmle_symm_fit.prior, x),
    legend = :bottomright,
    color = :purple,
    xlabel = L"\textrm{odds ratio } \mu",
    ylabel = "Distribution",
    label = "Estimated",
    xlim = (-3, 3),
    xticks = ([-log(10), -log(3), 0, log(3), log(10)], ["0.1", "0.33", "1", "3", "10"]),
)
plot!(
    support(npmle_symm_fit.prior),
    x -> log_OR_cdf(x),
    color = :grey,
    linestyle = :dash,
    label = "Empirical",
)
savefig("../output/trialverify_distributions.tex")

# Compute denoised log-odds ratios
postmeans = PosteriorMean.(Zs_all).(npmle_symm_fit.prior)
lower = exp.( quantile.(Empirikos.posterior.(Zs_all, Ref(npmle_symm_fit.prior)), 0.025))
upper = exp.( quantile.(Empirikos.posterior.(Zs_all, Ref(npmle_symm_fit.prior)), 0.975))

plot(
    abs.(log_ORs),
    abs.(postmeans),
    seriestype = :scatter,
    alpha = 0.4,
    markerstrokealpha = 0,
    markersize = 1,
    xlim = (-0.1, 3.45),
    ylim = (-0.1, 3.45),
    label = "",
    xlabel = "Sample odds ratio",
    ylabel = "Denoised odds ratio",
    yticks = ([0, log(3), log(10), log(30)], ["1", "3", "10", "30"]),
    xticks = ([0, log(3), log(10), log(30)], ["1", "3", "10", "30"]),
)
plot!([0, 3.45], [0, 3.45], linestyle = :dot, color = :grey, alpha = 0.7, label = "")
savefig("../output/denoise_odds_ratios.tex")

# Evaluation

for (i, postmean, lower, upper, p) in zip(_idx, postmeans, lower, upper, Ps_raw)
    entries[i]["postmean"] = postmean
    entries[i]["lower"] = lower
    entries[i]["upper"] = upper
    entries[i]["p"] = p
end

first = _idx[1]
print(first, entries[first])
print(postmeans[1], Ps_raw[1], Zs_all[1])

result_lines = JSON.json.(entries)

open("reference_set.txt", "w") do io
    for line in result_lines
        println(io, line)
    end
end
