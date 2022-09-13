using Pkg
Pkg.activate(".")

using CSV
using StatsBase
using DataFrames
using Plots
using LaTeXStrings
using Empirikos
using Hypatia
using ShapeConstrainedStats

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

tbl = CSV.File("trialverify_full.csv") |> DataFrame

Zs_raw =
    Empirikos.NonCentralHypergeometricSample.(
        tbl.X1,
        tbl.N1 .+ tbl.X1,
        tbl.N2 .+ tbl.X2,
        tbl.X1 .+ tbl.X2,
    )
_idx = findall(getfield.(Zs_raw, :Z) .> 0)
Zs_all = Zs_raw[_idx]

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
savefig("trialverify_distributions.tex")

# Compute denoised log-odds ratios
postmeans = PosteriorMean.(Zs_all).(npmle_symm_fit.prior)

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
savefig("denoise_odds_ratios.tex")

# Evaluation

tbl_rows = tbl[_idx, :]
tbl_rows.postmeans = copy(postmeans)

CSV.write("with_means.csv", tbl_rows)