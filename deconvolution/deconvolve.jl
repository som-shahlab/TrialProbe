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
    :default;
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
_idx = findall(getfield.(Zs_raw, :n) .> 0)
Zs_all = Zs_raw[_idx]

# Deconvolve
npmle_symm = NPMLE(
    convexclass = Empirikos.SymmetricDiscretePriorClass(; support = 0:0.01:8),
    solver = Hypatia.Optimizer,
)
npmle_symm_fit = fit(npmle_symm, Zs_all)

log_ORs = log.(Empirikos.odds_ratio.(Zs_all; offset = 0.5))
log_OR_cdf = ecdf([-log_ORs; log_ORs])

# Plot deconvolved distribution
plot(
    support(npmle_symm_fit.prior),
    x -> cdf(npmle_symm_fit.prior, x),
    legend = :bottomright,
    color = :purple,
    xlabel = L"\textrm{odds ratio } \mu",
    ylabel = "Distribution",
    label = "Deconvolved",
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
# savefig("trialverify_distributions.pdf")

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
# savefig("denoise_odds_ratios.pdf")

# Evaluation

tbl_rows = tbl[_idx, :]
tbl_rows.postmeans = copy(postmeans)

CSV.write("with_means.csv", tbl_rows)

sort!(tbl_rows, order(:postmeans, by = abs, rev = true))
tbl_rows.unadj_z = tbl_rows.unadjusted_Z ./ tbl_rows.unadjusted_se
tbl_rows.ipw_z = tbl_rows.ipw_Z ./ tbl_rows.ipw_se
sorted_postmeans = copy(tbl_rows.postmeans)


## Pointwise
true_discoveries_unadjusted_pointwise =
    (abs.(tbl_rows.unadj_z) .> 1.96) .&
    (sign.(tbl_rows.unadj_z) .== sign.(tbl_rows.postmeans))
false_discoveries_unadjusted_pointwise =
    (abs.(tbl_rows.unadj_z) .> 1.96) .&
    (sign.(tbl_rows.unadj_z) .!= sign.(tbl_rows.postmeans))
any_discovery_unadjusted_pointwise = abs.(tbl_rows.unadj_z) .> 1.96

true_discoveries_ipw_pointwise =
    (abs.(tbl_rows.ipw_z) .> 1.96) .& (sign.(tbl_rows.ipw_z) .== sign.(tbl_rows.postmeans))
false_discoveries_ipw_pointwise =
    (abs.(tbl_rows.ipw_z) .> 1.96) .& (sign.(tbl_rows.ipw_z) .!= sign.(tbl_rows.postmeans))
any_discovery_ipw_pointwise = abs.(tbl_rows.ipw_z) .> 1.96


## Cumulative
true_discoveries_unadj = cumsum(true_discoveries_unadjusted_pointwise)
false_discoveries_unadj = cumsum(false_discoveries_unadjusted_pointwise)
total_discoveries_unadj = false_discoveries_unadj .+ true_discoveries_unadj
sign_rate_unadj = false_discoveries_unadj ./ total_discoveries_unadj
sign_rate_unadj_isotonic =
    predict(fit(IsotonicRegression, sorted_postmeans, sign_rate_unadj))


true_discoveries_ipw = cumsum(true_discoveries_ipw_pointwise)
false_discoveries_ipw = cumsum(false_discoveries_ipw_pointwise)
total_discoveries_ipw = false_discoveries_ipw .+ true_discoveries_ipw
sign_rate_ipw = false_discoveries_ipw ./ total_discoveries_ipw
sign_rate_ipw_isotonic = predict(fit(IsotonicRegression, sorted_postmeans, sign_rate_ipw))

top_ns = 100:1000

plot(
    exp.(abs.(sorted_postmeans))[top_ns],
    [total_discoveries_unadj[top_ns] false_discoveries_unadj[top_ns]],
    legend = :topright,
    linestyle = [:solid :dash],
    color = :darkorange,
    label = ["Unadjusted Cox (Total)" "Unadjusted Cox (Discordant sign)"],
    xlabel = "Denoised odds ratio",
    ylabel = "Discoveries",
)

plot!(
    exp.(abs.(sorted_postmeans)[top_ns]),
    [total_discoveries_ipw[top_ns] false_discoveries_ipw[top_ns]],
    legend = :topright,
    linestyle = [:solid :dash],
    label = ["Matched Cox (Total)" "Matched Cox (Discordant sign)"],
    color = :purple,
)
savefig("trialverify_discoveries.pdf")

plot(
    exp.(abs.(sorted_postmeans))[top_ns],
    [sign_rate_unadj_isotonic[top_ns] sign_rate_ipw_isotonic[top_ns]],
    label = ["Unadjusted Cox" "Matched Cox"],
    color = [:orange :purple],
    xlabel = "Denoised odds ratio",
    ylabel = "Discordant sign rate",
    legend = :topright,
    ylim = (0, 0.5),
)
savefig("trialverify_sign_rate.pdf")
