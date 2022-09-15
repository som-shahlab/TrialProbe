
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

tbl_rows = CSV.File("with_means.csv") |> DataFrame


sort!(tbl_rows, order(:postmeans, by = abs, rev = true))
tbl_rows.unadj_z = tbl_rows.unadjusted_Z ./ tbl_rows.unadjusted_se
tbl_rows.ipw_z = tbl_rows.ipw_Z ./ tbl_rows.ipw_se
tbl_rows.match_z = tbl_rows.match_Z ./ tbl_rows.match_se
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

true_discoveries_match_pointwise =
    (abs.(tbl_rows.match_z) .> 1.96) .& (sign.(tbl_rows.match_z) .== sign.(tbl_rows.postmeans))
false_discoveries_match_pointwise =
    (abs.(tbl_rows.match_z) .> 1.96) .& (sign.(tbl_rows.match_z) .!= sign.(tbl_rows.postmeans))
any_discovery_match_pointwise = abs.(tbl_rows.match_z) .> 1.96


## Cumulative
true_discoveries_unadj = cumsum(true_discoveries_unadjusted_pointwise)
false_discoveries_unadj = cumsum(false_discoveries_unadjusted_pointwise)
total_discoveries_unadj = false_discoveries_unadj .+ true_discoveries_unadj
sign_rate_unadj = 1 .- (false_discoveries_unadj ./ total_discoveries_unadj)
sign_rate_unadj_isotonic = predict(fit(IsotonicRegression, sorted_postmeans, sign_rate_unadj))
true_frac_unadj = true_discoveries_unadj ./ (1:length(true_discoveries_unadj))


true_discoveries_ipw = cumsum(true_discoveries_ipw_pointwise)
false_discoveries_ipw = cumsum(false_discoveries_ipw_pointwise)
total_discoveries_ipw = false_discoveries_ipw .+ true_discoveries_ipw
sign_rate_ipw = 1 .- (false_discoveries_ipw ./ total_discoveries_ipw)
sign_rate_ipw_isotonic = predict(fit(IsotonicRegression, sorted_postmeans, sign_rate_ipw))
true_frac_ipw = true_discoveries_ipw ./ (1:length(true_discoveries_ipw))

true_discoveries_match = cumsum(true_discoveries_match_pointwise)
false_discoveries_match = cumsum(false_discoveries_match_pointwise)
total_discoveries_match = false_discoveries_match .+ true_discoveries_match
sign_rate_match = 1 .- (false_discoveries_match ./ total_discoveries_match)
sign_rate_match_isotonic = predict(fit(IsotonicRegression, sorted_postmeans, sign_rate_match))
true_frac_match = true_discoveries_match ./ (1:length(true_discoveries_match))

index = findfirst(x -> x == 0, sorted_postmeans)
print("Got", sorted_postmeans[index], 
    [sign_rate_unadj[index] sign_rate_match[index] sign_rate_ipw[index]],
    [true_frac_unadj[index] true_frac_match[index] true_frac_ipw[index]],
)

top_ns = vcat(100:1000, 1000:10:length(sorted_postmeans))

plot(
    exp.(abs.(sorted_postmeans))[top_ns],
    [sign_rate_unadj[top_ns] sign_rate_match[top_ns] sign_rate_ipw[top_ns]],
    label = ["Unadjusted Cox" "Propensity Score Matched Cox" "Inverse Propensity Score Weighted Cox"],
    color = [:orange :purple :green],
    xlabel = "Reference Set Odds Ratio Threshold",
    ylabel = "Fraction With Correct Sign",
    legend = :topright,
    xticks = 1:0.25:2.25,
    ylim = (0.5, 1),
    xflip= true,
)
savefig("trialverify_sign_rate.tex")

print(exp.(abs.(sorted_postmeans))[100])
print("\nwat\n")
print(exp.(abs.(sorted_postmeans))[1000])

plot(
    exp.(abs.(sorted_postmeans))[top_ns],
    [true_frac_unadj[top_ns] true_frac_match[top_ns] true_frac_ipw[top_ns]],
    label = ["Unadjusted Cox" "Propensity Score Matched Cox" "Inverse Propensity Score Weighted Cox"],
    color = [:orange :purple :green],
    xlabel = "Reference Set Odds Ratio Threshold",
    ylabel = "Fraction Recovered",
    legend = :topright,
    xflip= true,
    xticks = 1:0.25:2.25,
)
savefig("trialverify_found.tex")
