
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


#tbl_rows.metric = copy(-log.(tbl_rows.p_values))
#tbl_rows.metric = copy(-tbl_rows.p_values)
#tbl_rows.metric = copy(-tbl_rows.postmeans)
tbl_rows.metric = copy(exp.(-tbl_rows.postmeans))

sort!(tbl_rows, order(:metric, rev = true))
tbl_rows.unadj_z = tbl_rows.unadjusted_Z ./ tbl_rows.unadjusted_se
tbl_rows.ipw_z = tbl_rows.ipw_Z ./ tbl_rows.ipw_se
tbl_rows.match_z = tbl_rows.match_Z ./ tbl_rows.match_se
sorted_postmeans = copy(tbl_rows.metric)

print(sorted_postmeans)

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
true_frac_unadj = true_discoveries_unadj ./ (1:length(true_discoveries_unadj))


true_discoveries_ipw = cumsum(true_discoveries_ipw_pointwise)
false_discoveries_ipw = cumsum(false_discoveries_ipw_pointwise)
total_discoveries_ipw = false_discoveries_ipw .+ true_discoveries_ipw
sign_rate_ipw = 1 .- (false_discoveries_ipw ./ total_discoveries_ipw)
true_frac_ipw = true_discoveries_ipw ./ (1:length(true_discoveries_ipw))

true_discoveries_match = cumsum(true_discoveries_match_pointwise)
false_discoveries_match = cumsum(false_discoveries_match_pointwise)
total_discoveries_match = false_discoveries_match .+ true_discoveries_match
sign_rate_match = 1 .- (false_discoveries_match ./ total_discoveries_match)
true_frac_match = true_discoveries_match ./ (1:length(true_discoveries_match))


current = 1
x = Vector{Int}()
x_pos = Vector{Float64}()

start = sorted_postmeans[trunc(Int, length(sorted_postmeans) * 0.01)]
stop = sorted_postmeans[length(sorted_postmeans)]

for i in range(start, stop=stop, length=1000)
    global current
    print(i, ' ',current, ' ' , sorted_postmeans[current],'\n')
    while (sorted_postmeans[current] > i && current < length(sorted_postmeans))
        current += 1
    end

    if (length(x) == 0 || x[length(x)] != current)
        push!(x, current)
        push!(x_pos, i)
    end
    
end

x_pos = range(0, 1, length=101)
x = Int.(trunc.(x_pos .* length(sorted_postmeans)))
x[1] = 1

print(x_pos)
print(x)

print(length(x), " what ", length(x_pos), " sure ")

plot(
    x_pos,
    [sign_rate_unadj[x] sign_rate_match[x] sign_rate_ipw[x]],
    label = ["Unadjusted Cox" "PSM Cox" "IPSW Cox"],
    color = [:orange :purple :green],
    xlabel = "Reference Set Odds Ratio Threshold",
    ylabel = "Fraction With Concordant Sign",
    legend = :topright,
    ylim = (0.5, 1),
    xflip= true,
)
savefig("trialverify_sign_rate.tex")

print(exp.(abs.(sorted_postmeans))[100])
print("\nwat\n")
print(exp.(abs.(sorted_postmeans))[1000])

plot(
    x_pos,
    [true_frac_unadj[x] true_frac_match[x] true_frac_ipw[x]],
    label = ["Unadjusted Cox" "PSM Cox" "IPSW Cox"],
    color = [:orange :purple :green],
    xlabel = "Reference Set Odds Ratio Threshold",
    ylabel = "Fraction Recovered",
    legend = :topright,
    xflip= true,
    ylim = (0, 1),
)
savefig("trialverify_found.tex")
