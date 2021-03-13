using Random
using StatsBase
using Plots
using Distributions
using StatsPlots
using Plots.Measures

sigmoidal(x) =  exp(x)/(exp(x)+0.6)
#%%

plot(Normal(-2,0.4), fill=(0, .5,"#63C5DA"),linecolor ="#241571", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
plot!(Normal(2,0.4), fill=(0, .5,"#63C5DA"),linecolor ="#241571", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
plot!(Normal(-2,1), fill=(0, .5,"#EFB261"),linecolor ="#E66C2C", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
plot!(Normal(2,1), fill=(0, .5,"#EFB261"),linecolor ="#E66C2C", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
p = twinx()
plot!(p,sigmoidal, -6:0.01:6, legend = false, grid =false, dpi = 800,linecolor ="Red", linewidth = 2,box = true,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
savefig("Figure5A.svg")
#%%

plot(Normal(-1,0.5), fill=(0, .5,"#63C5DA"),linecolor ="#241571", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
plot!(Normal(-0.5,0.5), fill=(0, .5,"#80b280"),linecolor ="#338333", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
plot!(Normal(0,0.5), fill=(0, .5,"#EFB261"),linecolor ="#E66C2C", legend = false, grid =false, dpi = 800,left_margin = 15mm, right_margin = 15mm,fmt =:svg)
p = twinx()
plot!(p,sigmoidal, -4:0.01:4, legend = false, grid =false, dpi = 800,linecolor ="Red", linewidth = 2,box = true,left_margin = 15mm, right_margin = 15mm,fmt =:svg)

savefig("FigureS5B.svg")
