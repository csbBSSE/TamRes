source("simFuncs.R")
source("plotFuncs.R")

# Fig 4B----
simulateFull("core_ESHRMR.topo", paramIndex = 87, workers = 8, noiseLevel = 0) # set workers according to your system specs
dfTraj <- read.csv("core_ESHRMR/87_simDat.csv")
trajectoriesnIC(dfTraj, "Fig4B")

# Fig 4C----
# This was run multiple times to get the right plots. The function starts the simulation from 
# HS initial conditions so that there are better chances of capturing the switching aspects

simSiwtching("core_ESHRMR.topo", ParID = 87, scaledNoise = T, workers = 8)
dfTraj <- read.csv("core_ESHRMR/87_simDat.csv")
ids <- unique(df$Id)
sapply(ids, function(x){
    trajectoriesSwitching(dfTraj, x)
})

# For cases where we see switching, the correspoding time frame can be provided to 
# the above function using the arguments tMin and tMax to get a zoomed version.

# Fig 4E----
simulateFull("core_ESHSMR.topo", paramIndex = 81, workers = 8) # set workers according to your system specs
dfTraj <- read.csv("core_ESHSMR/81_simDat.csv")
trajectoriesnIC(dfTraj, "Fig4E")
landScape(dfTraj, "Fig4E")

# Fig 4F----
simulateFull("core_ESHRMR.topo", paramIndex = 17, workers = 8) # set workers according to your system specs
dfTraj <- read.csv("core_ESHRMR/17_simDat.csv")
trajectoriesnIC(dfTraj, "Fig4E")
landScape(dfTraj, "Fig4E")


#Fig S4A
simulateFull("core_ESHSMR.topo", paramIndex = 66, workers = 8) # set workers according to your system specs
dfTraj <- read.csv("core_ESHSMR/66_simDat.csv")
trajectoriesnIC(dfTraj, "Fig4E")
simulateFull("core_ESHSMR.topo", paramIndex = 8, workers = 8) # set workers according to your system specs
dfTraj <- read.csv("core_ESHSMR/8_simDat.csv")
trajectoriesnIC(dfTraj, "Fig4E")


#Fig S4B
simulateFull("core_ESHRMR.topo", paramIndex = 23, workers = 8) # set workers according to your system specs
dfTraj <- read.csv("core_ESHRMR/23_simDat.csv")
trajectoriesnIC(dfTraj, "FigS4B1")
simulateFull("core_ESHRMR.topo", paramIndex = 87, workers = 8) # set workers according to your system specs
dfTraj <- read.csv("core_ESHRMR/87_simDat.csv")
trajectoriesnIC(dfTraj, "FigS4B2")