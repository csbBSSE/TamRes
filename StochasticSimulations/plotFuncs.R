library(tidyverse)
options(stringsAsFactors = F)
library(plotly)
source("plot_theme.R")

landScape <- function(dfTraj, nam, pThresh = -3)
{
    kd <- with(dfTraj, MASS::kde2d(EMTScore,TamRes,n = 200))
    names(kd) <- c("EMTScore", "TamRes", "Probability")
    kd$Probability <- log10(kd$Probability)
    kd$Probability[kd$Probability < pThresh] <- NA
    colorAxis <- kd$Probability
    colorAxis[kd$EMTScore < -1, ] <- -1
    colorAxis[kd$EMTScore > 1, ] <- 1
    colorAxis[kd$EMTScore <=1& kd$EMTScore >=-1] <- 0
    colorAxis <- apply(colorAxis, 2, function(x){
        as.character(x)
    })
    
    colorAxis2 <- kd$Probability
    colorAxis2[kd$TamRes < 0, ] <- -1
    colorAxis2[kd$TamRes >= 0, ] <- 1
    colorAxis2 <- apply(colorAxis2, 2, function(x){
        as.character(x)
    })
    tryCatch({
        p <- plot_ly(x = kd$EMTScore, y = kd$TamRes, z = -kd$Probability) %>% 
            add_surface(surfacecolor = colorAxis, 
                        contours = list(
                            z = list(
                                show=TRUE,
                                usecolormap=TRUE,
                                highlightcolor="#ff0000",
                                project=list(z=TRUE)
                            )
                        )) %>%
            layout(
                scene = list(
                    camera=list(
                        eye = list(x=1.1, y=-1, z=0.64)
                    )
                )
            )
        
        orca(p, paste0(nam, "_EMTCol.png"))
    }, error = function(e){conditionMessage(e)})
}

# Generated trajectory plots for multiple initial conditions
trajectoriesnIC <- function(dfTraj, nam)
{
    dfInit <- dfTraj %>% filter(Time == max(dfTraj$Time))
    dfInit$State <- "Epithelial"
    dfInit$State[dfInit$EMTScore > 1] <- "Mesenchymal"
    dfInit$State[dfInit$EMTScore >= -1 & dfInit$EMTScore <= 1] <- "Hybrid"
    dfTraj <- merge(dfTraj, 
                    dfInit %>% select(Id, State), by = "Id",all = T)
    dfNew <- dfTraj %>% select(Id, Time, EMTScore, TamRes, State) %>% 
        gather(key = "Metric", value = "Score", - Time, -Id, -State)
    p <- ggplot()
    for (i in dfNew$Id %>% unique)
    {
        p <- p + geom_line(data = dfNew %>% filter(Id == i), mapping = aes(x = Time, y = Score, color = State))
    }
    p <- p + facet_wrap(~Metric) + theme_Publication() + 
        theme(legend.position = "top")
    ggsave(plot = p, filename = paste0(nam, "_trajecs.png"))
}

#Generates trajectory plots for switching.
trajectoriesSwitching <- function(dfTraj, id, tMin = 0, tMax = NULL)
{
    d <- dfTra %>% filter(Id == id) %>% select(Time, EMTScore, TamSens) %>% 
        gather(key = "Metric", value = "Score", -Time)
    p <- ggplot(d, aes(x = Time, y = Score, color = Metric)) + geom_line() + theme_Publication() + 
        theme(legend.position = "top")
    ggsave(paste0(id,"Normal.png"))
    print(p)
    if (!is.null(tMax))
    {
        p <- ggplot(d %>% filter(Time >=tMin, Time<=tMax), aes(x = Time, y = Score, color = Metric)) + 
            geom_line() + theme_Publication() + 
            theme(legend.position = "top")
        ggsave(paste0(id, "Zoom.png"))
    }
}





