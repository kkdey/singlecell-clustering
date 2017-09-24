
###############  Zero Inflated Poisson regression problem (Ecology)  #####################


require(ggplot2)
require(pscl)
require(boot)

zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv")
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})

summary(zinb)

ggplot(zinb, aes(count)) + geom_histogram() + scale_x_log10()

summary(m1 <- zeroinfl(count ~ child + camper | persons + camper, data = zinb))
