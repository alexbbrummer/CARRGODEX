#### plotting regression on terminal tumor cell number and cell index.

setwd("/Users/abrummer/odrive/Dropbox/Research/city_of_hope/car_T/data/pbt_second_round/")
flow <- read.csv("flow_ci_correlation_data.csv")[,1:6]

str(flow)
flow <- flow[-which(flow$PBT == "PBT138 mid"),]
flow$PBT <- as.factor(flow$PBT)
flow$Plated.Effector <- as.factor(flow$Plated.Effector)
flow$Dex.microgramsperml. <- as.factor(flow$Dex.microgramsperml.)


library(ggplot2)

png(filename = "/Users/abrummer/Desktop/cell_count_cell_index_regression.png", width = 8, height = 6, units = "in", res = 300)
  ggplot(data = flow, aes(x = Residual.CI, y = Residual.Tumor.Count, 
                          color = PBT))+
    geom_point(aes(shape = Dex.microgramsperml.), alpha = 0.5, size = 5)+
    scale_color_manual(name = "Cell Line", 
                         labels = c("PBT030", "PBT128", "PBT138"),
                         values = c("#1b9e77", "#d95f02", "#7570b3"))+
    scale_shape_discrete(name = expression(paste("Dex. dose ", "(", mu, "g/mL)")),
                         labels = c(0, expression(10^-4),expression(10^-3),
                                    expression(10^-2),expression(10^-1),1))+
    xlab("xCelligence cell index values")+
    ylab("Flow cytometry tumor cell count")+
    geom_smooth(data = flow[-which(flow$Plated.Effector == 0),], se = FALSE, method = "lm", 
                formula = y~x, aes(color = PBT))+
    geom_smooth(data = flow, se = FALSE, method = "lm", 
                formula = y~x, aes(color = PBT), linetype = 2)+
    geom_smooth(data = flow, se = FALSE, method = "lm", formula = y~x, color = "black", linetype = 2)+
    geom_smooth(data = flow[-which(flow$Plated.Effector == 0),], se = FALSE, method = "lm", 
                formula = y~x, color = "black")+
    geom_point(data = flow[which(flow$Plated.Effector == 0),], color = "black")+
    theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title = element_text(size=14))
dev.off()


flow <- flow[-which(flow$Plated.Effector == 0),]
# flow <- flow[-which(flow$Plated.Effector < 6),]

pbt138high_fit <- lm(formula = Residual.Tumor.Count~Residual.CI,data = flow, subset = PBT == "PBT138 high" & Plated.Effector != 0)
summary(pbt138high_fit)
pbt138high_fit <- lm(formula = Residual.Tumor.Count~Residual.CI,data = flow, subset = PBT == "PBT138 high")
summary(pbt138high_fit)

pbt128_fit <- lm(formula = Residual.Tumor.Count~Residual.CI,data = flow, subset = PBT == "PBT128+" & Plated.Effector != 0)
summary(pbt128_fit)
pbt128_fit <- lm(formula = Residual.Tumor.Count~Residual.CI,data = flow, subset = PBT == "PBT128+")
summary(pbt128_fit)

pbt030_fit <- lm(formula = Residual.Tumor.Count~Residual.CI,data = flow, subset = PBT == "PBT030" & Plated.Effector != 0)
summary(pbt030_fit)
pbt030_fit <- lm(formula = Residual.Tumor.Count~Residual.CI,data = flow, subset = PBT == "PBT030")
summary(pbt030_fit)

pbt_fit <- lm(formula = Residual.Tumor.Count~Residual.CI, data = flow, subset = Plated.Effector != 0)
summary(pbt_fit)
pbt_fit <- lm(formula = Residual.Tumor.Count~Residual.CI, data = flow)
summary(pbt_fit)

