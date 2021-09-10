### Analyzing CARRGO-Dex model parameters

library(ggplot2)
library(plyr)

setwd('/Users/abrummer/odrive/Dropbox/Research/city_of_hope/car_T/car_T_notebooks/third_round_exports/')

pbt_fits <- read.csv('pbt_1_2_4_fit_report_file_expdex_timeshift.csv')
pbt_fits$dex_.ug.ml. <- factor(pbt_fits$dex_.ug.ml.)
pbt_fits$data_column <- as.factor(pbt_fits$data_column)
pbt_fits$t_cell_start.ci. <- as.factor(pbt_fits$t_cell_start.ci.)

names(pbt_fits)

parameter_means <- names(pbt_fits)[seq(from = 10, to = 22, by = 2)]
parameter_sds <- names(pbt_fits)[seq(from = 11, to = 23, by = 2)]

## Practice filtering for subset of data from which to build barplot.
pbt_fits[which(pbt_fits$tumor_line == "PBT4" & pbt_fits$t_cell_start.ci. == 0.25), c(1,5,10:22)][order(pbt_fits[which(pbt_fits$tumor_line == "PBT4" & pbt_fits$t_cell_start.ci. == 0.25), c("dex_.ug.ml.")]),]

## Below code gets us our barplot for a given variable for all tumor lines and dex treatments.  Currently splitting between t cell treatment and no t cell treatment to match with graphs.

## t_cell == 500 cases are for full trials (growth and death)
## t_cell == 0 cases are for only the growth phase of trials.
parameter_means_fancy <- c(expression("Growth rate"~~rho~~(CI/day)), expression("Carrying capacity"~~Kappa~~(CI^-1)), expression("Dex. on tumor cells"~~c[0]~~(day^-1)), expression("CAR T-cell killing rate"~~kappa[1]~~(CI^-1*day^-1)), expression("CAR T-cell prolif./exhst."~~kappa[2]~~(CI^-1*day^-1)), expression("CAR T-cell death rate"~~theta~~(day^-1)), expression("Dex. on CAR T-cell"~~c[3]~~(day^-1)))


### Here we perform some conversions of constants to correct for working entirely in 
# the cell index space.

## kappa_1 and kappa_2 multiply by 5000.  Multiply all c_0 and c_3 by -1 (the model will convert those to additions).

# pbt_fits$kappac1_val = pbt_fits$kappac1_val*5000
# pbt_fits$kappac1_stderr = pbt_fits$kappac1_stderr*5000
# pbt_fits$c0_val = -1*pbt_fits$c0_val
# pbt_fits$c3_val = -1*pbt_fits$c3_val

# write.csv(x = pbt_fits, file = "pbt_fit_report_file_expdex_timeshift_CARTCIrescaled.csv")

## Now, convert all hours to days.
names(pbt_fits)

pbt_fits[,c("rho_val", "rho_stderr", "c0_val", "c0_stderr", "kappac1_val", "kappac1_stderr", "kappac2_val", "kappac2_stderr", "theta_val", "theta_stderr", "c3_val", "c3_stderr")] <- pbt_fits[,c("rho_val", "rho_stderr", "c0_val", "c0_stderr", "kappac1_val", "kappac1_stderr", "kappac2_val", "kappac2_stderr", "theta_val", "theta_stderr", "c3_val", "c3_stderr")]*24

## And, invert kappainv

pbt_fits$kappainv_val <- 1/pbt_fits$kappainv_val

# ymin <- list(c(0,0,0), c(0,0,0), c())
# ymax <- list(c(10,25,20), c(10, 10, 1), c(0,0,0), c(), )


tumor_line_tech <- c(rep("PBT030", 24), rep("PBT128+", 24), rep("PBT138high",24))

pbt_fits$tumor_line_tech <- tumor_line_tech

## Filter zero t-cell start rows and redefine column for levels.
pbt_fits <- pbt_fits[which(pbt_fits$t_cell_start.ci. != 0),]
pbt_fits$t_cell_start.ci. <- factor(as.character(pbt_fits$t_cell_start.ci.))

for(i in 1:length(parameter_means)){
  for(j in 1:length(levels(pbt_fits$tumor_line))){
    png(filename = paste("/Users/alexwork/Desktop/carrgo_fitting_second_round_boxplots/CARRGO_fitting_parameters_second_round_", strsplit(parameter_means[i], split = "_")[[1]][1], "_", levels(pbt_fits$tumor_line)[j], ".png", sep = ""), width = 5.5, height = 4, units = "in", res = 300, pointsize = 20)
    g <- ggplot(data = subset(pbt_fits,subset = (tumor_line == levels(pbt_fits$tumor_line)[j])), 
                aes_string(x = "t_cell_start.ci.", y = parameter_means[i], 
                           fill = "dex_.ug.ml."))+
      geomf_bar(stat = "identity", position = position_dodge())+
      # geom_errorbar(aes_string(ymin = paste(c(parameter_means[i],
      #                                         parameter_sds[i]), collapse='-'), 
      #                          ymax = paste(c(parameter_means[i],
      #                                         parameter_sds[i]), collapse='+')),
      #               width=0.8, position = position_dodge(0.9))+
      labs(x = "", y = "", title = parameter_means_fancy[i], fill = expression("Dex. "*(g/ml)))+
      labs(x="Initial CAR-T cell (CI)", y=parameter_means_fancy[i], fill = expression("Dex. "*(g/ml)))+
      scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02"))+
      theme(axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            title = element_text(size = 15))
      # ylim(c(0,25))
    print(g)
    dev.off()
  }
}


###############################################################################
###############################################################################

## Modifying datatable for custom boxplots.

library(tidyr)

# First graph for pbt1

# first subset table to columns of interest
brpt_tb_1 <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "rho_val", "c0_val", "theta_val", "c3_val")]

# pivot table for plotting format
brpt_tb_1 <- pivot_longer(data = brpt_tb_1, cols = c("rho_val", "c0_val", "theta_val", "c3_val"), names_to = "var_name_val", values_to = "value")

# subsample table by dex and car t values
brpt_tb_1 = brpt_tb_1[which((brpt_tb_1$tumor_line == "PBT1" & (brpt_tb_1$dex_.ug.ml. == 0 | brpt_tb_1$dex_.ug.ml. == 1e-6) & (brpt_tb_1$t_cell_start.ci. == 0 | brpt_tb_1$t_cell_start.ci. == 0.25) | (brpt_tb_1$tumor_line == "PBT1" & brpt_tb_1$t_cell_start.ci. == 0.125 & brpt_tb_1$dex_.ug.ml. == 1e-10))),]

# View(brpt_tb_1)
treatment <- c(rep("control", 4), rep("cart", 4), rep("dex", 4), rep("cartdex1", 4), rep("cartdex2", 4))
brpt_tb_1 <- data.frame(brpt_tb_1, "treatment" = treatment)

brpt_tb_1$var_name_val <- factor(x = brpt_tb_1$var_name_val, levels = c("rho_val", "c0_val", "theta_val", "c3_val"))
brpt_tb_1$treatment <- factor(x = brpt_tb_1$treatment, levels = c("control", "cart", "dex", "cartdex1", "cartdex2"))

png(filename = "/Users/abrummer/Desktop/pbt030_inverseday_parameters.png", width = 3.5, height = 3.5, units = "in", res = 300)
ggplot(data = brpt_tb_1, aes(fill = treatment, y = value, x = var_name_val))+
  geom_bar(position = "dodge", stat = "identity", alpha = 0.75)+
  scale_fill_manual(values = c('black', '#1b9e77', 'darkgrey', '#e6ab02', '#d95f02'), name = "", labels = c("Control", "CAR T only (0.25CI)", "Dex. only (1e-06g/ml)", "CAR T (0.25CI) + \n Dex. (1e-06g/ml)", "CAR T (0.125CI) + \n Dex. (1e-10g/ml)"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 16))+
  scale_x_discrete(name = "Variables", labels=c("rho_val" = expression(rho), "c0_val" = expression(c[0]), "theta_val" = expression(theta), "c3_val" = expression(c[3])))+
  scale_y_continuous(name = expression(paste("Value (", day^-1, ")")))+
  
  annotate(geom = "text", x = 2.95, y = 4, label = paste(expression(10^-12)), col = 'darkgrey', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 3, xend = 3, y = 3., yend = 0, colour = "darkgrey", size=0.5, arrow=arrow(length = unit(0.2, "cm")))+
  
  annotate(geom = "text", x = 3.25, y = 4, label = paste(expression(10^-3)), col = '#e6ab02', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 3.2, xend = 3.2, y = 3., yend = 0, colour = '#e6ab02', size=0.5, arrow=arrow(length = unit(0.2, "cm")))
    
dev.off()
  
# Second graph for pbt 1.

# first subset table to columns of interest
brpt_tb_1 <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "kappac1_val", "kappac2_val")]

# pivot table for plotting format
brpt_tb_1 <- pivot_longer(data = brpt_tb_1, cols = c("kappac1_val", "kappac2_val"), names_to = "var_name_val", values_to = "value")

# subsample table by dex and car t values
brpt_tb_1 = brpt_tb_1[which((brpt_tb_1$tumor_line == "PBT1" & (brpt_tb_1$dex_.ug.ml. == 0 | brpt_tb_1$dex_.ug.ml. == 1e-6) & (brpt_tb_1$t_cell_start.ci. == 0 | brpt_tb_1$t_cell_start.ci. == 0.25) | (brpt_tb_1$tumor_line == "PBT1" & brpt_tb_1$t_cell_start.ci. == 0.125 & brpt_tb_1$dex_.ug.ml. == 1e-10))),]

# View(brpt_tb_1)
treatment <- c(rep("control", 2), rep("cart", 2), rep("dex", 2), rep("cartdex1", 2), rep("cartdex2", 2))
brpt_tb_1 <- data.frame(brpt_tb_1, "treatment" = treatment)

brpt_tb_1$var_name_val <- factor(x = brpt_tb_1$var_name_val, levels = c("kappac1_val", "kappac2_val"))
brpt_tb_1$treatment <- factor(x = brpt_tb_1$treatment, levels = c("control", "cart", "dex", "cartdex1", "cartdex2"))

png(filename = "/Users/abrummer/Desktop/pbt030_inverseCIday_parameters.png", width = 3.8, height = 4, units = "in", res = 300)
ggplot(data = brpt_tb_1, aes(fill = treatment, y = value, x = var_name_val))+
  geom_bar(position = "dodge", stat = "identity", alpha = 0.75,
           width = 0.5)+
  scale_fill_manual(values = c('black', '#1b9e77', 'darkgrey', '#e6ab02', '#d95f02'), name = "", labels = c("Control", "CAR T only (0.25CI)", "Dex. only (1e-06g/ml)", "CAR T (0.25CI) + \n Dex. (1e-06g/ml)", "CAR T (0.125CI) + \n Dex. (1e-10g/ml)"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 16))+
  scale_x_discrete(name = "Variables", labels=c("kappac1_val" = expression(kappa[1]), "kappac2_val" = expression(kappa[2])))+
  scale_y_continuous(name = expression(paste("Value (", CI^-1, day^-1, ")")))
dev.off()



# Then produce same graphs for pbt 2 and pbt 4

### First graph for pbt2

# first subset table to columns of interest
brpt_tb_1 <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "rho_val", "c0_val", "theta_val", "c3_val")]

# pivot table for plotting format
brpt_tb_1 <- pivot_longer(data = brpt_tb_1, cols = c("rho_val", "c0_val", "theta_val", "c3_val"), names_to = "var_name_val", values_to = "value")

# subsample table by dex and car t values
brpt_tb_1 = brpt_tb_1[which((brpt_tb_1$tumor_line == "PBT2" & (brpt_tb_1$dex_.ug.ml. == 0 | brpt_tb_1$dex_.ug.ml. == 1e-9) & (brpt_tb_1$t_cell_start.ci. == 0 | brpt_tb_1$t_cell_start.ci. == 0.25) | (brpt_tb_1$tumor_line == "PBT2" & brpt_tb_1$t_cell_start.ci. == 0.25 & brpt_tb_1$dex_.ug.ml. == 1e-7))),]

# View(brpt_tb_1)
treatment <- c(rep("control", 4), rep("cart", 4), rep("dex", 4), rep("cartdex1", 4), rep("cartdex2", 4))
brpt_tb_1 <- data.frame(brpt_tb_1, "treatment" = treatment)

brpt_tb_1$var_name_val <- factor(x = brpt_tb_1$var_name_val, levels = c("rho_val", "c0_val", "theta_val", "c3_val"))
brpt_tb_1$treatment <- factor(x = brpt_tb_1$treatment, levels = c("control", "cart", "dex", "cartdex1", "cartdex2"))

png(filename = "/Users/abrummer/Desktop/pbt128_inverseday_parameters.png", width = 3.5, height = 3.5, units = "in", res = 300)
ggplot(data = brpt_tb_1, aes(fill = treatment, y = value, x = var_name_val))+
  geom_bar(position = "dodge", stat = "identity", alpha = 0.75)+
  scale_fill_manual(values = c('black', '#1b9e77', 'darkgrey', '#7570b3', '#66a61e'), name = "", labels = c("Control", "CAR T only (0.25CI)", "Dex. only (1e-09g/ml)", "CAR T (0.25CI) + \n Dex. (1e-09g/ml)", "CAR T (0.25CI) + \n Dex. (1e-7g/ml)"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 16))+
  scale_x_discrete(name = "Variables", labels=c("rho_val" = expression(rho), "c0_val" = expression(c[0]), "theta_val" = expression(theta), "c3_val" = expression(c[3])))+
  scale_y_continuous(name = expression(paste("Value (", day^-1, ")")))
dev.off()


# Second graph for pbt 2.

# first subset table to columns of interest
brpt_tb_1 <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "kappac1_val", "kappac2_val")]

# pivot table for plotting format
brpt_tb_1 <- pivot_longer(data = brpt_tb_1, cols = c("kappac1_val", "kappac2_val"), names_to = "var_name_val", values_to = "value")

# subsample table by dex and car t values
brpt_tb_1 = brpt_tb_1[which((brpt_tb_1$tumor_line == "PBT2" & (brpt_tb_1$dex_.ug.ml. == 0 | brpt_tb_1$dex_.ug.ml. == 1e-9) & (brpt_tb_1$t_cell_start.ci. == 0 | brpt_tb_1$t_cell_start.ci. == 0.25) | (brpt_tb_1$tumor_line == "PBT2" & brpt_tb_1$t_cell_start.ci. == 0.25 & brpt_tb_1$dex_.ug.ml. == 1e-7))),]

# View(brpt_tb_1)
treatment <- c(rep("control", 2), rep("cart", 2), rep("dex", 2), rep("cartdex1", 2), rep("cartdex2", 2))
brpt_tb_1 <- data.frame(brpt_tb_1, "treatment" = treatment)

brpt_tb_1$var_name_val <- factor(x = brpt_tb_1$var_name_val, levels = c("kappac1_val", "kappac2_val"))
brpt_tb_1$treatment <- factor(x = brpt_tb_1$treatment, levels = c("control", "cart", "dex", "cartdex1", "cartdex2"))

png(filename = "/Users/abrummer/Desktop/pbt128_inverseCIday_parameters.png", width = 3.8, height = 4, units = "in", res = 300)
ggplot(data = brpt_tb_1, aes(fill = treatment, y = value, x = var_name_val))+
  geom_bar(position = "dodge", stat = "identity", alpha = 0.75,
           width = 0.5)+
  scale_fill_manual(values = c('black', '#1b9e77', 'darkgrey', '#7570b3', '#66a61e'), name = "", labels = c("Control", "CAR T only (0.25CI)", "Dex. only (1e-07g/ml)", "CAR T (0.25CI) + \n Dex. (1e-09g/ml)", "CAR T (0.25CI) + \n Dex. (1e-06g/ml)"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 16))+
  scale_x_discrete(name = "Variables", labels=c("kappac1_val" = expression(kappa[1]), "kappac2_val" = expression(kappa[2])))+
  scale_y_continuous(name = expression(paste("Value (", CI^-1, day^-1, ")")))
dev.off()


### First graph for pbt4

# first subset table to columns of interest
brpt_tb_1 <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "rho_val", "c0_val", "theta_val", "c3_val")]

# pivot table for plotting format
brpt_tb_1 <- pivot_longer(data = brpt_tb_1, cols = c("rho_val", "c0_val", "theta_val", "c3_val"), names_to = "var_name_val", values_to = "value")

# subsample table by dex and car t values
brpt_tb_1 = brpt_tb_1[which((brpt_tb_1$tumor_line == "PBT4" & (brpt_tb_1$dex_.ug.ml. == 0 | brpt_tb_1$dex_.ug.ml. == 1e-9) & (brpt_tb_1$t_cell_start.ci. == 0 | brpt_tb_1$t_cell_start.ci. == 0.05) | (brpt_tb_1$tumor_line == "PBT4" & brpt_tb_1$t_cell_start.ci. == 0.05 & brpt_tb_1$dex_.ug.ml. == 1e-6))),]

# View(brpt_tb_1)
treatment <- c(rep("control", 4), rep("cart", 4), rep("dex", 4), rep("cartdex1", 4), rep("cartdex2", 4))
brpt_tb_1 <- data.frame(brpt_tb_1, "treatment" = treatment)

brpt_tb_1$var_name_val <- factor(x = brpt_tb_1$var_name_val, levels = c("rho_val", "c0_val", "theta_val", "c3_val"))
brpt_tb_1$treatment <- factor(x = brpt_tb_1$treatment, levels = c("control", "cart", "dex", "cartdex1", "cartdex2"))

png(filename = "/Users/abrummer/Desktop/pbt138_inverseday_parameters.png", width = 3.5, height = 3.5, units = "in", res = 300)
ggplot(data = brpt_tb_1, aes(fill = treatment, y = value, x = var_name_val))+
  geom_bar(position = "dodge", stat = "identity", alpha = 0.75)+
  scale_fill_manual(values = c('black', '#1b9e77', 'grey', '#7570b3', '#e6ab02'), name = "", labels = c("Control", "CAR T only (0.05CI)", "Dex. only (1e-09g/ml)", "CAR T (0.05CI) + \n Dex. (1e-09g/ml)", "CAR T (0.05CI) + \n Dex. (1e-06g/ml)"))+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 16))+
  scale_x_discrete(name = "Variables", labels=c("rho_val" = expression(rho), "c0_val" = expression(c[0]), "theta_val" = expression(theta), "c3_val" = expression(c[3])))+
  
  annotate(geom = "text", x = 2.95, y = 3, label = paste(expression(10^-4)), col = 'darkgrey', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 3, xend = 3, y = 2.5, yend = 0, colour = "darkgrey", size=0.5, arrow=arrow(length = unit(0.2, "cm")))+
  
  annotate(geom = "text", x = 3.25, y = 3, label = paste(expression(10^-4)), col = '#7570b3', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 3.2, xend = 3.2, y = 2.5, yend = 0, colour = '#7570b3', size=0.5, arrow=arrow(length = unit(0.2, "cm")))+
  
  scale_y_continuous(name = expression(paste("Value (", day^-1, ")")))
dev.off()

### Second graph for pbt 4.

# first subset table to columns of interest
brpt_tb_1 <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "kappac1_val", "kappac2_val")]

# pivot table for plotting format
brpt_tb_1 <- pivot_longer(data = brpt_tb_1, cols = c("kappac1_val", "kappac2_val"), names_to = "var_name_val", values_to = "value")

# subsample table by dex and car t values
brpt_tb_1 = brpt_tb_1[which((brpt_tb_1$tumor_line == "PBT4" & (brpt_tb_1$dex_.ug.ml. == 0 | brpt_tb_1$dex_.ug.ml. == 1e-9) & (brpt_tb_1$t_cell_start.ci. == 0 | brpt_tb_1$t_cell_start.ci. == 0.05) | (brpt_tb_1$tumor_line == "PBT4" & brpt_tb_1$t_cell_start.ci. == 0.05 & brpt_tb_1$dex_.ug.ml. == 1e-6))),]

# View(brpt_tb_1)
treatment <- c(rep("control", 2), rep("cart", 2), rep("dex", 2), rep("cartdex1", 2), rep("cartdex2", 2))
brpt_tb_1 <- data.frame(brpt_tb_1, "treatment" = treatment)

brpt_tb_1$var_name_val <- factor(x = brpt_tb_1$var_name_val, levels = c("kappac1_val", "kappac2_val"))
brpt_tb_1$treatment <- factor(x = brpt_tb_1$treatment, levels = c("control", "cart", "dex", "cartdex1", "cartdex2"))

png(filename = "/Users/abrummer/Desktop/pbt138_inverseCIday_parameters.png", width = 3.8, height = 4, units = "in", res = 300)
ggplot(data = brpt_tb_1, aes(fill = treatment, y = value, x = var_name_val))+
  geom_bar(position = "dodge", stat = "identity", alpha = 0.75, width = 0.5)+
  theme_bw()+
  theme(legend.position = "none", axis.title = element_text(size = 14), axis.text = element_text(size = 16))+
  scale_fill_manual(values = c('black', '#1b9e77', 'grey', '#7570b3', '#e6ab02'), name = "", labels = c("Control", "CAR T only (0.05CI)", "Dex. only (1e-09g/ml)", "CAR T (0.05CI) + \n Dex. (1e-09g/ml)", "CAR T (0.05CI) + \n Dex. (1e-06g/ml)"))+
  scale_x_discrete(name = "Variables", labels=c("kappac1_val" = expression(kappa[1]), "kappac2_val" = expression(kappa[2])))+
  scale_y_continuous(name = expression(paste("Value (", CI^-1, day^-1, ")")), labels = c("0", "25", "50", "75", "100"), limits = c(0,100))+
  
  annotate(geom = "text", x = 1.9, y = 50, label = paste(expression(10^-1)), col = '#1b9e77', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 1.9, xend = 1.9, y = 42, yend = 0, colour = '#1b9e77', size=0.5, arrow=arrow(length = unit(0.2, "cm")))+
  
  annotate(geom = "text", x = 2.075, y = 50, label = paste(expression(10^-1)), col = '#7570b3', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 2.1, xend = 2.1, y = 42, yend = 0, colour = '#7570b3', size=0.5, arrow=arrow(length = unit(0.2, "cm")))+
  
  annotate(geom = "text", x = 2.25, y = 50, label = paste(expression(10^-1)), col = '#e6ab02', size = 5, angle = 45, parse = TRUE)+
  annotate("segment", x = 2.2, xend = 2.2, y = 42, yend = 0, colour = '#e6ab02', size=0.5, arrow=arrow(length = unit(0.2, "cm")))

dev.off()



###########################################################################
###############  Building graph of theta/kappa2 for death/life comp.
###########################################################################


# first subset table to columns of interest
dead_alive <- pbt_fits[,c("tumor_line", "dex_.ug.ml.", "t_cell_start.ci.", "theta.kappac2")]
dead_alive <- dead_alive[which(dead_alive$theta.kappac2 != "#DIV/0!"),]
dead_alive$theta.kappac2 <- as.numeric(as.character(dead_alive$theta.kappac2))
# da_df <- gather(data = dead_alive, key = "series_type", value = "value", c(1,3,4))
test_paste <- paste(dead_alive$tumor_line, dead_alive$t_cell_start.ci., sep = "_")
dead_alive <- data.frame(dead_alive, test_paste)

# dead_alive$dex_.ug.ml. <- as.numeric(as.character(dead_alive$dex_.ug.ml.))

## Replacing extraneous values with NA's for better graphing clarity.
dead_alive$theta.kappac2[which(dead_alive$test_paste == "PBT2_0.05")] <- NA
dead_alive$theta.kappac2[which(dead_alive$test_paste == "PBT1_0.05")] <- NA
dead_alive$theta.kappac2[which(dead_alive$test_paste == "PBT4_0.25")] <- NA
dead_alive$theta.kappac2[which(dead_alive$test_paste == "PBT2_0.125")] <- NA
dead_alive$theta.kappac2[which(dead_alive$test_paste == "PBT4_0.125")] <- NA

png(filename = "/Users/abrummer/Desktop/tumor_prolif_death_thetakappa2.png", width = 6, height = 4, units = "in", res = 300)
ggplot(data = dead_alive, aes(x = dex_.ug.ml., y = theta.kappac2, shape = t_cell_start.ci., col = tumor_line, group = test_paste))+
  annotate(geom = "rect", xmin = 0, xmax = 7, ymin = 0.4, 
           ymax = 1.5,
           fill = "lightgrey", colour = "lightgrey", alpha = 0.4)+
  geom_point(size = 4, alpha = 0.75)+
  geom_line(alpha = 0.75)+
  scale_color_discrete(name = "Tumor Line", labels = c("PBT030", "PBT128+", "PBT138"))+
  scale_shape_discrete(name = "Initial CAR \nT-cell E:T", labels = c("1:20", "1:8", "1:4"))+
  scale_x_discrete(name = expression(paste("Initial dexamethasone concentration ", "(",mu,"g/mL)")), labels = c("0", expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(1)))+
  scale_y_continuous(name = expression(paste("Predicted final tumor cell population ", theta/kappa[2], " (CI)")),
                     #breaks = c(0, 0.4, 1, 2, 3, 4),
                     #minor_breaks = c(0.5, 1.5, 2.5, 3.5),
                     limits = c(-0.25, 1.5))+
  theme_bw()+
  annotate("text", x = 3.5, y = 1.25, label = "TUMOR CELL PROGRESSION", size = 5, col = "black")+
  annotate("text", x = 3.5, y = -0.15, label = "TUMOR CELL DEATH", size = 5, col = "black")
dev.off()
  # ylim(c(-0.5, 4))
  
