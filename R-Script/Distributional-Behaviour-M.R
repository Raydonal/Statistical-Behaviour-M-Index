#----------------------------------------------------------
# Statistical properties of an unassisted image quality 
# index for SAR imagery
# Raydonal Ospina  - 2018
# raydonal@de.ufpe.br
#
# Session working to source Location file 
# location in my machine
#setwd("~/R-script")

#----------------------------------------------------------
# Load my scripts
source("cullen.R")

# Require and install this packages
require(multcomp)
require(sandwich)
require(DescTools)
require(ggridges)
require(ggpubr)
require(robust)
require(DandEFA)
require(extrafont)
library(plyr)
#----------------------------------------------------------


#----------------------------------------------------------
# Load the data sets
#----------------------------------------------------------

# Review the path for data file location
#----------------------------------------------------------
# H0
data.0 <-  read.table("../../Data/H0.txt", header = TRUE,sep = "")
data.0$Filter <-  "H0"
data.0 <- data.0 %>% dplyr::select(-Replica)
#----------------------------------------------------------

#----------------------------------------------------------
# H0 with filters
# lee
data.1 <- read.table("../../Data/H0_150_150_E-Lee.txt",      header = TRUE,sep = "")
data.1$Filter <-  "E-LEE" 

# Fans
data.2 <- read.table("../../Data/H0_150_150_FANS.txt",      header = TRUE,sep = "")
data.2$Filter <-  "FANS" 

# Pbp
data.3 <- read.table("../../Data/H0_150_150_PPB.txt",      header = TRUE,sep = "")
data.3$Filter <- "PPB" 

# Srad
data.4 <- read.table("../../Data/H0_150_150_SRAD.txt",      header = TRUE,sep = "")
data.4$Filter <-  "SRAD" 
#----------------------------------------------------------


#----------------------------------------------------------
# H1 with filters
# Degrade

# lee
dataD.1 <- read.table("../../Data/PhantomDegradedH1E-Lee.txt",      header = TRUE,sep = "")
dataD.1$Filter <-  "E-LEE" 

# Fans
dataD.2 <- read.table("../../Data/PhantomDegradedH1FANS.txt",      header = TRUE,sep = "")
dataD.2$Filter <-  "FANS" 
# data.2 =   subset(data.2,  (data.2$ENL==1 | data.2$ENL==4) & data.2$Mask!=20)

# Pbp
dataD.3 <- read.table("../../Data/PhantomDegradedH1PPB.txt",      header = TRUE,sep = "")
dataD.3$Filter <- "PPB" 

# Srad
dataD.4 <- read.table("../../Data/PhantomDegradedH1SRAD.txt",      header = TRUE,sep = "")
dataD.4$Filter <-  "SRAD" 
#----------------------------------------------------------

#----------------------------------------------------------
# H1 with filters
# Steps

# lee
dataS.1 <- read.table("../../Data/PhantomStepE-Lee.txt",      header = TRUE,sep = "")
dataS.1$Filter <-  "E-LEE" 

# Fans
dataS.2 <- read.table("../../Data/PhantomStepFANS.txt",      header = TRUE,sep = "")
dataS.2$Filter <-  "FANS" 

# Pbp
dataS.3 <- read.table("../../Data/PhantomStepPPB.txt",      header = TRUE,sep = "")
dataS.3$Filter <- "PPB" 

# Srad
dataS.4 <- read.table("../../Data/PhantomStepSRAD.txt",      header = TRUE,sep = "")
dataS.4$Filter <-  "SRAD" 
#----------------------------------------------------------

#----------------------------------------------------------
# Join data sets

## All data H0 with H0 Filter - Stability of M 
data.H0  <- data.0 %>% bind_rows(data.1, data.2, data.3, data.4)

## All data H0 with Degrade
data.D  <- data.0 %>% bind_rows(dataD.1, dataD.2, dataD.3, dataD.4 )

## All data H0 with Step
data.S  <- data.0 %>% bind_rows(dataS.1, dataS.2, dataS.3, dataS.4 )
#----------------------------------------------------------

#----------------------------------------------------------
# Rename to evaluate each eperiment
#----------------------------------------------------------


#----------------------------------------------------------
# Experiment 1
# H0 

data <- data.0
#----------------------------------------------------------

#----------------------------------------------------------
# Experiment 2
# H0 with filters
data <- data.H0
#----------------------------------------------------------

#----------------------------------------------------------
# Experimen 3
# Experiments
data <- data.S
#----------------------------------------------------------

#----------------------------------------------------------
# Experiment 4
# Degrade
data <- data.D
#----------------------------------------------------------

#----------------------------------------------------------
# Remove NaN
data =   subset(data,  data$M_estimator!="NaN")
#----------------------------------------------------------

#----------------------------------------------------------
# reshape the data for ggplot2
data$ENL = as.factor(data$ENL)
data$Mask       = as.factor(data$Mask)
data$Tolerance  = as.factor(data$Tolerance)
data$Homogeneous_Areas  = as.numeric(data$Homogeneous_Areas)
data$Filter  = factor(data$Filter, levels = c("H0", "FANS", "PPB", "SRAD", "E-LEE"))
#----------------------------------------------------------

#----------------------------------------------------------
# to use tible and generate tables
data=as.tibble(data)
#----------------------------------------------------------

#----------------------------------------------------------
# mean values
mu <- ddply(data, c("ENL", "Mask", "Tolerance", "Filter"), 
            summarise, grp.mean=mean(M_estimator))
#----------------------------------------------------------

# palette colors
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7") #  "#F0E442", "#0072B2") 

#----------------------------------------------------------
# Plot with Tails

# ggplot - density - M - local
p <-  ggplot(data, aes(y=Filter,  x=M_estimator)) +
  facet_grid(ENL~Tolerance+Mask, labeller = label_parsed, switch='y')
p <-  p+geom_density_ridges(aes(fill=Filter, color=Filter), 
                            alpha=0.4,  scale=3,  size = 0.85, bandwidth=0.05
                            )+
  scale_y_discrete(expand = expand_scale(add = c(0.2, 3.3)))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Filter),
             linetype="dashed")
p <- p+scale_color_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)
p <- p <- p+labs(y = " ", x=TeX('$\\mathit{M}$'), title="")
p <- p +  theme(
  legend.key = element_rect(colour = "white",
                            fill = "white"),
  legend.key.size = unit(1.1, "cm"),
  legend.text = element_text(face = "bold",
                             size=15),
  legend.title = element_text(size=15, face="bold"),
  panel.grid.major = element_line(colour = "gray",
                                  linetype = "dotted"),
  panel.background = element_rect(fill = "white",
                                  colour="black"),
  strip.text.x = element_text(size=12,
                              hjust=0.5,
                              vjust=0.5,
                              face="bold", lineheight = 0.5),
  strip.background = element_rect(colour="black", fill="gray98"),
  axis.text=element_text(size=15, face="bold", colour="gray24"),
    axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
  axis.title=element_text(size=15,face="bold", vjust=0.5),
  plot.title = element_text(size = 13, colour = "black", face="bold", hjust=0.5),
  legend.position = "right", #c(0.1, 0.85), #"none", #  "right", #
  legend.direction="vertical"
)

# name o experiment plot
ggsave("name-experiment.png", p, device="png",  width=14, height = 8)
# #----------------------------------------------------------


# #----------------------------------------------------------
# Run robust rgression for anova

# intercept only
fit.0 = lmRob(M_estimator ~ 1, data=data) 

# + ENL
fit.1= update(fit.0, .~. + ENL, data=data)

# + Mask
fit.2= update(fit.1, .~. + Mask, data=data)

# + Tolerance
fit.3= update(fit.2, .~. + Tolerance, data=data)

# # + Filter
fit.4= update(fit.3, .~. + Filter, data=data)
# #----------------------------------------------------------


# #----------------------------------------------------------
# comparing the models...(Experiment 1)
m <- anova(fit.0, fit.1,fit.2, fit.3)

# comparing the models...(Experiment 2, 3 and 4)
m <- anova(fit.0, fit.1,fit.2, fit.3, fit.4)
# #----------------------------------------------------------


# #----------------------------------------------------------
# Variance (Experiment 1)
(fit.1$r.squared)*100
(fit.2$r.squared)*100
(fit.3$r.squared)*100

(fit.4$r.squared)*100 # + (Experiment 2,3 4)
# #----------------------------------------------------------

# #---------------------------------------------------------
# Run ANOVA

# last model (Experiment 1)
post.model = aov(fit.3)
# #---------------------------------------------------------

# #---------------------------------------------------------
# last model (Experiment 2, 3, 4)
post.model = aov(fit.4)
# #---------------------------------------------------------

# #---------------------------------------------------------
# For intereaction models
post.model = aov(lmRob(M_estimator ~ ENL * Mask * Filter, data = data))
# #---------------------------------------------------------

# #---------------------------------------------------------
# Post Hoc Tests  - Tukey's HSD
PostHocTest(post.model, method = "hsd")
# #---------------------------------------------------------

# #---------------------------------------------------------
# Shapiro Wilk Test - Example
# PPB Filter ENL =1 and ENL=4
ho.exem.1=data$M_estimator[data$ENL=="1" & data$Tolerance=="5" & 
                             data$Mask=="15" & data$Filter=="PPB"]

ho.exem.2=data$M_estimator[data$ENL=="4" & data$Tolerance=="5" & 
                             data$Mask=="15" & data$Filter=="PPB"]

shapiro.test(ho.exem.1)
shapiro.test(ho.exem.2)
# #---------------------------------------------------------


# #---------------------------------------------------------
# Cullen-Frey plots
# PPB Filter ENL =1 and ENL=4
pdf("CullenFreyGraph-H0-1.pdf", width=14, height = 12, pointsize =26)
cullen(ho.exem.1, boot = 1000, boot.col = "#E69F00", obs.pch = 16, obs.col = "red")
dev.off()

pdf("CullenFreyGraph-H0-2.pdf", width=14, height = 12, pointsize =26)
cullen(ho.exem.2, boot = 1000, boot.col = "#E69F00", obs.pch = 16, obs.col = "red")
dev.off()

# #---------------------------------------------------------





