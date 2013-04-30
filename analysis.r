# Analysis of CDI Data
setwd("/Volumes/Data/elofgren/Documents/Code/FecalTransplant/")
# Import Packages
library(ggplot2)
library(gridExtra)
library(plyr)

##########################################
## Import Many Data Sets and Merge Them ##
##########################################

# Baseline Data - Import
baseline <- read.csv("baseline.csv")
baseline$Intervention <- "Baseline"

# Antibiotics - Import, Merge and Cleanup

ABXOnly20 <- read.csv("abxonlyresults20.csv")
ABXOnly20$Intervention <- "Antibiotics"

ABXOnly40 <- read.csv("abxonlyresults40.csv")
ABXOnly40$Intervention <- "Antibiotics"

ABXOnly60 <- read.csv("abxonlyresults60.csv")
ABXOnly60$Intervention <- "Antibiotics"

ABXOnly80 <- read.csv("abxonlyresults80.csv")
ABXOnly80$Intervention <- "Antibiotics"

ABXOnly100 <- read.csv("abxonlyresults100.csv")
ABXOnly100$Intervention <- "Antibiotics"

ABX_Master <- rbind(ABXOnly20,ABXOnly40,ABXOnly60,ABXOnly80,ABXOnly100)
rm(ABXOnly20,ABXOnly40,ABXOnly60,ABXOnly80,ABXOnly100)

# PPI and Antibiotics - Import, Merge and Cleanup
ABXPPI20 <- read.csv("abxppiresults20.csv")
ABXPPI20$Intervention <- "Antibiotics/PPIs"

ABXPPI40 <- read.csv("abxppiresults40.csv")
ABXPPI40$Intervention <- "Antibiotics/PPIs"

ABXPPI60 <- read.csv("abxppiresults60.csv")
ABXPPI60$Intervention <- "Antibiotics/PPIs"

ABXPPI80 <- read.csv("abxppiresults80.csv")
ABXPPI80$Intervention <- "Antibiotics/PPIs"

ABXPPI100 <- read.csv("abxppiresults100.csv")
ABXPPI100$Intervention <- "Antibiotics/PPIs"

ABXPPI_Master <- rbind(ABXPPI20,ABXPPI40,ABXPPI60,ABXPPI80,ABXPPI100)
rm(ABXPPI20,ABXPPI40,ABXPPI60,ABXPPI80,ABXPPI100)

# Active CDI Cases Only - Import, Merge and Cleanup
CDIOnly20 <- read.csv("cdionlyresults20.csv")
CDIOnly20$Intervention <- "CDI Cases"

CDIOnly40 <- read.csv("cdionlyresults40.csv")
CDIOnly40$Intervention <- "CDI Cases"

CDIOnly60 <- read.csv("cdionlyresults60.csv")
CDIOnly60$Intervention <- "CDI Cases"

CDIOnly80 <- read.csv("cdionlyresults80.csv")
CDIOnly80$Intervention <- "CDI Cases"

CDIOnly100 <- read.csv("cdionlyresults100.csv")
CDIOnly100$Intervention <- "CDI Cases"

CDI_Master <- rbind(CDIOnly20,CDIOnly40,CDIOnly60,CDIOnly80,CDIOnly100)
rm(CDIOnly20,CDIOnly40,CDIOnly60,CDIOnly80,CDIOnly100)

# Combined Interventions - Import, Merge and Cleanup 
combined20 <- read.csv("combinedresults20.csv")
combined20$Intervention <- "Combined"

combined40 <- read.csv("combinedresults40.csv")
combined40$Intervention <- "Combined"

combined60 <- read.csv("combinedresults60.csv")
combined60$Intervention <- "Combined"

combined80 <- read.csv("combinedresults80.csv")
combined80$Intervention <- "Combined"

combined100 <- read.csv("combinedresults100.csv")
combined100$Intervention <- "Combined"

combined_Master <- rbind(combined20,combined40,combined60,combined80,combined100)
rm(combined20,combined40,combined60,combined80,combined100)
#####

# Merge All
FT_Master <- rbind(baseline,CDI_Master,ABX_Master, ABXPPI_Master, combined_Master)

# Append Baseline data to other data sets
ABXPPI_Master <- rbind(baseline,ABXPPI_Master)
ABX_Master <- rbind(baseline,ABX_Master)
CDI_Master <- rbind(baseline,CDI_Master)
combined_Master <- rbind(baseline,combined_Master)

###################
# Actual Analysis #
###################
a <- ggplot(CDI_Master, aes(factor(Level),Recur)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Recurrent Cases")
a2 <- ggplot(CDI_Master, aes(factor(Level),Incident)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Incident Cases")
CDIviz1 <- grid.arrange(a,a2,ncol=2)

b <- ggplot(ABX_Master, aes(factor(Level),Recur)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Recurrent Cases")
b2 <- ggplot(ABX_Master, aes(factor(Level),Incident)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Incident Cases")
CDIviz2 <- grid.arrange(b,b2,ncol=2)

c <- ggplot(ABXPPI_Master, aes(factor(Level),Recur)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Recurrent Cases")
c2 <- ggplot(ABXPPI_Master, aes(factor(Level),Incident)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Incident Cases")
CDIviz3 <- grid.arrange(c,c2,ncol=2)

d <- ggplot(combined_Master, aes(factor(Level),Recur)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Recurrent Cases")
d2 <- ggplot(combined_Master, aes(factor(Level),Incident)) + theme_bw(base_size=16) + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Proportion Treated")+ylab("Incident Cases")
CDIviz4 <- grid.arrange(d,d2,ncol=2)

FT_Master$Intervention2 <- factor(FT_Master$Intervention, c("Baseline","CDI Cases","Antibiotics","Antibiotics/PPIs","Combined"))
f <- ggplot(FT_Master, aes(factor(Intervention2),Treated)) + theme_bw(base_size=16)
f + geom_jitter(alpha=0.10) + geom_boxplot(outlier.size=0.0, alpha=0.20) +xlab("Intervention")+ylab("Treatments Administered")

#############################
# Distribution Calculations #
# Non-parametric Tests      #
#############################
# CDI Only
quantile(CDI_Master$Recur[CDI_Master$Level==0])
quantile(CDI_Master$Recur[CDI_Master$Level==20])
quantile(CDI_Master$Recur[CDI_Master$Level==40])
quantile(CDI_Master$Recur[CDI_Master$Level==60])
quantile(CDI_Master$Recur[CDI_Master$Level==80])
quantile(CDI_Master$Recur[CDI_Master$Level==100])
kruskal.test(Recur ~ Level, data=CDI_Master)

quantile(CDI_Master$Incident[CDI_Master$Level==0])
quantile(CDI_Master$Incident[CDI_Master$Level==20])
quantile(CDI_Master$Incident[CDI_Master$Level==40])
quantile(CDI_Master$Incident[CDI_Master$Level==60])
quantile(CDI_Master$Incident[CDI_Master$Level==80])
quantile(CDI_Master$Incident[CDI_Master$Level==100])
kruskal.test(Incident ~ Level, data=CDI_Master)

# Antibiotics Only
quantile(ABX_Master$Recur[CDI_Master$Level==0])
quantile(ABX_Master$Recur[CDI_Master$Level==20])
quantile(ABX_Master$Recur[CDI_Master$Level==40])
quantile(ABX_Master$Recur[CDI_Master$Level==60])
quantile(ABX_Master$Recur[CDI_Master$Level==80])
quantile(ABX_Master$Recur[CDI_Master$Level==100])
kruskal.test(Recur ~ Level, data=ABX_Master)

quantile(ABX_Master$Incident[CDI_Master$Level==0])
quantile(ABX_Master$Incident[CDI_Master$Level==20])
quantile(ABX_Master$Incident[CDI_Master$Level==40])
quantile(ABX_Master$Incident[CDI_Master$Level==60])
quantile(ABX_Master$Incident[CDI_Master$Level==80])
quantile(ABX_Master$Incident[CDI_Master$Level==100])
kruskal.test(Incident ~ Level, data=ABX_Master)

# Antibiotics & PPIs
quantile(ABXPPI_Master$Recur[CDI_Master$Level==0])
quantile(ABXPPI_Master$Recur[CDI_Master$Level==20])
quantile(ABXPPI_Master$Recur[CDI_Master$Level==40])
quantile(ABXPPI_Master$Recur[CDI_Master$Level==60])
quantile(ABXPPI_Master$Recur[CDI_Master$Level==80])
quantile(ABXPPI_Master$Recur[CDI_Master$Level==100])
kruskal.test(Recur ~ Level, data=ABXPPI_Master)

quantile(ABXPPI_Master$Incident[CDI_Master$Level==0])
quantile(ABXPPI_Master$Incident[CDI_Master$Level==20])
quantile(ABXPPI_Master$Incident[CDI_Master$Level==40])
quantile(ABXPPI_Master$Incident[CDI_Master$Level==60])
quantile(ABXPPI_Master$Incident[CDI_Master$Level==80])
quantile(ABXPPI_Master$Incident[CDI_Master$Level==100])
kruskal.test(Incident ~ Level, data=ABXPPI_Master)

# Combined
quantile(combined_Master$Recur[CDI_Master$Level==0])
quantile(combined_Master$Recur[CDI_Master$Level==20])
quantile(combined_Master$Recur[CDI_Master$Level==40])
quantile(combined_Master$Recur[CDI_Master$Level==60])
quantile(combined_Master$Recur[CDI_Master$Level==80])
quantile(combined_Master$Recur[CDI_Master$Level==100])
kruskal.test(Recur ~ Level, data=combined_Master)

quantile(combined_Master$Incident[CDI_Master$Level==0])
quantile(combined_Master$Incident[CDI_Master$Level==20])
quantile(combined_Master$Incident[CDI_Master$Level==40])
quantile(combined_Master$Incident[CDI_Master$Level==60])
quantile(combined_Master$Incident[CDI_Master$Level==80])
quantile(combined_Master$Incident[CDI_Master$Level==100])
kruskal.test(Incident ~ Level, data=combined_Master)