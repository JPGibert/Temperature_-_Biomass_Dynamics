
setwd("~/Desktop/JP/Papers_in_progress/JP_Tetra_Pheno_Temps")

## Load Needed packages
library("dplyr")
library("ggplot2")
library("ggExtra")
library("ggridges")
library("tidyr")
library("mgcv")
library("multispatialCCM")
library("RColorBrewer")
library("viridis")
library("lattice")
library("gridExtra")

##-------------------------------------------------------------------------------
## 1) DATA LOADING
data1 <- read.csv("full_dat_1.csv")
data2 <- read.csv("full_dat_2.csv")
data3 <- read.csv("full_dat_3.csv")
data4 <- read.csv("full_dat_4.csv")
data <- rbind(data1,data2, data3, data4)       
            
##-------------------------------------------------------------------------------
## 2) DENSITIES OVER TIME

# Summarize density data by treatment, day and jar
new_dat <- data %>%
  group_by(Treatment,Temperature, Fluctuations, Day, Jar) %>%
  summarize(Density = mean(Density, na.rm = TRUE))
print(new_dat, n=Inf)

# As before but only for means across Jars and coefficient of variation
means_new_dat <- data %>%
	group_by(Treatment, Temperature, Fluctuations, Day) %>%
	summarize(Abundance=mean(Density, na.rm=TRUE), CV=sd(Density, na.rm=TRUE)/mean(Density, na.rm=TRUE))

# As before but only for means and coefficient of variation across Treatments
treat_new_dat <- data %>%
	group_by(Treatment) %>%
	summarize(Abundance=mean(Density, na.rm=TRUE), CV=sd(Density, na.rm=TRUE)/mean(Density, na.rm=TRUE))

# Plot!
plot(Density~jitter(Day,2), data=new_dat, col=c("blue","orange")[as.factor(new_dat$Temperature)], pch=c(1,16)[as.factor(new_dat$Fluctuations)])
	


##-------------------------------------------------------------------------------
## 3) TRAITS OVER TIME

# Summarize trait data by treatment, day and jar
trait_dat <- data %>%
  group_by(Treatment,Temperature, Fluctuations,  Day, Jar) %>%
  summarize(Vol = mean(Volume, na.rm = TRUE), CV = sd(Volume, na.rm=TRUE)/mean(Volume, na.rm = TRUE),SD = sd(Volume, na.rm=TRUE), MIN=min(Volume, na.rm=TRUE), sumMass = sum(Volume*10^-12)) 

# We need to lower initial biomass as currently estimated from stock culture that had much larger density
trait_dat$sumMass[which(trait_dat$Day==1)] <- 3.57*10^-7

trait_dat <- trait_dat %>%
	group_by(Treatment,Temperature, Fluctuations,  Day, Jar) %>%
	mutate(Sel=c(NA,diff(Vol)))

# Summarize trait data by treatment, day
trait_mean_dat <- data %>%
  group_by(Treatment,Temperature, Fluctuations,  Day) %>%
  summarize(Vol = mean(Volume, na.rm = TRUE), CV = sd(Volume, na.rm=TRUE)/mean(Volume, na.rm = TRUE))

# Plot!
plot(log(trait_dat$Vol) ~ trait_dat$Day, col=c("black","grey","red","orange","blue","light blue")[as.factor(trait_dat$Treatment)], pch=16)

##-------------------------------------------------------------------------------
## 4) Effects of treatments on densities GAMM
	
 # Random effects model with Jar as random effect
new_dat$Temperature <- as.factor(new_dat$Temperature)
new_dat <- new_dat %>% mutate(logDensity=log(Density))
means_new_dat <- means_new_dat %>% mutate(logDensity=log(Abundance))
trait_dat$Fluctuations <- factor(new_dat$Fluctuations, levels=c("Stable","Fluct"))
trait_dat <- trait_dat %>% mutate(Mass=Vol*(10^-12),logMass=log(Mass), logVol=log(Vol))

	# For the conversion, we rescale Vol to convert to mass (multiply by 10^-12) and then rescale the biomass because density already accounts for flowcam inefficiency in counting cells but volume doesn't, so we use the fact that efficiency is roughly 70% to rescale the total biomass (divide by 0.28).
biomass_data <- merge(new_dat,trait_dat,by=c("Treatment","Temperature", "Fluctuations", "Day", "Jar")) %>%
	mutate(Bmass=sumMass/0.279855, logBmass=log(Bmass))

# CHOOSING BEST DENSITY MODEL (log vs non-log, is the equivalent of a good fit at first or a good fit at carrying capacity)

model.dens_1 <- gamm(logDensity ~ s(Day)+Temperature*Fluctuations,
			random=list(Jar=~1),
           data = new_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

model.dens_2 <- gamm(Density ~ s(Day, by=Fluctuations)+Fluctuations*Temperature,
			random=list(Jar=~1),
           data = new_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1), 
           na.action=na.omit)


# Are fluctuation interactions needed?
model.dens_3 <- gamm(logDensity ~ s(Day) + Temperature+Fluctuations,
			random=list(Jar=~1),
           data = new_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

anova(model.dens_1$lme,model.dens_2$lme) ## no "by=" better than with by=
anova(model.dens_1$lme,model.dens_3$lme) ## Comparable so we can drop interaction

summary(model.dens_1$gam)
summary(model.dens_3$gam) ## No effect of fluctuations but small positive effect of temperature


# CHOOSING BEST MASS MODEL
model.mass <- gamm(logMass ~ s(Day)+Temperature*Fluctuations,
			random=list(Jar=~1),
           data = trait_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

model.mass_1 <- gamm(logMass ~ s(Day)+Temperature+Fluctuations,
			random=list(Jar=~1),
           data = trait_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

model.mass_2 <- gamm(logMass ~ s(Day) + Fluctuations,
			random=list(Jar=~1),
           data = trait_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

model.mass_3 <- gamm(logMass ~ s(Day, by=Fluctuations)+Temperature*Fluctuations,
			random=list(Jar=~1),
           data = trait_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)


anova(model.mass$lme,model.mass_3$lme) # Without "by=" much better model.
anova(model.mass$lme,model.mass_1$lme) # Models are very similar, may drop interaction
anova(model.mass_1$lme,model.mass_2$lme) # Models very similar so almost no effect of temperature

summary(model.mass_1$gam) # Best model has additive effects of temp and fluctuations


# CHOOSING BEST BIOMASS MODEL
model.BM <- gamm(logBmass ~ s(Day)+Temperature*Fluctuations,
			random=list(Jar=~1),
           data = biomass_data, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

model.BM_1 <- gamm(logBmass ~ s(Day)+Temperature+Fluctuations,
			random=list(Jar=~1),
           data = biomass_data, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

anova(model.BM$lme,model.BM_1$lme) # No difference,m so can drop interaction

summary(model.BM_1$gam)

## Prep data for model predictions

##1) DENSITY GAMM Predictions
density_22_stable <- data.frame(logDensity=seq(min(new_dat$logDensity),max(new_dat$logDensity),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Stable",100),Temperature=rep(22,100))
density_22_fluct <- data.frame(logDensity=seq(min(new_dat$logDensity),max(new_dat$logDensity),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Fluct",100),Temperature=rep(22,100))
density_25_stable <- data.frame(logDensity=seq(min(new_dat$logDensity),max(new_dat$logDensity),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Stable",100),Temperature=rep(25,100))
density_25_fluct <- data.frame(logDensity=seq(min(new_dat$logDensity),max(new_dat$logDensity),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Fluct",100),Temperature=rep(25,100))

density_22_stable$Pred <- predict(model.dens_3$gam, density_22_stable)
density_22_fluct$Pred <- predict(model.dens_3$gam, density_22_fluct)
density_25_stable$Pred <- predict(model.dens_3$gam, density_25_stable)
density_25_fluct$Pred <- predict(model.dens_3$gam,density_25_fluct)

##3) TRAIT GAMM Predictions
Vol_22_stable <- data.frame(logMass=seq(min(trait_dat$logMass),max(trait_dat$logMass),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Stable",100),Temperature=rep(22,100))
Vol_22_fluct <- data.frame(logMass=seq(min(trait_dat$logMass),max(trait_dat$logMass),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Fluct",100),Temperature=rep(22,100))
Vol_25_stable <- data.frame(logMass=seq(min(trait_dat$logMass),max(trait_dat$logMass),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Stable",100),Temperature=rep(25,100))
Vol_25_fluct <- data.frame(logMass=seq(min(trait_dat$logMass),max(trait_dat$logMass),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Fluct",100),Temperature=rep(25,100))

Vol_22_stable$Pred <- predict(model.mass_1$gam, Vol_22_stable)
Vol_22_fluct$Pred <- predict(model.mass_1$gam, Vol_22_fluct)
Vol_25_stable$Pred <- predict(model.mass_1$gam, Vol_25_stable)
Vol_25_fluct$Pred <- predict(model.mass_1$gam,Vol_25_fluct)

## 4) BIOMASS GAMM Predictions
Bmass_22_stable <- data.frame(logBmass=seq(min(biomass_data$logBmass, na.rm=TRUE),max(biomass_data$logBmass, na.rm=TRUE),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Stable",100),Temperature=rep(22,100))
Bmass_22_fluct <- data.frame(logBmass=seq(min(biomass_data$logBmass, na.rm=TRUE),max(biomass_data$logBmass, na.rm=TRUE),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Fluct",100),Temperature=rep(22,100))
Bmass_25_stable <- data.frame(logBmass=seq(min(biomass_data$logBmass, na.rm=TRUE),max(biomass_data$logBmass, na.rm=TRUE),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Stable",100),Temperature=rep(25,100))
Bmass_25_fluct <- data.frame(logBmass=seq(min(biomass_data$logBmass, na.rm=TRUE),max(biomass_data$logBmass, na.rm=TRUE),length.out=100), Day=seq(1,15,length.out=100),Fluctuations=rep("Fluct",100),Temperature=rep(25,100))

Bmass_22_stable$Pred <- predict(model.BM_1$gam, Bmass_22_stable)
Bmass_22_fluct$Pred <- predict(model.BM_1$gam, Bmass_22_fluct)
Bmass_25_stable$Pred <- predict(model.BM_1$gam, Bmass_25_stable)
Bmass_25_fluct$Pred <- predict(model.BM_1$gam,Bmass_25_fluct)


###############################-------------------------------------------------------
## Plots


biomass_data <- merge(new_dat,trait_dat,by=c("Treatment","Temperature", "Fluctuations", "Day", "Jar")) %>%
	mutate(Bmass=sumMass/0.279855, logBmass=log(Bmass)) %>%
	mutate(unique_jar=paste(biomass_data$Treatment,biomass_data$Jar, sep="_")) %>%
	mutate(Treatment=factor(biomass_data$Treatment, levels = c("22_Stable","22_Flux","25_Stable","25_Flux")), Fluctuations=factor(biomass_data$Fluctuations, levels = c("Stable", "Fluct")))
  	
mean_biomass_data <- biomass_data %>% 	
  	group_by(Treatment, Day, Temperature, Fluctuations) %>%
  	summarize(Bmass = mean(Bmass, na.rm = TRUE), Mass = mean(Mass, na.rm=TRUE), Density = mean(Density, na.rm=TRUE))

## Fig 1

g_1 <- ggplot(biomass_data, aes(Day-1,Bmass, col=unique_jar)) +
	geom_line() +
	geom_point() +
	geom_hline(yintercept=0, linetype="dashed", col="grey", size=1) +
	scale_color_viridis(discrete=TRUE, option="viridis", alpha=1/2) +
	#theme_bw() +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	scale_y_continuous(limits=c(0, 1.05*10^-4),breaks=seq(0,10^-4,2.5*10^-5), labels=c(0,0.25,5.00,7.50,1.00)) +
 	labs(x="", y=paste("Biomass (","\U03bc","g/mL)", sep="")) +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )

g_2 <- ggplot(biomass_data, aes(Day-1,Density, col=unique_jar)) +
	geom_line() +
	geom_point() +
	geom_hline(yintercept=0, linetype="dashed", col="grey", size=1) +
	scale_color_viridis(discrete=TRUE, option="viridis", alpha=1/2) +
	#theme_bw() +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	scale_y_continuous(limits=c(0, 12*10^3),breaks=seq(0,12*10^3,3*10^3), labels=c(0,3,6,9,12)) +
 	labs(x="", y="Density (ind/mL)") +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )
	
g_3 <- ggplot(biomass_data, aes(Day-1,Mass, col=unique_jar)) +
	geom_line() +
	geom_point() +
	geom_hline(yintercept=0, linetype="dashed", col="grey", size=1) +
	scale_color_viridis(discrete=TRUE, option="viridis", alpha=1/2) +
	#theme_bw() +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	scale_y_continuous(limits=c(0.25*10^-8, 2*10^-8),breaks=seq(0,2*10^-8,0.5*10^-8), labels=c(0,0.5,1,1.5,2)) +
 	labs(x="", y=paste("Mass (","\U03bc","g)", sep="")) +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )

# Add some jitter
#jitter <- position_jitter(width = 0.2, height = 0)
g_4 <- ggplot(mean_biomass_data, aes(Day-1,Bmass, col=Temperature)) +
	geom_point(aes(col=Temperature, shape=Fluctuations), size=3) +
	scale_shape_manual(values=c(16, 1)) +
	geom_line(aes(linetype=Fluctuations)) +
	#geom_point(position= jitter, aes(shape= Temperature, colour=Treatment),size=2, alpha=1/2) +
	geom_line(Bmass_22_stable, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[1], size=1) +
	geom_line(Bmass_22_fluct, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[1], linetype="dashed", size=1) +
	geom_line(Bmass_25_stable, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[4], size=1) +
	geom_line(Bmass_25_fluct, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[4], linetype="dashed", size=1) +
	geom_hline(yintercept=0, linetype="dashed", col="grey", size=1) +
	scale_colour_manual(values=
		alpha(rev(c("#EF8A62","#67A9CF")),0.6)) +
	#theme_bw() +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	scale_y_continuous(limits=c(0, 0.75*10^-4),breaks=seq(0,10^-4,2.5*10^-5), labels=c(0,0.25,5.00,7.50,1.00)) +
 	labs(x="Time (days)", y=paste("Biomass (","\U03bc","g/mL)", sep="")) +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )

g_5 <- ggplot(mean_biomass_data, aes(Day-1,Density, col=Temperature)) +
	geom_point(aes(col=Temperature, shape=Fluctuations), size=3) +
	scale_shape_manual(values=c(16, 1)) +
	geom_line(aes(linetype=Fluctuations)) +
	#geom_point(position= jitter,aes(shape= Temperature, colour=Treatment),size=2, alpha=1/2) +
	geom_line(density_22_stable, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[1], size=1) +
	geom_line(density_22_fluct, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[1], linetype="dashed", size=1) +
	geom_line(density_25_stable, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[4], size=1) +
	geom_line(density_25_fluct, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[4], linetype="dashed", size=1) +
	geom_hline(yintercept=0, linetype="dashed", col="grey", size=1) +
	scale_colour_manual(values= alpha(rev(c("#EF8A62","#67A9CF")),0.6)) +
	#theme_bw() +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	scale_y_continuous(limits=c(0, 9*10^3),breaks=seq(0,12*10^3,3*10^3), labels=c(0,3,6,9,12)) +
 	labs(x="Time (days)", y="Density (ind/mL)") +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )
plot(Mass~Day,data=mean_biomass_data)
lines(Pred~Day,data=Vol_22_stable)

g_6 <- ggplot(mean_biomass_data, aes(Day-1,Mass, col=Temperature)) +
	geom_point(aes(col=Temperature, shape=Fluctuations), size=3) +
	scale_shape_manual(values=c(16, 1)) +
	geom_line(aes(linetype=Fluctuations)) +
	geom_line(Vol_22_stable, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[1], size=1) +
	geom_line(Vol_22_fluct, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[1], linetype="dashed", size=1) +
	geom_line(Vol_25_stable, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[4], size=1) +
	geom_line(Vol_25_fluct, mapping=aes(Day-1,exp(Pred)), color=rev(brewer.pal(n = 4, name = "RdBu"))[4], linetype="dashed", size=1) +
	geom_hline(yintercept=0, linetype="dashed", col="grey", size=1) +
	scale_colour_manual(values= alpha(rev(c("#EF8A62","#67A9CF")),0.6)) +
	#theme_bw() +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	scale_y_continuous(limits=c(0.25*10^-8, 1.7*10^-8),breaks=seq(0,2*10^-8,0.5*10^-8), labels=c(0,0.5,1,1.5,2)) +
 	labs(x="Time (days)", y=paste("Mass (","\U03bc","g)", sep="")) +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )   

Bmass_22_stable
biomass_data %>%
	group_by(Treatment, Jar) %>%
	mutate(detrended=)

day_22_stable <- rep(0,14)
day_25_stable <- rep(0,14)
day_22_fluct <- rep(0,14)
day_25_fluct <- rep(0,14)
for(i in 1:14){
	day_22_stable[i] <- which.min(abs(i-density_22_stable$Day))
	day_22_stable[i] <- which.min(abs(i-density_25_stable$Day))
	day_22_stable[i] <- which.min(abs(i-density_22_fluct$Day))
	day_22_stable[i] <- which.min(abs(i-density_25_fluct$Day))
	}

## De-Trending	
model.dens <- gamm(logDensity ~ s(Day),
			random=list(Jar=~1),
           data = new_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)
model.mass <- gamm(logMass ~ s(Day),
			random=list(Jar=~1),
           data = trait_dat, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)
model.bmass <- gamm(logBmass ~ s(Day),
			random=list(Jar=~1),
           data = biomass_data, correlation = corARMA(form = ~ Day|Jar/Treatment, p = 1),
           na.action=na.omit)

## Predictions
dens_pred <- data.frame(Day=seq(1,15,1))
mass_pred <- data.frame(Day=seq(1,15,1))
bmass_pred <- data.frame(Day=seq(1,15,1))
dens_pred$pred_dens <- predict(model.dens$gam, dens_pred)
mass_pred$pred_mass <- predict(model.mass$gam, mass_pred)
bmass_pred$pred_bmass <- predict(model.bmass$gam, bmass_pred)

biomass_dat <- 
	merge(merge(merge(biomass_data, dens_pred, by="Day"),mass_pred,by="Day"),bmass_pred) %>%
	mutate(bmass_unt=logBmass-pred_bmass, mass_unt=logMass-pred_mass, dens_unt=logDensity-pred_dens)

g_7 <- ggplot(biomass_dat, aes(Day,bmass_unt, color=unique_jar)) + 
	geom_point(size=10^-6) +
	geom_line() +
	scale_colour_manual(values= alpha(rev(rep(c("#EF8A62","#67A9CF"),each=12)),0.3)) +
	geom_hline(yintercept=mean(biomass_dat$bmass_unt[which(biomass_dat$Temperature==25)],na.rm=TRUE), linetype="dashed", col="#CA0020", size=1) +
	geom_hline(yintercept=mean(biomass_dat$bmass_unt[which(biomass_dat$Temperature==22)],na.rm=TRUE)
, linetype="dashed", col="#0571B0", size=1) +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	#scale_y_continuous(limits=c(0, 0.75*10^-4),breaks=seq(0,10^-4,2.5*10^-5), labels=c(0,0.25,5.00,7.50,1.00)) +
 	labs(x="Time (days)", y=paste("Biomass (","\U03bc","g/mL)", sep="")) +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )
	
g_8 <- ggplot(biomass_dat, aes(Day,dens_unt, color=unique_jar)) + 
	geom_point(size=10^-6) +
	geom_line() +
	scale_colour_manual(values= alpha(rev(rep(c("#EF8A62","#67A9CF"),each=12)),0.3)) +
	geom_hline(yintercept=mean(biomass_dat$dens_unt[which(biomass_dat$Temperature==25)],na.rm=TRUE), linetype="dashed", col="#CA0020", size=1) +
	geom_hline(yintercept=mean(biomass_dat$dens_unt[which(biomass_dat$Temperature==22)],na.rm=TRUE)
, linetype="dashed", col="#0571B0", size=1) +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	#scale_y_continuous(limits=c(0, 0.75*10^-4),breaks=seq(0,10^-4,2.5*10^-5), labels=c(0,0.25,5.00,7.50,1.00)) +
 	labs(x="Time (days)", y="Density (ind/mL)") +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )

g_9 <- ggplot(biomass_dat, aes(Day,mass_unt, color=unique_jar)) + 
	geom_point(size=10^-6) +
	geom_line() +
	scale_colour_manual(values= alpha(rev(rep(c("#EF8A62","#67A9CF"),each=12)),0.3)) +
	geom_hline(yintercept=mean(biomass_dat$mass_unt[which(biomass_dat$Temperature==25)],na.rm=TRUE), linetype="dashed", col="#CA0020", size=1) +
	geom_hline(yintercept=mean(biomass_dat$mass_unt[which(biomass_dat$Temperature==22)],na.rm=TRUE)
, linetype="dashed", col="#0571B0", size=1) +
	scale_x_continuous(limits=c(-0.1,15), name="",breaks=c(0,5,10,15)) +
 	#scale_y_continuous(limits=c(0, 0.75*10^-4),breaks=seq(0,10^-4,2.5*10^-5), labels=c(0,0.25,5.00,7.50,1.00)) +
 	labs(x="Time (days)", y=paste("Mass (","\U03bc","g)", sep="")) +
  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
  	axis.title=element_text(size=16),
  	axis.line = element_line(color='black', size=0.8),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.length=unit(-.15, "cm"),
    legend.position="none",
    axis.text.x = element_text(size=14, angle=0, vjust=-0.2),
    axis.text.y = element_text(size=14, angle=0, hjust=-0.2)
    )
    
g_10 <- ggMarginal(g_7, groupColour = TRUE, groupFill = TRUE, margins = "y")
g_11 <- ggMarginal(g_8, groupColour = TRUE, groupFill = TRUE, margins = "y")
g_12 <- ggMarginal(g_9, groupColour = TRUE, groupFill = TRUE, margins = "y")
	

#pdf("~filename.pdf", width = 14, height = 8.7) # Open a new pdf file
grid.arrange(g_1, g_2, g_3, g_4,g_5,g_6,g_10,g_11,g_12, nrow=3) # Write the grid.arrange in the file
#dev.off()


##-------------------------------------------------------------------------------
## 5) MULTISPATIAL CCM

## PREP data for Multispatial

Ab_TS_22Stable <- rep(0,15*6+6)
Ab_TS_25Stable <- rep(0,15*6+6)
Ab_TS_22Flux <- rep(0,15*6+6)
Ab_TS_25Flux <- rep(0,15*6+6)

Tr_TS_22Stable <- rep(0,15*6+6)
Tr_TS_25Stable <- rep(0,15*6+6)
Tr_TS_22Flux <- rep(0,15*6+6)
Tr_TS_25Flux <- rep(0,15*6+6)

for(i in 1:6){
	
	# DENSITY DATA
	Dat_22Stable <- new_dat %>% 
	filter(Treatment=="22_Stable" & Jar==i)
	Dat_22Flux <- new_dat %>% 
	filter(Treatment=="22_Flux" & Jar==i)
	Dat_25Stable <- new_dat %>% 
	filter(Treatment=="25_Stable" & Jar==i)
	Dat_25Flux <- new_dat %>% 
	filter(Treatment=="25_Flux" & Jar==i) 	
	Ab_TS_22Stable[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dat_22Stable$Density
	Ab_TS_22Flux[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dat_22Flux$Density
	Ab_TS_25Stable[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dat_25Stable$Density
	Ab_TS_25Flux[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dat_25Flux$Density
	if(i == 6){
		Ab_TS_22Stable <- Ab_TS_22Stable[-1]
		Ab_TS_22Flux <- Ab_TS_22Flux[-1]
		Ab_TS_25Stable <- Ab_TS_25Stable[-1]
		Ab_TS_25Flux <- Ab_TS_25Flux[-1]
	}
	
	# TRAIT DATA
	Dattr_22Stable <- trait_dat %>% 
	filter(Treatment=="22_Stable" & Jar==i)
	Dattr_22Flux <- trait_dat %>% 
	filter(Treatment=="22_Flux" & Jar==i)
	Dattr_25Stable <- trait_dat %>% 
	filter(Treatment=="25_Stable" & Jar==i)
	Dattr_25Flux <- trait_dat %>% 
	filter(Treatment=="25_Flux" & Jar==i) 	
	Tr_TS_22Stable[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dattr_22Stable$Vol
	Tr_TS_22Flux[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dattr_22Flux$Vol
	Tr_TS_25Stable[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dattr_25Stable$Vol
	Tr_TS_25Flux[(i+1+15*(i-1)):(i+15+15*(i-1))] <- Dattr_25Flux$Vol
	if(i == 6){
		Tr_TS_22Stable <- Tr_TS_22Stable[-1]
		Tr_TS_22Flux <- Tr_TS_22Flux[-1]
		Tr_TS_25Stable <- Tr_TS_25Stable[-1]
		Tr_TS_25Flux <- Tr_TS_25Flux[-1]
	}
}
	Ab_TS_22Stable[which(Ab_TS_22Stable==0)] <- NA
	Ab_TS_22Flux[which(Ab_TS_22Flux==0)] <- NA
	Ab_TS_25Stable[which(Ab_TS_25Stable==0)] <- NA
	Ab_TS_25Flux[which(Ab_TS_25Flux==0)] <- NA

	Tr_TS_22Stable[which(Tr_TS_22Stable==0)] <- NA
	Tr_TS_22Flux[which(Tr_TS_22Flux==0)] <- NA
	Tr_TS_25Stable[which(Tr_TS_25Stable==0)] <- NA
	Tr_TS_25Flux[which(Tr_TS_25Flux==0)] <- NA

## FIND EMBEDDING DIMENSION

#Calculate optimal E
maxE<-12 #Maximum E to test
#Matrix for storing output
Emat<-matrix(nrow=maxE-1, ncol=8); colnames(Emat)<-c("22Stable", "22Flux", "25Stable", "25Flux","TR22Stable", "TR22Flux", "TR25Stable", "TR25Flux")
#Loop over potential E values and calculate predictive ability
#of each process for its own dynamics
for(E in 2:maxE) {
#Uses defaults of looking forward one prediction step (predstep)
#And using time lag intervals of one time step (tau)
	Emat[E-1,"22Stable"]<-SSR_pred_boot(A=Ab_TS_22Stable, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"22Flux"]<-SSR_pred_boot(A=Ab_TS_22Flux, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"25Stable"]<-SSR_pred_boot(A=Ab_TS_25Stable, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"25Flux"]<-SSR_pred_boot(A=Ab_TS_25Flux, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"TR22Stable"]<-SSR_pred_boot(A=Tr_TS_22Stable, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"TR22Flux"]<-SSR_pred_boot(A=Tr_TS_22Flux, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"TR25Stable"]<-SSR_pred_boot(A=Tr_TS_25Stable, E=E, predstep=1, tau=1)$rho
	Emat[E-1,"TR25Flux"]<-SSR_pred_boot(A=Tr_TS_25Flux, E=E, predstep=1, tau=1)$rho
}

#Look at plots to find E for each process at which
#predictive ability rho is maximized
matplot(2:maxE, Emat, type="l", col=1:8, lty=1:2,xlab="E", ylab="rho", lwd=2)
legend("bottomleft", c("22Stable", "22Flux", "25Stable", "25Flux","TR22Stable", "TR22Flux", "TR25Stable", "TR25Flux"), lty=1:2, col=1:8, lwd=2, bty="n") 
## Not great FOR DENSITY, BETTER FOR TRAITS


## Say we choose, for simplicity:
E_22Stable <- 5
E_22Flux <- 3
E_25Stable <- 11
E_25Flux <- 12

E_TR22Stable <- 5
E_TR22Flux <- 4
E_TR25Stable <- 2
E_TR25Flux <- 2


E_max <- cbind(as.data.frame(Emat),E=seq(2,12,1)) %>%
	pivot_longer(1:8, names_to = "Treatment", values_to = "value") %>%
	group_by(Treatment) %>%
	filter(value==max(value))

## 22Stable
CCM_boot_A_22S<-CCM_boot(Ab_TS_22Stable, Tr_TS_22Stable, E_max$E[which(E_max$Treatment=="22Stable")], tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B_22S<-CCM_boot(Tr_TS_22Stable, Ab_TS_22Stable, E_max$E[which(E_max$Treatment=="TR22Stable")], tau=1, iterations=800)
## 22Fluct
CCM_boot_A_22F<-CCM_boot(Ab_TS_22Flux, Tr_TS_22Flux, E_max$E[which(E_max$Treatment=="22Flux")], tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B_22F<-CCM_boot(Tr_TS_22Flux, Ab_TS_22Flux, E_max$E[which(E_max$Treatment=="TR22Flux")], tau=1, iterations=800)
# 25Stable
CCM_boot_A_25S<-CCM_boot(Ab_TS_25Stable, Tr_TS_25Stable,E_max$E[which(E_max$Treatment=="25Stable")], tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B_25S<-CCM_boot(Tr_TS_25Stable, Ab_TS_25Stable, E_max$E[which(E_max$Treatment=="TR25Stable")], tau=1, iterations=800)
# 25Fluct
CCM_boot_A_25F<-CCM_boot(Ab_TS_25Flux, Tr_TS_25Flux, E_max$E[which(E_max$Treatment=="25Flux")], tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B_25F<-CCM_boot(Tr_TS_25Flux, Ab_TS_25Flux, E_max$E[which(E_max$Treatment=="TR25Flux")], tau=1, iterations=800)



## Extract values
Stabl22_N <- data.frame(rho=CCM_boot_A_22S$rho,Treatment=rep("Stabl22_N", length(CCM_boot_A_22S$rho)), Fluctuations=rep("Stable",length(CCM_boot_A_22S$rho)), Temperature=rep(22, length(CCM_boot_A_22S$rho)),Direction=rep("N", length(CCM_boot_A_22S$rho)))

Stabl22_M <- data.frame(rho=CCM_boot_B_22S$rho,Treatment=rep("Stabl22_M", length(CCM_boot_B_22S$rho)), Fluctuations=rep("Stable", length(CCM_boot_B_22S$rho)), Temperature=rep(22, length(CCM_boot_B_22S$rho)),Direction=rep("M", length(CCM_boot_B_22S$rho)))

Fluct22_N <- data.frame(rho=CCM_boot_A_22F$rho,Treatment=rep("Fluct22_N", length(CCM_boot_A_22F$rho)), Fluctuations=rep("Fluct", length(CCM_boot_A_22F$rho)), Temperature=rep(22, length(CCM_boot_A_22F$rho)),Direction=rep("N", length(CCM_boot_A_22F$rho)))

Fluct22_M <- data.frame(rho=CCM_boot_B_22F$rho,Treatment=rep("Fluct22_M", length(CCM_boot_B_22F$rho)), Fluctuations=rep("Fluct", length(CCM_boot_B_22F$rho)), Temperature=rep(22, length(CCM_boot_B_22F$rho)),Direction=rep("M", length(CCM_boot_B_22F$rho)))

Stabl25_N <- data.frame(rho=CCM_boot_A_25S$rho,Treatment=rep("Stabl25_N", length(CCM_boot_A_25S$rho)), Fluctuations=rep("Stable", length(CCM_boot_A_25S$rho)), Temperature=rep(25, length(CCM_boot_A_25S$rho)),Direction=rep("N", length(CCM_boot_A_25S$rho)))

Stabl25_M <- data.frame(rho=CCM_boot_B_25S$rho,Treatment=rep("Stabl25_M", length(CCM_boot_B_25S$rho)), Fluctuations=rep("Stable", length(CCM_boot_B_25S$rho)), Temperature=rep(25, length(CCM_boot_B_25S$rho)),Direction=rep("M", length(CCM_boot_B_25S$rho)))

Fluct25_N <- data.frame(rho=CCM_boot_A_25F$rho,Treatment=rep("Fluct25_N", length(CCM_boot_A_25F$rho)), Fluctuations=rep("Fluct", length(CCM_boot_A_25F$rho)), Temperature=rep(25, length(CCM_boot_A_25F$rho)),Direction=rep("N", length(CCM_boot_A_25F$rho)))

Fluct25_M <- data.frame(rho=CCM_boot_B_25F$rho,Treatment=rep("Fluct25_M", length(CCM_boot_B_25F$rho)), Fluctuations=rep("Fluct", length(CCM_boot_B_25F$rho)), Temperature=rep(25, length(CCM_boot_B_25F$rho)),Direction=rep("M", length(CCM_boot_B_25F$rho)))

##Make dataset
CCM_data <- rbind(Stabl22_N,Stabl22_M,Fluct22_N,Fluct22_M,Stabl25_N,Stabl25_M,Fluct25_N,Fluct25_M)
##Define factors
CCM_data$Treatment <- factor(CCM_data$Treatment,levels=c("Stabl22_N","Stabl22_M","Fluct22_N","Fluct22_M","Stabl25_N","Stabl25_M","Fluct25_N","Fluct25_M"))


CCM_N <- CCM_data %>%
	filter(Direction=="N") 
CCM_M <- CCM_data %>%
	filter(Direction=="M") 


plot(rho~Treatment, data=CCM_data)
plot(rho~Treatment, data=CCM_N)


## Change order of levels in Fluctuations
CCM_N$Fluctuations <- factor(CCM_N$Fluctuations, levels=c("Stable", "Fluct"))
CCM_M$Fluctuations <- factor(CCM_M$Fluctuations, levels=c("Stable", "Fluct"))

summary(lm(rho~as.factor(Temperature)*Fluctuations,data=CCM_N)) 
summary(lm(rho~as.factor(Temperature)*Fluctuations,data=CCM_M)) 	
	
library("gridExtra")	
	
g_N <- CCM_N %>%
		ggplot(aes(Treatment,rho, col=factor(Temperature))) +
			geom_boxplot(outlier.colour="white",aes(linetype=Fluctuations)) +
			geom_jitter(aes(shape=Fluctuations), position=position_jitter(width=0.3)) +
			scale_shape_manual(values=c(16, 1)) +
			scale_colour_manual(values=alpha(c("#0571B0","#CA0020"),0.6)) +
		 	scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.2)) +
		 	labs(x="Treatment", y="Cross Mapping Skill ( )") +
		  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
			  	axis.title=element_text(size=16),
			  	axis.line = element_line(color='black', size=0.8),
			    plot.background = element_blank(),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    panel.border = element_blank(),
			    axis.ticks.length=unit(-.15, "cm"),
			    legend.position="none",
			    axis.text.x = element_blank(),
			    axis.text.y = element_text(size=14, angle=0, margin = unit(c(0,0.3,0,0), "cm"))
			    #axis.title.x=element_blank(vjust=-0.2)	
			)


g_M <- CCM_M %>%
		ggplot(aes(Treatment,rho, col=factor(Temperature))) +
			geom_boxplot(outlier.colour="white",aes(linetype=Fluctuations)) +
			geom_jitter(aes(shape=Fluctuations), position=position_jitter(width=0.3)) +
			scale_shape_manual(values=c(16, 1)) +
			scale_colour_manual(values=alpha(c("#0571B0","#CA0020"),0.6)) +
		 	scale_y_continuous(limits=c(0, 1),breaks=seq(0,1,0.2)) +
		 	labs(x="Treatment", y="Cross Mapping Skill ( )") +
		  	theme(panel.background = element_rect(fill = "transparent",colour =NA),
			  	axis.title=element_text(size=16),
			  	axis.line = element_line(color='black', size=0.8),
			    plot.background = element_blank(),
			    panel.grid.major = element_blank(),
			    panel.grid.minor = element_blank(),
			    panel.border = element_blank(),
			    axis.ticks.length=unit(-.15, "cm"),
			    legend.position="none",
			    axis.text.x = element_blank(),
			    axis.text.y = element_text(size=14, angle=0, margin = unit(c(0,0.3,0,0), "cm"))
			    #axis.title.x=element_blank(vjust=-0.2)	
			)


#pdf("~/Desktop/JP/Papers_in_progress/JP_Tetra_Pheno_Temps/Figures/Fig_4/_Fig_4.pdf", width = 4, height = 12) # Open a new pdf file
grid.arrange(g_N,g_M, nrow=1) # Write the grid.arrange in the file
#dev.off()

## CCM Appendix
### FROM CCM MANUAL
#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)

pdf("~/Desktop/JP/Papers_in_progress/JP_Tetra_Pheno_Temps/Figures/Appendix/Fig_S2/Fig_S2.pdf", width = 8,height = 8)
par(mfrow=c(2,2))
## 22STABLE
# Does DENSITY "cause" TRAITS?
CCM_boot_A<-CCM_boot(Ab_TS_22Stable, Tr_TS_22Stable, E_22Stable, tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B<-CCM_boot(Tr_TS_22Stable, Ab_TS_22Stable, E_TR22Stable, tau=1, iterations=800)
#Test for significant causal signal
#See R function for details
(CCM_significance_test<-ccmtest(CCM_boot_A,
CCM_boot_B))

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
xlab="L", ylab="rho", main="22 Stable")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col=2)


## 22FLUX
# Does DENSITY "cause" TRAITS?
CCM_boot_A<-CCM_boot(Ab_TS_22Flux, Tr_TS_22Flux, E_22Flux, tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B<-CCM_boot(Tr_TS_22Flux, Ab_TS_22Flux, E_TR22Flux, tau=1, iterations=800)
#Test for significant causal signal
#See R function for details
(CCM_significance_test<-ccmtest(CCM_boot_A,
CCM_boot_B))

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
xlab="L", ylab="rho", main="22 Fluct")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col=2)


## 25STABLE
# Does DENSITY "cause" TRAITS?
CCM_boot_A<-CCM_boot(Ab_TS_25Stable, Tr_TS_25Stable, E_25Stable, tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B<-CCM_boot(Tr_TS_25Stable, Ab_TS_25Stable, E_TR25Stable, tau=1, iterations=800)
#Test for significant causal signal
#See R function for details
(CCM_significance_test<-ccmtest(CCM_boot_A,
CCM_boot_B))

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
xlab="L", ylab="rho", main="25 Stable")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col=2)


## 25FLUX
# Does DENSITY "cause" TRAITS?
CCM_boot_A<-CCM_boot(Ab_TS_25Flux, Tr_TS_25Flux, E_25Flux, tau=1, iterations=800)
# Does TRAITS "cause" DENSITY?
CCM_boot_B<-CCM_boot(Tr_TS_25Flux, Ab_TS_25Flux, E_TR25Flux, tau=1, iterations=800)
#Test for significant causal signal
#See R function for details
(CCM_significance_test<-ccmtest(CCM_boot_A,
CCM_boot_B))

#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot DENSITY "cause" TRAITS
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
xlab="L", ylab="rho", main="25 Fluct")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
CCM_boot_A$rho+CCM_boot_A$sdevrho),
lty=3, col=1)
#Plot TRAITS "cause" DENSITY?
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
CCM_boot_B$rho+CCM_boot_B$sdevrho),
lty=3, col=2)
#legend("topleft",
#c("Density causes Traits", "Traits causes density"),
#lty=c(1,2), col=c(1,2), lwd=2, bty="n")

dev.off()

    



## THE END
