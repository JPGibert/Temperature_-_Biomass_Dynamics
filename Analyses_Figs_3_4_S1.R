### Clear cache
rm(list=ls())

#
### 1) Load packages -----
packages = c("tidyverse", "mgcv", "data.table", "multcompView", "ggpubr") 

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#
### 2) Data analysis -----

## Set working directory (this directory must contain downloaded raw data files)
data_dir<-"~/"


## Import raw data | Calculate individual masses | Set treatment factor levels
read_func<-function(file){
  data<-fread(paste0(data_dir, file)) %>%
    select(Day, Treatment, Jar, Volume, Density)
  if(file=="full_dat_2.csv"){data[-1]} else {data}
}

rawdata<-lapply(c("full_dat_1.csv", "full_dat_2.csv", "full_dat_3.csv", "full_dat_4.csv"), read_func) %>%
  rbindlist(use.names=T, fill=T) %>%
  drop_na() %>%
  separate(Treatment, c("Temperature", "Fluctuations"), sep=fixed("_"), remove=F) %>%
  mutate(Day=Day-1, Mass_g=Volume*10^-12,
         Treatment=factor(Treatment, levels=c("22_Stable", "22_Flux", "25_Stable", "25_Flux")),
         Temperature=factor(Temperature, levels=c("22", "25")),
         Fluctuations=factor(Fluctuations, levels=c("Stable", "Flux")))


## Summarize data within replicates (Jars)
# Part 1: summarize sample data
jar_data_1<-rawdata %>%
  group_by(Day, Treatment, Temperature, Fluctuations, Jar, Density) %>%
  summarize(Count=n(), Biomass_g_Sample=sum(Mass_g, na.rm=T), Volume=mean(Volume, na.rm=T), Mass_g=mean(Mass_g, na.rm=T)) %>%
  ungroup() %>%
  mutate(BiomassDensity_g=Mass_g*Density) 

# Part 2: perform regression to quantify differences between the number of individuals sampled and actual density estimates
fit<-lm(Count ~ 0 + Density, data=filter(jar_data_1, Day!=0))
summary(fit)

# Part 3: use regression slope to linearly transform Biomass according to sample-density differences
jar_data_2<-jar_data_1 %>%
  mutate(Count_2=ifelse(Day==0, 10, Count/coef(fit)[[1]]), Biomass_g=ifelse(Day==0, BiomassDensity_g, Biomass_g_Sample/coef(fit)[[1]]))


## Calculate population intrinsic growth rates
linear_r_estimates<-jar_data_2 %>%
  select(Day, Treatment, Temperature, Fluctuations, Jar, Count_2) %>%
  filter(Day %in% c(0,1)) %>%
  pivot_wider(names_from=Day, names_prefix="Day", values_from=Count_2) %>%
  mutate(rC=log(Day1/Day0))


## Calculate contributions of Density (d ln(N)/dt) and Mass (d ln(m)/dt) to changes in Biomass over time
jar_data<-jar_data_2 %>%
  mutate(across(c(Volume, Mass_g, Density, Count_2, Biomass_g, BiomassDensity_g), ~log(.x))) %>%
  group_by(Treatment, Temperature, Fluctuations, Jar) %>%
  arrange(Day) %>%
  mutate(dM=(Mass_g-lag(Mass_g)),
         dN=(Density-lag(Density)),
         dC=(Count_2-lag(Count_2)),
         dB=(Biomass_g-lag(Biomass_g)),
         dBD=(BiomassDensity_g-lag(BiomassDensity_g))) %>%
  ungroup %>%
  mutate(dB_dM=dB-dM, dB_dC=dB-dC, dBD_dM=dBD-dM, dBD_dN=dBD-dN)

jar_diff_data<-jar_data %>%
  select(Day:Jar, dM:dBD_dN) %>%
  pivot_longer(dM:dBD_dN, names_to="Var", values_to="Val")  %>%
  mutate(Var=factor(Var, levels=c("dN", "dC", "dM", "dB", "dBD", "dB_dC", "dB_dM", "dBD_dN", "dBD_dM")))

# Fit GAMMs
gam_diff_data<-filter(jar_diff_data, Var %in% c("dN", "dM")) %>%
  mutate(Treatment=factor(Treatment, levels=c('22_Stable', '22_Flux', '25_Stable', '25_Flux')))

gam_fit_diffs<-gamm(Val ~ s(Day) + Temperature + Var,
                    random=list(Jar=~1),
                    data = gam_diff_data, correlation = corARMA(form = ~Day|Jar/Treatment/Var, p = 1),
                    na.action=na.omit)
summary(gam_fit_diffs$gam)

gam_fits<-expand.grid(Day=seq(0, 14, .1), Temperature=unique(gam_diff_data$Temperature), Fluctuations=unique(gam_diff_data$Fluctuations), Var=unique(gam_diff_data$Var)) %>%
  mutate(Treatment=paste(Temperature, Fluctuations, sep="_"),
         Val=predict(gam_fit_diffs$gam, .)) %>%
  mutate(Treatment=factor(Treatment, levels=c('22_Stable', '22_Flux', '25_Stable', '25_Flux')))


#
### 3) Plot results -----

## Set base graphics theme
plot_theme<-theme(plot.background=element_blank(), panel.background=element_blank(), panel.border=element_rect(color="black", fill=NA), panel.grid=element_blank(),
                  axis.text=element_text(size=12), axis.title=element_text(size=14))
color_palette<-c("#e41a1c", "#ff7f00", "#4daf4a", "#377eb8", "#a65628", "#984ea3", "#f781bf", "#984ea3", "#f781bf")


## Figure 3: Boxplots of dynamics summaries -----

# Peak Biomass
temp_df<-filter(jar_data, Day %in% c(3,4,5)) %>%
  mutate(id=row_number())

fit_pB<-lm(exp(Biomass_g)~Temperature*Fluctuations, data=temp_df, na.action=na.omit)
summary(fit_pB)

fit_b<-aov(exp(Biomass_g)~Treatment, data=temp_df, na.action=na.omit)
summary(fit_b)
fit_tuk<-TukeyHSD(fit_b)
tuk_grps<-multcompLetters(fit_tuk$Treatment[,4])
tuk_grps_df<-data.frame(tuk_grp=tuk_grps$Letters) %>% rownames_to_column(var="Treatment")

tuk_grps_g_df<-temp_df %>%
  group_by(Treatment, Temperature, Fluctuations) %>%
  summarize(y_val=max(exp(Biomass_g))+0.000005) %>%
  left_join(tuk_grps_df)

g1<-
  ggplot() +
  geom_boxplot(data=temp_df, aes(Treatment, exp(Biomass_g), color=Temperature, fill=Temperature, linetype=Fluctuations), alpha=0.5) +
  geom_text(data=tuk_grps_g_df, aes(x=Treatment, y=y_val, group=, label=tuk_grp), size=5) +
  scale_color_manual(values=color_palette[c(4,1)]) +
  scale_fill_manual(values=color_palette[c(4,1)]) +
  labs(x="Treatment", y="Peak biomass"~(g~mL^-1)) +
  plot_theme +
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        aspect.ratio=1.2)

# Intrinsic Growth Rate
fit_r<-lm(rC~Temperature*Fluctuations, data=linear_r_estimates, na.action=na.omit)
summary(fit_r)

fit_r<-aov(rC~Treatment, data=linear_r_estimates, na.action=na.omit)
summary(fit_r)
fit_tuk<-TukeyHSD(fit_r)
tuk_grps<-multcompLetters(fit_tuk$Treatment[,4])
tuk_grps_df<-data.frame(tuk_grp=tuk_grps$Letters) %>% rownames_to_column(var="Treatment")

tuk_grps_g_df<-linear_r_estimates %>%
  group_by(Treatment, Temperature, Fluctuations) %>%
  summarize(y_val=max(rC)+0.1) %>%
  left_join(tuk_grps_df)

g2<-
  ggplot() +
  geom_boxplot(data=linear_r_estimates, aes(Treatment, rC, color=Temperature, fill=Temperature, linetype=Fluctuations), alpha=0.5) +
  geom_text(data=tuk_grps_g_df, aes(x=Treatment, y=y_val, group=, label=tuk_grp), size=5) +
  scale_color_manual(values=color_palette[c(4,1)]) +
  scale_fill_manual(values=color_palette[c(4,1)]) +
  labs(x="Treatment", y=Intrinsic~growth~rate~(ind~ind^-1~mL^-1)) + 
  plot_theme +
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        aspect.ratio=1.2)

# Carrying Capacity
temp_df<-filter(jar_data, Day>=13)

fit_k<-lm(exp(Density)~Temperature*Fluctuations, data=temp_df, na.action=na.omit)
summary(fit_k)

fit_k<-aov(exp(Density)~Treatment, data=temp_df, na.action=na.omit)
summary(fit_k)
fit_tuk<-TukeyHSD(fit_k)
tuk_grps<-multcompLetters(fit_tuk$Treatment[,4])
tuk_grps_df<-data.frame(tuk_grp=tuk_grps$Letters) %>% rownames_to_column(var="Treatment")

tuk_grps_g_df<-temp_df %>%
  group_by(Treatment, Temperature, Fluctuations) %>%
  summarize(y_val=max(exp(Density))+150) %>%
  left_join(tuk_grps_df)

g3<-
  ggplot() +
  geom_boxplot(data=temp_df, aes(Treatment, exp(Density), color=Temperature, fill=Temperature, linetype=Fluctuations), alpha=0.5) +
  geom_text(data=tuk_grps_g_df, aes(x=Treatment, y=y_val, group=, label=tuk_grp), size=5) +
  scale_color_manual(values=color_palette[c(4,1)]) +
  scale_fill_manual(values=color_palette[c(4,1)]) +
  labs(x="Treatment", y='Carrying capacity'~(ind~mL^-1)) +
  plot_theme +
  theme(legend.position="none", 
        axis.text.x=element_blank(), axis.title.x=element_blank(),
        aspect.ratio=1.2)

# Export combined figure panels
svg(file=paste0(data_dir, "f3.svg"), width=10, height=4, bg=NA) 
ggarrange(g1, g2, g3, 
          nrow=1, ncol=3, align="hv")
dev.off()


## Figure 4: Contributions of Density and Mass to Biomass dynamics -----

# GAMM fits to Density time series
svg(file=paste0(data_dir, "f4ab.svg"), bg=NA, width=3, height=3.5)
ggplot() +
  geom_hline(yintercept=0, size=0.1, linetype=1) +
  geom_line(data=filter(jar_diff_data, Var %in% c("dN", "dM")), aes(Day, Val, color=Var, linetype=Fluctuations, group=interaction(Jar, Var, Fluctuations)), alpha=0.2) +
  geom_line(data=gam_fits, aes(Day, Val, color=Var), size=1.2, lineend="round") +
  scale_color_manual(values=color_palette[c(3,6)], labels=c("dN/dt", "dM/dt", "dB/dt")) +
  labs(y="Value", color="", linetype="", size="") +
  scale_x_continuous(limits=c(1,14), breaks=seq(0,14,2)) +
  ylim(-1,4.5) +
  facet_wrap(~Temperature, ncol=1) +
  plot_theme +
  theme(legend.key=element_blank(), 
        legend.position="none",
        axis.title=element_blank(), aspect.ratio=1/2)
dev.off()

# Density distributions for Day<=2
svg(file=paste0(data_dir, "f4cd.svg"), bg=NA, width=3, height=3.5)
ggplot() +
  geom_density(data=filter(jar_diff_data, Var %in% c("dN", "dM"), Day<=2), aes(Val, ..density.., fill=Var), color=NA, alpha=0.65) + 
  geom_vline(xintercept=0, size=0.1) +
  scale_color_manual(values=color_palette[c(3,6)]) +
  scale_fill_manual(values=color_palette[c(3,6)]) +
  xlim(-1,4.5) +
  coord_flip() +
  facet_wrap(~Temperature, ncol=1) +
  plot_theme +
  theme(legend.position="none", axis.title=element_blank(), aspect.ratio=2)
dev.off()

# Density distributions for Day>2
svg(file=paste0(data_dir, "4ef.svg"), bg=NA, width=3, height=3.5)
ggplot() +
  geom_density(data=filter(jar_diff_data, Var %in% c("dN", "dM"), Day>2), aes(Val, ..density.., fill=Var), color=NA, alpha=0.65) +
  geom_vline(xintercept=0, size=0.1) +
  scale_color_manual(values=color_palette[c(3,6)]) +
  scale_fill_manual(values=color_palette[c(3,6)]) +
  coord_flip() +
  facet_wrap(~Temperature, ncol=1) +
  plot_theme +
  theme(legend.position="none", axis.title=element_blank(), aspect.ratio=2)
dev.off()

# Density boxplots across treatments for Day<=2 and Day>2
diff_aov_data<-
  jar_diff_data %>%
  mutate(Period=ifelse(Day<=2, "<=2", ">2"), Treatment=paste(Temperature, Var, Period, sep="_")) %>%
  filter(Var %in% c("dN", "dM"))

fit_diff<-aov(Val~Treatment, data=diff_aov_data, na.action=na.omit)
summary(fit_diff)
fit_tuk<-TukeyHSD(fit_diff)
fit_tuk

tuk_grps<-multcompLetters(fit_tuk$Treatment[,4])
tuk_grps_df<-data.frame(tuk_grp=tuk_grps$Letters) %>% rownames_to_column(var="Treatment")
tuk_grps_df

temp_data<-jar_diff_data %>%
  filter(Var %in% c("dN", "dM")) %>%
  mutate(Period=ifelse(Day<=2, "<=2", ">2"), Treatment=paste(Temperature, Var, Period, sep="_")) %>%
  arrange(Period, Var, Temperature)
temp_data$Treatment<-factor(temp_data$Treatment, levels=unique(temp_data$Treatment))

tuk_grps_g_df<-temp_data %>%
  group_by(Treatment) %>%
  summarize(y_val=max(Val, na.rm=T)+0.3) %>%
  left_join(tuk_grps_df)
  

svg(file=paste0(data_dir, "f4g.svg"), bg=NA, width=3, height=3.1)
ggplot() +
  geom_hline(yintercept=0, size=0.1) +
  geom_boxplot(data=temp_data, aes(Treatment, Val, color=Var, fill=Var), alpha=0.5) + 
  geom_text(data=tuk_grps_g_df, aes(x=Treatment, y=y_val, label=tuk_grp), size=5) +
  scale_color_manual(values=color_palette[c(3,6)]) +
  scale_fill_manual(values=color_palette[c(3,6)]) +
  ylim(-1,4.6) +
  plot_theme +
  theme(axis.title=element_blank(), axis.text.x=element_blank(),  
        strip.text=element_blank(), strip.background=element_blank(),
        legend.key=element_blank(), 
        legend.position="none")
dev.off()


## Figure S1: Sample counts vs. Density estimates -----
temp_data<-rawdata %>%
  filter(Day!=0) %>%
  group_by(Day, Treatment, Temperature, Fluctuations, Jar, Density) %>%
  summarize(Count=n(), Biomass_g=sum(Mass_g, na.rm=T), Mass_g=mean(Mass_g, na.rm=T)) %>%
  ungroup()

fit<-lm(Count ~ 0 + Density, data=temp_data)
summary(fit)

temp_data_2<-temp_data %>%
  mutate(Count_2=Count/coef(fit)[[1]], Biomass_g_2=Biomass_g/coef(fit)[[1]])

fit_2<-lm(Count_2 ~ 0 + Density, data=temp_data_2)
summary(fit_2)

svg(file=paste0(data_dir, "fs1.svg"), width=5, height=5, bg=NA) 
ggplot(data=temp_data_2) +
  geom_point(aes(Density, Count), color="black", alpha=0.5) +
  geom_point(aes(Density, Count_2), color="blue", alpha=0.5) +
  geom_abline(intercept=0, slope=coef(fit)[[1]], color="black") +
  geom_abline(intercept=0, slope=1, color="blue") +
  lims(x=c(0, 12000), y=c(0, 12000)) +
  theme(plot.background=element_blank(), panel.background=element_blank(), panel.border=element_rect(color="black", fill=NA), 
        axis.text=element_text(size=12), axis.title=element_text(size=14), aspect.ratio=1)
dev.off()


#
### END -----