###-----------------------------------------------------------------------------
# Title: Seasonal effects of farmer-managed livestock grazing exclusions on bird 
#        communities in Burkina Faso
# Author: Ian Quintas
# Date: 23.05.2024
# Journal: Ecological Applications
# R version: 4.2.1
# Revised by Gabriel Marcacci (15.05.2025) / Pius Korner (20.05.2025)
###-----------------------------------------------------------------------------
setwd("C:/Users/gm/OneDrive - Vogelwarte/Documents/Burkina/Manuscript/Ecological Applications/zenodo")

# Theme 
theme_newtree <- function () { 
  theme(panel.border = element_rect(fill = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background  = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 18, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 12, colour = "black"),
        legend.key = element_rect(fill = "white")
        
  )
}

# Libraries
library(tidyverse) 
library(dplyr) 
library(ggplot2) 
library(gridExtra) 
library(rstanarm)
library(broom.mixed)  
library(spdep)
library(cowplot)
library(geosphere)

#### Effect of grazing exclusion on bird richness-------------------------------

### Load and check data
by_survey <- read.csv("by_survey.csv", header = T, sep = ";") 
str(by_survey)
# habitat and season must be factors
by_survey$habitat <- factor(by_survey$habitat)
by_survey$season <- factor(by_survey$season)

### Overview of the data 

# mean bird rich per survey
mean(by_survey$bird_rich) #22.78 bird species per survey
sd(by_survey$bird_rich) #sd of bird richness per survey = 4.24

# mean bird rich per survey depending on habitat
tapply(by_survey$bird_rich,by_survey$habitat, mean)

# mean bird rich per survey depending on habitat and season
tapply(by_survey$bird_rich, list(by_survey$habitat, by_survey$season), mean) 

# SD
tapply(by_survey$bird_rich, list(by_survey$habitat, by_survey$season), sd)



##### Bird richness model

# check correlation between variables
by_survey %>%
  dplyr::select(tot_strata, plant_rich_tot, Shannon_vertical, herb_cov, 
                herb_height, diff.max.NDVI.500, nb_houses) %>%
  cor()


### model
set.seed(123)

# with only habitat * season
bird_rich <- stan_glmer(bird_rich ~ habitat * season + (1|landscape_id),
                        data = by_survey,
                        family = "poisson")

# full model
bird_rich_full <- stan_glmer(bird_rich ~ habitat * season + scale(age_enclos) + 
                               scale(diff.max.NDVI.500) + scale(herb_cov) + 
                               scale(tot_strata) + scale(nb_houses) +
                               scale(Shannon_vertical) + (1|landscape_id),
                             data = by_survey,
                             family = "poisson")

## check model

# trace plots
plot(bird_rich_full, plotfun = "trace") # good

# posterior predictive check
pp_check(bird_rich_full) # good

# Rhat & ESS
summary(bird_rich_full) # good
as.data.frame(summary(bird_rich_full, pars = c("alpha", "beta"), probs = c(0.05, 0.95)))
ranef(bird_rich_full) # random effect

## spatial autocorrelation
residuals <- residuals(bird_rich_full)
coords <- by_survey[, c("dec_long", "dec_lat")]
# create neighbor list
nb <- knn2nb(knearneigh(coords, k = 4))
# create spatial weights
lw <- nb2listw(nb, style = "W")
moran_test <- moran.test(residuals, lw)
print(moran_test) # p = 0.54 => no spatial autocorrelation

# plot residuals
coords$residuals <- residuals
ggplot(coords, aes(x = dec_long, y = dec_lat, color = residuals)) +
  geom_point(size = 5) +
  coord_cartesian(ylim = c(12.25, 12.75)) +
  scale_color_gradient2() +
  theme_minimal()


## summarize results
mod_rich_results <- as.matrix(bird_rich_full)
colnames(mod_rich_results)
mean_rich <- round(apply(mod_rich_results[,c(1:12,40)], 2, mean), 2)
lower_rich <- round(apply(mod_rich_results[,c(1:12,40)], 2, quantile, probs = 0.025), 2)
upper_rich <- round(apply(mod_rich_results[,c(1:12,40)], 2, quantile, probs = 0.975), 2)
paste0(mean_rich, " [", lower_rich, ", ", upper_rich, "]")


## plot estimates
tidy_model_full <- broom.mixed::tidy(bird_rich_full, effects = "fixed", conf.int = TRUE)
print(tidy_model_full)

# remove intercept for plotting
tidy_model_full_no_intercept <- tidy_model_full %>%
  filter(term != "(Intercept)") %>%
  mutate(significant = ifelse(conf.low > 0 | conf.high < 0, "*", "")) 

# plot
effect_plot_full_no_intercept <- ggplot(tidy_model_full_no_intercept, aes(x = term, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_text(aes(label = significant), nudge_x = 0.2, size = 5, vjust = -1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  labs(x = "Fixed Effects", y = "Estimate with 95% CI") +
  theme_minimal() +
  theme_newtree() +
  coord_flip()  

# display the plot
print(effect_plot_full_no_intercept)



### Compute predictions for bird richness

newdat <- expand.grid(habitat = factor(levels(by_survey$habitat), 
                                       levels = levels(by_survey$habitat)),
                      season = factor(levels(by_survey$season),
                                      levels = levels(by_survey$season)))

fitmat <- posterior_epred(bird_rich, newdata = newdat, re.form = NA)

newdat$fit <- apply(fitmat, 2, mean)
newdat$upr <- apply(fitmat, 2, quantile, probs = 0.975)
newdat$lwr <- apply(fitmat, 2, quantile, probs = 0.025)

## calculate median % change with CI

# exclos vs woody (dry)
quantile((fitmat[,1] - fitmat[,3]) / fitmat[,3], probs = 0.5) * 100 # 18.92
quantile((fitmat[,1] - fitmat[,3]) / fitmat[,3], probs = 0.975) * 100 # 32.85
quantile((fitmat[,1] - fitmat[,3]) / fitmat[,3], probs = 0.025) * 100 # 6.38

# exclos vs woody (wet)
quantile((fitmat[,4] - fitmat[,6]) / fitmat[,6], probs = 0.5) * 100 # 19.30
quantile((fitmat[,4] - fitmat[,6]) / fitmat[,6], probs = 0.975) * 100 # 33.92
quantile((fitmat[,4] - fitmat[,6]) / fitmat[,6], probs = 0.025) * 100 # 6.85

# exclos vs open (dry)
quantile((fitmat[,1] - fitmat[,2]) / fitmat[,2], probs = 0.5) * 100 # 29.28
quantile((fitmat[,1] - fitmat[,2]) / fitmat[,2], probs = 0.975) * 100 # 44.53
quantile((fitmat[,1] - fitmat[,2]) / fitmat[,2], probs = 0.025) * 100 # 15.88

# exclos vs open (wet)
quantile((fitmat[,4] - fitmat[,5]) / fitmat[,5], probs = 0.5) * 100 # 13.06
quantile((fitmat[,4] - fitmat[,5]) / fitmat[,5], probs = 0.975) * 100 # 26.52
quantile((fitmat[,4] - fitmat[,5]) / fitmat[,5], probs = 0.025) * 100 # 0.91


### Point range plot

# add levels to interaction between habitat and season with the levels we want
newdat <- newdat %>%
  mutate(interaction_var = factor(interaction(habitat, season, sep = "_"),
                                  levels = c("enclos_dry", "enclos_post rainy", "open_dry", 
                                             "open_post rainy", "woody_dry", "woody_post rainy")))
newdat$habitat <- as.character(newdat$habitat)
newdat$habitat[newdat$habitat == "enclos"] <- "exclos"
newdat$habitat <- factor(newdat$habitat, levels = c("exclos","open","woody"))


## Figure 2

ggplot(newdat) +  
  geom_pointrange(aes(x = habitat, y = fit, ymin = lwr, ymax = upr, col = interaction_var),
                  size = 10, fatten = 0.5, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("enclos_dry" = "chartreuse4", "enclos_post rainy" = "#A3E786", "open_dry" = "darkorange", 
                                "open_post rainy" = "#FFDED1", "woody_dry" = "dodgerblue2","woody_post rainy" = "#C9D9FF"),
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Habitat") +
  ylab("Bird richness") +
  theme_newtree() +
  guides(color = guide_legend(override.aes = list(size = .8),
                              title = "Habitat and Season")) +  
  theme(legend.position = "bottom",  
        legend.justification = "center",  
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(10,10,10,10),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 15))




##### Species-specific responses to grazing exclusion---------------------------

### Import data with > 20 occurrences and corresponding 0s (when the species is absent)
t.dat0 <- read.csv("t.dat0.csv", header = T, row.names = "X") 
t.dat0$season <- factor(t.dat0$season, levels = c("dry","post rainy"))
t.dat0$habitat <- factor(t.dat0$habitat, levels = c("open","woody","enclos"))
t.dat0$species_sci <- factor(t.dat0$species_sci)
t.sel <- levels(t.dat0$species_sci)

#remove control, t.dat0 should have 5055 rows
t.dat0 <- subset(t.dat0, landscape_id != "26C_07")

### one model with all species (!!! can take hours to run !!!)
#mod1 <- stan_glmer(cbind(Npres, Nabs.interv) ~ habitat * season + (habitat * season|species_sci) + 
#                     (1|landscape_id), 
#                   data = t.dat0, family = binomial, iter = 6000)

## check model
load("mod1.RData")
pp_check(mod1) # quite bad, maybe overdispersion and zero-inflation
plot(mod1, plotfun = "trace")
summary(mod1, probs = c(0.025,0.975))

# check for zero-inflation
pp <- posterior_predict(mod1)
pp0 <- apply(pp, 1, FUN = function(x) sum(x == 0))
hist(pp0, xlim = range(c(pp0, sum(t.dat0$Npres == 0))))
abline(v = sum(t.dat0$Npres == 0), col = "red")     # quite some zero-inflation.
par(mfrow = c(4, 4), mar = c(1, 1, 1, 1))
plot(table(t.dat0$Npres), xlim = c(0,30),xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for(i in 1:15) {
  plot(table(pp[i,]), xlim = c(0, 30),xaxt = "n", yaxt = "n", xlab = "", ylab = "")
}
# the first plot = raw data: more zeros (first column) compared to the 14 replicate data shown

t.dat0$ol <- 1:nrow(t.dat0)   # add ol = observation-level random factor

# refit the model with an observation-level random intercept to account for 
# overdispersion and zero-inflation

# Sys.time()
#mod2 <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season + (habitat * season|species_sci) +
#                   (1|landscape_id) + (1|ol),
#                data = t.dat0, family = binomial, iter = 6000)
#save(mod2, file = "mod2.RData")
load("mod2.RData")

pp_check(mod2)   # now looks perfect!
summary(mod2)    # all Rhat are good


## summarize results
mod_species_results <- as.matrix(mod2)
colnames(mod_species_results) # matrix too large need to find position of the other random factors
which(colnames(mod_species_results) %in% c("Sigma[landscape_id:(Intercept),(Intercept)]",
                                           "Sigma[ol:(Intercept),(Intercept)]"))
colnames(mod_species_results[,c(1:6, 5419, 5441)]) # correct:5419 = SD landscape_id, 5441 = SD observation level random factor
mean_occ <- round(apply(mod_species_results[,c(1:6, 5419, 5441)], 2, mean), 2)    
lower_occ <- round(apply(mod_species_results[,c(1:6, 5419, 5441)], 2, quantile, probs = 0.025), 2)
upper_occ <- round(apply(mod_species_results[,c(1:6, 5419, 5441)], 2, quantile, probs = 0.975), 2)
paste0(mean_occ, " [", lower_occ, ", ", upper_occ, "]")


### Make predictions

# expand the dataset to make predictions
newdat2 <- expand.grid(habitat = factor(levels(t.dat0$habitat), levels = levels(t.dat0$habitat)),
                       season = factor(levels(t.dat0$season), levels = levels(t.dat0$season)),
                       species_sci = factor(levels(t.dat0$species_sci), levels = levels(t.dat0$species_sci)))

# we can improve that by using the landscape_id level with the smallest effect:
tt <- ranef(mod2)$landscape_id
t.landscape_id_min <- rownames(tt)[abs(tt) == min(abs(tt))][1]  # (the [1] at the end only for the unlikely case that two levels are equally minimal)
newdat2$landscape_id <- t.landscape_id_min
# same with ol
tt <- ranef(mod2)$ol
t.ol.min <- names(tt)[abs(tt == min(abs(tt)))][1]
newdat2$ol <- t.ol.min


# make predictions
fitmat <- posterior_epred(mod2, newdata = newdat2)

# calculate the quantiles of the predictions:
newdat2$lwr <- apply(fitmat, 2, quantile, probs = 0.025)
newdat2$upr <- apply(fitmat, 2, quantile, probs = 0.975)
newdat2$fit <- apply(fitmat, 2, mean)
newdat2$x <- as.numeric(factor(newdat2$habitat))
newdat2$col <- as.numeric(factor(newdat2$species_sci))

# which change is significant?: we get that info out of fitmat:
str(fitmat)   # the 330 columns in fitmat correspond to the 330 rows in newdat2

# for 1 species we have to compare the difference of occurrence separately
# dry season exclos vs open // dry season exclos vs woody // wet season exclos vs open // wet season exclos vs woody

# first we start with the comparison between exclos and open
tt <- apply(fitmat, 1, FUN = function(x) tapply(x,
                                            paste(newdat2$species_sci, newdat2$season),
                                            FUN = function(y) y[levels(newdat2$habitat) == "enclos"]-
                                              y[levels(newdat2$habitat) == "open"]))
t.lwr <- apply(tt, 1, quantile, probs = 0.025)  # 2.5% quantile for the difference, per species and season
t.upr <- apply(tt, 1, quantile, probs = 0.975)  # 97.5% quantile for the difference, per species and season
t.sign <- rownames(tt)[0 <= t.lwr | 0 >= t.upr]
t.fit <- apply(tt, 1, quantile, probs = 0.5)

ci_sign_dry <- as.data.frame(cbind(t.lwr, t.upr))
ci_sign_dry <- ci_sign_dry[0 <= t.lwr | 0 >= t.upr,]


### Fig 3 (exclos vs open)

tt_df <- data.frame(
  species_season = rownames(tt),
  lower_bound = t.lwr,
  upper_bound = t.upr,
  median_effect = t.fit
)

# create species and season columns separately
split_columns <- str_split(tt_df$species_season, "\\s", n = 3)
# extracting species and season from the split result
tt_df$genre <- sapply(split_columns, function(x) x[[1]])
tt_df$sp <- sapply(split_columns, function(x) x[[2]])
tt_df$season <- sapply(split_columns, function(x) x[[3]])

tt_df$species <- paste(tt_df$genre,tt_df$sp, sep = " ")

tt_df$species_season <- factor(tt_df$species_season)
tt_df$species_season <- fct_rev(factor(tt_df$species_season)) # reorder species

# create color column
tt_df$color <- ifelse(tt_df$lower_bound > 0 & tt_df$upper_bound > 0 & grepl("dry", tt_df$species_season), "chartreuse4",
                      ifelse(tt_df$lower_bound > 0 & tt_df$upper_bound > 0 & grepl("post rainy", tt_df$species_season), "chartreuse3",
                             ifelse(tt_df$lower_bound < 0 & tt_df$upper_bound < 0 & grepl("dry", tt_df$species_season), "orangered3",
                                    ifelse(tt_df$lower_bound < 0 & tt_df$upper_bound < 0 & grepl("post rainy", tt_df$species_season), "tomato",
                                           ifelse(grepl("dry", tt_df$species_season), "darkgrey",
                                                  "lightgrey")))))



## plot

## better to make one plot per season
tt_df <- tt_df %>% mutate(species = as.factor(species))
tt_df_dry <- tt_df %>% filter(season == "dry")
tt_df_wet <- tt_df %>% filter(season == "post rainy")

# dry season
enclosVSopen_dry <- ggplot(tt_df_dry, aes(x = median_effect, y = fct_rev(species), color = color)) +
  geom_errorbarh(aes(xmin = lower_bound, xmax = upper_bound), height = 0) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of exclos-open (logit-scale)",  
       y = "",  
       title = "Dry season",
       fill = "Classification") +  
  theme_minimal() +
  theme_newtree() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        legend.position.inside = c(0.8, 0.9)) +  
  scale_x_continuous(limits = range(c(t.lwr, t.upr))) +  
  scale_color_manual(values = c("chartreuse4" = "chartreuse4",
                                "orangered3" = "orangered3", 
                                "darkgrey" = "darkgrey"),
                     labels = c("chartreuse4" = "Positive",
                                "orangered3" = "Negative",
                                "darkgrey" = "Unclear")) + 
  guides(color = guide_legend(title = "Classification"))


# wet season
enclosVSopen_wet <- ggplot(tt_df_wet, aes(x = median_effect, y = fct_rev(species), color = color)) +
  geom_errorbarh(aes(xmin = lower_bound, xmax = upper_bound), height = 0) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "",  
       y = "",  
       title = "Wet season",
       fill = "Classification") +  
  theme_minimal() +
  theme_newtree() +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        legend.position.inside = c(0.8, 0.9)) +  
  scale_x_continuous(limits = range(c(t.lwr, t.upr))) +  
  scale_color_manual(values = c("chartreuse3" = "chartreuse3",
                                "tomato" = "tomato",
                                "lightgrey" ="lightgrey"),
                     labels = c("chartreuse3" = "Positive",
                                "tomato" = "Negative",
                                "lightgrey" = "Unclear")) + 
  guides(color = guide_legend(title = "Classification"))




# align the two plots
combined_plot_enclosVSopen <- plot_grid(enclosVSopen_dry, enclosVSopen_wet, 
                                        align = "hv", 
                                        axis = "tb",  
                                        nrow = 1)  

# display the combined plot
print(combined_plot_enclosVSopen)

#ggsave("treeplot_enclosVS_open2.pdf", plot = combined_plot_enclosVSopen, width = 14, height = 10)


### Fig 4 (exclos vs woody)

tt <- apply(fitmat, 1, FUN = function(x) tapply(x,
                                            paste(newdat2$species_sci, newdat2$season),
                                            FUN = function(y) y[levels(newdat2$habitat) == "enclos"]-
                                              y[levels(newdat2$habitat) == "woody"]))

t.lwr <- apply(tt,1,quantile,probs=0.025)  # 2.5% quantile for the difference, per species and season
t.upr <- apply(tt,1,quantile,probs=0.975)  # 97.5% quantile for the difference, per species and season
t.sign <- rownames(tt)[0<=t.lwr | 0>=t.upr]
t.fit <- apply(tt,1,quantile,probs=0.5)

ci_sign_rainy <- as.data.frame(cbind(t.lwr,t.upr))
ci_sign_rainy<- ci_sign_rainy[0<=t.lwr | 0>=t.upr,]

tt_df <- data.frame(
  species_season = rownames(tt),
  lower_bound = t.lwr,
  upper_bound = t.upr,
  median_effect = t.fit
)

# create species and season columns separately
split_columns <- str_split(tt_df$species_season, "\\s", n = 3)

# extracting species and season from the split result
tt_df$genre <- sapply(split_columns, function(x) x[[1]])
tt_df$sp <- sapply(split_columns, function(x) x[[2]])
tt_df$season <- sapply(split_columns, function(x) x[[3]])
tt_df$species <- paste(tt_df$genre,tt_df$sp, sep = " ")
tt_df$species_season <- factor(tt_df$species_season)
tt_df$species_season <- fct_rev(factor(tt_df$species_season))
tt_df$species <- factor(tt_df$species)

# create color column
tt_df$color <- ifelse(tt_df$lower_bound > 0 & tt_df$upper_bound > 0 & grepl("dry", tt_df$species_season), "chartreuse4",
                      ifelse(tt_df$lower_bound > 0 & tt_df$upper_bound > 0 & grepl("post rainy", tt_df$species_season), "chartreuse3",
                             ifelse(tt_df$lower_bound < 0 & tt_df$upper_bound < 0 & grepl("dry", tt_df$species_season), "orangered3",
                                    ifelse(tt_df$lower_bound < 0 & tt_df$upper_bound < 0 & grepl("post rainy", tt_df$species_season), "tomato",
                                           ifelse(grepl("dry", tt_df$species_season), "darkgrey",
                                                  "lightgrey")))))

tt_df_dry <- tt_df %>% filter(season == "dry")
tt_df_wet <- tt_df %>% filter(season == "post rainy")


## plot

# dry season
enclosVSwoody_dry <- ggplot(tt_df_dry, aes(x = median_effect, y = fct_rev(species), color = color)) +
  geom_errorbarh(aes(xmin = lower_bound, xmax = upper_bound), height = 0) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect of exclos-woody (logit-scale)",  
       y = "",  
       title = "Dry season",
       fill = "Classification") + 
  theme_minimal() +
  theme_newtree() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        legend.position.inside = c(0.2, 0.9)) +  
  scale_x_continuous(limits = range(c(t.lwr, t.upr))) + 
  scale_color_manual(values = c("chartreuse4" = "chartreuse4",
                                "orangered3" = "orangered3", 
                                "darkgrey" = "darkgrey"),
                     labels = c("chartreuse4" = "Positive",
                                "orangered" = "Negative",
                                "darkgrey" = "Unclear")) + 
  guides(color = guide_legend(title = "Classification"))

# wet season
enclosVSwoody_wet <- ggplot(tt_df_wet, aes(x = median_effect, y = fct_rev(species), color = color)) +
  geom_errorbarh(aes(xmin = lower_bound, xmax = upper_bound), height = 0) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "",  
       y = "",  
       title = "Wet season",
       fill = "Classification") +  
  theme_minimal() +
  theme_newtree() +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(),
        legend.position.inside = c(0.8, 0.9)) +  
  scale_x_continuous(limits = range(c(t.lwr, t.upr))) +  
  scale_color_manual(values = c("chartreuse3" = "chartreuse3",
                                "tomato" = "tomato",
                                "lightgrey" ="lightgrey"),
                     labels = c("chartreuse3" = "Positive",
                                "tomato" = "Negative",
                                "lightgrey" = "Unclear")) + 
  guides(color = guide_legend(title = "Classification"))


# adjust plot margins if necessary
enclosVwoody_dry <- enclosVSwoody_dry + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
enclosVwoody_wet <- enclosVSwoody_wet + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# align the two plots
combined_plot_enclosVSwoody <- plot_grid(enclosVSwoody_dry, enclosVSwoody_wet, 
                                         align = "hv",  
                                         axis = "tb",  
                                         nrow = 1)  

# display the combined plot
print(combined_plot_enclosVSwoody)

#ggsave("treeplot_enclosVSwoody2.pdf", plot = combined_plot_enclosVSwoody, width = 14, height = 10)


##### Guild-specific responses to grazing exclusion-----------------------------

### load and prepare data
# ecological traits dataset for all recorded species
all_species <- read.csv("all_species.csv", header = T, sep = ";")

# merge the 2 df
t.dat0 <- t.dat0 %>%
  left_join(all_species, by = "species_sci")

ix <- !duplicated(t.dat0$species_sci)
table(t.dat0$migration[ix],t.dat0$habitat_pref[ix])
table(t.dat0$migration[ix],t.dat0$diet[ix])
table(t.dat0$habitat_pref[ix],t.dat0$diet[ix])

# too many diets so we merge some
t.dat0$diet2 <- ifelse(t.dat0$diet %in% c("Vertivore", "Nectarivore", "Omnivore"), "Other", t.dat0$diet)

t.dat0$diet2 <- factor(t.dat0$diet2, levels = c("Frugivore", "Granivore", "Invertivore", "Other"))
t.dat0$migration  <- factor(t.dat0$migration, levels = c("resident","migrant"))
t.dat0$habitat_pref <- factor(t.dat0$habitat_pref, levels = c("open","intermediate","forest"))
t.dat0$species_sci <- factor(t.dat0$species_sci) 


##### Models (1 model for each guild)

### Model diet

#mod_diet <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * diet2 +
#                         (1|landscape_id), 
#                       data = t.dat0, family = binomial, iter = 6000)
load("mod_diet.RData")
pp_check(mod_diet) # looks terrible
plot(mod1, plotfun = "trace")
summary(mod_diet, probs = c(0.025,0.975))
posterior_interval(mod_diet, prob = 0.95)

# refit model with ol random factor
#mod_diet2 <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * diet2 +
 #                        (1|landscape_id) + (1|ol), 
#                       data = t.dat0, family = binomial, iter = 6000)


#save(mod_diet2, file = "mod_diet2.RData")
load("mod_diet2.RData")

pp_check(mod_diet2)   # looks perfect
summary(mod_diet2)    # all Rhat are good


# summarise results
mod_diet_results <- as.matrix(mod_diet2)
colnames(mod_diet_results)
which(colnames(mod_diet_results) %in% c("Sigma[landscape_id:(Intercept),(Intercept)]",
                                           "Sigma[ol:(Intercept),(Intercept)]"))
Predictor_diet <- colnames(mod_diet_results[, c(1:24, 5107, 5108)])  # ol: 5107, landscape:5108
mean_diet <- round(apply(mod_diet_results[, c(1:24, 5107, 5108)], 2, mean), 2)
lower_diet <- round(apply(mod_diet_results[,c(1:24, 5107, 5108)], 2, quantile, probs = 0.025), 2)
upper_diet <- round(apply(mod_species_results[,c(1:24, 5107, 5108)], 2, quantile, probs = 0.975), 2)
CrI_diet <- c(paste0("[", lower_diet, ", ", upper_diet, "]"))
summary_diet <- as.data.frame(cbind(Predictor_diet, mean_diet, CrI_diet))
#write.csv(summary_diet, "summary_diet.csv")


### Model migration
#mod_migr <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * migration
#                       + (1|landscape_id), 
#                       data = t.dat0, family = binomial, iter = 6000)

load("mod_migr.RData")
pp_check(mod_migr) # looks terrible
summary(mod_migr, probs = c(0.025,0.975))
posterior_interval(mod_migr, prob = 0.95)


# refit model with ol random factor
#mod_migr2 <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * migration
#                       + (1|landscape_id) + (1|ol), 
#                       data = t.dat0, family = binomial, iter = 6000)

#save(mod_migr2, file = "mod_migr2.RData")
load("mod_migr2.RData")

pp_check(mod_migr2)   # looks perfect
summary(mod_migr2)    # all Rhat are good


# summarise results
mod_migr_results <- as.matrix(mod_migr2)
colnames(mod_migr_results)
which(colnames(mod_migr_results) %in% c("Sigma[landscape_id:(Intercept),(Intercept)]",
                                        "Sigma[ol:(Intercept),(Intercept)]"))
Predictor_migr <- colnames(mod_migr_results[, c(1:12, 5095, 5096)])  # ol: 5095, landscape:5096
mean_migr <- round(apply(mod_migr_results[, c(1:12, 5095, 5096)], 2, mean), 2)
lower_migr <- round(apply(mod_migr_results[,c(1:12, 5095, 5096)], 2, quantile, probs = 0.025), 2)
upper_migr <- round(apply(mod_species_results[,c(1:12, 5095, 5096)], 2, quantile, probs = 0.975), 2)
CrI_migr <- c(paste0("[", lower_migr, ", ", upper_migr, "]"))
summary_migr <- as.data.frame(cbind(Predictor_migr, mean_migr, CrI_migr))
#write.csv(summary_migr, "summary_migr.csv")


### Model habitat
# here, Euplecte franciscanus is an open habitat species it needs tall herb veg
# so almost only present in exclosures: we model with and without it
#With euplectes
#mod_habitat <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * habitat_pref
#                          + (1|landscape_id), 
#                          data = t.dat0, family = binomial, iter = 6000)
load("mod_habitat.RData")
pp_check(mod_habitat) # looks terrible
summary(mod_habitat, probs = c(0.025,0.975))
posterior_interval(mod_habitat, prob = 0.95)

# same without Euplectes (because it is an open habitat species that was only present in grazing exclusions due to specific habitat requirements)
# We removed Euplectes to see if it changes sth for the results BUT ONLY WITH HABITAT model
t.dat0_no_euplectes <- t.dat0[t.dat0$species_sci !="Euplectes franciscanus",]

#mod_habitat_no_euplectes <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * habitat_pref
#                                       + (1|landscape_id), 
#                                       data = t.dat0_no_euplectes, family = binomial, iter = 6000)
load("mod_habitat_no_euplectes.RData")
pp_check(mod_habitat_no_euplectes) # looks terrible
summary(mod_habitat_no_euplectes, probs = c(0.025,0.975))
posterior_interval(mod_habitat_no_euplectes, prob = 0.95)


# refit model with ol random factor
#mod_habitat_no_euplectes2 <- stan_glmer(cbind(Npres,Nabs.interv) ~ habitat * season * habitat_pref
#                                        + (1|landscape_id) + (1|ol), 
#                                        data = t.dat0_no_euplectes, family = binomial, iter = 6000)


#save(mod_habitat_no_euplectes2, file = "mod_habitat_no_euplectes2.RData")
load("mod_habitat_no_euplectes2.RData")

pp_check(mod_habitat_no_euplectes2)   # looks perfect
summary(mod_habitat_no_euplectes2)    # all R_hat are good


# summarise results
mod_hab_results <- as.matrix(mod_habitat_no_euplectes2)
colnames(mod_hab_results)
which(colnames(mod_hab_results) %in% c("Sigma[landscape_id:(Intercept),(Intercept)]",
                                        "Sigma[ol:(Intercept),(Intercept)]"))
Predictor_hab <- colnames(mod_hab_results[, c(1:18, 5014, 5015)])  # ol: 5014, landscape:5015
mean_hab <- round(apply(mod_hab_results[, c(1:18, 5014, 5015)], 2, mean), 2)
lower_hab <- round(apply(mod_hab_results[,c(1:18, 5014, 5015)], 2, quantile, probs = 0.025), 2)
upper_hab <- round(apply(mod_species_results[,c(1:18, 5014, 5015)], 2, quantile, probs = 0.975), 2)
CrI_hab <- c(paste0("[", lower_hab, ", ", upper_hab, "]"))
summary_hab <- as.data.frame(cbind(Predictor_hab, mean_hab, CrI_hab))
#write.csv(summary_hab, "summary_hab.csv")


### plots and predictions

newdat_diet <- expand.grid(habitat = factor(levels(t.dat0$habitat), levels = levels(t.dat0$habitat)),
                           season = factor(levels(t.dat0$season), levels = levels(t.dat0$season)),
                           diet2 = factor(levels(t.dat0$diet2), levels = levels(t.dat0$diet2)))

newdat_migr <- expand.grid(habitat = factor(levels(t.dat0$habitat), levels = levels(t.dat0$habitat)),
                           season = factor(levels(t.dat0$season), levels = levels(t.dat0$season)),
                           migration = factor(levels(t.dat0$migration), levels = levels(t.dat0$migration)))

newdat_habitat <- expand.grid(habitat = factor(levels(t.dat0$habitat), levels = levels(t.dat0$habitat)),
                              season = factor(levels(t.dat0$season), levels = levels(t.dat0$season)),
                              habitat_pref = factor(levels(t.dat0$habitat_pref), levels = levels(t.dat0$habitat_pref)))

newdat_habitat_no_euplectes <- expand.grid(habitat = factor(levels(t.dat0_no_euplectes$habitat), levels = levels(t.dat0_no_euplectes$habitat)),
                                           season = factor(levels(t.dat0_no_euplectes$season), levels = levels(t.dat0_no_euplectes$season)),
                                           habitat_pref = factor(levels(t.dat0_no_euplectes$habitat_pref), levels = levels(t.dat0_no_euplectes$habitat_pref)))

## first for diet

posterior_interval(mod_diet2, prob = 0.95)

#add random effect to newdat
tt <- ranef(mod_diet2)$landscape_id
t.landscape_id_min <- rownames(tt)[abs(tt) == min(abs(tt))][1]  
newdat_diet$landscape_id <- t.landscape_id_min
tt <- ranef(mod_diet2)$ol
t.ol.min <- rownames(tt)[abs(tt) == min(abs(tt))][1]
newdat_diet$ol <- t.ol.min

# predict the values
fitmat <- posterior_epred(mod_diet2, newdata = newdat_diet, re.form = NA)

# calculate quantiles:
newdat_diet$lwr <- apply(fitmat, 2, quantile, probs = 0.025)
newdat_diet$upr <- apply(fitmat, 2, quantile, probs = 0.975)
newdat_diet$fit <- apply(fitmat, 2, mean)

newdat_diet <- newdat_diet %>%
  mutate(interaction_var = factor(interaction(habitat, season, sep = "_"),
                                  levels = c("enclos_dry", "enclos_post rainy", "open_dry", 
                                             "open_post rainy", "woody_dry", "woody_post rainy")))


# Plot for diet
plot_diet <- ggplot(newdat_diet) +  
  geom_pointrange(aes(x = diet2, y = fit, ymin = lwr, ymax = upr, col = interaction_var),
                  size = 5, fatten = 0.5, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("woody_dry" = "dodgerblue2", "open_dry" = "darkorange", "enclos_dry" = "chartreuse4", 
                                "woody_post rainy" = "#C9D9FF", "open_post rainy" = "#FFDED1", "enclos_post rainy" = "#A3E786"), 
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Diet") +
  ylab("Occurence probability") +
  theme_newtree() +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = .6),
                              title = "Habitat and Season", title.position = "top", title.hjust = .5)) +  
  theme(legend.position = "bottom",  # Position legend at bottom right
        legend.justification = "center",  # Justify the legend to the bottom right
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(10,10,10,10),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 15),
        plot.margin = margin(1,2,1,1, "cm"))


# repeat the process for migration
posterior_interval(mod_migr2, prob = 0.95)

tt <- ranef(mod_migr2)$landscape_id
t.landscape_id_min <- rownames(tt)[abs(tt) == min(abs(tt))][1]  
newdat_migr$landscape_id <- t.landscape_id_min
tt <- ranef(mod_migr2)$ol
t.ol.min <- rownames(tt[abs(tt) == min(abs(tt))])[1]
newdat_migr$ol <- t.ol.min

fitmat <- posterior_epred(mod_migr2, newdata = newdat_migr, re.form = NA)

newdat_migr$lwr <- apply(fitmat, 2, quantile, probs = 0.025)
newdat_migr$upr <- apply(fitmat, 2, quantile, probs = 0.975)
newdat_migr$fit <- apply(fitmat, 2, mean)

newdat_migr <- newdat_migr %>%
  mutate(interaction_var = factor(interaction(habitat, season, sep = "_"),
                                  levels = c("enclos_dry", "enclos_post rainy", "open_dry", 
                                             "open_post rainy", "woody_dry", "woody_post rainy")))

# plot for migration
plot_migr <- ggplot(newdat_migr) +  
  geom_pointrange(aes(x = migration, y = fit, ymin = lwr, ymax = upr, col = interaction_var),
                  size = 5, fatten = 0.5, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("woody_dry" = "dodgerblue2", "open_dry" = "darkorange", "enclos_dry" = "chartreuse4", 
                                "woody_post rainy" = "#C9D9FF", "open_post rainy" = "#FFDED1", "enclos_post rainy" = "#A3E786"), 
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Migration") +
  ylab("Occurence probability") +
  theme_newtree() +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = .6),
                              title = "Habitat and Season", title.position = "top", title.hjust = .5)) +  
  theme(legend.position = "bottom",  # Position legend at bottom right
        legend.justification = "center",  # Justify the legend to the bottom right
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(10,10,10,10),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 15),
        plot.margin = margin(1,2,1,1, "cm"))

# And finally for habitat preference
posterior_interval(mod_habitat, prob = 0.95)

tt <- ranef(mod_habitat)$landscape_id
t.landscape_id_min <- rownames(tt)[abs(tt) == min(abs(tt))][1]  
newdat_habitat$landscape_id <- t.landscape_id_min

fitmat <- posterior_epred(mod_habitat, newdata=newdat_habitat, re.form = NA)

newdat_habitat$lwr <- apply(fitmat,2,quantile,probs=0.025)
newdat_habitat$upr <- apply(fitmat,2,quantile,probs=0.975)
newdat_habitat$fit <- apply(fitmat,2,mean)
newdat_habitat <- newdat_habitat %>%
  mutate(interaction_var = factor(interaction(habitat, season, sep = "_"),
                                  levels = c("enclos_dry", "enclos_post rainy", "open_dry", 
                                             "open_post rainy", "woody_dry", "woody_post rainy")))

# plot for habitat
ggplot(newdat_habitat) +  
  geom_pointrange(aes(x = habitat_pref, y = fit, ymin = lwr, ymax = upr, col = interaction_var),
                  size = 5, fatten = 0.5, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("woody_dry" = "dodgerblue2", "open_dry" = "darkorange", "enclos_dry" = "chartreuse4", 
                                "woody_post rainy" = "#C9D9FF", "open_post rainy" = "#FFDED1", "enclos_post rainy" = "#A3E786"), 
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Preferred habitat") +
  ylab("Occurence probability") +
  theme_newtree() +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = .6),
                              title = "Habitat and Season", title.position = "top", title.hjust = .5)) +  
  theme(legend.position = "bottom",  # Position legend at bottom right
        legend.justification = "center",  # Justify the legend to the bottom right
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(10,10,10,10),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 15),
        plot.margin = margin(1,2,1,1, "cm"))




# and with habitat no euplectes
posterior_interval(mod_habitat_no_euplectes2, prob = 0.95)

tt <- ranef(mod_habitat_no_euplectes2)$landscape_id
t.landscape_id_min <- rownames(tt)[abs(tt) == min(abs(tt))][1]  
newdat_habitat_no_euplectes$landscape_id <- t.landscape_id_min
tt <- ranef(mod_habitat_no_euplectes2)$ol
t.ol.min <- rownames(tt)[abs(tt) == min(abs(tt))][1]
newdat_habitat_no_euplectes$ol <- t.ol.min

fitmat <- posterior_epred(mod_habitat_no_euplectes2, newdata = newdat_habitat_no_euplectes, re.form = NA)

newdat_habitat_no_euplectes$lwr <- apply(fitmat, 2, quantile, probs = 0.025)
newdat_habitat_no_euplectes$upr <- apply(fitmat, 2, quantile, probs = 0.975)
newdat_habitat_no_euplectes$fit <- apply(fitmat, 2, mean)
newdat_habitat_no_euplectes <- newdat_habitat_no_euplectes %>%
  mutate(interaction_var = factor(interaction(habitat, season, sep = "_"),
                                  levels = c("enclos_dry", "enclos_post rainy", "open_dry", 
                                             "open_post rainy", "woody_dry", "woody_post rainy")))

# plot for habitat (without euplectes)
plot_habitat <- ggplot(newdat_habitat_no_euplectes) +  
  geom_pointrange(aes(x = habitat_pref, y = fit, ymin = lwr, ymax = upr, col = interaction_var),
                  size = 5, fatten = 0.5, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("woody_dry" = "dodgerblue2", "open_dry" = "darkorange", "enclos_dry" = "chartreuse4", 
                                "woody_post rainy" = "#C9D9FF", "open_post rainy" = "#FFDED1", "enclos_post rainy" = "#A3E786"), 
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  xlab("Preferred habitat") +
  ylab("Occurence probability") +
  theme_newtree() +
  theme(legend.key.size = unit(0.5, "lines")) +
  guides(color = guide_legend(override.aes = list(size = .6),
                              title = "Habitat and Season", title.position = "top", title.hjust = .5)) +  
  theme(legend.position = "bottom",  # Position legend at bottom right
        legend.justification = "center",  # Justify the legend to the bottom right
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(10,10,10,10),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size = 15),
        plot.margin = margin(1,2,1,1, "cm"))

## combine plots (Fig 5)

#pdf("guilds.pdf", width = 10, height = 22)
#grid.arrange(plot_habitat, plot_migr, plot_diet, nrow = 3)
#dev.off()


### Supplementary figures-------------------------------------------------------


### Herbaceous vegetation

by_survey <- by_survey %>% 
  mutate(interaction_var = paste(habitat, season, sep = "_"))

#pdf("figure_sup_herb.pdf", height = 5, width = 10)

grid.arrange(

ggplot(by_survey) +  
  geom_boxplot(aes(x = habitat, y = herb_cov, col = interaction_var)) +
  scale_color_manual(values = c("enclos_dry" = "chartreuse4", "enclos_post rainy" = "#A3E786", "open_dry" = "darkorange", 
                                "open_post rainy" = "#FFDED1", "woody_dry" = "dodgerblue2","woody_post rainy" = "#C9D9FF"),
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Habitat") +
  ylab("Herbaceous vegetation coverage (%)") +
  theme_newtree(), 


ggplot(by_survey) +  
  geom_boxplot(aes(x = habitat, y = herb_height, col = interaction_var)) +
  scale_color_manual(values = c("enclos_dry" = "chartreuse4", "enclos_post rainy" = "#A3E786", "open_dry" = "darkorange", 
                                "open_post rainy" = "#FFDED1", "woody_dry" = "dodgerblue2","woody_post rainy" = "#C9D9FF"),
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Habitat") +
  ylab("Height of herbaceous vegetation (cm)") +
  theme_newtree() ,

nrow = 1)

#dev.off()


### Trees
col_hab <- c("chartreuse4", "darkorange", "dodgerblue")

#pdf("figure_sup_tree.pdf", width = 15, height = 5)

grid.arrange(

ggplot(by_survey) +  
  geom_boxplot(aes(x = habitat, y = plant_rich_tot, col = habitat)) +
  scale_color_manual(values = col_hab) +
  xlab("Habitat") +
  ylab("Tree richness") +
  theme_newtree() ,


ggplot(by_survey) +  
  geom_boxplot(aes(x = habitat, y = tot_strata, col = habitat)) +
  scale_color_manual(values = col_hab) +
  xlab("Habitat") +
  ylab("Number of trees") +
  theme_newtree() ,


ggplot(by_survey) +  
  geom_boxplot(aes(x = habitat, y = Shannon_vertical, col = habitat)) +
  scale_color_manual(values = col_hab) +
  xlab("Habitat") +
  ylab("Vegetation heterogeneity (Shannon)") +
  theme_newtree() ,

nrow = 1)

#dev.off()


### NDVI

#pdf("figure_sup_NDVI.pdf", width = 10, height = 5)

grid.arrange(
  
ggplot(by_survey) +  
  geom_boxplot(aes(x = habitat, y = max.NDVI.50, col = interaction_var)) +
  scale_color_manual(values = c("enclos_dry" = "chartreuse4", "enclos_post rainy" = "#A3E786", "open_dry" = "darkorange", 
                                "open_post rainy" = "#FFDED1", "woody_dry" = "dodgerblue2","woody_post rainy" = "#C9D9FF"),
                     labels = c("Exclos (Dry)", "Exclos (Wet)", 'Open (Dry)', "Open (Wet)", "Woody (Dry)", "Woody (Wet)")) +
  xlab("Habitat") +
  ylab("max NDVI 50") +
  theme_newtree() ,


ggplot(by_survey) +  
  geom_boxplot(aes(x = season, y = max.NDVI.500, col = season)) +
  scale_color_manual(values = c("darkgrey", "grey"),
                     labels = c("Dry", "Wet")) +
  xlab("Habitat") +
  ylab("max NDVI 500") +
  theme_newtree() ,

nrow = 1)

#dev.off()
