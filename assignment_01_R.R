###===============================
## BINF6210 Assignment 1
## Kenneth Gamueda
## 2025-10-12

### === COMMENT CODES ===========
# Coding explanations (# XXX)
# Justification for a section of code (### XXX)
# Solutions/results/interpretations (#==> XXX)
# Conclusions (# ***)
# Performance check for code (#PC)

### === PACKAGES USED ===========
library("tidyverse")
library("tidyr")
library("ggplot2")
library("dplyr")
library("vegan")
library("cowplot")
library("geodata")
library("terra")
library("rayshader")
library("png")
library("grid")

# 1 === INPUTTING AND FILTERING DATA ========
data_BOLD_Araneae <- read_tsv("../data/BOLD_results_Araneae.tsv")

names(data_BOLD_Araneae) 
#==> sampleid will be used as measure of the number of samples, BIN will be used as proxy for species and be the dependent variable, and elevation will be the independent variable

df_Araneae <- data_BOLD_Araneae %>% 
  select("processid", "sampleid", "bin_uri", "country/ocean", "elev") %>%
  filter(!is.na(elev) & !is.na(bin_uri) & elev > 0)
#count of records reduced from 205455 to 122027
#filtering criteria more specified to remove records with no elevation values or elevation values below ground (0) and missing BIN assignment
any(is.na(df_Araneae$elev), is.na(df_Araneae$bin_uri), df_Araneae < 0) #PC

# Conducting preliminary analysis of elevation values
summary(df_Araneae$elev)
hist(df_Araneae$elev)
#==> Histogram shows a right-skew for global data

sort(table(df_Araneae$"country/ocean"), decreasing = TRUE)
#==> Test the hypothesis within Canada as it has the most records

df_Araneae_CA <- df_Araneae %>% 
  filter(`country/ocean` == "Canada")
#==> Number of records for went from 122027 to 56424

hist(df_Araneae_CA$elev) 

#by adding prior filters the distribution is more even. Consider constructing elevational bins prior to removing records of low frequency, especially when using an arbitrary cut-off value. Retaining the full elevation range would be preferable for reproducibility and transparency, especially since the banding step already controls for sparse sampling at extreme elevations.

# # Filter elevations with low number of records
# df_Araneae_CA <- df_Araneae_CA %>% 
#   filter (elev > 0 & elev < 2500) ### Removing records with too low of a frequency
# #==> Number of records went from 58207 to 58008
# any(df_Araneae_CA$elev < 0) #PC
# any(df_Araneae_CA$elev > 2500) #PC
# 
# # Creating elevation bands to reduce noise from uneven sampling
# df_Araneae_CA <- df_Araneae_CA %>%
#   mutate(elev_band_CA = cut(elev, 
#                             breaks = seq(0, 2500, by = 100), 
#                             include.lowest = T, 
#                             right = F, 
#                             dig.lab = 5)) 
# ### 25 elevation bands were created to ensure that each band contained enough samples but was also sensitive to any smaller-scale trends
# table(df_Araneae_CA$elev_band_CA) #PC
# 
# midpoints <- (head(seq(0, 2500, by = 100), n = -1) + tail(seq(0, 2500, by = 100), n = -1)) / 2

#Use zero as min bound and determine upper bound
min_elev <- 0
max_elev <- ceiling(max(df_Araneae_CA$elev, na.rm = TRUE))

bin_width <- 100

# generate the sequence of breaks into width = 100
breaks <- seq(min_elev, max_elev + bin_width, by = bin_width)

#Create elevation bands
df_Araneae_CA <- df_Araneae_CA %>%
  mutate(elev_band_CA = cut(elev,
                            breaks = breaks,
                            include.lowest = TRUE,
                            right = FALSE,
                            dig.lab = 5))

#Compute midpoints for each bin
midpoints <- (head(breaks, -1) + tail(breaks, -1)) / 2

df_Araneae_CA <- df_Araneae_CA %>%
  mutate(elev_midpoint_CA = midpoints[as.numeric(elev_band_CA)]) ### Calculating midpoints for modelling and plotting
table(df_Araneae_CA$elev_midpoint_CA) #PC


#remove variables as you are working to keep environment tidy
rm(data_BOLD_Araneae, min_elev, max_elev, bin_width, breaks, midpoints)

# 2 === CHECKING FOR SAMPLING BIAS/EFFORT ========
# Before plotting BIN richness (proxy for species richness) against elevation, need to check if sampling bias plays a role (i.e. sample sizes should be relatively even at different elevations)
df_sub_Araneae_CA <- df_Araneae_CA %>%
  group_by(elev_midpoint_CA, elev_band_CA) %>%
  summarize(samples = n_distinct(sampleid), BIN_richness = n_distinct(bin_uri)) ### To check for sampling effort, need to know the number of samples per elevation band
# names(df_sub_Araneae_CA) #PC
# 
# # Quick visual check of number of samples in each elevation band
# effort_Araneae_plot_1 <- ggplot(df_sub_Araneae_CA, aes(x = elev_band_CA, y = samples)) + 
#   geom_col(fill = "#61497b") +
#   labs(title = "Sampling effort of Araneae across elevation in Canada", 
#        x = "Elevation (m)", 
#        y = "Number of samples") +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   theme_light()
# effort_Araneae_plot_1
# #==> Large sample sizes for low elevation band (between 0 and 100), but very low sample sizes for higher elevation bands
# 
# # Need to create models to check if sampling effort plays a role in predicting BIN richness
# lm_effort_Araneae_CA <- lm(samples ~ elev_midpoint_CA, df_sub_Araneae_CA)
# summary(lm_effort_Araneae_CA)
# plot(lm_effort_Araneae_CA)
# #==> p < 0.001, but R-squared = 0.43 and residuals show a non-linear relationship
# 
# pm_effort_Araneae_CA <- lm(samples ~ poly(elev_midpoint_CA, degree = 5), df_sub_Araneae_CA) ### Creating a model with polynomial regression to better reflect non-linear relationship
# summary(pm_effort_Araneae_CA)
# plot(pm_effort_Araneae_CA)
# #==> p < 0.001 and R-squared = 0.94
# 
# anova(lm_effort_Araneae_CA, pm_effort_Araneae_CA)
# #==> p < 0.001, so fit is much better with the polynomial regression model
# 
# effort_Araneae_plot_2 <- ggplot(df_sub_Araneae_CA, 
#        aes(x = elev_midpoint_CA, y = samples)) +
#   geom_point(size = 2.5) +
#   stat_smooth(method = lm, 
#               formula = y ~ poly(x, degree = 5), 
#               colour = "#df4b68") + 
#   labs(title = "Sampling effort of Araneae across elevation in Canada", 
#        x = "Elevation (m)", 
#        y = "Number of samples") +
#   theme_minimal()
# effort_Araneae_plot_2
# #==> There are much less samples at higher elevations, so need to account for sampling effort and adjust BIN richness
# 
# # Is there also a relationship between the number of samples and BIN richness?
# summary(lm(BIN_richness ~ poly(samples, degree = 2), df_sub_Araneae_CA))
# plot(lm(BIN_richness ~ poly(samples, degree = 2), df_sub_Araneae_CA))
# #==> p < 0.001 and R-squared = 0.94
# 
# effort_Araneae_plot_3 <- ggplot(df_sub_Araneae_CA, 
#        aes(x = samples, y = BIN_richness)) +
#   geom_point(size = 2.5) +
#   stat_smooth(method = lm, 
#               se = FALSE,
#               formula = y ~ poly(x, degree = 2), 
#               colour = "#df4b68") +
#   labs(title = "Number of samples vs BIN richness of Araneae in Canada", 
#        x = "Number of samples", 
#        y = "BIN richness") +
#   theme_minimal()
# effort_Araneae_plot_3
# #==> Richness also increases with sampling effort; thus another reason to account for sampling effort

#visual check for sampling bias
ggplot(df_sub_Araneae_CA, aes(x = elev_midpoint_CA, y = samples)) +
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "Sampling effort across elevation",
    x = "Elevation (m)",
    y = "Number of samples"
  ) +
  theme_minimal()

#correlation test to confirm
cor.test(df_sub_Araneae_CA$samples, df_sub_Araneae_CA$elev_midpoint_CA)
#==>p < 0.001

#check relationship between samples and BIN richness
ggplot(df_sub_Araneae_CA, aes(x = samples, y = BIN_richness)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "BIN richness vs sampling effort",
    x = "Number of samples",
    y = "BIN richness"
  ) +
  theme_minimal()

#based on the trend between samples and BIN richness, test log-transformed samples 
summary(lm(BIN_richness ~ log(samples), df_sub_Araneae_CA))
#==>p < 0.001

#*** Account for sampling effort using two methods: 1) including log-transformed sampling effort as a covariate in regression model 2) rarefying richness

# 3 === SAMPLING EFFORT AS COVARIATE ========
# # Creating new model with sampling effort as a covariate
# model_Araneae_CA <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 5) + samples, df_sub_Araneae_CA)
# summary(model_Araneae_CA)
# #==> p < 0.001 and R-squared = 0.98, so elevation is a statistically significant predictor of richness after correcting for sampling effort

#choose polynomial degree by comparing AIC scores
model1 <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 1), data = df_sub_Araneae_CA)
model2 <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 2), data = df_sub_Araneae_CA)
model3 <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 3), data = df_sub_Araneae_CA)
model4 <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 4), data = df_sub_Araneae_CA)
model5 <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 5), data = df_sub_Araneae_CA)

# Compare AIC
AIC(model1, model2, model3, model4, model5)
#==> lowest AIC score = model 5
#==> delta AIC < 2 between model 4 and 5 
#==> simpler model is preferred as it achieves a comparable fit with less complexity

rm(model1, model2, model3, model4, model5)

#model BIN richness as a function of elevation using polynomial degree 4 and log-transformed samples covariate
model_Araneae_CA <- lm(BIN_richness ~ poly(elev_midpoint_CA, degree = 4) + log(samples), df_sub_Araneae_CA)
summary(model_Araneae_CA)

# Since sampling effort is now included as a covariate, can now plot corrected richness with elevation
sampling_effort <- mean(df_sub_Araneae_CA$samples) ### Quantify sampling effort as the mean of sample of numbers per elevation band
df_sub_Araneae_CA$corr_BIN_richness <- predict(object = model_Araneae_CA, 
                                               newdata = data.frame (samples = sampling_effort, 
                                                                     elev_midpoint_CA = df_sub_Araneae_CA$elev_midpoint_CA))

Araneae_CA_corr_plot <- ggplot(df_sub_Araneae_CA, 
                               aes(x = elev_midpoint_CA, y = corr_BIN_richness)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(y = corr_BIN_richness),
              colour = "#df4b68") +
  labs(title = "BIN richness of Araneae across elevation in Canada", 
       x = "Elevation (m)", 
       y = "Corrected BIN richness") +
  theme_light()
Araneae_CA_corr_plot
##==> After sampling effort has been accounted for, the overall trend is a decrease in species richness as elevation increases
#==> After sampling effort has been accounted for, the overall trend is mostly decreasing, with some increase after 2000m

rm(df_sub_Araneae_CA, model_Araneae_CA, Araneae_CA_corr_plot, sampling_effort)

# 4 === RAREFYING RICHNESS ========
# Calculating rarefied richness to account for sampling bias
matrix_comm_Araneae_CA <- df_Araneae_CA %>%
  group_by(elev_band_CA, bin_uri) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = bin_uri, 
              values_from = n, 
              values_fill = 0) %>%
  column_to_rownames(var = "elev_band_CA") ### vegan requires data to be in the form of a matrix; rows are elevation bands (i.e. sites) and BINs are columns

raremin <- min(rowSums(matrix_comm_Araneae_CA))
#==> minimum sample size is 46; sufficient sample size for rarefaction
rare_BIN_richness <- rarefy(matrix_comm_Araneae_CA, raremin)

df_rare_Araneae_CA <- data.frame(elev_band_CA = names(rare_BIN_richness), 
                                 rare_BIN_richness = as.numeric(rare_BIN_richness))
df_rare_Araneae_CA <- df_sub_Araneae_CA %>%
  left_join(df_rare_Araneae_CA, by = "elev_band_CA")

Araneae_CA_rare_plot <- ggplot(df_rare_Araneae_CA, 
                               aes(x = elev_midpoint_CA, y = rare_BIN_richness)) +
  geom_point(size = 2.5) +
  stat_smooth(method = lm, 
              formula = y ~ poly(x, 2), 
              colour = "#df4b68") +
  labs(title = "BIN richness of Araneae across elevation in Canada", 
       x = "Elevation (m)", 
       y = "Rarefied BIN richness") + 
  theme_light()
Araneae_CA_rare_plot
#==> After richness has been rarefied, overall trend is a decrease in species richness as elevation increases

plot_grid(Araneae_CA_corr_plot, Araneae_CA_rare_plot, labels = c("A", "B"))
#==> Rarefied richness has a bigger confidence interval likely due to the low sample size used during rarefaction

# 5 === EXPLORING TRENDS IN BIRDS ========
# The workflow above will be repeated for the Class Aves to test if the elevational diversity gradient exists in other taxa
data_BOLD_Aves <- read_tsv("../data/BOLD_results_Aves.tsv")
df_Aves <- data_BOLD_Aves %>% 
  select("processid", "sampleid", "bin_uri", "country/ocean", "elev") %>%
  filter(!is.na(elev) & !is.na(bin_uri) & elev > 0)
#count of records reduced from 67047 to 8956
#filtering criteria more specified to remove records with no elevation values or elevation values below ground (0) and missing BIN assignment
any(is.na(df_Aves$elev), is.na(df_Aves$bin_uri), df_Araneae < 0) #PC

# Conducting preliminary analysis of elevation values
summary(df_Aves$elev)
hist(df_Aves$elev)
#==> bigger range but also shows same right-skew in values

sort(table(df_Aves$"country/ocean"), decreasing = TRUE)
#==> Test the hypothesis within Peru as it has the most records

df_Aves_PE <- df_Aves %>% filter(`country/ocean` == "Peru")
hist(df_Aves_PE$elev) 
# df_Aves_PE <- df_Aves_PE %>% filter (elev < 4000) ### Removing records with too low of a frequency
# #==> Number of records went from 2474 to 2445
# any(df_Aves_PE$elev > 4000) #PC
# 
# df_Aves_PE <- df_Aves_PE %>%
#   mutate(elev_band_PE = cut(elev, 
#                             breaks = seq(0, 4000, by = 150), 
#                             include.lowest = T, 
#                             right = F, 
#                             dig.lab = 5)) 
# ### From previous analysis with Araneae, ~25 bands were sufficient
# table(df_Aves_PE$elev_band_PE) #PC
# 
# midpoints_Aves <- (head(seq(0, 4000, by = 150), n = -1) + tail(seq(0, 4000, by = 150), n = -1)) / 2

#Use zero as min bound and determine upper bound
min_elev <- 0 
max_elev <- ceiling(max(df_Aves_PE$elev, na.rm = TRUE))

#consider doubling the Araneae bin width for Aves given the lower overall record count and elevation bin count distribution
#Generate the sequence of breaks using bin width 200
bin_width <- 200

breaks <- seq(min_elev, max_elev + bin_width, by = bin_width)

#Create elevation bins and organize data into bins
df_Aves_PE <- df_Aves_PE %>%
  mutate(elev_band_PE = cut(elev,
                            breaks = breaks,
                            include.lowest = TRUE,
                            right = FALSE,
                            dig.lab = 5))

#Compute midpoints for each bin
midpoints_Aves <- (head(breaks, -1) + tail(breaks, -1)) / 2

df_Aves_PE <- df_Aves_PE %>%
  mutate(elev_midpoint_PE = midpoints_Aves[as.numeric(elev_band_PE)]) 
table(df_Aves_PE$elev_midpoint_PE) #PC

#remove variables as you are working to keep environment tidy
rm(data_BOLD_Aves, min_elev, max_elev, bin_width, breaks, midpoints_Aves)

df_sub_Aves_PE <- df_Aves_PE %>%
  group_by(elev_midpoint_PE, elev_band_PE) %>%
  summarize(samples = n_distinct(sampleid), BIN_richness = n_distinct(bin_uri)) 

# # Quick check of number of samples in each elevation band
# ggplot(df_sub_Aves_PE, 
#        aes(x = elev_band_PE, y = samples)) + 
#   geom_col(fill = "#61497b") +
#   labs(title = "Sampling effort of Aves across elevation in Peru", 
#        x = "Elevation (m)", 
#        y = "Number of samples") +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   theme_light()
# #==> Different distribution compared to Araneae; looks to be trimodal
# 
# ggplot(df_sub_Aves_PE, aes(x = elev_midpoint_PE, weight = samples)) +
#   geom_density(fill = "#61497b") + 
#   labs(title = "Density plot of Aves samples across elevation", x = "Elevation (m)", y = "Density of samples")
# #==> Density plot shows a bimodal structure

#visual check for sampling bias
ggplot(df_sub_Aves_PE, aes(x = elev_midpoint_PE, y = samples)) +
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "Sampling effort across elevation",
    x = "Elevation (m)",
    y = "Number of samples"
  ) +
  theme_minimal()

#correlation test to confirm
cor.test(df_sub_Aves_PE$samples, df_sub_Aves_PE$elev_midpoint_PE)
#==>p < 0.05

#check relationship between samples and BIN richness
ggplot(df_sub_Aves_PE, aes(x = samples, y = BIN_richness)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "BIN richness vs sampling effort",
    x = "Number of samples",
    y = "BIN richness"
  ) +
  theme_minimal()

#based on the trend between samples and BIN richness, test log-transformed samples
summary(lm(BIN_richness ~ log(samples), df_sub_Aves_PE))
#==>p < 0.001

#ensure to apply same methods for sampling bias check for Aves as Araneae for standardization

# # Since some elevation bands have considerably more samples than others, also need to account for sampling effort 
# model_Aves_PE <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 8) + samples, df_sub_Aves_PE) ### Degree 8 offered highest R-squared value without overfitting
# summary(model_Aves_PE)
# #==> p < 0.001 and R-squared = 0.88, so elevation is a statistically significant predictor of richness after correcting for sampling effort

#choose polynomial degree by comparing AIC scores
model1 <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 1), data = df_sub_Aves_PE)
model2 <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 2), data = df_sub_Aves_PE)
model3 <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 3), data = df_sub_Aves_PE)
model4 <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 4), data = df_sub_Aves_PE)
model5 <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 5), data = df_sub_Aves_PE)

# Compare AIC
AIC(model1, model2, model3, model4, model5)
#==> lowest AIC score = model 5
#==> delta AIC < 2 between model 4 and 5 
#==> simpler model is preferred as it achieves a comparable fit with less complexity

rm(model1, model2, model3, model4, model5)

#model BIN richness as a function of elevation using polynomial degree 4 and log-transformed samples covariate
model_Aves_PE <- lm(BIN_richness ~ poly(elev_midpoint_PE, degree = 4) + log(samples), df_sub_Aves_PE)
summary(model_Aves_PE)

# Since sampling effort is now included as a covariate, can now plot corrected richness with elevation
sampling_effort_Aves <- mean(df_sub_Aves_PE$samples) ### Quantify sampling effort as the mean of sample of numbers per elevation band
df_sub_Aves_PE$corr_BIN_richness <- predict(object = model_Aves_PE, 
                                            newdata = data.frame(samples = sampling_effort_Aves, 
                                                                 elev_midpoint_PE = df_sub_Aves_PE$elev_midpoint_PE))

Aves_PE_corr_plot <- ggplot(df_sub_Aves_PE, 
                            aes(x = elev_midpoint_PE, y = corr_BIN_richness)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(y = corr_BIN_richness), 
              colour = "#df4b68") +
  labs(title = "BIN richness of Aves across elevation in Peru", 
       x = "Elevation (m)", 
       y = "Corrected BIN richness") +
  theme_light()
Aves_PE_corr_plot
##==> Unimodal distribution, which is another type of trend observed in elevational diversity gradients
#==> mostly unimodal but one one spike due to extreme elevation (5000m)

rm(df_sub_Aves_PE, model_Aves_PE, Aves_PE_corr_plot, sampling_effort_Aves)

# Rarefying richness in Aves
matrix_comm_Aves_PE <- df_Aves_PE %>%
  group_by(elev_band_PE, bin_uri) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = bin_uri, 
              values_from = n, 
              values_fill = 0) %>%
  column_to_rownames(var = "elev_band_PE") ### vegan requires data to be in the form of a matrix; rows are elevation bands (i.e. sites) and BINs are columnns

raremin <- min(rowSums(matrix_comm_Aves_PE))
#==> Sample size to be used for rarefaction is only 4, so rarefied richness will likely be inaccurate; continuing analysis but not expecting accurate results
rare_BIN_richness <- rarefy(matrix_comm_Aves_PE, raremin)

df_rare_Aves_PE <- data.frame(elev_band_PE = names(rare_BIN_richness), 
                              rare_BIN_richness = as.numeric(rare_BIN_richness))
df_rare_Aves_PE <- df_sub_Aves_PE %>%
  left_join(df_rare_Aves_PE, by = "elev_band_PE")

pm_rare_Aves_PE <- lm(rare_BIN_richness ~ poly(elev_midpoint_PE, degree = 3), df_rare_Aves_PE) ### Degree 3 had highest R-squared without overftting
summary(pm_rare_Aves_PE)
#==> p > 0.05 and R-squared = 0.06; as expected, sample size for rarefaction was too small

# Plotting the relationship with rarefied richness to see it visually
Aves_PE_rare_plot <- ggplot(df_rare_Aves_PE, 
                            aes(x = elev_midpoint_PE, y = rare_BIN_richness)) +
  geom_point(size = 2.5) +
  stat_smooth(method = lm, 
              formula = y ~ poly(x, degree = 3), 
              colour = "#df4b68") +
  labs(title = "BIN richness of Aves across elevation in Peru", 
       x = "Elevation (m)", 
       y = "Rarefied BIN richness") + 
  theme_light()
Aves_PE_rare_plot
#==> Trend does not match that of other models, most likely due to the low sample size used during rarefaction

# *** For populations with a low minimum sample size, regression model with covariate will be used to examine the elevational diversity gradient

# Avian studies of elevational diversity gradients sometimes set minimum value of elevation as 1000m
Aves_PE_1000_plot <- ggplot(filter(df_sub_Aves_PE, elev_midpoint_PE > 1000), 
                            aes(x = elev_midpoint_PE, y = corr_BIN_richness)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(y = corr_BIN_richness), 
              colour = "#df4b68") +
  labs(title = "BIN richness of Aves across elevation in Peru", 
       x = "Elevation (m)", 
       y = "Corrected BIN richness") +
  theme_light()
Aves_PE_1000_plot
#==> Clear unimodal distribution, which better fits in line with what was observed in previous research

plot_grid(Aves_PE_corr_plot, Aves_PE_1000_plot, labels = c("A", "B"))
# *** Depending on the elevation cutoffs (which is related to where samples are taken from), the overall trend observed may differ

# 6 === EXPLORING RELATIONSHIPS ACROSS COUNTRIES ========
# For different taxa in the same country, what pattern of elevational diversity gradient is observed?
sort(table(df_Araneae$"country/ocean"), decreasing = TRUE)
sort(table(df_Aves$"country/ocean"), decreasing = TRUE)
#==> For both taxa, Argentina has sufficient sample size (1651 for Araneae, 1992 for Aves)

df_Araneae_AR <- df_Araneae %>% 
  filter(`country/ocean` == "Argentina")
df_Aves_AR <- df_Aves %>% 
  filter(`country/ocean` == "Argentina")

# To simplify analysis, polynomial regression with sampling efffort as the covariate will be used for both taxa; previous workflow will be repeated for the new subsets

## 6.A === Repeating workflow for Araneae in Argentina ========
summary(df_Araneae_AR$elev)
hist(df_Araneae_AR$elev)
#==> binomial distribution of samples

# df_Araneae_AR <- df_Araneae_AR %>% 
#   filter (elev > 0 & elev < 1600) ### Removing records with too low of a frequency
# #==> Number of records went from 1651 to 1606
# any(df_Araneae_AR$elev < 0) #PC
# any(df_Araneae_AR$elev > 1600) #PC
# 
# df_Araneae_AR <- df_Araneae_AR %>%
#   mutate(elev_band_AR = cut(elev, 
#                             breaks = seq(0, 1600, by = 100), 
#                             include.lowest = T, 
#                             right = F, 
#                             dig.lab = 5)) 
# table(df_Araneae_AR$elev_band_AR) #PC
# 
# midpoints_Araneae_AR <- (head(seq(0, 1600, by = 100), n = -1) + tail(seq(0, 1600, by = 100), n = -1)) / 2

min_elev <- 0
max_elev <- ceiling(max(df_Araneae_AR$elev, na.rm = TRUE))

bin_width <- 100

# generate the sequence of breaks into width = 100
breaks <- seq(min_elev, max_elev + bin_width, by = bin_width)

#Create elevation bands
df_Araneae_AR <- df_Araneae_AR %>%
  mutate(elev_band_AR = cut(elev,
                            breaks = breaks,
                            include.lowest = TRUE,
                            right = FALSE,
                            dig.lab = 5))

#Compute midpoints for each bin
midpoints_Araneae_AR <- (head(breaks, -1) + tail(breaks, -1)) / 2

df_Araneae_AR <- df_Araneae_AR %>%
  mutate(elev_midpoint_AR = midpoints_Araneae_AR[as.numeric(elev_band_AR)]) 
table(df_Araneae_AR$elev_midpoint_AR) #PC

rm(min_elev, max_elev, breaks, bin_width, midpoints_Araneae_AR)

df_sub_Araneae_AR <- df_Araneae_AR %>%
  group_by(elev_midpoint_AR, elev_band_AR) %>%
  summarize(samples = n_distinct(sampleid), BIN_richness = n_distinct(bin_uri)) 

# # Quick check of number of samples in each elevation band
# ggplot(df_sub_Araneae_AR, 
#        aes(x = elev_band_AR, y = samples)) + 
#   geom_col(fill = "#61497b") +
#   labs(title = "Sampling effort of Araneae across elevation in Argentina", 
#        x = "Elevation (m)", 
#        y = "Number of samples") +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   theme_light()

#visual check for sampling bias
ggplot(df_sub_Araneae_AR, aes(x = elev_midpoint_AR, y = samples)) +
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "Sampling effort across elevation",
    x = "Elevation (m)",
    y = "Number of samples"
  ) +
  theme_minimal()

#correlation test to confirm
cor.test(df_sub_Araneae_AR$samples, df_sub_Araneae_AR$elev_midpoint_AR)
#==>p < 0.05

#check relationship between samples and BIN richness
ggplot(df_sub_Araneae_AR, aes(x = samples, y = BIN_richness)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "BIN richness vs sampling effort",
    x = "Number of samples",
    y = "BIN richness"
  ) +
  theme_minimal()

#based on the trend between samples and BIN richness, test log-transformed samples 
summary(lm(BIN_richness ~ log(samples), df_sub_Araneae_AR))
#==>p < 0.001

# model_Araneae_AR <- lm(BIN_richness ~ elev_midpoint_AR + samples, df_sub_Araneae_AR)
# summary(model_Araneae_AR)
# #==> p < 0.001 and R-squared = 0.94, so elevation is a statistically significant predictor of richness after correcting for sampling effort

#choose polynomial degree by comparing AIC scores
model1 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 1), data = df_sub_Araneae_AR)
model2 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 2), data = df_sub_Araneae_AR)
model3 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 3), data = df_sub_Araneae_AR)
model4 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 4), data = df_sub_Araneae_AR)
model5 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 5), data = df_sub_Araneae_AR)
model6 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 6), data = df_sub_Araneae_AR)
model7 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 7), data = df_sub_Araneae_AR)
model8 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 8), data = df_sub_Araneae_AR)
model9 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 9), data = df_sub_Araneae_AR)
model10 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 10), data = df_sub_Araneae_AR)
model11 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 11), data = df_sub_Araneae_AR)
model12 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 12), data = df_sub_Araneae_AR)

# Compare AIC
AIC(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12)
#==> lowest AIC score = model 10

rm(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12)

#model BIN richness as a function of elevation using polynomial degree 4 and log-transformed samples covariate
model_Araneae_AR <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 10) + log(samples), df_sub_Araneae_AR)
summary(model_Araneae_AR)

# Since sampling effort is now included as a covariate, can now plot corrected richness with elevation
sampling_effort <- mean(df_sub_Araneae_AR$samples) ### Quantify sampling effort as the mean of sample of numbers per elevation band
df_sub_Araneae_AR$corr_BIN_richness <- predict(object = model_Araneae_AR, 
                                                  newdata = data.frame(samples = sampling_effort, 
                                                                       elev_midpoint_AR = df_sub_Araneae_AR$elev_midpoint_AR))

Araneae_AR_corr_plot <- ggplot(df_sub_Araneae_AR, 
                               aes(x = elev_midpoint_AR, y = corr_BIN_richness)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(y = corr_BIN_richness),
              colour = "#df4b68") +
  labs(title = "BIN richness of Araneae across elevation in Argentina", 
       x = "Elevation (m)", 
       y = "Corrected BIN richness") +
  theme_light()
Araneae_AR_corr_plot
##==> BIN richness decreases linearly with increasing elevation
#==> BIN richness displays exponential decline and saturation after 1000m

rm(df_sub_Araneae_AR, model_Araneae_AR, Araneae_AR_corr_plot, sampling_effort)

## 6.B === Repeating workflow for Aves in Argentina ========
summary(df_Aves_AR$elev)
hist(df_Aves_AR$elev)
#==> right-skewed distribution of samples

#remove biologically impossible elevations / outliers
df_Aves_AR <- df_Aves_AR %>% 
  filter (elev > 0 & elev < 7000) ### Removing records with too low of a frequency
#==> Number of records went from 1949 to 1940
any(df_Aves_AR$elev < 0) #PC
any(df_Aves_AR$elev > 7000) #PC
#view new histogram
hist(df_Aves_AR$elev)

# df_Aves_AR <- df_Aves_AR %>% 
#   filter (elev > 0 & elev < 4000) ### Removing records with too low of a frequency
# #==> Number of records went from 1992 to 1948
# any(df_Aves_AR$elev < 0) #PC
# any(df_Aves_AR$elev > 4000) #PC
# 
# df_Aves_AR <- df_Aves_AR %>%
#   mutate(elev_band_AR = cut(elev, 
#                             breaks = seq(0, 4000, by = 250), 
#                             include.lowest = T, 
#                             right = F, 
#                             dig.lab = 5)) 
# table(df_Aves_AR$elev_band_AR) #PC
# 
# midpoints_Aves_AR <- (head(seq(0, 4000, by = 250), n = -1) + tail(seq(0, 4000, by = 250), n = -1)) / 2

min_elev <- 0
max_elev <- ceiling(max(df_Aves_AR$elev, na.rm = TRUE))

bin_width <- 100

# generate the sequence of breaks into width = 100
breaks <- seq(min_elev, max_elev + bin_width, by = bin_width)

#Create elevation bands
df_Aves_AR <- df_Aves_AR %>%
  mutate(elev_band_AR = cut(elev,
                            breaks = breaks,
                            include.lowest = TRUE,
                            right = FALSE,
                            dig.lab = 5))

#Compute midpoints for each bin
midpoints_Aves_AR <- (head(breaks, -1) + tail(breaks, -1)) / 2

df_Aves_AR <- df_Aves_AR %>%
  mutate(elev_midpoint_AR = midpoints_Aves_AR[as.numeric(elev_band_AR)]) 
table(df_Aves_AR$elev_midpoint_AR) #PC

rm(data_BOLD_Aves, min_elev, max_elev, breaks, bin_width, midpoints_Aves_AR)

df_sub_Aves_AR <- df_Aves_AR %>%
  group_by(elev_midpoint_AR, elev_band_AR) %>%
  summarize(samples = n_distinct(sampleid), BIN_richness = n_distinct(bin_uri)) 

#visual check for sampling bias
ggplot(df_sub_Aves_AR, aes(x = elev_midpoint_AR, y = samples)) +
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "Sampling effort across elevation",
    x = "Elevation (m)",
    y = "Number of samples"
  ) +
  theme_minimal()

#correlation test to confirm
cor.test(df_sub_Aves_AR$samples, df_sub_Aves_AR$elev_midpoint_AR)
#==>p < 0.05

#check relationship between samples and BIN richness
ggplot(df_sub_Aves_AR, aes(x = samples, y = BIN_richness)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "#df4b68") +
  labs(
    title = "BIN richness vs sampling effort",
    x = "Number of samples",
    y = "BIN richness"
  ) +
  theme_minimal()

#based on the trend between samples and BIN richness, test log-transformed samples 
summary(lm(BIN_richness ~ log(samples), df_sub_Aves_AR))
#==>p < 0.001

# model_Aves_AR <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 2) + samples, df_sub_Aves_AR)
# summary(model_Aves_AR)
# #==> p < 0.001 and R-squared = 0.97, so elevation is a statistically significant predictor of richness after correcting for sampling effort

#choose polynomial degree by comparing AIC scores
model1 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 1), data = df_sub_Aves_AR)
model2 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 2), data = df_sub_Aves_AR)
model3 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 3), data = df_sub_Aves_AR)
model4 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 4), data = df_sub_Aves_AR)
model5 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 5), data = df_sub_Aves_AR)
model6 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 6), data = df_sub_Aves_AR)
model7 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 7), data = df_sub_Aves_AR)
model8 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 8), data = df_sub_Aves_AR)
model9 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 9), data = df_sub_Aves_AR)
model10 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 10), data = df_sub_Aves_AR)
model11 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 11), data = df_sub_Aves_AR)
model12 <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 12), data = df_sub_Aves_AR)

# Compare AIC
AIC(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12)
#==> lowest AIC score = model 9

rm(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10, model11, model12)

#model BIN richness as a function of elevation using polynomial degree  and log-transformed samples covariate
model_Aves_AR <- lm(BIN_richness ~ poly(elev_midpoint_AR, degree = 9) + log(samples), df_sub_Aves_AR)
summary(model_Aves_AR)

# Since sampling effort is now included as a covariate, can now plot corrected richness with elevation
sampling_effort <- mean(df_sub_Aves_AR$samples) ### Quantify sampling effort as the mean of sample of numbers per elevation band
df_sub_Aves_AR$corr_BIN_richness <- predict(object = model_Aves_AR, 
                                            newdata = data.frame(samples = sampling_effort, 
                                                                 elev_midpoint_AR = df_sub_Aves_AR$elev_midpoint_AR))

Aves_AR_corr_plot <- ggplot(df_sub_Aves_AR, 
                            aes(x = elev_midpoint_AR, y = corr_BIN_richness)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(y = corr_BIN_richness), 
              colour = "#df4b68") +
  labs(title = "BIN richness of Aves across elevation in Argentina", 
       x = "Elevation (m)", 
       y = "Corrected BIN richness") +
  theme_light()
Aves_AR_corr_plot
##==> BIN richness also decreases with increasing elevation, but relationship is non-linear
#==> BIN richness shows decline and saturation after 1000m

plot_grid(Araneae_AR_corr_plot, Aves_AR_corr_plot, labels = c("A", "B"))
# *** For different taxa within the same country, the elevational diversity gradient appears to follow the same trend (i.e. monotonic decrease in richness as elevation increases)

rm(df_sub_Aves_AR, model_Aves_AR, Aves_AR_corr_plot, sampling_effort)

# 7 === FIGURE GENERATION ========
# Figure of graphs generated to check for sampling effort
summary_plots <- plot_grid(effort_Araneae_plot_1, 
                           effort_Araneae_plot_3,
                           labels = c("A", "B"))
ggsave("../figs/summary_plots.png",
       plot = summary_plots,
       width = 11.8,
       height = 5.74,
       units = "in",
       dpi = 300)

# Figure of BIN richness of spiders and birds
taxa_richness_plots <- plot_grid(Araneae_CA_corr_plot, 
                           Araneae_CA_rare_plot,
                           Aves_PE_corr_plot,
                           Aves_PE_1000_plot,
                           labels = c("A", "B", "C", "D"))
ggsave("../figs/taxa_richness_plots.png", 
       plot = taxa_richness_plots,
       width = 11.8,
       height = 5.74,
       units = "in",
       dpi = 300)

# Figure of BIN richness of spiders and birds in same country
region_richness_plots <- plot_grid(Araneae_AR_corr_plot, 
                                   Aves_AR_corr_plot,
                                   labels = c("A", "B"))
ggsave("../figs/region_richness_plots.png", 
       plot = region_richness_plots, 
       width = 11.8,
       height = 5.74,
       units = "in",
       dpi = 300)

# 8 === VISUALIZING GENERAL TREND WITH 3D MAP ========
# Exploratory project investigating how topographic maps could be used to visualize trends in elevational diversity gradients
# As an example, the trend of decreasing BIN richness will be visualized in a 3D projection

# Obtaining Digital Elevation Model of Mt. Robson (chosen as an arbitrary representation of a mountain with a max elevation of 3954 m)
DEM_mtRobson <- elevation_3s(lon = -119.15639, 
                             lat = 53.11056, 
                             path = "../data")
plot(DEM_mtRobson) #PC

lat <- 53.11056
lon <- -119.15639
boundary_box <- c(lon -0.1, lon + 0.1, lat - 0.1, lat + 0.1)
DEM_extent <- ext(boundary_box) ### Need a boundary box to restrict dimensions of DEM
DEM_cropped <- crop(DEM_mtRobson, DEM_extent)
plot(DEM_cropped) #PC

# Converting DEM to matrix so it can be used in rayshader package
DEM_values <- values(DEM_cropped, mat = FALSE)
matrix_DEM <- as.matrix(DEM_cropped) 
storage.mode(matrix_DEM) <- "double"
matrix_DEM <- matrix(DEM_values,
                     nrow = nrow(DEM_cropped), 
                     ncol = ncol(DEM_cropped), 
                     byrow = TRUE)

# Generating 3D topographic map of Mt. Robson
matrix_DEM %>%
  sphere_shade(texture = "imhof1") %>%
  plot_3d(matrix_DEM, 
          zscale = 30, 
          solid = TRUE) 

# Was unable to make spheres3d function work with the projection, so saving projection as 2D image
# render_snapshot("../figs/mountain.png") ### image is already saved in figs folder

img_mtRobson <- readPNG("../figs/mountain.png")
graphics.off() ### Clearing Plots pane
grid.raster(img_mtRobson)


x_pixels <- c(1250, 850, 1280) ### Arbitrary x-coordinates chosen to match low, mid, and high elevations
y_pixels <- c(1000, 630, 170) 
point_sizes <- c(50, 30, 20) ### Size of sphere reflects BIN richness
point_colors <- c("red", "orange", "yellow")
grid.points(x = x_pixels, 
            y = y_pixels, 
            pch = 16, 
            size = unit(point_sizes, "mm"), 
            gp = gpar(col = point_colors, cex = 2))
#==> Placement of points is not consistent across different display screens as they use coordinates instead of scaling to the dimensions of the image
# Image saved in figs folder as 3D_visualization_richness

# *** With further exploration of different packages, may be able to create a 3D visualization where spheres are plotted at different elevations to show the elevational diversity gradient




