# Dataset generation to illustrate models from the paper:
# "Regression approaches for modelling genotype-environment interaction and making predictions into a target population of environments"
# using simulated data.
# 2025, University of Hohenheim
# Corresponding author:
# Prof. Hans-Peter Piepho, piepho@uni-hohenheim.de
# Code written by:
# Maksym Hrachov, PhD student, maksym.hrachov@gmail.com

## THIS FILE SHOULD BE CALLED FROM WITHIN "Example models.r"
library(dplyr)
library(tidyverse)
# Define levels
genotypes <- paste0("G", 1:10)
locations <- paste0("L", 1:5)
years <- paste0("Y", 1:4)

# Create full factorial design
df <- expand.grid(G = genotypes, L = locations, Y = years)

# Set seed for reproducibility
set.seed(341323)

# Simulate random effects based on previous information https://doi.org/10.1007/s00122-023-04260-x (Table 5, Monsoon rice)
G_eff     <- rnorm(10, 0, sqrt(0.262))
GY_eff    <- rnorm(10*4, 0, sqrt(0.017))
GL_eff    <- rnorm(10*5, 0, sqrt(0.046))
GLY_eff   <- rnorm(nrow(df), 0, sqrt(0.270))
Y_eff     <- rnorm(4, 0, sqrt(0.017))
L_eff     <- rnorm(5, 0, sqrt(0.134))
LY_eff    <- rnorm(5*4, 0, sqrt(0.491))

# Add random effects to data frame
df <- df %>%
  mutate(
    G_eff   = G_eff[match(G, genotypes)],
    Y_eff   = Y_eff[match(Y, years)],
    L_eff   = L_eff[match(L, locations)],
    GY_eff  = GY_eff[as.integer(factor(paste(G,Y)))],
    GL_eff  = GL_eff[as.integer(factor(paste(G,L)))],
    LY_eff  = LY_eff[as.integer(factor(paste(L,Y)))],
    GLY_eff = GLY_eff
  ) %>%
  mutate(
    YIELD = G_eff + Y_eff + L_eff + GY_eff + GL_eff + LY_eff + GLY_eff
  ) %>% 
  group_by(L,Y) %>% 
  mutate(LY_yield_means = mean(YIELD)) %>% 
  ungroup()
hist(df$YIELD)

### simulate 5 covariates with different level of relatedness

EC1 <- rnorm(20, 0, sqrt(0.111)) + unique(df$LY_yield_means)*0.1
EC2 <- rnorm(20, 0, sqrt(0.222)) + EC1 * 0.2 + unique(df$LY_yield_means)*0.1
EC3 <- rnorm(20, 0, sqrt(0.333)) + EC2 * 0.3 - EC1 * 0.5 + unique(df$LY_yield_means)*0.1
EC4 <- rnorm(20, 0, sqrt(0.555)) + EC3 * 0.8 + EC2 * 0.8 + EC1 * 0.8 + unique(df$LY_yield_means)*0.1
EC5 <- rnorm(20, 0, sqrt(0.777)) - EC4 * 0.8 - EC3 * 0.3 - EC2 * 0.2 - EC1 * 0.7 - unique(df$LY_yield_means)*0.4


df_covariates <- expand.grid(L = locations, Y = years) %>% 
  mutate(EC1 = EC1,
         EC2 = EC2,
         EC3 = EC3,
         EC4 = EC4,
         EC5 = EC5)


df2 <- left_join(df, df_covariates, by = c("L", "Y"))

cor(df2 %>% select(YIELD, EC1:EC5))

rm(list = setdiff(ls(), "df2"))
