# Models from the paper:
# "Regression approaches for modelling genotype-environment interaction and making predictions into a target population of environments".
# With simulated data.
# 2025, University of Hohenheim
# Corresponding author:
# Prof. Hans-Peter Piepho, piepho@uni-hohenheim.de
# Code written by:
# Maksym Hrachov, PhD student, maksym.hrachov@gmail.com

# source the right file with the script for the dataset simulation
source("Dataset simulation.r")

library(asreml)
# Increase the total number of iterations and force the model to do a few extra iterations
asreml.options(maxit = 50, extra = 2)

data <- df2 %>% mutate(x0 = 0, x00 = 0, z0 = 1,
                       ENV = as.factor(paste(L, Y, sep = "_")))

## Baseline model without fixed main effect for EC:
# Every further model can also be fitted without the main effect for EC.
# This will not affect predictions much.
# However, the variance of prediction, as outlined in the paper, needs the main EC effect.
baseline <- asreml(
  fixed = YIELD ~ 1,
  random = ~ G + L + Y + Y:L +L:G + Y:G,
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)
# summary(baseline)$varcomp

## Baseline model with fixed main effect for EC
baseline_EC <- asreml(
  fixed = YIELD ~ EC1 + EC2 + EC3 + EC4 + EC5,
  random = ~ G + L + Y + Y:L +L:G + Y:G,
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)
# summary(baseline_EC)$varcomp

## Reduced Rank Regression (RRR) of rank one
# This model in ASReml requires a dummy variable of x=0 for each level of RR.
# Here, rr(1) means x0=0 has to be fit so that the model converges.
# It seems to be something special about the str(), not rr() itself.
# Refer to the manual if you need more info on str() and rr().

rr1 <- asreml(
  fixed = YIELD ~ EC1 + EC2 + EC3 + EC4 + EC5,
  random = ~str(~ G + G:EC1 + G:EC2+ G:EC3 + G:EC4 + G:EC5 + G:x0, vmodel = ~ rr(6,1):id(G)) + L + Y + Y:L +L:G + Y:G,
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)
# summary(rr1)$varcomp

## Reduced Rank Regression (RRR) of rank two
# For rr2 mode put the genotype on the last place using the dummy z0=1 variable.
# One of the loadings has to be fixed at 0 for rr(2), and this will be the first loading in ASReml.
# Thus, we want to avoid it to be the genotype main effect related loading
rr2 <- asreml(
  fixed = YIELD ~ EC1 + EC2 + EC3 + EC4 + EC5,
  random = ~str(~ G:EC1 + G:EC2+ G:EC3 + G:EC4 + G:EC5 + G:z0 + G:x0 + G:x00, vmodel = ~ rr(6,2):id(G)) + L + Y + Y:L +L:G + Y:G,
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)
# summary(rr1)$varcomp

## Unstructured regression
# The USR model needs many more iterations to converge. 
# Please use the while() loop to get to stable LogLik values.
# Keep in mind that it can take significant time (a few minutes with 6-digit precision)
# If you inspect the "us" object, "converge = FALSE".
# This model will take a lot of time to converge with real-world data in case the number of covariates is large.
# us_uptd$loglik > 5000 exemplifies how to deal with sudden LogLik jumps that can happen in practice due to overfitting
us <- asreml(
  fixed = YIELD ~ EC1 + EC2 + EC3 + EC4 + EC5,
  random = ~str(~ G + G:EC1 + G:EC2+ G:EC3 + G:EC4 + G:EC5 , vmodel = ~ us(6):id(G)) + L + Y + Y:L +L:G + Y:G,
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)
convergence <- TRUE
while(convergence){
  try(us_uptd <- update(us))
  if( us$loglik > us_uptd$loglik | us_uptd$loglik > 5000 ){
    rm("us_uptd", "h1")
    print("Loop exit criteria were met")
    break()
  } else if(exists("h1") &&  round(us_uptd$loglik,6) == round(h1$loglik,6)){
    us <- us_uptd
    rm("us_uptd", "h1")
    print("Convergence was reached")
    break()
  } else {
    h1 <- us
    us <- us_uptd
    rm("us_uptd")
  }
}
# summary(rr1)$varcomp

## Extended Finlay-Wilkinson regression (FW1-US)
# See "Extending Finlay–Wilkinson regression with environmental covariates" by Piepho and Blancon (2023)
# The model from respective paper is slightly modified. Here, it includes random genotype-location and genotype-year deviations.
# The model needs to be fitted in two steps. 
# In the simulated example the first step in both FW1-US and FW2-US has a very flat convergence curve.
# This should not be a problem for extracting the synthetic covariate(s).
FW_1_rr_step <- asreml(
  fixed = YIELD ~ G + L + Y + Y:L,
  random = ~ L:G + Y:G + str(~ G:EC1 + G:EC2+ G:EC3 + G:EC4 + G:EC5 + G:x0, vmodel = ~ rr(5,1):id(G)),
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)

lambda <- summary(FW_1_rr_step)$varcomp[grep(x = rownames(summary(FW_1_rr_step)$varcomp),pattern = "fa",fixed = T),"component"]
# Computing z1 and merging
C <- data %>% select(EC1:EC5) %>% unique
z1 <- t(lambda %*% t(C))

z1_env <- data.frame(z1 = z1, ENV = unique(data$ENV))

data_augmented_fw1 <- merge(data, z1_env, by = "ENV")

FW_1_us_step <- asreml(
  fixed = YIELD ~ z1,
  random = ~str(~ G + G:z1, vmodel = ~ us(2):id(G)) + L + Y + Y:L +L:G + Y:G,
  data = data_augmented_fw1,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)
# summary(FW_1_us_step)$varcomp

## Extended Finlay-Wilkinson regression (FW2-US)
# With 3 components for the unstructured variance-covariance this model need many more iterations for the simulated data.
# The LogLik reaches a plateau after 10000 iterations and improves only slightly. 
FW_2_rr_step <- asreml(
  fixed = YIELD ~ G + L + Y + Y:L,
  random = ~ L:G + Y:G + str(~ G:EC1 + G:EC2+ G:EC3 + G:EC4 + G:EC5 + G:x0 + G:x00, vmodel = ~ rr(5,2):id(G)),
  data = data,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic")
)

lambdas <- summary(FW_2_rr_step)$varcomp[grep(x = rownames(summary(FW_2_rr_step)$varcomp),pattern = "fa",fixed = T),"component"]
lambda1 = lambdas[1:5]
lambda2 = lambdas[6:10]
# Computing z1 and z2 and merging
C <- data %>% select(EC1:EC5) %>% unique
z1 <- t(lambda1 %*% t(C))
z2 <- t(lambda2 %*% t(C))

z_env <- data.frame(z1 = z1, z2 = z2, ENV = unique(data$ENV))

data_augmented_fw2 <- merge(data, z_env, by = "ENV")

FW_2_us_step <- asreml(
  fixed = YIELD ~ z1 + z2,
  random = ~str(~ G + G:z1 + G:z2, vmodel = ~ us(3):id(G)) + L + Y + Y:L +L:G + Y:G,
  data = data_augmented_fw2,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic"),
  maxit = 10000
)
# summary(FW_1_us_step)$varcomp

library(asremlPlus)
infoCriteria(baseline_EC)
infoCriteria(rr1)
infoCriteria(rr2)
infoCriteria(FW_1_us_step)
infoCriteria(FW_2_us_step)


## How to get variance of prediction

# Unfortunately, ASRemlR version 4.2.0.302 does not provide complete C-inverse matrix to extract variance-covariance of coefficients.
# Hence, it has to be built by hand.
# The large function was put in another file, and it is recommended to go through it in order to check if it does what it has for the model that you specified.
# The functions in the sourced file are adapted solely for the purpose of the paper and were not tested for the ASReml output of ther configurations.
source("Support functions for the variance.r")

# As for prediction scenario, let's consider an unknown enviroment, just like in the paper (corresponds to LYLO CV)
year <- "Y1"
location <- "L1"
vars <- c("z1", "z2")
multi_df1 <- data_augmented_fw2 %>% 
  mutate(across(all_of(vars), ~if_else(Y == year, NA, .))) %>% 
  select(all_of(vars),L,Y) %>%
  unique()

multivar <- asreml(fixed = cbind(z1,z2) ~  trait + L:trait,
                   random = ~  Y:us(trait),
                   residual = ~ units:us(trait), # in theory corresponds to L:Y:us(trait), because L:Y would specify the units
                   data = multi_df1)

multivar.predict <- predict(multivar, classify = "Y:L:trait", 
                            levels = list(Y = c(paste(year)), L = c(unique(multi_df1$L)), trait = c("z1", "z2")), 
                            ignore = c("Y:trait", "units:trait"),
                            maxit = 2,
                            vcov = T)

plug.in <-  multivar.predict$pvals %>% select(Y,L,trait, predicted.value) %>% 
  pivot_wider(., names_from = "trait", values_from = "predicted.value")

data_fw2_for_prediction <- data_augmented_fw2 %>%
  mutate(YIELD = if_else(Y == year | L == location, NA, YIELD),
         across(c(all_of(vars)), ~if_else(Y == year, NA, .))
         ) %>% 
  left_join(., plug.in, by = c("Y","L")) %>% 
  mutate(z1 = coalesce(z1.x, z1.y), z2 = coalesce(z2.x, z2.y)) %>% select(!c(z1.x, z1.y, z2.x, z2.y))

asreml.options(design = TRUE)

model <- asreml(
  fixed = YIELD ~ z1 + z2,
  random = ~str(~ G + z1:G + z2:G, vmodel = ~ us(3):id(G)) + L + Y + Y:L +L:G + Y:G,
  data = data_fw2_for_prediction,
  na.action = na.method(x = "include"),
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1), # see the note below
  maxit = 10000
)
# Note: For this dataset ASRemlR will internally set the residual variance to 1, and the variance components in the "model" object will be scaled accordingly
# if we call `summary(model)$varcomp`, then we get variance components on the original scale.
# The functions I have written work with by default with the "model" object, and not with the summary().
# Also, in a two-stage analysis that was implemented in the paper the option
# family = asr_gaussian(dispersion=1)
# is required to fit the variances from the 1st stage correctly.

xi_prime <- extract_xi_for_unseen_location(multivar.predict, new_location = location)

sigma_x_prime <- get_sigma_x_prime(multivar,  year = "Y", variables = vars)

gammas <- get_gammas_from_the_main_model(model = model, type = "us", slopes = c("z1:G", "z2:G"), intercept = "G")
# Eq. 3,4,5 show that elements like gamma can be dissected into an intercept (e.g. z1) and a deviation from it (e.g. z1:G),
# thus the "correction" means that the intercept value is added on top of all related deviations
gammas_corected <- gammas$gamma_prime  
for(i in 1:nrow(gammas$gamma_prime)){
  gammas_corected[i,] <- gammas$gamma_prime[i,] + t(model$coefficients$fixed)
}
gammas$gamma_prime <- gammas_corected

R.variances <- as.vector(rep(model[["vparameters"]]["units!R"], 
                 length(data_fw2_for_prediction$YIELD[!is.na(data_fw2_for_prediction$YIELD)]))) # if you look into the "model" object's R.param part, the units were set to 1 and variance components are scaled accordingly

C_inverse <- get_C_inv(model, R.variances, model_terms = "G+z1:G+z2:G", n_terms = 3, fixed.effects = vars)

combinations <- expand.grid( vars = c("G_", paste0(vars, ":G_")), pars = paste0(unique(data_fw2_for_prediction$G))  )
column_names <- c(paste0(combinations$vars, combinations$pars))

vars_int <- c("(Intercept)",vars)
C_subset <- as.matrix(C_inverse)[c(vars_int, column_names), c(vars_int, column_names)]
X_yn <- do.call(rbind, rep(list(diag(rep(1, length(vars_int)))), length(unique(data_augmented_fw2$G))))
Z_yn <- bdiag(rep(list(diag(rep(1,length(vars_int)))),length(unique(data_augmented_fw2$G))))

# the operations below now combines the variance assosiated with the intercepts and the deviations
sandwich <- cbind(X_yn, Z_yn)

var_gamma_matrix <- sandwich %*% C_subset %*% t(sandwich)
dimnames(var_gamma_matrix)[[1]] <- column_names
dimnames(var_gamma_matrix)[[2]] <- column_names


names <- rownames(gammas$gamma_prime)
var_gamma <- get_var_gamma(var_gamma_matrix, slopes = c("z1", "z2"), intercept = "G", names)

# the part below is calculating the parts of the total prediciton variance of each genotype
est.var <- c() # Eq.44
for(i in 1:length(names)){
  est.var[[i]] <- as.matrix(as.data.frame(gammas[["gamma_prime"]][i, ])) %*% xi_prime[["var_xi_BLUE"]] %*% t(as.matrix(gammas[["gamma_prime"]][i, ])) +
    as.matrix(xi_prime[["xi_prime_BLUE"]]) %*% var_gamma[[i]] %*% t(as.matrix(xi_prime[["xi_prime_BLUE"]])) -
    sum(diag(var_gamma[[i]] %*% xi_prime[["var_xi_BLUE"]]))
}

phi_i_x_prime <- c() # Eq.45
for(i in 1:length(names)){
  phi_i_x_prime[[i]] <- as.matrix(as.data.frame(gammas[["gamma_prime"]][i, ])) %*% as.matrix(sigma_x_prime[["sigma_x_prime"]]) %*% t(as.matrix(gammas[["gamma_prime"]][i, ])) -
    sum(diag(var_gamma[[i]] %*% sigma_x_prime[["sigma_x_prime"]]))
}

non_vi_r <- summary(model)$varcomp[grep(x = rownames(summary(model)$varcomp),pattern = "!",fixed = T),"component"]
vi_r <- sum(setdiff(summary(model)$varcomp$component, non_vi_r)) # related to Eq. 46-47

estimated_variance_of_each_genotype <- c()
for(i in 1:length(names)){
  estimated_variance_of_each_genotype[[i]] <- est.var[[i]] + phi_i_x_prime[[i]] + vi_r
}

## for the pairwise differences (section 3.6)
I <- length(unique(data_fw2_for_prediction$G))
ones <- c(rep(1, I))
J <- ones %*% t(ones)
K <- 1/I * ones %*% t(ones)
P <- 2*(1/(I - 1)) * (diag(rep(1,I))-K)

gamma_prime_row <- as.matrix(c(t(gammas$gamma_prime)))


mean_vd1 <- sum(diag( P %x% ( xi_prime[["var_xi_BLUE"]] + sigma_x_prime[["sigma_x_prime"]] ) %*% gamma_prime_row %*% t(gamma_prime_row) ) )
mean_vd2 <- sum(diag( var_gamma_matrix %*% ( P %x% (t(xi_prime[["xi_prime_BLUE"]]) %*% xi_prime[["xi_prime_BLUE"]]) ) ))
mean_vd3 <- sum(diag( (var_gamma_matrix) %*%
                        (as.matrix( P %x% ( xi_prime[["var_xi_BLUE"]] + sigma_x_prime[["sigma_x_prime"]] ) ))    ))

mean_v_r <- summary(model)$varcomp %>% mutate(name = rownames(.)) %>% 
  filter(name %in% c("L:G", "Y:G", "Y:L:G")) %>% select(component) %>% sum()

estimated_variance_of_the_pairwise_difference <- mean_vd1 + mean_vd2 - mean_vd3 + 2*mean_v_r


## Compare with the conventional ASRemlR
predict.normal <- predict.asreml(model, classify = "G:L:Y:z1:z2", maxit = 1)$pvals %>% 
  filter(L == location, Y == year) %>% 
  mutate(var.pred = std.error^2)


# estimated_variance_of_the_pairwise_difference
# estimated_variance_of_each_genotype

# The var.pred from ASReml is the same for each genotype, and hand-coded est.var is not.
# This reflects the contribution of genotype-specific slopes to overall uncertainty.
