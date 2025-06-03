# An unchanged CV code file that was used to produce results for FW2-US model in:
# "Regression approaches for modelling genotype-environment interaction and making predictions into a target population of environments"
# using simulated data.
# Corresponding author:
# Prof. Hans-Peter Piepho, University of Hohenheim, piepho@uni-hohenheim.de
# Code written by:
# Maksym Hrachov, University of Hohenheim, PhD student, maksym.hrachov@gmail.com

# THIS FILE IS NOT MEANT TO BE EXECUTED
# BECAUSE THE DATASET USED IN THE PUBLICATION IS NOT FREELY AVAILABLE RIGHT NOW

# THIS FILE IS HERE TO SHOW HOW THE CV WAS DONE


# Initiation ----
# file_number <- "YYY"
# 
# save_here <- paste0("XXX", file_number, " XXX")
# 
# create_output = T # make TRUE to overwrite the saved files
# 
# if (!dir.exists(save_here) & create_output){
#   dir.create(save_here)
# }
# setwd(save_here)
# libraries 
# source("XXX/XXX Support functions.r")
# source("XXX/XXX The function for variance.r")
# library(asreml)
# library(asremlPlus)
# library(tidyverse)
# library(tictoc)
# library(readxl)
# 
# asreml.options(maxit=100, extra=5, pworkspace="1gb", workspace="1gb", design=TRUE)
# tic()
# 
# get_dataset1 <- read.csv("XXX")
# 
# prep1 <- get_dataset1 %>% 
#   mutate(across(c("G", "Y", "L"), factor)) %>% 
#   mutate(x0 = 0) %>% 
#   mutate(wt = 1/(std.error^2))
# 
# soil.vars <- c("SOL_K1", "SILT1", "SAND1", "SOL_ALB1", "CLAY2") # to exclude them later
# 
# prep6 <- prep1 %>% 
#   mutate(ENV_BRRI = as.factor(paste(Y,L, sep="_")))
# 
# C2 <- prep1 %>% 
#   mutate(ENV = paste(Y,L, sep="_")) %>% 
#   select(!any_of(c("Y", "L", "G", "predicted.value", "std.error", "wt", "X", "Trial.year", "Release.year", "G_names", soil.vars))) %>% 
#   mutate(across(c("ENV"), factor)) %>% 
#   mutate(x0 = 0, x00 = 0) %>% unique()
# 
# # IF THIS MODEL PRODUCES AN ERROR, DOWNLOAD R >4.4.1 
# mod1 <- asreml(
#   fixed = predicted.value ~ G + L + Y + Y:L,
#   random = ~ L:G + Y:G + Y:L:G + str(~mbf(Co):G, vmodel = ~ rr(40,2):id(G)),
#   mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C2")),
#   data = prep6,
#   na.action = na.method(x = "include"),
#   weights = wt,
#   wald = list(denDF="algebraic"),
#   family = asr_gaussian(dispersion=1)
# )
# for(i in 1:20){
#   mod1 <- update(mod1)
# }
# lambdas=summary(mod1)$varcomp[grep(x = rownames(summary(mod1)$varcomp),pattern = "fa",fixed = T),"component"]
# lambda1 = lambdas[1:40]
# lambda2 = lambdas[41:80]
# C <- C2 %>%  select(!any_of(c("x0", "x00", "ENV")))
# z1 <- t(lambda1 %*% t(C))
# z2 <- t(lambda2 %*% t(C))
# 
# z_env <- data.frame(z1 = z1, z2 = z2, ENV_BRRI = unique(prep6$ENV_BRRI))
# 
# prep6_1 <- merge(prep6, z_env, by = "ENV_BRRI")
# 
# vars <- z_env %>% select(z1, z2) %>% names() # here only names of EC are stored
# 
# prep2 <- prep6_1 %>% select(ENV_BRRI, Y, L, G, predicted.value, std.error, wt, all_of(vars))
# 
# values_list <- c()
# year = unique(prep2$Y)[19]
# location = unique(prep2$L)[2]
# asreml.options(Csparse = ~ G + z1:G + z2:G)
# for (year in unique(prep2$Y)){
#   
#   multi_df1 <- prep2 %>% 
#   mutate(across(all_of(vars), ~if_else(Y == year, NA, .))) %>% 
#   select(all_of(vars),L,Y) %>%
#   unique()
#   
#   sink(file=paste0(save_here, '/', file_number, ' ASRemlR messages.txt'), append = TRUE)
#   multivar <- asreml(fixed = cbind(z1,z2) ~  trait + L:trait,
#                      random = ~  Y:us(trait),
#                      residual = ~ units:us(trait), # in theory corresponds to L:Y:us(trait), because L:Y would specify the units
#                      data = multi_df1)
#   
#   multivar.predict <- predict(multivar, classify = "Y:L:trait", 
#                               levels = list(Y = c(paste(year)), L = c(unique(multi_df1$L)), trait = c("z1", "z2")), ignore = c("Y:trait", "units:trait"), vcov = T)
#   sink()
#   
#  for (location in unique(prep2$L)){
#           
#    plug.in <-  multivar.predict$pvals %>% select(Y,L,trait, predicted.value) %>% 
#      pivot_wider(., names_from = "trait", values_from = "predicted.value")
# 
#     prep3 <- prep2 %>%
#       mutate(predicted.value = if_else(Y == year | L == location, NA, predicted.value),
#              across(c(wt, std.error), ~if_else(Y == year | L == location, NA, .)),
#              across(c(all_of(vars)), ~if_else(Y == year, NA, .)),
#              wt = if_else(is.na(wt), 1, wt)) %>% left_join(., plug.in, by = c("Y","L")) %>% 
#       mutate(z1 = coalesce(z1.x, z1.y), z2 = coalesce(z2.x, z2.y)) %>% select(!c(z1.x, z1.y, z2.x, z2.y))
#     
#     sink(file=paste0(save_here, '/', file_number, ' ASRemlR messages.txt'), append = TRUE)
#     
#     model <- asreml(
#       fixed = predicted.value ~ z1 + z2,
#       random = ~ str(~G + z1:G + z2:G, ~us(3):id(G)) + L + Y + Y:L + L:G + Y:G + Y:L:G,
#       data = prep3,
#       na.action = na.method(x = "include"),
#       weights = wt,
#       wald = list(denDF="algebraic"),
#       family = asr_gaussian(dispersion=1)
#     )
#     sink()
#     
#     xi_prime <- extract_xi_for_unseen_location(multivar.predict, new_location = location)
#     
#     sigma_x_prime <- get_sigma_x_prime(multivar,  year = "Y", variables = vars)
#     
#     gammas <- get_gammas_from_the_main_model(model = model, type = "us", slopes = c("z1:G", "z2:G"), intercept = "G")
#     gammas_corected <- gammas$gamma_prime
#     for(i in 1:nrow(gammas$gamma_prime)){
#       gammas_corected[i,] <- gammas$gamma_prime[i,] + t(model$coefficients$fixed)
#     }
#     gammas$gamma_prime <- gammas_corected
#     
#     variances <- prep3$std.error^2
# 
#     C_inverse <- get_C_inv(model, variances, model_terms = "G+z1:G+z2:G", n_terms = 3, fixed.effects = vars)
#     
#     combinations <- expand.grid( vars = c("G_", paste0(vars, ":G_")), pars = paste0(unique(prep3$G))  )
#     column_names <- c(paste0(combinations$vars, combinations$pars))
# 
#     vars_int <- c("(Intercept)",vars)
#     C_subset <- as.matrix(C_inverse)[c(vars_int, column_names), c(vars_int, column_names)]
#     X_yn <- do.call(rbind, rep(list(diag(rep(1, length(vars_int)))), length(unique(prep1$G))))
#     Z_yn <- bdiag(rep(list(diag(rep(1,length(vars_int)))),length(unique(prep1$G))))
#     
#     sandwich <- cbind(X_yn, Z_yn)
#     
#     var_gamma_matrix <- sandwich %*% C_subset %*% t(sandwich)
#     dimnames(var_gamma_matrix)[[1]] <- column_names
#     dimnames(var_gamma_matrix)[[2]] <- column_names
#     
#     
#     names <- rownames(gammas$gamma_prime)
#     var_gamma <- get_var_gamma(var_gamma_matrix, slopes = c("z1", "z2"), intercept = "G", names)
#     
#     est.var <- c()
#     for(i in 1:length(names)){
#       est.var[[i]] <- as.matrix(as.data.frame(gammas[["gamma_prime"]][i, ])) %*% xi_prime[["var_xi_BLUE"]] %*% t(as.matrix(gammas[["gamma_prime"]][i, ])) +
#         as.matrix(xi_prime[["xi_prime_BLUE"]]) %*% var_gamma[[i]] %*% t(as.matrix(xi_prime[["xi_prime_BLUE"]])) -
#           sum(diag(var_gamma[[i]] %*% xi_prime[["var_xi_BLUE"]]))
#     }
#     
#     vi_i_x_prime <- c()
#     for(i in 1:length(names)){
#     vi_i_x_prime[[i]] <- as.matrix(as.data.frame(gammas[["gamma_prime"]][i, ])) %*% as.matrix(sigma_x_prime[["sigma_x_prime"]]) %*% t(as.matrix(gammas[["gamma_prime"]][i, ])) -
#       sum(diag(var_gamma[[i]] %*% sigma_x_prime[["sigma_x_prime"]]))
#     }
#     
#     non_vi_r <- summary(model)$varcomp[grep(x = rownames(summary(model)$varcomp),pattern = "!",fixed = T),"component"]
#     vi_r <- sum(setdiff(summary(model)$varcomp$component, non_vi_r))
#     
#     ###########################################
#     
#     I <- length(unique(prep3$G))
#     ones <- c(rep(1, I))
#     J <- ones %*% t(ones)
#     K <- 1/I * ones %*% t(ones)
#     P <- 2*(1/(I - 1)) * (diag(rep(1,I))-K)
#     
#     gamma_prime_row <- as.matrix(c(t(gammas$gamma_prime)))
# 
#      
#     mean_vd1 <- sum(diag( P %x% ( xi_prime[["var_xi_BLUE"]] + sigma_x_prime[["sigma_x_prime"]] ) %*% gamma_prime_row %*% t(gamma_prime_row) ) )
#     mean_vd2 <- sum(diag( var_gamma_matrix %*% ( P %x% (t(xi_prime[["xi_prime_BLUE"]]) %*% xi_prime[["xi_prime_BLUE"]]) ) ))
#     mean_vd3 <- sum(diag( (var_gamma_matrix) %*%
#                           (as.matrix( P %x% ( xi_prime[["var_xi_BLUE"]] + sigma_x_prime[["sigma_x_prime"]] ) ))    ))
#     
#     mean_v_r <- summary(model)$varcomp %>% mutate(name = rownames(.)) %>% 
#       filter(name %in% c("L:G", "Y:G", "Y:L:G")) %>% select(component) %>% sum()
#     
#     total_mean_vd <- mean_vd1 + mean_vd2 - mean_vd3 + 2*mean_v_r
# 
#     ###########################################
#     
#     deleted_rows <- prep2 %>%
#       filter(Y == year & L == location) %>%
#       select(G, L, Y) %>%
#       summarise(across(everything(), ~list(unique(.))))
#      
#     EC_TPE_average <- xi_prime[["xi_prime_BLUE"]] %>%
#       as.data.frame() %>% 
#       select(all_of(vars)) %>%
#       summarise(across(everything(), ~list(mean(.))))
#     
#     deleted_rows2 <- cbind(deleted_rows, EC_TPE_average)
#     
#     my_levels <- setNames(lapply(deleted_rows2, unlist), names(deleted_rows2))
#     
#     sink(file=paste0(save_here, '/', file_number, ' ASRemlR messages.txt'), append = TRUE)
#     values_list[[paste(year)]][[paste(location)]] <- predict(model,
#                                                              classify = paste(names(my_levels), collapse = ":"),
#                                                              levels = my_levels, maxit = 2)
#     sink()
# 
#     vector_1 <- prep2 %>%
#       filter(Y == year & L == location) %>% 
#       select(G, predicted.value) 
#     
#     vector2 <- values_list[[paste(year)]][[paste(location)]][["pvals"]]
#     
#     intersection <- semi_join(vector2, vector_1, by = c("G"))
#     
#     values_list[[paste(year)]][[paste(location)]][["corr"]] <- cor(intersection$predicted.value, vector_1$predicted.value)
#     values_list[[paste(year)]][[paste(location)]][["msepd"]] <- get_MSEPD(observedValues = vector_1$predicted.value, 
#                                                                           predictedValues = intersection$predicted.value)
#     values_list[[paste(year)]][[paste(location)]][["mspe"]] <- get_MSPE(observedValues = vector_1$predicted.value, 
#                                                                         predictedValues = intersection$predicted.value)
#     
#     values_list[[paste(year)]][[paste(location)]][["var.pred.diff"]] <- total_mean_vd
#     
#     for(i in 1:length(names)){
#       values_list[[paste(year)]][[paste(location)]][["var.pred.list"]][[i]] <- est.var[[i]] + vi_i_x_prime[[i]] + vi_r
#     }
#         
#     values_list[[paste(year)]][[paste(location)]][["var.pred"]] <- sumirize_var.pred(values_list[[paste(year)]][[paste(location)]][["var.pred.list"]], vector_1)
#     print(paste(year, location, "r2:", round(values_list[[paste(year)]][[paste(location)]][["corr"]],3)))
#     rm("model", "intersection")
#   }
#   
#   if(create_output){
#     save.image(file=paste0(save_here, '/', file_number, ' Environment.RData'))
#     saveRDS(values_list, file=paste0(save_here, '/', file_number, ' values_list.rds'))
#   }
# }
# 
# time <- toc()
# if(create_output){
#   save.image(file=paste0(save_here, '/', file_number, ' Environment.RData'))
# }
