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
# prep0 <- get_dataset1 %>% 
#   mutate(across(c("G", "Y", "L"), factor)) %>% 
#   mutate(wt = 1/(std.error^2)) %>% 
#   mutate(ENV_BRRI = as.factor(paste(Y,L, sep="_")))

# stock_vars <- prep0 %>% 
#   select(!any_of(c("wt", "ENV_BRRI", colnames(.)[1:8]))) %>% 
#   names()

# rescaled_vars <- prep0 %>% 
#   select(L,Y,all_of(stock_vars)) %>% 
#   distinct() %>% 
#   mutate(across(all_of(stock_vars), ~scale(.)[,1])) %>% 
#   data.frame()

# prep1 <- prep0 %>% 
#   select(-all_of(stock_vars)) %>% 
#   left_join(., rescaled_vars, by = c("L", "Y"))

# C1 <- prep1 %>% 
#   rename(ENV = ENV_BRRI) %>% 
#   mutate(x0 = 0, x00 = 0) %>% 
#   select(ENV, all_of(stock_vars), x0, x00) %>% 
#   unique() %>% 
#   mutate(across(all_of(stock_vars), ~scale(.)[,1]))

# random1 <- as.formula(paste0("~L:G + Y:G + Y:L:G + str(~mbf(Co):G, vmodel = ~ rr(", length(stock_vars),
#                              ",2):id(G))"))

# RRR_observed <- asreml(
#   fixed = predicted.value ~ G + L + Y + Y:L,
#   random = random1,
#   mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C1")),
#   data = prep1,
#   na.action = na.method(x = "include"),
#   weights = wt,
#   wald = list(denDF="algebraic"),
#   family = asr_gaussian(dispersion=1)
# )

# RRR_starting <- asreml(
#   fixed = predicted.value ~ G + L + Y + Y:L,
#   random = random1,
#   mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C1")),
#   data = prep1,
#   na.action = na.method(x = "include"),
#   weights = wt,
#   wald = list(denDF="algebraic"),
#   family = asr_gaussian(dispersion=1),
#   G.param = RRR_observed$G.param,
#   start.values = T
# )
# start_values <- RRR_starting$vparameters.table %>% 
#   mutate(Constraint = if_else(Component %in% c("L:G", "Y:G", "Y:L:G"), "P", Constraint))

# rm("C1", "RRR_observed", "RRR_starting")

# values_list <- c()
# year = unique(prep1$Y)[1]
# location = unique(prep1$L)[1]

# for (year in unique(prep1$Y)){
#  for (location in unique(prep1$L)){
          
#    prep1_na <- prep1 %>%
#      mutate(predicted.value = if_else(Y == year | L == location, NA, predicted.value),
#             across(c(wt, std.error), ~if_else(Y == year | L == location, NA, .)),
#             across(c(all_of(stock_vars)), ~if_else(Y == year, NA, .)),
#             wt = if_else(is.na(wt), 1, wt)) 
   
#    rescaled_vars <- prep1_na %>%
#      select(L,Y,all_of(stock_vars)) %>%
#      distinct() %>%
#      mutate(across(all_of(stock_vars), ~scale(., )[,1])) %>%
#      data.frame()
   
#    prep1_1r <- prep1_na %>% 
#      select(-all_of(stock_vars)) %>%
#      left_join(., rescaled_vars, by = c("L", "Y"))
   
#    plug.in <- all_univar_BLUPs[[year]] %>% select(Y,L, predicted.value, EC) %>% 
#      pivot_wider(names_from = EC, values_from = predicted.value)
   
#    prep1.1 <- prep1_1r %>% 
#      left_join(plug.in, by = c("Y", "L"), suffix = c(".xtab", ".ytab")) %>%
#      mutate(ENV_BRRI = as.factor(paste(Y,L, sep="_")), 
#             Y = as.factor(Y)) %>% 
#      mutate(across(ends_with(".xtab"), 
#                    .fns = ~ coalesce(., get(sub(".xtab", ".ytab", cur_column()))),
#                    .names = "{sub('.xtab', '', .col)}")) %>%
#      select(-ends_with(".xtab"), -ends_with(".ytab")) %>% 
#      select(G,L,Y, predicted.value, std.error, wt, ENV_BRRI, all_of(stock_vars)) %>%
#      arrange(Y, L, G)
    
#     C2 <- prep1.1 %>% 
#       rename(ENV = ENV_BRRI) %>% 
#       select(ENV, all_of(stock_vars)) %>% 
#       mutate(across(c("ENV"), factor)) %>% 
#       mutate(x0 = 0, x00 = 0) %>% unique()
    
#     sink(file=tempfile(), append = TRUE)
#     try(mod1 <- asreml(
#       fixed = predicted.value ~ G + L + Y + Y:L,
#       random = random1,
#       mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C2")),
#       data = prep1.1,
#       na.action = na.method(x = "include"),
#       weights = wt,
#       wald = list(denDF="algebraic"),
#       family = asr_gaussian(dispersion=1),
#       trace = F
#     ), silent = F)
#     sink()
    
#     if(!exists("mod1")){
#       sink(file=tempfile(), append = TRUE)
#       try(mod1 <- asreml(
#         fixed = predicted.value ~ G + L + Y + Y:L,
#         random = random1,
#         mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C2")),
#         data = prep1.1,
#         na.action = na.method(x = "include"),
#         weights = wt,
#         wald = list(denDF="algebraic"),
#         family = asr_gaussian(dispersion=1),
#         G.param = start_values,
#         trace = F
#       )
#       )
#       sink()
#       if(!exists("mod1")){stop("Starting values didn't help. Inspect the fold")} else {
#         print("INFO: Starting values were successfully used.")
#       }
#     }
    
#     convergence <- TRUE
#     while(convergence){
#       sink(file=tempfile(), append = TRUE)
#       try(model_uptd <- update(mod1, maxit = 5, trace = F))
#       sink()
#       if( mod1$loglik > model_uptd$loglik | model_uptd$loglik > 5000 ){
#         rm("model_uptd", "h1")
#         print("Loop exit criteria were met")
#         break()
#       } else if(exists("h1") &&  round(model_uptd$loglik,5) == round(h1$loglik,5)){
#         mod1 <- model_uptd
#         rm("model_uptd", "h1")
#         print("Convergence was reached")
#         break()
#       } else {
#         h1 <- mod1
#         mod1 <- model_uptd
#         rm("model_uptd")
#       }
#     }
    
#     lambdas=summary(mod1)$varcomp[grep(x = rownames(summary(mod1)$varcomp),pattern = "fa",fixed = T),"component"]
#     lambda1 = lambdas[1:length(stock_vars)]
#     lambda2 = lambdas[c(length(stock_vars)+1):c(length(stock_vars)*2)]
#     # Computing z1 and merging
#     C <- C2 %>%  select(!any_of(c("x0", "x00", "ENV")))
#     z1 <- t(lambda1 %*% t(C))
#     z2 <- t(lambda2 %*% t(C))
    
#     z_env <- data.frame(z1 = z1, z2 = z2, ENV_BRRI = unique(prep1.1$ENV_BRRI))
    
#     prep1.1_1 <- merge(prep1.1, z_env, by = "ENV_BRRI")
    
#     rm("mod1")

#     vars <- c("z1", "z2") # here only names of EC are stored
    
#     prep3 <- prep1.1_1 %>% 
#       select(-all_of(stock_vars)) %>% 
#       mutate(predicted.value = if_else(Y == year | L == location, NA, predicted.value),
#              across(c(wt), ~if_else(Y == year | L == location, NA, .)),
#              wt = if_else(is.na(wt), 1, wt))
    
#     fixed_part <- paste0("predicted.value ~ ", paste(vars, collapse = " + "), collapse = " ")
#     fixed_formula <- formula(paste(fixed_part))

#     model <- asreml(
#       fixed = fixed_formula,
#       random = ~ str(~G + z1:G + z2:G, ~us(3):id(G)) + L + Y + Y:L + L:G + Y:G + Y:L:G,
#       data = prep3,
#       na.action = na.method(x = "include"),
#       weights = wt,
#       wald = list(denDF="algebraic"),
#       family = asr_gaussian(dispersion=1),
#       trace = F
#     )

#     multi_df1 <- prep1.1_1 %>% 
#       mutate(across(all_of(vars), ~if_else(Y == year, NA, .))) %>% 
#       select(all_of(vars),L,Y) %>%
#       distinct()
    
#     # see part on bivariate models here https://tomhouslay.com/wp-content/uploads/2017/02/indivvar_mv_tutorial_asreml.pdf
#     # for mor information on how to fit multivariate models
#     # trait is an internal ASReml keyword which specifies that we have a multi-trait situation
    
#     multivar <- asreml(fixed = cbind(z1,z2) ~  trait + L:trait,
#                       random = ~  Y:us(trait),
#                       residual = ~ units:us(trait), # in theory corresponds to L:Y:us(trait), because L:Y would specify the units
#                       data = multi_df1, workspace = "0.2gb", trace = F)
    
#     multivar.predict <- predict(multivar, classify = "Y:L:trait", 
#                                 levels = list(Y = c(paste(year)), L = c(unique(multi_df1$L)), trait = c("z1", "z2")), ignore = c("Y:trait", "units:trait"), vcov = T)
   
    
#     if (any(round(c(multivar.predict$pvals %>% filter(L==location, Y == year))$predicted.value,6) != round(c(prep3 %>% 
#                                                                            filter(L==location, Y == year) %>% 
#                                                                            select(z1,z2) %>% distinct() %>% unlist()),6))) {
#       stop("Predicted Z values != estimated means.")
#     }
    
#     ####### variance of prediction ########
    
#     xi_prime <- extract_xi_for_unseen_location(multivar.predict, new_location = location)
    
#     sigma_x_prime <- get_sigma_x_prime(multivar,  year = "Y", variables = vars)
      
#     gammas <- get_gammas_from_the_main_model(model = model, type = "us", slopes = c("z1:G", "z2:G"), intercept = "G")
#     gammas_corected <- gammas$gamma_prime
#     for(i in 1:nrow(gammas$gamma_prime)){
#       gammas_corected[i,] <- gammas$gamma_prime[i,] + t(model$coefficients$fixed)
#     }
#     gammas$gamma_prime <- gammas_corected
    
#     variances <- prep3$std.error^2
#     #variances[which(is.na(variances))] <- 1
    
#     C_inverse <- get_C_inv(model, variances, model_terms = "G+z1:G+z2:G", n_terms = 3, fixed.effects = vars)
    
#     combinations <- expand.grid( vars = c("G_", paste0(vars, ":G_")), pars = paste0(unique(prep3$G))  )
#     column_names <- c(paste0(combinations$vars, combinations$pars))

#     vars_int <- c("(Intercept)",vars)
#     C_subset <- as.matrix(C_inverse)[c(vars_int, column_names), c(vars_int, column_names)]
#     X_yn <- do.call(rbind, rep(list(diag(rep(1, length(vars_int)))), length(unique(prep1$G))))
#     Z_yn <- bdiag(rep(list(diag(rep(1,length(vars_int)))),length(unique(prep1$G))))
    
#     sandwich <- cbind(X_yn, Z_yn)
    
#     var_gamma_matrix <- sandwich %*% C_subset %*% t(sandwich)
#     dimnames(var_gamma_matrix)[[1]] <- column_names
#     dimnames(var_gamma_matrix)[[2]] <- column_names
    
    
#     names <- rownames(gammas$gamma_prime)
#     var_gamma <- get_var_gamma(var_gamma_matrix, slopes = c("z1", "z2"), intercept = "G", names)
    
#     est.var <- c()
#     for(i in 1:length(names)){
#       est.var[[i]] <- as.matrix(as.data.frame(gammas[["gamma_prime"]][i, ])) %*% xi_prime[["var_xi_BLUE"]] %*% t(as.matrix(gammas[["gamma_prime"]][i, ])) +
#         as.matrix(xi_prime[["xi_prime_BLUE"]]) %*% var_gamma[[i]] %*% t(as.matrix(xi_prime[["xi_prime_BLUE"]])) -
#           sum(diag(var_gamma[[i]] %*% xi_prime[["var_xi_BLUE"]]))
#     }
    
#     vi_i_x_prime <- c()
#     for(i in 1:length(names)){
#     vi_i_x_prime[[i]] <- as.matrix(as.data.frame(gammas[["gamma_prime"]][i, ])) %*% as.matrix(sigma_x_prime[["sigma_x_prime"]]) %*% t(as.matrix(gammas[["gamma_prime"]][i, ])) -
#       sum(diag(var_gamma[[i]] %*% sigma_x_prime[["sigma_x_prime"]]))
#     }
    
#     vi_r_select <- summary(model)$varcomp %>%
#     tibble::rownames_to_column("items") %>%
#     dplyr::filter(!grepl("!", items, fixed = TRUE))

#     vi_r <- vi_r_select$component %>% sum()
    
#     ## pairwise differences
    
#     I <- length(unique(prep3$G))
#     ones <- c(rep(1, I))
#     J <- ones %*% t(ones)
#     K <- 1/I * ones %*% t(ones)
#     P <- 2*(1/(I - 1)) * (diag(rep(1,I))-K)
    
#     gamma_prime_row <- as.matrix(c(t(gammas$gamma_prime)))

     
#     mean_vd1 <- sum(diag( P %x% ( xi_prime[["var_xi_BLUE"]] + sigma_x_prime[["sigma_x_prime"]] ) %*% gamma_prime_row %*% t(gamma_prime_row) ) )
#     mean_vd2 <- sum(diag( var_gamma_matrix %*% ( P %x% (t(xi_prime[["xi_prime_BLUE"]]) %*% xi_prime[["xi_prime_BLUE"]]) ) ))
#     mean_vd3 <- sum(diag( (var_gamma_matrix) %*%
#                           (as.matrix( P %x% ( xi_prime[["var_xi_BLUE"]] + sigma_x_prime[["sigma_x_prime"]] ) ))    ))
    
#     mean_v_r <- summary(model)$varcomp %>% mutate(name = rownames(.)) %>% 
#       filter(name %in% c("L:G", "Y:G", "Y:L:G")) %>% select(component) %>% sum()
    
#     total_mean_vd <- mean_vd1 + mean_vd2 - mean_vd3 + 2*mean_v_r

#     ###########################################  
#     # prepare data for prediction
    
#     deleted_rows <- prep3 %>%
#       filter(Y == year & L == location) %>%
#       select(G, L, Y) %>%
#       summarise(across(everything(), ~list(unique(.))))
     
 #    EC_values_prediction <- xi_prime[["xi_prime_BLUE"]] %>%
 #      as.data.frame() %>% 
 #      select(all_of(vars)) %>%
 #      summarise(across(everything(), ~list(.)))
    
 #    deleted_rows2 <- cbind(deleted_rows, EC_values_prediction)
    
#     my_levels <- setNames(lapply(deleted_rows2, unlist), names(deleted_rows2))
    
#     values_list[[paste(year)]][[paste(location)]] <- predict(model,
#                                                              classify = paste(names(my_levels), collapse = ":"),
#                                                              levels = my_levels, maxit = 2, trace = F)
    
#     # prepare data for precision values calculation
#     vector_1 <- prep1 %>%
#       filter(Y == year & L == location) %>% 
#       select(G, predicted.value) 
    
#     vector2 <- values_list[[paste(year)]][[paste(location)]][["pvals"]]
    
#     intersection <- semi_join(vector2, vector_1, by = c("G"))
    
#     values_list[[paste(year)]][[paste(location)]][["corr"]] <- cor(intersection$predicted.value, vector_1$predicted.value)
#     values_list[[paste(year)]][[paste(location)]][["msepd"]] <- get_MSEPD(observedValues = vector_1$predicted.value, 
#                                                                           predictedValues = intersection$predicted.value)
#     values_list[[paste(year)]][[paste(location)]][["mspe"]] <- get_MSPE(observedValues = vector_1$predicted.value, 
#                                                                         predictedValues = intersection$predicted.value)
    
#     values_list[[paste(year)]][[paste(location)]][["var.pred.diff"]] <- total_mean_vd
    
#     for(i in 1:length(names)){
#       values_list[[paste(year)]][[paste(location)]][["var.pred.list"]][[i]] <- est.var[[i]] + vi_i_x_prime[[i]] + vi_r
#     }
        
#     values_list[[paste(year)]][[paste(location)]][["var.pred"]] <- sumirize_var.pred(values_list[[paste(year)]][[paste(location)]][["var.pred.list"]], vector_1)
#     print(paste(year, location, "r2:", round(values_list[[paste(year)]][[paste(location)]][["corr"]],3)))
#     rm("model", "intersection")
#   }
  
#   if(create_output){
#     save.image(file=paste0(save_here, '/', file_number, ' Environment.RData'))
#     saveRDS(values_list, file=paste0(save_here, '/', file_number, ' values_list.rds'))
#   }
# }

# time <- toc()
# if(create_output){
#   save.image(file=paste0(save_here, '/', file_number, ' Environment.RData'))
# }
