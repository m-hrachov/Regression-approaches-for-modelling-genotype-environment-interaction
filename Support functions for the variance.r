# Be aware when working with small response values that this function can produce wrong results due to a rounding error
# (i.e. I had a problem that ASREml-based C inverse had values 25 times larger. When the responce was multiplied by 10 and the variance adjusted, it produced a correct result)
#4  
get_C_inv <- function(model, variances, X = NULL, Z = NULL, G = NULL, 
                        model_terms = "G+z1:G+z2:G", n_terms = 3, type = "us", G_rr = NULL, 
                        fixed.effects = NULL){
    
    fixed.effects <- c("(Intercept)", fixed.effects) %>% rev(.) # ASReml as of 4.2.0 has reversed direction for the fixed effects in the design matrix
    R <- diag(variances[!is.na(variances)])  # alphabetic order of L Y G
    solve_R <- solve(R)
    
    if(is.null(Z)){
      Z <- model$design[,!dimnames(model$design)[[2]] %in% fixed.effects]
      NA_count <- 0
      try(NA_count <- Z %>% as.matrix() %>% as.data.frame() %>%  select(starts_with("mv_")) %>% ncol())
      if(NA_count != 0){
        Z <- Z[!apply(Z, 1, function(row) any(row == -1)), ] %>% as.matrix() %>% as.data.frame() %>% select(-starts_with("mv_")) %>% as.matrix()
      }
      if(type == "rr"){
        Z <- Z[!apply(Z, 1, function(row) any(row == -1)), ] %>% as.matrix() %>% as.data.frame() %>% select(-starts_with("mbf(Co)_x0")) %>% as.matrix()
        
      }
    } else {stop("This function doesn't support custom Z. Under development.")}
    
    if( is.null(X)){
      if(NA_count != 0){
        X <- model$design[!apply(model$design, 1, function(row) any(row == -1)), ] %>% as.matrix() %>% as.data.frame() %>% 
          select(-starts_with("mv_")) %>% 
          select(all_of(fixed.effects)) %>% 
          as.matrix()
        
      } else {
        X <- model$design[,dimnames(model$design)[[2]] %in% fixed.effects]
      }
    } else {stop("This function doesn't support custom X. Under development.")}
    
    G_l <- diag(rep(model$G.param$L$variance$initial,model$G.param$L$variance$size))
    G_y <- diag(rep(model$G.param$Y$variance$initial,model$G.param$Y$variance$size))
    if(type == "us"){
      if(n_terms == 1){
        G_us <- diag(rep(model$G.param[[model_terms]]$variance$initial[1],model$G.param[[model_terms]]$variance$size))
      } else if(n_terms == 2){
        G_g <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[1],model$G.param[[model_terms]]$variance$size/n_terms))
        G_z1 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[3],model$G.param[[model_terms]]$variance$size/n_terms))
        G_gz1 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[2],(model$G.param[[model_terms]]$variance$size/n_terms)))
        G_us <- cbind(rbind(G_g, G_gz1), rbind(G_gz1, G_z1)) 
      } else if(n_terms == 3){
        G_g <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[1],model$G.param[[model_terms]]$variance$size/n_terms))
        G_z1 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[3],model$G.param[[model_terms]]$variance$size/n_terms))
        G_z2 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[6],model$G.param[[model_terms]]$variance$size/n_terms))
        G_gz1 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[2],(model$G.param[[model_terms]]$variance$size/n_terms)))
        G_gz2 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[4],(model$G.param[[model_terms]]$variance$size/n_terms)))
        G_z1z2 <- diag(rep(model$G.param[[model_terms]][[paste(n_terms)]]$initial[5],(model$G.param[[model_terms]]$variance$size/n_terms)))
        
        G_us <-  cbind(rbind(G_g, G_gz1), rbind(G_gz1, G_z1)) %>% rbind(., cbind(G_gz2,G_z1z2)) %>% cbind(. ,rbind(rbind(G_gz2, G_z1z2), G_z2))
        
      } else {stop("Can't adjust US structure to more than 2 covariates.")}
    } else if (type == "rr"){
      
      G_covariates <- as.matrix(expand_G(G_rr, model$G.param[[model_terms]]$variance$size/n_terms))
      
    } else {
      stop("Only unstructured and reduced rank variance-covariance are supported")
    }
    G_yl <- diag(rep(model$G.param$`Y:L`$variance$initial,model$G.param$`Y:L`$variance$size))
    G_lg <- diag(rep(model$G.param$`L:G`$variance$initial,model$G.param$`L:G`$variance$size))
    G_yg <- diag(rep(model$G.param$`Y:G`$variance$initial,model$G.param$`Y:G`$variance$size))
    
    G_ylg <- diag(rep(model$G.param$`Y:L:G`$variance$initial,model$G.param$`Y:L:G`$variance$size))
    if(type == "us"){
      G <- bdiag(G_l, G_y, G_us, G_yl, G_lg, G_yg, G_ylg)
    } else if(type == "rr"){
      G <- bdiag(G_l, G_y, G_yl, G_lg, G_yg, G_covariates, G_ylg)
    }
    
    try(C_11  <-  t(X) %*% solve_R %*% X, silent = TRUE)
    try(C_21  <-  t(Z) %*% solve_R %*% X, silent = TRUE)
    try(C_12  <-  t(X) %*% solve_R %*% Z, silent = TRUE)
    try(C_22  <-  t(Z) %*% solve_R %*% Z + solve(G), silent = TRUE)
    
    try(C_1 <- cbind(C_11,C_12), silent = TRUE)
    try(C_2 <- cbind(C_21,C_22), silent = TRUE)
    
    try(C <- rbind(C_1, C_2), silent = TRUE)
    
    try( C_inv <- solve(as.matrix(C)), silent = TRUE)
    # in case the above method fails, use the one from Henderson, "Application of linear models in animal breeding" (1984)
    if(exists("C_inv")){
      if(any(is.na(C_inv))){
        print("C_inv produced NA! Alternative derivation is used. Check the inputs if the alternative fails.")
        M_11 = t(X) %*% solve_R %*% X
        M_21 = G %*% t(Z) %*% solve_R %*% X
        M_12 = t(X) %*% solve_R %*% Z %*% G
        M_22 = G %*% t(Z) %*% solve_R %*% Z %*% G + G
        
        M_1 <- cbind(M_11,M_12)
        M_2 <- cbind(M_21,M_22)
        
        M <- rbind(M_1, M_2)
        M_inv <- MASS::ginv(as.matrix(M))
        
        C_inv <- bdiag(diag(rep(1, length(fixed.effects))), G) %*% M_inv %*% bdiag(diag(rep(1, length(fixed.effects))), G)
        colnames(C_inv) <- c(fixed.effects, colnames(Z))
        rownames(C_inv) <- c(fixed.effects, colnames(Z))
      } 
    } else {
      print("Standard derivation failed. Alternative derivation of C is used.")
      ### if G is singular
      M_11 = t(X) %*% solve_R %*% X
      M_21 = G %*% t(Z) %*% solve_R %*% X
      M_12 = t(X) %*% solve_R %*% Z %*% G
      M_22 = G %*% t(Z) %*% solve_R %*% Z %*% G + G
      
      M_1 <- cbind(M_11,M_12)
      M_2 <- cbind(M_21,M_22)
      
      M <- rbind(M_1, M_2)
      M_inv <- MASS::ginv(as.matrix(M))
      
      C_inv <- bdiag(diag(rep(1, length(fixed.effects))), G) %*% M_inv %*% bdiag(diag(rep(1, length(fixed.effects))), G)
      colnames(C_inv) <- c(fixed.effects, colnames(Z))
      rownames(C_inv) <- c(fixed.effects, colnames(Z))
    }
 
  return(C_inv)
  
}

expand_G <- function(G, block_size) {
  # for symmetric input matrix G!!!!!!
  # Get the dimensions of the original G matrix
  p <- nrow(G)

  # Initialize the resulting matrix with zeros
  resulting_G <- matrix(0, nrow = p * block_size, ncol = p * block_size)

  # Fill the resulting matrix based on the structure
  for (i in 1:p) {
    for (j in 1:p) {
      block <- diag(G[i, j], block_size)
      # Corrected indexing with proper parentheses
      resulting_G[((i - 1) * block_size + 1):(i * block_size), ((j - 1) * block_size + 1):(j * block_size)] <- block
      # Ensure symmetry by copying the block to the symmetric position if i != j
      if (i != j) {
        resulting_G[((j - 1) * block_size + 1):(j * block_size), ((i - 1) * block_size + 1):(i * block_size)] <- block
      }
    }
  }

  return(resulting_G)
}
#3
get_gammas_from_the_main_model <- function(model, type = c("us", "rr"), slopes = c("z1:G", "z2:G"), intercept = "G"){
  # assume, variances of gammas are the variances of individual coefficients, while covariances are given in 
  # $varcomp of the summary respectively, and are fixed values between any two genotypes.
  if(!type %in% c("us", "rr")){stop("The function supports only us() and rr() variance covaraince for intercepts and slopes!")}
  
  all_values <- cbind(as.data.frame(model[["coefficients"]][["random"]]), data.frame(var = model[["vcoeff"]][["random"]]))
  
  effects_list <- list()
  var_list <- list()
  
  i <-  0
  # Loop over each slope and extract the corresponding effects
  for (slope in c(intercept, slopes)) {
    if (slope == intercept){
      filtered_values <- all_values %>% 
        .[grep(x = rownames(.), pattern = paste0(slope, "_"), fixed = TRUE), c("effect", "var"), drop = FALSE] %>% 
        .[grep(x = rownames(.), pattern = paste0("^(?!.*:", slope, "|", slope, ":)", "|", slope, "^(?!.*::)"), perl = TRUE), c("effect", "var"), drop = FALSE]
    rownames(filtered_values) <- gsub("mbf\\(Co\\)_z0:", "", rownames(filtered_values))
      } else {
      filtered_values <- all_values %>% 
      .[grep(x = rownames(.), pattern = paste0(slope, "_"), fixed = TRUE), c("effect", "var"), drop = FALSE]
      rownames(filtered_values) <- gsub("mbf\\(Co\\)_", "", rownames(filtered_values))
    }
    # Store the effect in the list with the name equal to the slope
    effects_list[[slope]] <- filtered_values[,1, drop = FALSE]
    var_list[[slope]] <- filtered_values[,2]
    i <- i + 1
  }
  
  # Combine the effects into a data frame
  gamma_prime <- do.call(data.frame, effects_list)
  colnames(gamma_prime) <- names(effects_list)
  
  var_gamma <- do.call(data.frame, var_list)
  colnames(var_gamma) <- names(var_list)
  rownames(var_gamma) <- rownames(gamma_prime)
  
  both_values <- list(gamma_prime = gamma_prime, simple_var_gamma = var_gamma)
  
  return(both_values)
}

#1
extract_xi_for_unseen_location <- function(multivar.predict, new_location, location = "L"){
  
  xi_prime_blue <- multivar.predict$pvals %>% 
    select(c(Y, L, trait, predicted.value)) %>% 
    pivot_wider(names_from = trait, values_from = predicted.value) %>%
    as.data.frame() %>% 
    filter(L == new_location) %>% 
    mutate(intercept = 1,
           LY = paste(L,Y, sep = "_")) %>% 
    select(-c(L,Y)) %>% 
    column_to_rownames(., var = "LY") %>% 
    select(intercept, everything()) %>% as.matrix()
    
  
  var_xi_blue_all <- multivar.predict$vcov
  dimnames(var_xi_blue_all)[[1]] <- paste(multivar.predict$pvals[[location]], multivar.predict$pvals$trait, sep = "_")
  dimnames(var_xi_blue_all)[[2]] <- paste(multivar.predict$pvals[[location]], multivar.predict$pvals$trait, sep = "_")
  var_xi_blue_no_intercept <- var_xi_blue_all %>% 
    .[grep(x = rownames(.), pattern = paste0(new_location), fixed = TRUE), grep(x = rownames(.), pattern = paste0(new_location), fixed = TRUE)]
  var_xi_blue <- rbind(0, cbind(0, var_xi_blue_no_intercept))
  
  both <- list(xi_prime_BLUE = xi_prime_blue, var_xi_BLUE = var_xi_blue)
  return(both)
}

# reorder_symmetric_matrix <- function(I, n_slopes, subset_C) {
#   # Generate intercepts
#   intercepts <- 1:I
#   
#   # Generate slopes using a loop for n_slopes
#   slopes_list <- lapply(1:n_slopes, function(s) {
#     start <- I * s + 1
#     end <- I * (s + 1)
#     return(start:end)
#   })
#   
#   # Combine intercepts and slopes in the desired order
#   new_order <- as.vector(do.call(rbind, c(list(intercepts), slopes_list)))
#   
#   # Reorder the symmetric matrix subset_C based on new_order
#   subset_C_ordered <- subset_C[new_order, new_order]
#   
#   return(subset_C_ordered)
# }

#5
get_var_gamma <- function(subset_C, slopes = c("z1", "z2"), intercept = "G", geno_names, type = "us"){
  if(type == "us"){
    all_list <- list()
    for (name in geno_names){
      corrected_names <- c(name)
      for(slope in slopes){
        corrected_names <- c(corrected_names, paste(slope, name,sep = ":"))
      }
      one_entry <- as.matrix(subset_C)[corrected_names, corrected_names]
      all_list <- c(all_list, list(one_entry))
    }
    names(all_list) <- geno_names
    
  } else if(type == "rr")
    {
    all_list <- list()
    for (name in geno_names){
      corrected_names <- c(paste(intercept, name, sep = ":"))
      for(slope in slopes){
        corrected_names <- c(corrected_names, paste(slope, name,sep = ":"))
      }
      one_entry <- as.matrix(subset_C)[corrected_names, corrected_names]
      all_list <- c(all_list, list(one_entry))
    }
    names(all_list) <- geno_names
  }
  return(all_list)
}

#2
get_sigma_x_prime <- function(multivar, year = "Y", variables = c("z1", "z2"), reverse.triangle = F){
  df_sigma_x_Y <-  summary(multivar)$varcomp[grep(x = rownames(summary(multivar)$varcomp),pattern = year,fixed = T), "component", drop = F] %>% 
    mutate(term = rownames(.)) %>% 
    select(term, component)
  df_sigma_x_LY <-  summary(multivar)$varcomp[grep(x = rownames(summary(multivar)$varcomp),pattern = "units:trait!trait",fixed = T), "component", drop = F] %>% 
    mutate(term = rownames(.)) %>% 
    select(term, component)
  
  extract_covar <- function(dataframe, variables){
  
  df <- dataframe %>%
    mutate(stripped_terms = str_extract_all(term, paste0("(", paste(variables, collapse = "|"), ")"))) %>%
    mutate(stripped_terms = sapply(stripped_terms, function(x) paste(x, collapse = ":")))  %>%
    select(stripped_terms, component) %>% 
    separate_wider_delim(stripped_terms, names = c("var1", "var2"), delim = ":")
  
  combs <-  expand.grid(var1 = variables, var2 = variables)
  
  cov_matrix <- combs %>%
    left_join(df, by = c("var1", "var2")) %>%
    replace_na(list(component = 0)) %>%
    spread(var2, component) 
  
  rownames(cov_matrix) <- cov_matrix$var1
  cov_matrix <- cov_matrix %>% select(-var1) %>% as.matrix()
  if(reverse.triangle){
    cov_matrix[lower.tri(cov_matrix)] <- t(cov_matrix)[lower.tri(cov_matrix)]
  } else {
    cov_matrix[upper.tri(cov_matrix)] <- t(cov_matrix)[upper.tri(cov_matrix)]
  }
  dimnames(cov_matrix) <- list(variables, variables)
  
  return(cov_matrix)
  }
  
  sigma_x_Y <- rbind(intercept = 0, cbind(intercept = 0, extract_covar(df_sigma_x_Y, variables)))
  sigma_x_LY <- rbind(intercept = 0, cbind(intercept = 0, extract_covar(df_sigma_x_LY, variables)))
  
  sigma_x_prime <- sigma_x_Y + sigma_x_LY
  
  return(list(sigma_x_Y = sigma_x_Y, sigma_x_LY = sigma_x_LY, sigma_x_prime = sigma_x_prime))
}

# 
# extract_vi_i_x_prime <- function(gamma_prime, sigma_x_prime, var_gamma){
#   vi_x_prime <- c()
#   for(i in 1:length(rownames(gamma_prime))){
#     vi_x_prime[i] <- as.matrix(gamma_prime[i,]) %*% (sigma_x_prime) %*% t(as.matrix(gamma_prime[i,])) - sum(diag(var_gamma) %*% sigma_x_prime)
#   }
#   names(vi_x_prime) <- rownames(gamma_prime)
#   
#   return(vi_x_prime)
# }


# get_est.var_gamma_prime_xi_prime <- function(gamma_prime, var_xi, xi_prime, var_gamma_i_prime){
#   est.var_gamma_prime_xi_prime <- c()
#   for(i in 1:length(rownames(gamma_prime))){
#     est.var_gamma_prime_xi_prime[i] <- as.matrix(gamma_prime[i,]) %*% as.matrix(var_xi) %*% t(as.matrix(gamma_prime[i,])) +
#       as.matrix(xi_prime) %*% as.matrix(var_gamma_i_prime[i,]) %*% t(as.matrix(xi_prime)) - 
#       sum(diag(as.matrix(var_gamma_i_prime[i,]) %*% as.matrix(var_xi)))
#   }
#   names(est.var_gamma_prime_xi_prime) <- rownames(gamma_prime)
#   return(est.var_gamma_prime_xi_prime)
# }
# 
