####################################################################
## Supplementary script for the submission:
## "Regression approaches for modelling genotype-environment interaction and 
## making predictions into unseen environments"
##
## Script author: Maksym Hrachov
## Contact: maksym.hrachov@uni-hohenheim.de
###################################################################

# This script shows similarities (and differences) of C-inverse matrix obtained from ASReml-R and 
# approach from Henderson (1984, Chapter 5).
# The initial reason for this comparison is an observation that ASReml-R 4.2 returned block-diagonal matrix
# instead of a complete one for the part pertaining to genotype-specific covariate slopes. 
# We investigated if two matrices are equal for the non-zero part in both matrices.
# Short conclusion: they are almost identical for non-zero entries, and small differences are likely due to rounding error.
# We also examined if direct specification of ~str(G + z1:G + z2:G, ~us(3)) and ~str(mbf(X):G, ~us(3)) resulting into 
# the same C-iverse. Short answer: yes, they do, but small differences exist likely do to the rounding error during
# matrix operations.


#
# this function expands a dense variance-covariance into a block-diagonal G matrix
#
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


#
# This one is Henderson's method (1984, Chapter 5) for both non-singular and singular G matrices.
#
get_C_inv <- function(model, variances, X = NULL, Z = NULL, G = NULL, num_G = NULL,
                      model_terms = "G+z1:G+z2:G", n_terms = 3, type = "us", G_rr = NULL,
                      fixed.effects = NULL){
  
  if(is.null(num_G)){ num_G = model$G.param[[model_terms]]$variance$size/n_terms }
  X = NULL 
  Z = NULL
  
  fixed.effects <- c("(Intercept)", fixed.effects) %>% rev(.) # ASReml as of 4.2.0 has reversed direction for the fixed effects in the design matrix
  R <- diag(variances[!is.na(variances)]) 
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
    us_estimates_df <- model$G.param[[model_terms]][[paste(n_terms)]]$initial %>% 
      data.frame(.) %>% 
      tibble::rownames_to_column(., "component") %>% 
      rename(value = ".") %>% 
      separate(component, into = c("part1", "part2"), sep = "_", remove = TRUE) %>%
      separate(part2, into = c("row", "column"), sep = ":", remove = TRUE) %>% 
      select(-part1) %>% 
      mutate(across(all_of(c("row","column")), ~as.numeric(.)))
    
    us_dense_var_mat <- matrix(rep(NA, n_terms^2), ncol = n_terms)
    for (i in 1:nrow(us_estimates_df)){
      us_dense_var_mat[us_estimates_df$row[i], us_estimates_df$column[i]] <- us_estimates_df$value[i]
    }
    us_dense_var_mat[upper.tri(us_dense_var_mat)] <- t(us_dense_var_mat)[upper.tri(us_dense_var_mat)]
    
    G_g_and_cov <- as.matrix(expand_G(us_dense_var_mat, num_G))
    
  } else if (type == "rr"){
    
    G_g_and_cov <- as.matrix(expand_G(G_rr, num_G))
    
  } else {
    stop("Only unstructured and reduced rank variance-covariance are supported")
  }
  G_yl <- diag(rep(model$G.param$`Y:L`$variance$initial,model$G.param$`Y:L`$variance$size))
  G_lg <- diag(rep(model$G.param$`L:G`$variance$initial,model$G.param$`L:G`$variance$size))
  G_yg <- diag(rep(model$G.param$`Y:G`$variance$initial,model$G.param$`Y:G`$variance$size))
  
  G_ylg <- diag(rep(model$G.param$`Y:L:G`$variance$initial,model$G.param$`Y:L:G`$variance$size))
  
  tag_component <- function(nm) {
    nm <- trimws(nm)
    if (nm == "L")              return("G_l")
    if (grepl("^mbf\\(", nm))   return("G_g_and_cov")
    if (grepl("^z\\d+", nm)) return("G_g_and_cov")
    if (grepl("^G\\+z\\d+", nm)) return("G_g_and_cov")
    if (nm == "Y")              return("G_y")
    if (nm == "L:G")            return("G_lg")
    if (nm == "Y:G")            return("G_yg")
    if (nm == "Y:L")            return("G_yl")
    if (nm == "Y:L:G")          return("G_ylg")
    if (nm == "units!R")        return(NA_character_) # residual, not in G
    NA_character_
  }
  
  # 2) Determine the block order from first occurrence in vparameters
  vp_names <- names(model$vparameters)
  block_order <- vp_names %>% 
    map_chr(tag_component)  %>% 
    discard(is.na)  %>% 
    {\(x) x[!duplicated(x)]}()
  
  # 3) List your candidate G blocks (some may be NULL/absent)
  blocks <- list(
    G_l = G_l,
    G_g_and_cov = G_g_and_cov,
    G_y = G_y,
    G_lg = G_lg,
    G_yg = G_yg,
    G_yl = G_yl,
    G_ylg = G_ylg
  )
  
  # 4) Keep only blocks that exist and are in the detected order
  blocks_ordered <- blocks[intersect(block_order, names(blocks))]
  blocks_ordered <- blocks_ordered[!vapply(blocks_ordered, is.null, logical(1))]
  
  # 5) Build bdiag in the dynamic order
  G <- do.call(bdiag, blocks_ordered)
  message("Detected order from vparameters: ",
          paste(block_order, collapse = ", "))

  try(C_11  <-  t(X) %*% solve_R %*% X, silent = TRUE)
  try(C_21  <-  t(Z) %*% solve_R %*% X, silent = TRUE)
  try(C_12  <-  t(X) %*% solve_R %*% Z, silent = TRUE)
  try(C_22  <-  t(Z) %*% solve_R %*% Z + solve(G), silent = TRUE)
  
  try(C_1 <- cbind(C_11,C_12), silent = TRUE)
  try(C_2 <- cbind(C_21,C_22), silent = TRUE)
  
  try(C <- rbind(C_1, C_2), silent = TRUE)
  
  try( C_inv <- solve(as.matrix(C)), silent = TRUE)
  
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


# Libraries 
library(asreml)
library(asremlPlus)
library(tidyverse)
library(tictoc)
library(dae)
asreml.options(maxit=20, extra=2, pworkspace="2gb", workspace="2gb", design = T)
tic()

# this is a dummy dataset that mimics real data
load(file ="D:/bwSync/Outputs/03 BRRI/fork_v1/Proof that C_inv is fine/simulated_dataset.RData")

prep0 <- sim_dataset %>% 
  mutate(across(c("G", "Y", "L"), factor)) %>% 
  mutate(wt = 1/(std.error^2)) %>% 
  mutate(ENV_BRRI = as.factor(paste(Y,L, sep="_")))

stock_vars <- prep0 %>% 
  select(!any_of(c("wt", "ENV_BRRI", "L", "Y", "G", "predicted.value", "std.error", "Trial.year", "Release.year", "G_names"))) %>% 
  names()

rescaled_vars <- prep0 %>% 
  select(L,Y,all_of(stock_vars)) %>% 
  distinct() %>% 
  mutate(across(all_of(stock_vars), ~scale(.)[,1])) %>% 
  data.frame()

prep1 <- prep0 %>% 
  select(-all_of(stock_vars)) %>% 
  left_join(., rescaled_vars, by = c("L", "Y"))

C1 <- prep1 %>% 
  rename(ENV = ENV_BRRI) %>% 
  mutate(x0 = 0, x00 = 0) %>% 
  select(ENV, all_of(stock_vars), x0, x00) %>% 
  unique() %>% 
  mutate(across(all_of(stock_vars), ~scale(.)[,1]))

random1 <- as.formula(paste0("~L:G + Y:G + Y:L:G + str(~mbf(Co):G, vmodel = ~ rr(", length(stock_vars),
                             ",2):id(G))"))

RRR_observed <- asreml(
  fixed = predicted.value ~ G + L + Y + Y:L,
  random = random1,
  mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C1")),
  data = prep1,
  na.action = na.method(x = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1)
)

# here we create starting values to have faster convergence
RRR_starting <- asreml(
  fixed = predicted.value ~ G + L + Y + Y:L,
  random = random1,
  mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "C1")),
  data = prep1,
  na.action = na.method(x = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1),
  G.param = RRR_observed$G.param,
  start.values = T
)
start_values <- RRR_starting$vparameters.table %>% 
  mutate(Constraint = if_else(Component %in% c("L:G", "Y:G", "Y:L:G"), "P", Constraint))

rm("C1", "RRR_observed", "RRR_starting")

values_list <- c()
year = unique(prep1$Y)[1]
location = unique(prep1$L)[1]

# below we mimic what happens in LYLO CV described in the paper: leave one year and location out, and we keep correct scaling of covariates

prep1_na <- prep1 %>%
  mutate(predicted.value = if_else(Y == year | L == location, NA, predicted.value),
         across(c(wt, std.error), ~if_else(Y == year | L == location, NA, .)),
         across(c(all_of(stock_vars)), ~if_else(Y == year, NA, .)),
         wt = if_else(is.na(wt), 1, wt)) 

rescaled_stock_vars <- prep1_na %>%
  select(L,Y,all_of(stock_vars)) %>%
  distinct() %>%
  mutate(across(all_of(stock_vars), ~scale(., )[,1])) %>%
  data.frame()

prep1_1r <- prep1_na %>% 
  select(-all_of(stock_vars)) %>%
  left_join(., rescaled_stock_vars, by = c("L", "Y"))

plug.in <- prep1_1r %>% 
  select(L,Y,all_of(stock_vars)) %>% 
  distinct() %>% 
  group_by(L) %>% 
  summarise(across(-Y, ~mean(., na.rm = T))) %>% 
  mutate(Y = year) %>% 
  ungroup() %>% 
  data.frame() 

prep1.1 <- prep1_1r %>% 
  left_join(plug.in, by = c("Y", "L"), suffix = c(".xtab", ".ytab")) %>%
  mutate(ENV_BRRI = as.factor(paste(Y,L, sep="_")), 
         Y = as.factor(Y)) %>% 
  mutate(across(ends_with(".xtab"), 
                .fns = ~ coalesce(., get(sub(".xtab", ".ytab", cur_column()))),
                .names = "{sub('.xtab', '', .col)}")) %>%
  select(-ends_with(".xtab"), -ends_with(".ytab")) %>% 
  select(G,L,Y, predicted.value, std.error, wt, ENV_BRRI, all_of(stock_vars)) %>%
  arrange(Y, L, G)

X2 <- prep1.1 %>% 
  rename(ENV = ENV_BRRI) %>% 
  select(ENV, all_of(stock_vars)) %>% 
  mutate(across(c("ENV"), factor)) %>% 
  mutate(x0 = 0, x00 = 0) %>% unique()

mod1 <- asreml(
  fixed = predicted.value ~ G + L + Y + Y:L,
  random = random1,
  mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "X2")),
  data = prep1.1,
  na.action = na.method(x = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1)
  )

# this is the method we used to ensure some stable convergence of complex models
convergence <- TRUE
while(convergence){
  sink(file=tempfile(), append = TRUE)
  try(model_uptd <- update(mod1, maxit = 5, trace = F))
  sink()
  if( mod1$loglik > model_uptd$loglik | model_uptd$loglik > 5000 ){ # custom overfitting threshold (LogLik 5000 should not be observed in this CV for this data at all)
    rm("model_uptd", "h1")
    print("Loop exit criteria were met")
    break()
  } else if(exists("h1") &&  round(model_uptd$loglik,5) == round(h1$loglik,5)){ # we decided that this threshold makes sense for this data and this model
    mod1 <- model_uptd
    rm("model_uptd", "h1")
    print("Convergence was reached")
    break()
  } else {
    h1 <- mod1
    mod1 <- model_uptd
    rm("model_uptd")
  }
}

lambdas <- summary(mod1)$varcomp[grep(x = rownames(summary(mod1)$varcomp),pattern = "fa",fixed = T),"component"]
lambda1 <- lambdas[1:length(stock_vars)]
lambda2 <- lambdas[c(length(stock_vars)+1):c(length(stock_vars)*2)]
# Computing z1 and merging
X <- X2 %>%  select(!any_of(c("x0", "x00", "ENV")))
z1 <- t(lambda1 %*% t(X))
z2 <- t(lambda2 %*% t(X))

z_env <- data.frame(z1 = z1, z2 = z2, ENV_BRRI = unique(prep1.1$ENV_BRRI))

prep1.1_1 <- merge(prep1.1, z_env, by = "ENV_BRRI")

rm("mod1")

vars <- c("z1", "z2") # here only names of EC are stored

prep3 <- prep1.1_1 %>% 
  select(-all_of(stock_vars)) %>% 
  mutate(predicted.value = if_else(Y == year | L == location, NA, predicted.value),
         across(c(wt), ~if_else(Y == year | L == location, NA, .)),
         wt = if_else(is.na(wt), 1, wt))

# specify the model and C-spares (which is C-inverse from ASReml-R)

asreml.options(Csparse = ~ G + z1:G + z2:G)
FW_m1 <- asreml(
  fixed = predicted.value ~ z1 + z2,
  random = ~ str(~G + z1:G + z2:G, ~us(3):id(G)) + L + Y + Y:L + L:G + Y:G + Y:L:G,
  data = prep3,
  na.action = na.method(x = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1)
)

X_fw <- prep3 %>% 
  rename(ENV = ENV_BRRI) %>% 
  select(ENV, all_of(vars)) %>% 
  mutate(across(c("ENV"), factor)) %>% 
  mutate(z0 = 1) %>% unique() %>% 
  select(z0, z1, z2, ENV)

asreml.options(ai.sing = T, Csparse = ~ mbf(Co):G)
FW_m2 <- asreml(
  fixed = predicted.value ~ z1 + z2,
  random = ~ str(~mbf(Co):G, ~us(3):id(G)) + L + Y + Y:L + L:G + Y:G + Y:L:G,
  data = prep3,
  mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "X_fw")),
  na.action = na.method(x = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1)
)

fixed_part <- paste0("predicted.value ~ ", paste(stock_vars, collapse = " + "), collapse = " ")
fixed_formula <- formula(paste(fixed_part))
random1 <- as.formula(paste0("~ str(~mbf(Co):G, vmodel = ~ rr(", length(stock_vars)+1,
                             ",2):id(G)) + L + Y + Y:L + L:G + Y:G + Y:L:G"))

X_str2 <- prep1.1 %>% 
  mutate(ENV = ENV_BRRI,
         z0 = 1, x0 = 0, x00 = 0) %>% 
  select(ENV, all_of(stock_vars), z0, x0, x00) %>% 
  mutate(across(c("ENV"), factor)) %>% 
  unique()

asreml.options(ai.sing = T, Csparse = ~ mbf(Co):G)
RR_m3 <- asreml(
  fixed = fixed_formula,
  random = ~str(~mbf(Co):G, vmodel = ~rr(9, 2):id(G)) + L + Y + Y:L + L:G + 
    Y:G + Y:L:G,
  mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "X_str2")),
  data = prep1.1,
  na.action = na.method(x = "include", y = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1),
  maxit = 100
)
RR_m3 <- update(RR_m3)

X_str3 <- prep1.1 %>% 
  mutate(ENV = ENV_BRRI,
         z0 = 1) %>% 
  select(ENV, all_of(stock_vars), z0) %>% 
  mutate(across(c("ENV"), factor)) %>% 
  unique()

asreml.options(ai.sing = T, Csparse = ~ mbf(Co):G)
US_m4 <- asreml(
  fixed = fixed_formula,
  random = ~str(~mbf(Co):G, vmodel = ~us(9):id(G)) + L + Y + Y:L + L:G + 
    Y:G + Y:L:G,
  mbf= list(Co = list(key= c("ENV_BRRI","ENV"), cov = "X_str3")),
  data = prep1.1,
  na.action = na.method(x = "include", y = "include"),
  weights = wt,
  wald = list(denDF="algebraic"),
  family = asr_gaussian(dispersion=1),
  maxit = 100
)
US_m4 <- update(US_m4)

# the C_sparse from the two different uses of FW model will have different lengths (this means the rest of entries are 0).

sparse_FW_m1_original <- as.data.frame(FW_m1[["Csparse"]])

sparse_FW_m1_processed <- sparse_FW_m1_original %>%
  as_tibble(rownames = "id") %>%                        # 1. bring rownames into a column "id"
  mutate(
    id_clean = str_remove(id, "\\.[0-9]+$"),            # 2. strip trailing ".###"
    id_factor = as.integer(factor(id_clean)) - 1        # 3. convert to 0,1,2,... sequential
  ) %>%
  mutate(
    row    = case_when(id_factor == 1 ~ row + length(unique(prep1.1$G)),
                       id_factor == 2 ~ row + length(unique(prep1.1$G))*2,
                       TRUE ~ row),                     # 4. adjust "row"
    column = case_when(id_factor == 1 ~ column + length(unique(prep1.1$G)),
                       id_factor == 2 ~ column + length(unique(prep1.1$G))*2,
                       TRUE ~ column)                   # 5. adjust "column"
  ) %>% 
  as.data.frame
rownames(sparse_FW_m1_processed) <- sparse_FW_m1_processed$id
sparse_FW_m1 <- sparse_FW_m1_processed %>% select(row, column, Cij)
sparse_FW_m2 <- as.data.frame(FW_m2[["Csparse"]])
sparse_RR_m3 <- as.data.frame(RR_m3[["Csparse"]])
sparse_US_m4 <- as.data.frame(US_m4[["Csparse"]])
# extract my C_inv for each model
variances <- prep1.1$std.error^2

lambdas_rr <- summary(RR_m3)$varcomp[grep(x = rownames(summary(RR_m3)$varcomp),pattern = "fa",fixed = T),"component"]
variances_rr <- summary(RR_m3)$varcomp[grep(x = rownames(summary(RR_m3)$varcomp),pattern = "var",fixed = T),"component"]
G_rr <- matrix(lambdas_rr, ncol = 2) %*% t(matrix(lambdas_rr, ncol = 2)) 

C_FW_m1 <- get_C_inv(FW_m1, variances, model_terms = "G+z1:G+z2:G", 
                           n_terms = c(length(vars)+1), type = "us", fixed.effects = vars)

C_FW_m2 <- get_C_inv(FW_m2, variances, model_terms = "mbf(Co):G", 
                           n_terms = c(length(vars)+1), type = "us", fixed.effects = vars)

C_RR_m3 <- get_C_inv(RR_m3, variances, model_terms = "mbf(Co):G", 
                           n_terms = c(length(stock_vars)+3), G_rr = G_rr, type = "rr", fixed.effects = stock_vars)

C_US_m4 <- get_C_inv(US_m4, variances, model_terms = "mbf(Co):G", 
                         n_terms = c(length(stock_vars)+1), type = "us", fixed.effects = stock_vars)


sum(round(C_FW_m1,8) != round(C_FW_m2,8)) # differences after the 8th digit can be due to some numerical instability or not 100% identical convergence
# don't compare these with RR yet, because RR has different ordering of effects.
# This is just to show that both ways of fitting the model result to essentially the same estimates.

##################
dim_of_C_FW <- length(unique(prep3$G))*(length(vars)+1) # plus G main, plus x0 and x00

Cm_FW_m1 <- matrix(rep(0,dim_of_C_FW*dim_of_C_FW), ncol = dim_of_C_FW)
Cm_FW_m2 <- matrix(rep(0,dim_of_C_FW*dim_of_C_FW), ncol = dim_of_C_FW)

for (i in 1:nrow(FW_m1[["Csparse"]])){
  Cm_FW_m1[sparse_FW_m1[[1]][i], sparse_FW_m1[[2]][i]] <- sparse_FW_m1[[3]][i]
}
for (i in 1:nrow(FW_m1[["Csparse"]])){
  Cm_FW_m2[sparse_FW_m2[[1]][i], sparse_FW_m2[[2]][i]] <- sparse_FW_m2[[3]][i]
}

Cm_FW_m1[upper.tri(Cm_FW_m1)] <- t(Cm_FW_m1)[upper.tri(Cm_FW_m1)]
Cm_FW_m1 <- as.matrix(Cm_FW_m1)
names_Cm_FW_m1  <- FW_m1$design[,-1] %>% as.matrix() %>% as.data.frame() %>%  
  select(-starts_with("mv_")) %>%  select(starts_with("G"), starts_with("z")) %>% select(-any_of(c("z1","z2"))) %>% colnames()
dimnames(Cm_FW_m1) <- list(names_Cm_FW_m1,names_Cm_FW_m1)

Cm_FW_m2[upper.tri(Cm_FW_m2)] <- t(Cm_FW_m2)[upper.tri(Cm_FW_m2)]
Cm_FW_m2 <- as.matrix(Cm_FW_m2)
names_Cm_FW_m2  <- FW_m2$design[,-1] %>% as.matrix() %>% as.data.frame() %>%  
  select(-starts_with("mv_")) %>%  select(starts_with("mbf(Co)")) %>% colnames()
dimnames(Cm_FW_m2) <- list(names_Cm_FW_m2,names_Cm_FW_m2)

# the colnames of the two objects are aligned and are going sequentially
# colnames(Cm_FW_m1)
# colnames(Cm_FW_m2)

# and if we mask entries that are not estimated by the model without mbf we get:
sum(Cm_FW_m1 == 0)
sum(Cm_FW_m2 == 0)
mask <- ((Cm_FW_m1 != 0) & (Cm_FW_m2 != 0))       # the two matrices have different number of non-zero entries
rel <- 1 - (Cm_FW_m1[mask] / Cm_FW_m2[mask])      # relative difference on valid cells
vals <- abs(rel)
# mean.diffs[[year]][[location]]   <- mean(vals, na.rm = TRUE) 
mean(vals, na.rm = TRUE) # very good alignment
# but also a conclusion: two approaches produce a good overlap, but both contain parts of output that are absent in the other implementation



############### now compare the RR vs US variance-covariances ###########

combinations_RR_m3 <- expand.grid( vars = paste0("mbf(Co)_", c("z0", stock_vars)), pars = paste0("G_", unique(prep3$G))  )
column_names_RR_m3 <- paste0(combinations_RR_m3$vars, ":", combinations_RR_m3$pars)

dim_of_C_RR <- length(unique(prep3$G))*(length(stock_vars)+3)
Cm_RR_m3 <- matrix(rep(0,dim_of_C_RR*dim_of_C_RR), ncol = dim_of_C_RR)
for (i in 1:nrow(RR_m3[["Csparse"]])){
  Cm_RR_m3[sparse_RR_m3[[1]][i], sparse_RR_m3[[2]][i]] <- sparse_RR_m3[[3]][i]
}
Cm_RR_m3[upper.tri(Cm_RR_m3)] <- t(Cm_RR_m3)[upper.tri(Cm_RR_m3)]
Cm_RR_m3 <- as.matrix(Cm_RR_m3)
names_Cm_RR_m3  <- RR_m3$design[,-1] %>% as.matrix() %>% as.data.frame() %>%  
  select(-starts_with("mv_")) %>%  select(starts_with("mbf(Co)")) %>% colnames()
dimnames(Cm_RR_m3) <- list(names_Cm_RR_m3,names_Cm_RR_m3)

Cm_RR_m3_s <- Cm_RR_m3[column_names_RR_m3,column_names_RR_m3]
C_RR_m3_s <- C_RR_m3[column_names_RR_m3,column_names_RR_m3] %>% as.matrix()


dim_of_C_US <- length(unique(prep3$G))*(length(stock_vars)+1)
Cm_US_m4 <- matrix(rep(0,dim_of_C_US*dim_of_C_US), ncol = dim_of_C_US)
for (i in 1:nrow(US_m4[["Csparse"]])){
  Cm_US_m4[sparse_US_m4[[1]][i], sparse_US_m4[[2]][i]] <- sparse_US_m4[[3]][i]
}
Cm_US_m4[upper.tri(Cm_US_m4)] <- t(Cm_US_m4)[upper.tri(Cm_US_m4)]
Cm_US_m4 <- as.matrix(Cm_US_m4)
names_Cm_US_m4  <- US_m4$design[,-1] %>% as.matrix() %>% as.data.frame() %>%  
  select(-starts_with("mv_")) %>%  select(starts_with("mbf(Co)")) %>% colnames()
dimnames(Cm_US_m4) <- list(names_Cm_US_m4,names_Cm_US_m4)

Cm_US_m4_s <- Cm_US_m4[column_names_RR_m3,column_names_RR_m3]
C_US_m4_s <- C_US_m4[column_names_RR_m3,column_names_RR_m3]

# we created C-inverse matrices from Asreml model (Cm) and from hand-coded method (C).
# we perform the same comparison, but between model-based and hand-coded (rather than between two model-based):
rm("mask", "rel", "vals")

sum(C_US_m4_s == 0)
sum(Cm_US_m4_s == 0)
mask <- ((C_US_m4_s != 0) & (Cm_US_m4_s != 0)) # compare only entries that are non-zero in both cases
rel <- 1 - (C_US_m4_s[mask] / Cm_US_m4_s[mask])     # relative difference on valid cells
vals <- abs(rel)
mean(vals)  # the value shows average percent difference between the two matrices

#### now compare RR vs RR hand-coded matrices:

rm("mask", "rel", "vals")

mask <- ((Cm_RR_m3_s != 0) & (C_RR_m3_s != 0)) # compare only entries that are non-zero in both cases
rel <- 1 - (Cm_RR_m3_s[mask] / C_RR_m3_s[mask])     # relative difference on valid cells
vals <- abs(rel)
mean(vals) 

# Conclusion: looks fine.

