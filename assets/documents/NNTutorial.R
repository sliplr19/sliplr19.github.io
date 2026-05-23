### 3-node network


### Functions

#### Generate network function

netgen3 <- function(n, meanA, sdA, meanB, sdB, meanC, sdC,
                    beta_con, beta_ab, beta_lin, beta_non, 
                    func = c("quad", "inter", "linear")){
  A <- rnorm(n, meanA, sdA) 
  if(func == "quad"){
    B <- beta_ab*A + rnorm(n, meanB, sdB)
    C <- beta_non*A^2 + beta_lin*A + beta_con*B + rnorm(n, meanC, sdC)
  }
  else if(func == "inter"){
    B <- beta_ab*A + rnorm(n, meanB, sdB)
    C <- beta_non*A*B + beta_lin*A + beta_con*B + rnorm(n, meanC, sdC)
  }
  else if(func == "linear"){
    B <- beta_ab*A + rnorm(n, meanB, sdB)
    C <- beta_lin*A + beta_con*B + rnorm(n, meanC, sdC) 
  }
  else if(func == "log"){
    B <- beta_ab*A + rnorm(n, meanB, sdB)
    C <- log(abs(A)) + beta_lin*A + beta_con*B + rnorm(n, meanC, sdC)
  }
  dat <- data.frame("A" = A, "B" = B, "C" = C)
  return(dat)
}


#### Residualize function

get_residuals <- function(dat, A, C, B_variables){
  require(mgcv)
  dat <- data.frame(dat)
  colnames(dat) <- c("A", "B", "C")
  Bs <- paste0("s(", B_variables, ")", collapse = " + ")
  Aform <- as.formula(paste("A ~", Bs))
  Cform <- as.formula(paste("C ~", Bs))
  A_fit <- gam(Aform, data=dat, method="REML")
  C_fit <- gam(Cform, data=dat, method="REML")
  
  A_resid <- resid(A_fit)
  C_resid <- resid(C_fit)
  
  A_resid_lm <- resid(lm(A_resid ~ C_resid))
  
  newdat <- data.frame("A" = A_resid_lm, 
                       "B"= dat$B,
                       "C" = C_resid)
  
  return(data.frame(newdat))
}


#### Conditional mutual information function

CondInfoN <- function(dat, n_perm = 1000) {
  require(infotheo)
  
  nbins <- min(floor(sqrt(nrow(dat))), apply(dat, 2, function(x) length(unique(x))) - 1)
  nbins <- max(2, min(nbins))
  dat <- discretize(dat, "globalequalwidth", nbins)
  dim_cor <- ncol(dat)
  k <- 1:dim_cor
  
  cor_mat <- matrix(0, dim_cor, dim_cor)
  pval_mat <- matrix(NA, dim_cor, dim_cor)
  colnames(cor_mat) <- colnames(dat)
  rownames(cor_mat) <- colnames(dat)
  colnames(pval_mat) <- colnames(dat)
  rownames(pval_mat) <- colnames(dat)
  
  # Compute observed CMI
  for (i in 1:dim_cor) {
    for (j in 1:dim_cor) {
      k_0 <- setdiff(k, c(i, j))
      cor_mat[j, i] <- condinformation(dat[, i], dat[, j], dat[, k_0], method = "sg")
    }
  }
  
  set.seed(0619)
  # Store permutation statistics
  for (i in 1:dim_cor) {
    for (j in 1:dim_cor) {
      k_0 <- setdiff(k, c(i, j))
      I_obs <- cor_mat[j, i]
      perm_vals <- numeric(n_perm)
      
      for (p in 1:n_perm) {
        xi_perm <- sample(dat[, i])
        perm_vals[p] <- condinformation(xi_perm, dat[, j], dat[, k_0], method = "sg")
      }
      
      pval_mat[j, i] <- (sum(perm_vals >= I_obs) + 1) / (n_perm + 1)
    }
  }
  
  return(as.data.frame(pval_mat))
}


### Set parameters

n = 500
meanA = 0
 sdA = 1
 meanB = 0
 sdB = 1
 meanC = 0
 sdC = 1
 beta_con = 1
 beta_ab = 1
 beta_lin = 1
 beta_non = 1
 func = "quad"
group = "Residual"
 nperm = 1000
 
### Set up
 
 cond_results <- data.frame(
                            "n" = n, "meanA" = meanA,
                            "sdA" = sdA, "meanB" = meanB, "sdB" = sdB,
                            "meanC" = meanC, "sdC" = sdC, "beta_lin" = beta_lin,
                            "beta_non" = beta_non, "beta_con" = beta_con,
                            "beta_ab" = beta_ab, "func" = func, "group" = group)
 
 require(stats)
 require(energy)
 
 res_temp <- cond_results
 
 
### Generate data
 
 dat <- data.frame(netgen3(n, meanA, sdA, meanB, sdB, meanC, sdC,
                           beta_con, beta_ab, beta_lin, beta_non, func))
 
### Residualize
 
 dat <- data.frame(get_residuals(data.frame(dat), "A", "C", "B"))
 
 
### Get p-values
 
 
 res_pdpvAC = dcor.test(dat$A, dat$C, R = 1000)$p.value
 res_sppvAC = cor.test(dat$A, dat$C, method = "spearman")$p.value
 res_prpvAC = cor.test(dat$A, dat$C, method = "pearson")$p.value
 
 
 res_pdpvBC = dcor.test(dat$B, dat$C, R = 1000)$p.value
 res_sppvBC = cor.test(dat$B, dat$C, method = "spearman")$p.value
 res_prpvBC = cor.test(dat$B, dat$C, method = "pearson")$p.value
 
 res_CMI <- CondInfoN(dat, n_perm = 1000)
 res_cmi_full <- data.frame(
   ACres_cmi_pv = res_CMI["C", "A"],       
   BCres_cmi_pv = res_CMI["C", "B"]
 )
 res_temp <- cbind(res_temp, res_pdpvAC, res_sppvAC, res_prpvAC,
                   res_pdpvBC, res_sppvBC, res_prpvBC,
                   res_cmi_full)
 
### Output
 
res_temp