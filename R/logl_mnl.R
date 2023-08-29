#' @keywords internal
# Evaluate loglikelihood for intercept only Multinomial logit (MNL) model., bassed on P. Rossi 2004
# 本质就是最大似然估计MLE
llmnl_int <- function(int, Obs, n_cat) {
  Obs <- c(Obs,c(1,2,1,2))
  n_Obs <- length(Obs)
  betas <- c(0, int)
  if(is.infinite(log(sum(exp(betas))))){
    alpha <- 500
  }
  return( sum(betas[Obs]) - log(sum(exp(betas))) * n_Obs)

}
#' @keywords internal
# Evaluate loglikelihood for intercept only Multinomial logit (MNL) model., bassed on P. Rossi 2004
llmnl_int2 <- function(int, Obs, n_cat) {
  #browser()
  Obs[[1]] <- c(Obs[[1]],c(1,2,1,2))
  for(i in 1:9){
  Obs[[2]][i][[1]] <- c(Obs[[2]][i][[1]],c(0,1,0,1))
  }
  n_Obs <- length(Obs[[1]])
  betas <- c(0, int)
  #betas1 <- Obs[[3]][2,][1]
  #betas2 <- Obs[[3]][3,][1]
  sum <- 0
  alpha <- 0
  for(i in 1:n_Obs){
    alpha <- alpha + betas[Obs[[1]]][i]
    for(j in 1:9){
      alpha <- alpha + Obs[[3]][j+1,][1] * Obs[[2]][[j]][i]
    }
    #alpha <- betas[Obs[[1]]][i] + betas1*Obs[[2]][[1]][i] + betas2*Obs[[2]][[2]][i]
    if(is.infinite(log(exp(alpha)+1))){
      alpha <- 500
    }
    sum <- sum + alpha - log(exp(alpha)+1)
  }
  return(sum)
}

#' @keywords internal
# Obtain fractional log likelihood for multinomial intercept only model, bassed on P. Rossi 2004
llmnl_int_frac <- function(int, Obs, n_cat, pooled_likel, w, wgt){
  return((1 - w) * llmnl_int(int = int, Obs = Obs, n_cat = n_cat) + w * wgt * pooled_likel)
}
