

mHMM <- function(s_data, gen, xx = NULL, xt, emiss_mu_int_bar, start_val, mcmc, return_path = FALSE, print_iter, show_progress = TRUE,
                 gamma_hyp_prior = NULL, emiss_hyp_prior = NULL, gamma_sampler = NULL, emiss_sampler = NULL){
# s_data id obs1 obs2 obs3 obs4
#  t=1   1   0    1    2    3
#  t=2   1   1    2    2    1
# gen m隐变量个数 n_dep因变量obs的个数 q_emiss每个标签的类别个数
# xx表示协变量 由多个matrix组成的 第一个预测TM 其他预测EM matrix每一列代表一个协变量

# 清除掉一些exception
  if(!missing(print_iter)){
    warning("The argument print_iter is deprecated; please use show_progress instead to show the progress of the algorithm.")
  }
  if(sum(objects(gen) %in% "m") != 1 | sum(objects(gen) %in% "n_dep") != 1 | sum(objects(gen) %in% "q_emiss") != 1){
    stop("The input argument gen should contain the elements m, n_dep and q_emiss.")
  }
  if(sum(sapply(s_data, is.factor)) > 0 ){
    stop("Your data contains factorial variables, which cannot be used as input in the function mHMM. All variables have to be numerical.")
  }


# 一、定义一些基本参数------
  n_dep			 <- gen$n_dep  # 因变量obs的个数4
  dep_labels <- colnames(s_data[,2:(n_dep+1)])  # 每一个因变量的标签
  id         <- unique(s_data[,1])  # unique id  1,2,3...10
  n_subj     <- length(id)  # subject个数 10个
  subj_data  <- rep(list(NULL), n_subj) #replicate复制 意思是有多少subject就复制多少次
  for(s in 1:n_subj){
    subj_data[[s]]$y <- as.matrix(s_data[s_data[,1] == id[s],][,-1], ncol = n_dep)
    #相当于把s_data换个格式存放起来 对应subject
  }

  ypooled    <- n_t <- NULL
  n_vary     <- numeric(n_subj) #生成一个10个0的向量
  m          <- gen$m #隐变量个数2
  q_emiss 		 <- gen$q_emiss #每个标签的类别数【3，2，3，2】
  if(length(q_emiss) != n_dep){
    stop("The lenght of q_emiss specifying the number of output categories for each of the number of dependent variables should equal the number of dependent variables specified in n_dep")
  }
  emiss_mhess   <- rep(list(NULL), n_dep) #构造EM s和类别l一一对应
  for(q in 1:n_dep){
    emiss_mhess[[q]] <- matrix(, (q_emiss[q] - 1) * m, (q_emiss[q] - 1))
    #构建n*n的空矩阵
  }

  for(s in 1:n_subj){
    ypooled   <- rbind(ypooled, subj_data[[s]]$y)
    n_t         <- dim(subj_data[[s]]$y)[1]
    n_vary[s] <- n_t
    subj_data[[s]]	<- c(subj_data[[s]], n_t = n_t, list(gamma_mhess = matrix(, (m - 1) * m, (m - 1)), emiss_mhess = emiss_mhess)) #c是合并向量的意思
  }
  n_total 		<- dim(ypooled)[1]

# 二、 covariates构建（没考虑t）xx变换还是叫xx -------------------------------------------------------------------------------------------------------------
  n_dep1 <- 1 + n_dep
  nx <- numeric(n_dep1)
  if (is.null(xx)){
    xx <- rep(list(matrix(1, ncol = 1, nrow = n_subj)), n_dep1)
    nx[] <- 1
  }
  else {
    #browser()
    if(!is.list(xx) | length(xx) != n_dep1){
      stop("If xx is specified, xx should be a list, with the number of elements equal to the number of dependent variables + 1")
    }

    for(i in 1:n_dep1){

      if (is.null(xx[[i]])){
        xx[[i]] <- matrix(1, ncol = 1, nrow = n_subj)
        nx[i] <- 1
      }

      else {
        nx[i] <- ncol(xx[[i]]) #存放xx第i个element的列数，以此判断有没有协变量
        if (sum(xx[[i]][,1] != 1)){
          stop("If xx is specified, the first column in each element of xx has to represent the intercept. That is, a column that only consists of the value 1")
        }

        if(nx[i] > 1){ #大于1表示存在协变量

          for(j in 2:nx[i]){ #从第i个element的第2列到最后一列 相当于遍历每一个cov变量

            if(is.factor(xx[[i]][,j])){
              stop("Factors currently cannot be used as covariates, see help file for alternatives")

            }
            if((length(unique(xx[[i]][,j])) == 2) & (sum(xx[[i]][,j] != 0 & xx[[i]][,j] !=1) > 0)){
              stop("Dichotomous covariates in xx need to be coded as 0 / 1 variables. That is, only conisting of the values 0 and 1")

            }
            if(length(unique(xx[[i]][,j])) > 2){
              xx[[i]][,j] <- xx[[i]][,j] - mean(xx[[i]][,j]) #完事儿就是把每一项减去了平均数

            }
          }
        }
      }
    }
  }

# 三、 MCMC算法参数/sampler 和 超参数piror 的初始化 -------------------------------------------------------------------------------------------------------------
  J 				<- mcmc$J
  burn_in			<- mcmc$burn_in

# browser()
# Initialize gamma sampler
  if(is.null(gamma_sampler)) {
    gamma_int_mle0  <- matrix(0, nrow = m, ncol = m - 1)
    gamma_scalar    <- 2.93 / sqrt(m - 1)
    gamma_w         <- .1
  } else {
    if (!is.mHMM_pdRW_gamma(gamma_sampler)){
      stop("The input object specified for gamma_sampler should be from the class mHMM_pdRW_gamma, obtained by using the function pd_RW_gamma")
    }
    if (gamma_sampler$m != m){
      stop("The number of states specified in m is not equal to the number of states specified when setting the proposal distribution of the RW Metropolis sampler on gamma using the function pd_RW_gamma")
    }
    gamma_int_mle0  <- gamma_sampler$gamma_int_mle0
    gamma_scalar    <- gamma_sampler$gamma_scalar
    gamma_w         <- gamma_sampler$gamma_w
  }

# Initialize emiss sampler
  if(is.null(emiss_sampler)){
    emiss_int_mle0 <- rep(list(NULL), n_dep)
    emiss_scalar	<- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
      emiss_int_mle0[[q]] <- matrix(0, ncol = q_emiss[q] - 1, nrow = m)
      emiss_scalar[[q]] 	<- 2.93 / sqrt(q_emiss[q] - 1)
    }
    emiss_w		<- rep(list(.1), n_dep)
  } else {
    if (!is.mHMM_pdRW_emiss(emiss_sampler)){
      stop("The input object specified for emiss_sampler should be from the class mHMM_pdRW_emiss, obtained by using the function pd_RW_emiss_cat")
    }
    if (emiss_sampler$gen$m != m){
      stop("The number of states specified in m is not equal to the number of states specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
    }
    if (emiss_sampler$gen$n_dep != n_dep){
      stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
    }
    if (sum(emiss_sampler$gen$q_emiss != q_emiss) > 0){
      stop("The number of number of observed categories for each of the dependent variable specified in q_emiss is not equal to q_emiss specified when setting the proposal distribution of the RW Metropolis sampler on the emission distribution(s) using the function pd_RW_emiss_cat.")
    }
    emiss_int_mle0	<- emiss_sampler$emiss_int_mle0
    emiss_scalar 	<- emiss_sampler$emiss_scalar
    emiss_w    		<- emiss_sampler$emiss_w
  }

# Initialize TM hyper prior
  if(is.null(gamma_hyp_prior)){
    gamma_mu0	  <- rep(list(matrix(0,nrow = nx[1], ncol = m - 1)), m)
    gamma_K0			<- diag(1, nx[1])
    gamma_nu			<- 3 + m - 1
    gamma_V			  <- gamma_nu * diag(m - 1)
  } else {
    if (!is.mHMM_prior_gamma(gamma_hyp_prior)){
      stop("The input object specified for gamma_hyp_prior should be from the class mHMM_prior_gamma, obtained by using the function prior_gamma.")
    }
    if (gamma_hyp_prior$m != m){
      stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution gamma using the function prior_gamma.")
    }
    if(is.null(gamma_hyp_prior$n_xx_gamma) & nx[1] > 1){
      stop("Covariates were specified to predict gamma, but no covariates were specified when creating the informative hyper-prior distribution on gamma using the function prior_gamma.")
    }
    if(!is.null(gamma_hyp_prior$n_xx_gamma)){
      if(gamma_hyp_prior$n_xx_gamma != nx[1]){
        stop("The number of covariates specified to predict gamma is not equal to the number of covariates specified when creating the informative hper-prior distribution on gamma using the function prior_gamma.")
      }
    }
    gamma_mu0			<- gamma_hyp_prior$gamma_mu0
    gamma_K0			<- gamma_hyp_prior$gamma_K0
    gamma_nu			<- gamma_hyp_prior$gamma_nu
    gamma_V			  <- gamma_hyp_prior$gamma_V
  }

# Initialize EM hyper prior
  if(is.null(emiss_hyp_prior)){
    emiss_mu0	  <- rep(list(vector("list", m)), n_dep)
    emiss_nu	    <- rep(list(NULL), n_dep)
    emiss_V	    <- rep(list(NULL), n_dep)
    emiss_K0     <- rep(list(NULL), n_dep)
    for(q in 1:n_dep){
      for(i in 1:m){
        emiss_mu0[[q]][[i]]		<- matrix(0, ncol = q_emiss[q] - 1, nrow = nx[1 + q])
      }
      emiss_nu[[q]]		<- 3 + q_emiss[q] - 1
      emiss_V[[q]]		  <- emiss_nu[[q]] * diag(q_emiss[q] - 1)
      emiss_K0[[q]]		<- diag(1, nx[1 + q])
    }
  } else {
    if (!is.mHMM_prior_emiss(emiss_hyp_prior)){
      stop("The input object specified for emiss_hyp_prior should be from the class mHMM_prior_emiss, obtained by using the function prior_emiss_cat.")
    }
    if (emiss_hyp_prior$gen$m != m){
      stop("The number of states specified in m is not equal to the number of states specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
    }
    if (emiss_hyp_prior$gen$n_dep != n_dep){
      stop("The number of dependent variables specified in n_dep is not equal to the number of dependent variables specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
    }
    if (sum(emiss_hyp_prior$gen$q_emiss != q_emiss) > 0){
      stop("The number of number of observed categories for each of the dependent variable specified in q_emiss is not equal to q_emiss specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
    }
    if(is.null(emiss_hyp_prior$n_xx_emiss) & sum(nx[2:n_dep1] > 1) > 0){
      stop("Covariates were specified to predict the emission distribution(s), but no covariates were specified when creating the informative hyper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
    }
    if(!is.null(emiss_hyp_prior$n_xx_emiss)){
      if(sum(emiss_hyp_prior$n_xx_emiss != nx[2:n_dep1]) > 0){
        stop("The number of covariates specified to predict the emission distribution(s) is not equal to the number of covariates specified when creating the informative hper-prior distribution on the emission distribution(s) using the function prior_emiss_cat.")
      }
    }
    emiss_mu0	 <- emiss_hyp_prior$emiss_mu0
    emiss_nu   <- emiss_hyp_prior$emiss_nu
    emiss_V	   <- emiss_hyp_prior$emiss_V
    emiss_K0	 <- emiss_hyp_prior$emiss_K0
  }


# 四、 MCMC中间/输出变量设置和初始化(其实就是给框框，赋空值)  -------------------------------------------------------------------------------------------------------------
  # 定义MCMC迭代中间变量
  # overall
  c <- llk <- numeric(1)
  sample_path <- lapply(n_vary, dif_matrix, cols = J)
  trans <- rep(list(vector("list", m)), n_subj)

  # gamma
  gamma_int_mle_pooled <- gamma_pooled_ll <- vector("list", m)
  gamma_c_int <- rep(list(matrix(, n_subj, (m-1))), m)
  gamma_mu_int_bar <- gamma_V_int <- vector("list", m)
  gamma_mu_prob_bar <- rep(list(numeric(m)), m)
  gamma_naccept <- matrix(0, n_subj, m)

  # emiss
  cond_y <- lapply(rep(n_dep, n_subj), nested_list, m = m)
  emiss_int_mle_pooled <- emiss_pooled_ll <- rep(list(vector("list", n_dep)), m)
  emiss_c_int <- rep(list(lapply(q_emiss - 1, dif_matrix, rows = n_subj)), m)
  emiss_V_int <- rep(list(vector("list", n_dep)), m)
  emiss_mu_prob_bar <- rep(list(lapply(q_emiss, dif_vector)), m)
  emiss_naccept <- rep(list(matrix(0, n_subj, m)), n_dep)
  #browser()
  emiss_xt1 <- rep(list(vector("list", n_subj)), m)
  emiss_xt2 <- rep(list(vector("list", n_subj)), m)
  emiss_xt3 <- rep(list(vector("list", n_subj)), m)
  emiss_xt4 <- rep(list(vector("list", n_subj)), m)
  emiss_xt5 <- rep(list(vector("list", n_subj)), m)
  emiss_xt6 <- rep(list(vector("list", n_subj)), m)
  emiss_xt7 <- rep(list(vector("list", n_subj)), m)
  emiss_xt8 <- rep(list(vector("list", n_subj)), m)
  emiss_xt9 <- rep(list(vector("list", n_subj)), m)
  # 定义MCMC输出变量

  # 1定义PD/PD_subj
  if(length(start_val) != n_dep + 1){
    stop("The number of elements in the list start_val should be equal to 1 + the number of dependent variables,
         and should not contain nested lists (i.e., lists within lists)")
  }
  PD 					  <- matrix(, nrow = J, ncol = sum(m * q_emiss) + m * m + 1)
  PD_emiss_names   <- paste("q", 1, "_emiss", rep(1:q_emiss[1], m), "_S", rep(1:m, each = q_emiss[1]), sep = "")
  if(n_dep > 1){
    for(q in 2:n_dep){
      PD_emiss_names <- c(PD_emiss_names, paste("q", q, "_emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = ""))
    }
  }
  colnames(PD) 	<- c(PD_emiss_names, paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = ""), "LL")

  PD[1, ((sum(m * q_emiss) + 1)) :((sum(m * q_emiss) + m * m))] <- unlist(sapply(start_val, t))[1:(m*m)]
  PD[1, 1:((sum(m * q_emiss)))] <- unlist(sapply(start_val, t))[(m*m + 1): (m*m + sum(m * q_emiss))]

  PD_subj				<- rep(list(PD), n_subj)
  #browser()


  # 2定义gamma/emiss_prob_bar gamma/emiss_int_bar
  gamma_prob_bar		<- matrix(, nrow = J, ncol = (m * m))
  colnames(gamma_prob_bar) <- paste("S", rep(1:m, each = m), "toS", rep(1:m, m), sep = "")
  gamma_prob_bar[1,] <- PD[1,(sum(m*q_emiss) + 1):(sum(m * q_emiss) + m * m)]

  emiss_prob_bar			<- lapply(q_emiss * m, dif_matrix, rows = J)
  names(emiss_prob_bar) <- dep_labels
  for(q in 1:n_dep){
    colnames(emiss_prob_bar[[q]]) <- paste("Emiss", rep(1:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q]), sep = "")
    start <- c(0, q_emiss * m)
    emiss_prob_bar[[q]][1,] <- PD[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))]
  }

  gamma_int_bar				<- matrix(, nrow = J, ncol = ((m-1) * m))
  colnames(gamma_int_bar) <- paste("int_S", rep(1:m, each = m-1), "toS", rep(2:m, m), sep = "")
  gamma_int_bar[1,] <- as.vector(t(prob_to_int(matrix(gamma_prob_bar[1,], byrow = TRUE, ncol = m, nrow = m))))
  if(nx[1] > 1){
    gamma_cov_bar				<- matrix(, nrow = J, ncol = ((m-1) * m) * (nx[1] - 1))
    colnames(gamma_cov_bar) <- paste( paste("cov", 1 : (nx[1] - 1), "_", sep = ""), "S", rep(1:m, each = (m-1) * (nx[1] - 1)), "toS", rep(2:m, m * (nx[1] - 1)), sep = "")
    gamma_cov_bar[1,] <- 0
  } else{
    gamma_cov_bar <- "No covariates where used to predict the transition probability matrix"
  }

  emiss_int_bar			<- lapply((q_emiss-1) * m, dif_matrix, rows = J)
  names(emiss_int_bar) <- dep_labels
  for(q in 1:n_dep){
    colnames(emiss_int_bar[[q]]) <-  paste("int_Emiss", rep(2:q_emiss[q], m), "_S", rep(1:m, each = q_emiss[q] - 1), sep = "")
    emiss_int_bar[[q]][1,] <- as.vector(prob_to_int(matrix(emiss_prob_bar[[q]][1,], byrow = TRUE, ncol = q_emiss[q], nrow = m)))
  }
  if(sum(nx[-1]) > n_dep){
    emiss_cov_bar			<- lapply((q_emiss-1) * m * (nx[-1] - 1 ), dif_matrix, rows = J)
    names(emiss_cov_bar) <- dep_labels
    for(q in 1:n_dep){
      if(nx[1 + q] > 1){
        colnames(emiss_cov_bar[[q]]) <-  paste( paste("cov", 1 : (nx[1 + q] - 1), "_", sep = ""), "emiss", rep(2:q_emiss[q], m * (nx[1 + q] - 1)), "_S", rep(1:m, each = (q_emiss[q] - 1) * (nx[1 + q] - 1)), sep = "")
        emiss_cov_bar[[q]][1,] <- 0
      } else {
        emiss_cov_bar[[q]] <- "No covariates where used to predict the emission probabilities for this outcome"
      }
    }
  } else{
    emiss_cov_bar <- "No covariates where used to predict the emission probabilities"
  }
  #browser()
  # 特定subeject的截距，10个matrix
  gamma_int_subj			<- rep(list(gamma_int_bar), n_subj)
  emiss_int_subj			<- rep(list(emiss_int_bar), n_subj)

  # Put starting values in place for fist run forward algorithm
  emiss_sep 			<- vector("list", n_dep)
  for(q in 1:n_dep){
    start <- c(0, q_emiss * m)
    emiss_sep[[q]] <- matrix(PD[1,(sum(start[1:q]) + 1):(sum(start[1:q]) + (m * q_emiss[q]))], byrow = TRUE, ncol = q_emiss[q], nrow = m)
  }
  emiss				  <- rep(list(emiss_sep), n_subj)
  gamma 			<- rep(list(matrix(PD[1,(sum(m*q_emiss) + 1):(sum(m * q_emiss) + m * m)], byrow = TRUE, ncol = m)), n_subj)
  delta 			<- rep(list(solve(t(diag(m) - gamma[[1]] + 1), rep(1, m))), n_subj)


# 五、 Start analysis -------------------------------------------------------------------
  # Run the MCMC algorithm
  itime <- proc.time()[3]
  if(show_progress == TRUE){
    cat("Progress of the Bayesian mHMM algorithm:", "\n")
    pb <- utils::txtProgressBar(min = 2, max = J, style = 3)
  }

  for (iter in 2 : J){
    emiss_mu_int_bar1 <- emiss_mu_int_bar
    
   # 1⃣️ 每一个subject第J次迭代的隐变量序列 sampled state sequence------
    for(s in 1:n_subj){
      # 前向算法，获得subject specific forward proababilities 和 log likelihood
#browser()
      forward				<- cat_mult_fw_r_to_cpp(x = subj_data[[s]]$y, m = m, emiss = emiss[[s]], gamma = gamma[[s]], n_dep = n_dep, delta=NULL)
      alpha         <- forward[[1]]
      c             <- max(forward[[2]][, subj_data[[s]]$n_t])
      llk           <- c + log(sum(exp(forward[[2]][, subj_data[[s]]$n_t] - c)))
      PD_subj[[s]][iter, sum(m * q_emiss) + m * m + 1] <- llk

      # Using the forward probabilites, sample the state sequence in a backward manner.
      # In addition, saves state transitions in trans, and conditional observations within states in cond_y
      trans[[s]]					                  <- vector("list", m)
      sample_path[[s]][n_vary[[s]], iter] 	<- sample(1:m, 1, prob = c(alpha[, n_vary[[s]]])) # n_vary[[s]]=900 alpha

      #这里还是空值
      for(t in (subj_data[[s]]$n_t - 1):1){
        sample_path[[s]][t,iter] 	              <- sample(1:m, 1, prob = (alpha[, t] * gamma[[s]][,sample_path[[s]][t + 1, iter]]))
        trans[[s]][[sample_path[[s]][t,iter]]]	<- c(trans[[s]][[sample_path[[s]][t,iter]]], sample_path[[s]][t + 1, iter])
      }

      for (i in 1:m){
        trans[[s]][[i]] <- c(trans[[s]][[i]], 1:m)
        trans[[s]][[i]] <- rev(trans[[s]][[i]])
        #上面两行看不懂在干嘛
        

        for(q in 1:n_dep){
          cond_y[[s]][[i]][[q]] <- subj_data[[s]]$y[sample_path[[s]][, iter] == i, q]
        }
        #browser()
      emiss_xt1[[i]][[s]] <- xt[[1]][,s][sample_path[[s]][, iter] == i]
      emiss_xt2[[i]][[s]] <- xt[[2]][,s][sample_path[[s]][, iter] == i]
      emiss_xt3[[i]][[s]] <- xt[[3]][,s][sample_path[[s]][, iter] == i]
      emiss_xt4[[i]][[s]] <- xt[[4]][,s][sample_path[[s]][, iter] == i]
      emiss_xt5[[i]][[s]] <- xt[[5]][,s][sample_path[[s]][, iter] == i]
      emiss_xt6[[i]][[s]] <- xt[[6]][,s][sample_path[[s]][, iter] == i]
      emiss_xt7[[i]][[s]] <- xt[[7]][,s][sample_path[[s]][, iter] == i]
      emiss_xt8[[i]][[s]] <- xt[[8]][,s][sample_path[[s]][, iter] == i]
      emiss_xt9[[i]][[s]] <- xt[[9]][,s][sample_path[[s]][, iter] == i]
      }
    }





  # 2⃣️ State specific （不是subject） (The remainder of the mcmc algorithm)------------
    for(i in 1:m){
  # 1得到最大似然估计的协方差矩阵和likelihood值at subject and group level 作为参数
    # group-level
    # group-level参数 gamma_mle_pooled (gamma_int_mle_pooled,gamma_pooled_ll)
      
      trans_pooled			  <- factor(c(unlist(sapply(trans, "[[", i)), c(1:m)))
      #browser()
      gamma_mle_pooled		<- optim(gamma_int_mle0[i,], llmnl_int, Obs = trans_pooled, #gamma_int_mle0是初始值 llmnl_int是函数
                                   n_cat = m, method = "BFGS", hessian = FALSE,
                                   control = list(fnscale = -1))
      gamma_int_mle_pooled[[i]]  <- gamma_mle_pooled$par
      gamma_pooled_ll[[i]]			<- gamma_mle_pooled$value

     # group-level参数 emiss_mle_pooled (emiss_int_mle_pooled,emiss_pooled_ll)
      for(q in 1:n_dep){
        cond_y_pooled					      <- numeric()
        emiss_xt1_pooled					      <- numeric()
        emiss_xt2_pooled					      <- numeric()
        emiss_xt3_pooled					      <- numeric()
        emiss_xt4_pooled					      <- numeric()
        emiss_xt5_pooled					      <- numeric()
        emiss_xt6_pooled					      <- numeric()
        emiss_xt7_pooled					      <- numeric()
        emiss_xt8_pooled					      <- numeric()
        emiss_xt9_pooled					      <- numeric()
        ### MOET OOK ECHT BETER KUNNEN, eerst # cond_y_pooled				<- unlist(sapply(cond_y, "[[", m))
        for(s in 1:n_subj){
          cond_y_pooled             <- c(cond_y_pooled, cond_y[[s]][[i]][[q]])
          emiss_xt1_pooled           <- c(emiss_xt1_pooled, emiss_xt1[[i]][[s]])
          emiss_xt2_pooled           <- c(emiss_xt2_pooled, emiss_xt2[[i]][[s]])
          emiss_xt3_pooled           <- c(emiss_xt3_pooled, emiss_xt3[[i]][[s]])
          emiss_xt4_pooled           <- c(emiss_xt4_pooled, emiss_xt4[[i]][[s]])
          emiss_xt5_pooled           <- c(emiss_xt5_pooled, emiss_xt5[[i]][[s]])
          emiss_xt6_pooled           <- c(emiss_xt6_pooled, emiss_xt6[[i]][[s]])
          emiss_xt7_pooled           <- c(emiss_xt7_pooled, emiss_xt7[[i]][[s]])
          emiss_xt8_pooled           <- c(emiss_xt8_pooled, emiss_xt8[[i]][[s]])
          emiss_xt9_pooled           <- c(emiss_xt9_pooled, emiss_xt9[[i]][[s]])
        }
        emiss_xt_pooled <- rep(list(0), 9)
        emiss_xt_pooled[[1]] <- emiss_xt1_pooled
        emiss_xt_pooled[[2]] <- emiss_xt2_pooled
        emiss_xt_pooled[[3]] <- emiss_xt3_pooled
        emiss_xt_pooled[[4]] <- emiss_xt4_pooled
        emiss_xt_pooled[[5]] <- emiss_xt5_pooled
        emiss_xt_pooled[[6]] <- emiss_xt6_pooled
        emiss_xt_pooled[[7]] <- emiss_xt7_pooled
        emiss_xt_pooled[[8]] <- emiss_xt8_pooled
        emiss_xt_pooled[[9]] <- emiss_xt9_pooled

        # emiss_mle_pooled		<- optim(emiss_int_mle0[[q]][i,], llmnl_int, Obs = c(cond_y_pooled, c(1:q_emiss[q])),
        #                              n_cat = q_emiss[q], method = "BFGS", hessian = FALSE,
        #                              control = list(fnscale = -1))
        #browser()
        obs_all <- rep(list(0), 3)
        obs_all[[1]] <- cond_y_pooled
        obs_all[[2]] <- emiss_xt_pooled
        obs_all[[3]] <- emiss_mu_int_bar1[[i]][[q]]
        #browser()
        #a <- llmnl_int2(int = obs_all[[3]][1,], Obs = obs_all, n_cat = q_emiss[q])
        emiss_mle_pooled		<- optim(emiss_int_mle0[[q]][i,], llmnl_int2, Obs = obs_all,
                                   n_cat = q_emiss[q], method = "Nelder-Mead", hessian = FALSE,
                                   control = list(fnscale = -1))
        emiss_int_mle_pooled[[i]][[q]]  <- emiss_mle_pooled$par
        emiss_pooled_ll[[i]][[q]]				<- emiss_mle_pooled$value
      }




    # subject-level
      for (s in 1:n_subj){
        wgt 				<- subj_data[[s]]$n_t / n_total

    # subject level参数 gamma_out->subj_data$gamma_mhess
        gamma_out					<- optim(gamma_int_mle_pooled[[i]], llmnl_int_frac, Obs = c(trans[[s]][[i]], c(1:m)),
                                 n_cat = m, pooled_likel = gamma_pooled_ll[[i]], w = gamma_w, wgt = wgt,
                                 method="BFGS", hessian = FALSE, control = list(fnscale = -1))
        if(gamma_out$convergence == 0){
          subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<-
            mnlHess_int(int = gamma_out$par, Obs = c(trans[[s]][[i]], c(1:m)), n_cat =  m)
        } else {
          subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ]	<- diag(m-1)
        }
        # if this is first iteration, use MLE for current values RW metropolis sampler
        if (iter == 2){
          gamma_c_int[[i]][s,]		<- gamma_out$par
        }

    # subject level参数 emiss_out->subj_data$emiss_mhess
        for(q in 1:n_dep){
          emiss_out				<- optim(emiss_int_mle_pooled[[i]][[q]], llmnl_int_frac, Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])),
                                  n_cat = q_emiss[q], pooled_likel = emiss_pooled_ll[[i]][[q]],
                                  w = emiss_w[[q]], wgt = wgt, method = "Nelder-Mead", hessian = FALSE, control = list(fnscale = -1))
          if(emiss_out$convergence == 0){
            subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]		<-
              mnlHess_int(int = subj_data[[s]]$emiss_int_mle[[q]][i,], Obs = c(cond_y[[s]][[i]][[q]], c(1:q_emiss[q])), n_cat =  q_emiss[q])
          } else {
            subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ]	<- diag(q_emiss[q] - 1)
          }
          # if this is first iteration, use MLE for current values RW metropolis sampler
          if (iter == 2){
            emiss_c_int[[i]][[q]][s,]	<- emiss_out$par
          }
        }
      }
#browser()

   # 2 用上述的参数 使用Gibbs sampler 对group-level的TM和EM进行采样
    #得到gamma_mu_int_bar和gamma_mu_prob_bar
     # gamma_mu0_n and gamma_mu_int_bar are matrices,
     # with the number of rows equal to the number of covariates, and ncol equal to number of intercepts estimated
#browser()
      gamma_mu0_n           <-   # 用矩阵乘法对10个subject数据加和，得到新的超先验均值
        solve(t(xx[[1]]) %*% xx[[1]] + gamma_K0) %*% #t是矩阵转置的意思 solve表示取倒数 %*%是矩阵乘积
        (t(xx[[1]]) %*% gamma_c_int[[i]] + gamma_K0 %*% gamma_mu0[[i]])

      gamma_V_n             <-   # 与subject单位的差值作为协方差矩阵
        gamma_V +
        t(gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) %*% (gamma_c_int[[i]] - xx[[1]] %*% gamma_mu0_n) +
        t(gamma_mu0_n - gamma_mu0[[i]]) %*% gamma_K0 %*% (gamma_mu0_n - gamma_mu0[[i]])

      gamma_V_int[[i]]      <-   # Wishart分布的随机生成，S表示逆比例矩阵，v表示自由度
        solve(rwish(S = solve(gamma_V_n), v = gamma_nu + n_subj)) # 可以理解为是在V_n的基础上采样迭代出了V_int

      gamma_mu_int_bar[[i]] <-   # 从随机生成的方差中提取截距
        gamma_mu0_n +
        solve(chol(t(xx[[1]]) %*% xx[[1]] + gamma_K0)) %*%
        matrix(rnorm((m - 1) * nx[1]), nrow = nx[1]) %*%
        t(solve(chol(solve(gamma_V_int[[i]]))))

      gamma_exp_int				  <-
        matrix(exp(c(0, gamma_mu_int_bar[[i]][1,] )), nrow  = 1)

      gamma_mu_prob_bar[[i]] 	<-  # 从int_bar到这里只是做计算而已 截距转概率
        gamma_exp_int / as.vector(gamma_exp_int %*% c(rep(1,(m))))

    #得到emiss_mu_int_bar和emiss_mu_prob_bar
#browser()
      for(q in 1:n_dep){
        emiss_mu0_n                 <-
          solve(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]]) %*%
          (t(xx[[1 + q]]) %*% emiss_c_int[[i]][[q]] +
             emiss_K0[[q]] %*% emiss_mu0[[q]][[i]])

        emiss_V_n                   <-
          emiss_V[[q]] +
          t(emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) %*%
            (emiss_c_int[[i]][[q]] - xx[[1 + q]] %*% emiss_mu0_n) +
          t(emiss_mu0_n - emiss_mu0[[q]][[i]]) %*% emiss_K0[[q]] %*%
            (emiss_mu0_n - emiss_mu0[[q]][[i]])

        emiss_V_int[[i]][[q]]       <-
          solve(rwish(S = solve(emiss_V_n), v = emiss_nu[[q]] + n_subj))

        emiss_mu_int_bar[[i]][[q]]	 <-
          emiss_mu0_n +
          solve(chol(t(xx[[1 + q]]) %*% xx[[1 + q]] + emiss_K0[[q]])) %*%
          matrix(rnorm((q_emiss[q] - 1) * nx[1 + q]), nrow = nx[1 + q]) %*%
          t(solve(chol(solve(emiss_V_int[[i]][[q]]))))

        emiss_exp_int				       <-
          matrix(exp(c(0, emiss_mu_int_bar[[i]][[q]][1, ])), nrow  = 1)

        emiss_mu_prob_bar[[i]][[q]] <-
          emiss_exp_int / as.vector(emiss_exp_int %*% c(rep(1, (q_emiss[q]))))
      }

##browser()

  # 3 用上述的参数 使用RW Metropolis sampler 对subject-level的TM和EM进行采样
    #得到gamma_int_subj
      for (s in 1:n_subj){

        gamma_candcov_comb 			<-
          chol2inv(
            chol(
              subj_data[[s]]$gamma_mhess[(1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1)), ] + chol2inv( chol( gamma_V_int[[i]] ) )     ))

        gamma_RWout					    <- mnl_RW_once(int1 = gamma_c_int[[i]][s,], Obs = trans[[s]][[i]], n_cat = m, mu_int_bar1 = c(t(gamma_mu_int_bar[[i]]) %*% xx[[1]][s,]), V_int1 = gamma_V_int[[i]], scalar = gamma_scalar, candcov1 = gamma_candcov_comb)
        gamma[[s]][i,]  	<- PD_subj[[s]][iter, c((sum(m * q_emiss) + 1 + (i - 1) * m):(sum(m * q_emiss) + (i - 1) * m + m))] <- gamma_RWout$prob + .0001
        gamma_naccept[s, i]			<- gamma_naccept[s, i] + gamma_RWout$accept
        gamma_c_int[[i]][s,]		<- gamma_RWout$draw_int
        gamma_int_subj[[s]][iter, (1 + (i - 1) * (m - 1)):((m - 1) + (i - 1) * (m - 1))] <- gamma_c_int[[i]][s,]
    #得到emiss_int_subj
        start <- c(0, q_emiss * m)

        for(q in 1:n_dep){
          emiss_xt <- rep(list(0), 9)
          emiss_xt[[1]] <- emiss_xt1[[i]][[s]]
          emiss_xt[[2]] <- emiss_xt2[[i]][[s]]
          emiss_xt[[3]] <- emiss_xt3[[i]][[s]]
          emiss_xt[[4]] <- emiss_xt4[[i]][[s]]
          emiss_xt[[5]] <- emiss_xt5[[i]][[s]]
          emiss_xt[[6]] <- emiss_xt6[[i]][[s]]
          emiss_xt[[7]] <- emiss_xt7[[i]][[s]]
          emiss_xt[[8]] <- emiss_xt8[[i]][[s]]
          emiss_xt[[9]] <- emiss_xt9[[i]][[s]]
          emiss_candcov_comb		     <- chol2inv(chol(subj_data[[s]]$emiss_mhess[[q]][(1 + (i - 1) * (q_emiss[q] - 1)):((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1)), ] + chol2inv(chol(emiss_V_int[[i]][[q]]))))
#browser()
          emiss_RWout				       <-
            mnl_RW_once2(int1 = emiss_c_int[[i]][[q]][s,],
                        Obs = cond_y[[s]][[i]][[q]],
                        emiss_xt = emiss_xt,
                        beta_cov = emiss_mu_int_bar1[[i]][[q]],
                        mu_int_bar1 = c(t(emiss_mu_int_bar[[i]][[q]]) %*% xx[[1 + q]][s,]),
                        n_cat = q_emiss[q],
                        V_int1 = emiss_V_int[[i]][[q]],
                        scalar = emiss_scalar[[q]],
                        candcov1 = emiss_candcov_comb)

          emiss[[s]][[q]][i,]		   <- PD_subj[[s]][iter, (sum(start[1:q]) + 1 + (i - 1) * q_emiss[q]):(sum(start[1:q]) + (i - 1) * q_emiss[q] + q_emiss[q])]<- emiss_RWout$prob + .0001

          emiss_naccept[[q]][s, i]	 <- emiss_naccept[[q]][s, i] + emiss_RWout$accept
          emiss_c_int[[i]][[q]][s,] <- emiss_RWout$draw_int
          emiss_int_subj[[s]][[q]][iter, (1 + (i - 1) * (q_emiss[q] - 1)) : ((q_emiss[q] - 1) + (i - 1) * (q_emiss[q] - 1))]	<- emiss_c_int[[i]][[q]][s,]
        }

        if(i == m){
          delta[[s]] 		<- solve(t(diag(m) - gamma[[s]] + 1), rep(1, m))
        }
      }
    }
  #browser()


  # 3⃣️ End of 1 MCMC iteration, save output values --------
    gamma_int_bar[iter, ]				   	<- unlist(lapply(gamma_mu_int_bar, "[",1,))
    if(nx[1] > 1){
      gamma_cov_bar[iter, ]      	<- unlist(lapply(gamma_mu_int_bar, "[",-1,))
    }
#browser()
    gamma_prob_bar[iter,]			<- unlist(gamma_mu_prob_bar)
    for(q in 1:n_dep){
      emiss_int_bar[[q]][iter, ]	<- as.vector(unlist(lapply(
        lapply(emiss_mu_int_bar, "[[", q), "[",1,)
      ))
      if(nx[1+q] > 1){
        emiss_cov_bar[[q]][iter, ]  <- as.vector(unlist(lapply(
          lapply(emiss_mu_int_bar, "[[", q), "[",-1,)
        ))
      }
      emiss_prob_bar[[q]][iter,]	<- as.vector(unlist(sapply(emiss_mu_prob_bar, "[[", q)))
    }
    if(show_progress == TRUE){
      utils::setTxtProgressBar(pb, iter)
    }
    #browser()

  }# 这个是迭代循环J的大括号

  browser()
  
  if(show_progress == TRUE){
     close(pb)
  }

# End of function, return output values --------
  ctime = proc.time()[3]
  message(paste("Total time elapsed (hh:mm:ss):", hms(ctime-itime)))
  if(return_path == TRUE){
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, emiss_int_bar = emiss_int_bar,
                emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                emiss_prob_bar = emiss_prob_bar, gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept,
                sample_path = sample_path)
  } else {
    out <- list(input = list(m = m, n_dep = n_dep, q_emiss = q_emiss, J = J,
                             burn_in = burn_in, n_subj = n_subj, n_vary = n_vary, dep_labels = dep_labels),
                PD_subj = PD_subj, gamma_int_subj = gamma_int_subj, emiss_int_subj = emiss_int_subj,
                gamma_int_bar = gamma_int_bar, gamma_cov_bar = gamma_cov_bar, emiss_int_bar = emiss_int_bar,
                emiss_cov_bar = emiss_cov_bar, gamma_prob_bar = gamma_prob_bar,
                emiss_prob_bar = emiss_prob_bar, gamma_naccept = gamma_naccept, emiss_naccept = emiss_naccept)
  }
  class(out) <- append(class(out), "mHMM")
  return(out)
  }


