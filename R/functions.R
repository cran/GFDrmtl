library(lpSolve)
library(GFDrmst)

# show which groups correspond to which factor values
formula2groups <- function(formula, event = "event", data, cens){
  formula2 <- paste0(formula, "*", event)
  dat <- model.frame(formula2, data)
  nf <- ncol(dat) - 2
  nadat <- names(dat)
  if (anyNA(data[, nadat])) {
    stop("Data contains NAs!")
  }
  names(dat) <- c("time", nadat[2:(1 + nf)], "event")
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
  }
  lev_names <- expand.grid(levels)
  if(nf > 1) lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),]

  event2 <- dat$event
  censored <- (event2 == cens)
  dat[!censored, "event"] <- as.numeric(droplevels(factor(event2[!(censored)])))
  dat[censored, "event"]  <- 0
  M <- max(dat[, "event"])

  lev_names <- lev_names[rep(seq_len(nrow(lev_names)), each = M),]

  lev_names <- cbind(1:nrow(lev_names), lev_names, levels(droplevels(factor(event2[!(censored)]))))
  colnames(lev_names) <- c("Subgroup", colnames(dat)[2:(nf + 1)], event)

  return(lev_names)
}


# calculate the p value
p.value <- function(data_mat, teststat){
  data_order <- t( apply(data_mat, 1, sort) )
  pvalues <- numeric(nrow(data_mat))
  B <- ncol(data_mat)

  for(l in 1:length(teststat)){
    if(teststat[l] < data_order[l,1]){
      pvalues[l] <- 1
    }else{
      beta <- mean(data_mat[l,] >= teststat[l])
      if(beta < 1){
        x <- round((1-beta)*B)
        quants <- data_order[,x]
      }else{
        quants <- -Inf
      }
      pvalues[l] <- mean(apply(data_mat > quants, 2, any))
    }
  }
  pvalues
}

# is there a solution mu for the hypothesis H %*% mu = c in (0,tau]^k
existence <- function(H, c, tau, k){
  # set c >= 0
  H[c < 0,] <- -H[c < 0,]
  c <- abs(c)

  r <- nrow(H)
  kM <- ncol(H)

  Objective.in <- c(rep(0,kM),rep(1,r))
  zero_mat_r <- matrix(0, nrow=r, ncol=kM)
  zero_mat_kM <- matrix(0, nrow=kM, ncol=r)
  zero_mat_k <- matrix(0, nrow=k, ncol=r)
  Const.mat <- rbind(cbind(H,diag(r)),
                     cbind(diag(kM),zero_mat_kM),
                     cbind(kronecker(diag(k),t(rep(1,kM/k))),zero_mat_k),
                     cbind(zero_mat_r,diag(r)))
  Const.rhs <- c(c, numeric(kM), rep(tau,k), numeric(r))
  Const.dir <- c(rep("=",r),rep(">",kM),rep("<=",k),rep(">=",r))

  Optimum <- lp(direction="min",Objective.in,Const.mat,Const.dir,Const.rhs)
  Optimum$objval == 0
}

# function for DA_hat, DN_hat and Y
estimators <- function(X, D, tau, M = NULL, times = NULL, var = TRUE){
  if(is.null(M)) M <- max(D)
  if(is.null(times)) times <- sort(unique(X[D > 0]))
  if(length(times) == 0){
    return(list(Yt = numeric(0), DNt = matrix(ncol=M, nrow=0),
                DNjt = lapply(1:M,function(x) matrix(ncol=0, nrow=0)),
                DA_hat = matrix(ncol=M, nrow=0), F_hat = matrix(ncol=M, nrow=0),
                fm = matrix(ncol=M, nrow=0), gm = matrix(ncol=M, nrow=0),
                sqrt_fac = numeric(0), ind = numeric(0), diffs = numeric(0)))
  }else{
    Yt <- sapply(times, function(t) sum(X >= t))
    #Nt <- matrix(sapply(1:M, function(m) sapply(times, function(t) sum(X[D == m] <= t))), ncol = M) # M cols, times rows
    DNt <- matrix(0, nrow = length(times), ncol = M)
    # Loop through each m and t combination
    for (m in 1:M) {
      # Find indices where D == m
      indices <- (D == m)
      # Count occurrences of each unique value of X within indices
      counts <- table(X[indices])
      # Update the result matrix with the counts
      for (t_index in 1:length(times)) {
        t <- times[t_index]
        DNt[t_index, m] <- ifelse(t %in% names(counts), counts[[as.character(t)]], 0)
      }
    }
    if(var){ DNjt <- lapply(1:M, function(m){
      if(any(D==m)){ input <- sapply(times, function(t) (X[D == m] == t)) }else{ input <- NA }
      matrix(input,nrow = length(times),ncol=sum(D==m),byrow = TRUE)
    })
    }
    A_hat  <- matrix(apply(ifelse(is.na(DNt / Yt),0,DNt / Yt),2,cumsum),ncol=M)
    DA_hat <- A_hat - rbind(0,A_hat[-length(times),,drop=FALSE])

    S_hat  <- cumprod(1 - rowSums(DA_hat))
    S_hat_minus  <- c(1,S_hat[-length(times)])
    F_hat  <- matrix(apply(DA_hat*S_hat_minus,2,cumsum), ncol = M)

    ind <- (times <= tau)
    if(!any(ind)){
      return(list(Yt = Yt, DNt = DNt, DNjt = DNjt, DA_hat = DA_hat, F_hat = F_hat, fm = matrix(ncol=M, nrow=0), gm = matrix(ncol=M, nrow=0),
                  sqrt_fac = numeric(0), ind = ind, diffs = numeric(0)))
    }else{
      timestau <- times[ind]
      diffs <- diff(c(timestau,tau))

      if(!var){
        return(list(Yt = Yt, DNt = DNt, DA_hat = DA_hat, F_hat = F_hat, ind = ind, diffs = diffs))
      }else{

        F_hat_diffs <- F_hat[ind,,drop=FALSE]*diffs
        IntF_hat <- matrix(apply(F_hat_diffs,2,function(x) rev(cumsum(rev(x)))),ncol=M)
        #IntF_hat <- t(sapply(timestau,
        #                     function(t){
        #                       ind2 <- (timestau >= t)
        #                       colSums(F_hat_diffs[ind2,,drop=FALSE])
        #                     }))
        fm <- (tau-timestau)*matrix(sapply(1:M, function(m) 1 - rowSums(F_hat[ind,-m,drop=FALSE])),ncol=M) - IntF_hat
        gm <- (tau-timestau)*F_hat[ind,,drop=FALSE] - IntF_hat

        sqrt_fac <- numeric(sum(ind))
        nenner <- (1-rowSums(DA_hat[ind,,drop=FALSE]))
        sqrt_fac[nenner > 0]  <- 1/nenner[nenner > 0]

        return(list(Yt = Yt, DNt = DNt, DNjt = DNjt, DA_hat = DA_hat, F_hat = F_hat, fm = fm, gm = gm,
                    sqrt_fac = sqrt_fac, ind = ind, diffs = diffs))
      }}}
}

# function for the RMTL
RMTL <- function(X, D, tau, M = NULL, var = TRUE){

  # number of competing risks
  if(is.null(M)) M <- max(D)
  # number of observations
  #n <- length(X)

  # values of Y, N,.. at t_1,...,t_m
  temp <- estimators(X, D, tau, M, var = var)
  F_hat <- temp$F_hat
  ind <- temp$ind
  diffs <- temp$diffs

  #DF_hat <- F_hat - rbind(0,F_hat[-length(times),])
  mu_hat <- colSums(F_hat[ind,,drop=FALSE] * diffs) ## the RMTL estimator

  if(var){
    Yt <- temp$Yt
    DNt <- temp$DNt
    DA_hat <- temp$DA_hat

    fm <- temp$fm
    gm <- temp$gm
    fac <- (temp$sqrt_fac)^2

    # now: compute the covariance
    Sigma_hat <- matrix(NA, ncol = M, nrow = M)
    sigma_hat <- array(dim = c(sum(ind),M,M))

    for(m in 1:M){ # without factor n
      sigma_hat[,m,m] <- cumsum(DA_hat[ind,m] * (1 - DA_hat[ind,m]) / Yt[ind])
      if(m < M){for(m2 in (m+1):M){
        sigma_hat[,m,m2] <- sigma_hat[,m2,m] <- - cumsum(DA_hat[ind,m2] * DA_hat[ind,m] / Yt[ind])
      }
      }
    }

    Dsigma_hat <- sigma_hat
    Dsigma_hat[-1,,] <- sigma_hat[-1,,] - sigma_hat[-sum(ind),,]


    for(m in 1:M){
      Sigma_hat[m,m] <- sum(fac * (Dsigma_hat[,m,m] *  fm[,m]^2 +
                                     2*rowSums(Dsigma_hat[,m,-m,drop=FALSE]) * gm[,m] * fm[,m] +
                                     rowSums(Dsigma_hat[,-m,-m,drop=FALSE]) * gm[,m]^2), na.rm = TRUE)

      if(m < M){
        for(m2 in (m+1):M){
          Sigma_hat[m,m2] <-
            Sigma_hat[m2,m] <- sum(fac * (Dsigma_hat[,m,m2] * fm[,m] * fm[,m2] +
                                            rowSums(Dsigma_hat[,m,-m2,drop=FALSE]) * gm[,m2] * fm[,m] +
                                            rowSums(Dsigma_hat[,m2,-m,drop=FALSE]) * gm[,m] * fm[,m2] +
                                            rowSums(Dsigma_hat[,-m,-m2,drop=FALSE]) * gm[,m] * gm[,m2]), na.rm = TRUE)

        }}
    }
    return(list("rmtl" = mu_hat, "var_rmtl" = Sigma_hat))
  }else{return(list("rmtl" = mu_hat))}
}

# function for multiple critical values
crit_values <- function(data_mat, alpha = 0.05){
  n <- ncol(data_mat)
  dimension <- nrow(data_mat)

  # First forget the connection of the coordinates, and sort the values per coordinate
  data_order <- t( apply(data_mat, 1, sort) )

  j <- min(c(floor(n*alpha/dimension), n-1))
  j_upper <- n
  #nonDupl <- which(!duplicated(t(data_order))) # instead j+1, nonDupl[nonDupl > j][1] can be used, because there is no change in between (faster for discrete data)
  if(j != j_upper - 1){
    while(mean(apply(data_mat > data_order[,n - j - 1],2,any)) <= alpha){ # check if j is the maximum
      j <- j+1 # increase j
      if(j == j_upper-1) break # check if j is the maximum
    }}

  data_order[,n-j] # critical values
}

# if a formula is given, we need to extract the time, status, group vector from the formula + data
formula2input <- function(formula, event = "event", data){
  formula2 <- paste0(formula, "*", event)
  dat <- model.frame(formula2, data)
  subject <- 1:nrow(dat)
  n_all <- length(subject)
  formula <- as.formula(formula)
  nf <- ncol(dat) - 2
  nadat <- names(dat)
  if (anyNA(data[, nadat])) {
    stop("Data contains NAs!")
  }
  names(dat) <- c("time", nadat[2:(1 + nf)], "event")
  dat2 <- data.frame(dat, subject = subject)
  nadat2 <- nadat[-c(1, nf + 2)]
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
  }
  lev_names <- expand.grid(levels)

  if (nf == 1) {
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    hypo_matrices <- list(diag(fl) - matrix(1/fl, ncol = fl,
                                            nrow = fl))
    dat2$group <- rep(1:length(n), n)
  }else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),
    ]
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    dat2$group <- rep(1:length(n), n)
  }
  return(dat2[,c("time", "event", "group")])
}

# return vector of number of factor levels
num_factor_levels <- function(formula, event = "event", data){
  formula2 <- paste0(formula, "*", event)
  dat <- model.frame(formula2, data)
  subject <- 1:nrow(dat)
  n_all <- length(subject)
  formula <- as.formula(formula)
  nf <- ncol(dat) - 2
  nadat <- names(dat)
  if (anyNA(data[, nadat])) {
    stop("Data contains NAs!")
  }
  names(dat) <- c("time", nadat[2:(1 + nf)], "event")
  dat2 <- data.frame(dat, subject = subject)
  nadat2 <- nadat[-c(1, nf + 2)]
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
  }
  return(fl)
}

# function for getting the global hypothesis matrix from a list of single hypothesis matrices
global_mat <- function(c_mat, kM) matrix(unlist(c_mat), ncol = kM, byrow=TRUE)

# function for the parts of the studentized permutation test statistic
test_stat_perm_bonf <- function(perm_values, tau, mu_hat, n_group, M, kM){
  group <- perm_values$group
  rmtl <- numeric(kM)
  var <- matrix(0, ncol=kM,nrow=kM)

  for(i in 1:n_group) {
    values2 <- perm_values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M = M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    var[((i-1)*M +1):(i*M), ((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
  }

  return(list(vec1 = (rmtl- mu_hat), var = var))
  #calculate the pooled BS test statistic
  #vec <- c_mat %*% (rmtl- mu_hat)
  #return( t(vec) %*% MASS::ginv(c_mat %*% Sigma_hat %*% t(c_mat)) %*% vec )
  # the factor n is eliminated by the other factors
}

# Multiple RMTL asymptotic tests
# input:
# time, status, group:    vectors of the same length
# OR
# formula, event, data:   formula, event variable name and data set

# hyp_mat:                partitionized matrix as list, i.e. list of length L with matrices as entries
#                         alternatively, "center", "Dunnett", "Tukey" as 'type' in ?contr_mat can be used
# hyp_vec:                partitionized vector as list (default = zero vector)
# M:                      number of competing events (default = number of observed events)
# tau:                    end point of time window
# stepwise:               Option for closed testing procedure. If TRUE possibly more hypotheses can be rejected.
# alpha:                  level of significance
# Nres:                   number of resampling iterations
# seed:                   seed

# output:
# Object of class "GFDrmst"
# method: character of used method
# test_stat: vector of test statistics
# p.value: vector of adjusted p values
# res: list with hypothesis matrices, estimators, confidence intervals if the
#       matrices are vectors,  test statistics, critical values, decisions
# alpha: level of significance
RMTL.asymptotic.test <- function(time = NULL, status = NULL, group = NULL,
                                 formula = NULL, event = NULL, data = NULL,
                                 hyp_mat, hyp_vec = NULL, M = NULL,
                                 tau, stepwise = FALSE, alpha = 0.05,
                                 Nres = 4999, seed = 1){
  set.seed(seed)

  if(is.null(time) | is.null(status) | is.null(group)){
    dat2 <- formula2input(formula, event, data)
    time   <- dat2$time
    status <- dat2$event
    group  <- dat2$group
  }

  k <- length(unique(group))
  group <- factor(group, labels = 1:k)
  status <- as.numeric(status)
  if(is.null(M)) M <- max(status)
  kM <- k*M
  rmtl <- numeric(kM)
  Sigma_hat <- matrix(0, ncol=kM,nrow=kM)

  # some possible values for hyp_mat
  if(hyp_mat %in% c("center", "Dunnett", "Tukey")){
    Y <- kronecker(GFDmcv::contr_mat(k, type = hyp_mat),diag(M))
    X <- data.frame(t(Y))
    colnames(X) <- rownames(Y)
    hyp_mat <- as.list(X)
    hyp_mat <- lapply(hyp_mat, t)
  }else{
    if(hyp_mat== "2by2"){
      if(k != 4) stop("More/less than 4 groups is not possible for 2-by-2 design.")
      nadat <- c("", "factor A", "factor B")
      if(!is.null(formula)){
        if(length(num_factor_levels(formula, event, data)) != 2 ||
           !all.equal(as.integer(c(2,2)), num_factor_levels(formula, event, data))){
          stop("Input data set does not have two factors with two levels each.")
        }
        formula2 <- paste0(formula, "*", event)
        dat <- model.frame(formula2, data)
        nadat <- names(dat)
      }
      HA <- cbind(diag(M), diag(M), -diag(M), -diag(M))
      HB <- cbind(diag(M), -diag(M), diag(M), -diag(M))
      HAB <- cbind(diag(M), -diag(M), -diag(M), diag(M))
      hyp_mat <- list(HA, HB, HAB)
      names(hyp_mat) <- c(paste("main effect of",nadat[2]),
                          paste("main effect of",nadat[3]), "interaction effect")

      }else{
      if(hyp_mat == "2by2 cause-wisely"){
        if(k != 4) stop("More/less than 4 groups is not possible for 2-by-2 design.")
        nadat <- c("", "factor A", "factor B")
        if(!is.null(formula)){
          if(length(num_factor_levels(formula, event, data)) != 2 ||
             !all.equal(as.integer(c(2,2)), num_factor_levels(formula, event, data))){
            stop("Input data set does not have two factors with two levels each.")
          }
          formula2 <- paste0(formula, "*", event)
          dat <- model.frame(formula2, data)
          nadat <- names(dat)
        }
        HA <- cbind(diag(M), diag(M), -diag(M), -diag(M))
        HB <- cbind(diag(M), -diag(M), diag(M), -diag(M))
        HAB <- cbind(diag(M), -diag(M), -diag(M), diag(M))
        Y <- rbind(HA, HB, HAB)
        X <- data.frame(t(Y))
        hyp_mat <- as.list(X)
        hyp_mat <- lapply(hyp_mat, t)
        names(hyp_mat) <- c(paste("main effect of",nadat[2],"on cause",1:M),
                            paste("main effect of",nadat[3],"on cause",1:M),
                            paste("interaction effect on cause",1:M))
      }else{
    if(is.null(hyp_vec)){
      warning("hyp_vec is chosen as zero vector since it is unspecified!")
    }}}
  }
  if(is.matrix(hyp_mat)) hyp_mat <- list(hyp_mat)
  L <- length(hyp_mat)
  if(is.null(hyp_vec)){
    hyp_vec <- lapply(hyp_mat, function(x) rep(0, nrow(x)))
    #warning("hyp_vec is chosen as zero vector since it is unspecified!")
  }
  if(!is.list(hyp_vec)) hyp_vec <- list(hyp_vec)

  # check whether all hypotheses have a solution in (0,tau]
  solution <- logical(length(hyp_mat))
  for(i in 1:length(hyp_mat)){
    solution[i] <- existence(hyp_mat[[i]], hyp_vec[[i]], tau, k)
  }
  if(any(!solution)){
    stop(paste0("Hypothesis ", which(!solution),
                " does not have a possible solution. "))
  }

  #time <- factor(time, exclude = c(NA, NaN))
  values <- data.frame(X=time, D=status, group=group)

  for(i in 1:k) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    Sigma_hat[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
  }
  teststats <- sapply(1:L, function(c_ind){
    c_mat <- hyp_mat[[c_ind]]
    vec <- c_mat %*% rmtl - hyp_vec[[c_ind]]
    out <- t(vec) %*% MASS::ginv(c_mat %*% Sigma_hat %*% t(c_mat)) %*% vec
    if(out == 0) out <- ifelse(all(vec == 0), 0, Inf)
    return(out)
  })
  attr(teststats, "names") <- paste0("W_n(H_", 1:L, ", c_", 1:L, ")")

  # if all matrices are row vectors
  if(all(sapply(hyp_mat, function(x) nrow(x)) == 1) ){
    # calculate the scaling matrix
    D <- diag(lapply(hyp_mat, function(vec){
      MASS::ginv(sqrt(vec %*% Sigma_hat %*% t(vec))) #using the g-inverse to avoid problems with zero variances
    }),nrow=length(hyp_mat))
    H <- global_mat(hyp_mat, kM)
    # equicoordinate normal quantiles
    sigma <- (D%*%H%*%Sigma_hat%*%t(H)%*%D)
    # since errors occur when any of diag(sigma) = 0, we do not consider these components
    my_index <- (diag(sigma) != 0)
    if(any(my_index)){
      quant <- (mvtnorm::qmvnorm(1 - alpha , tail="both.tails", sigma=sigma[my_index, my_index])$quantile)
    }else{ quant <- 0 }
    pvalues <- numeric(L)

    res <- matrix(list() , nrow = L, ncol = 8)
    colnames(res) <- c("hyp_matrix", "estimator", "lwr_conf", "upr_conf", "test_stat", "critical value", "adj_pvalue", "decision")

    for(l in 1:L){
      if(sum(my_index) > 0){
        bound <- rep(sqrt(teststats[l]), sum(my_index))
        pvalues[l] <- (1-mvtnorm::pmvnorm(lower = -bound, upper = bound, sigma=sigma[my_index, my_index]))
        # calculate the simultaneous confidence intervals
        conf.int <- c(hyp_mat[[l]]%*%rmtl) + sqrt(c(hyp_mat[[l]]%*% Sigma_hat %*% t(hyp_mat[[l]]))) * quant * t(c(-1,1))
        if(!(my_index[l])) conf.int <- c(-Inf, Inf)
        # determine the smallest and largest possible value
        conf.int[1] <- max(conf.int[1], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, 0, tau))
        conf.int[2] <- min(conf.int[2], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, tau, 0))
        attr(conf.int, "conf.level") <- 1-alpha
      }else{
        pvalues[l] <- as.numeric(hyp_mat[[l]]%*%rmtl == hyp_vec[[l]])
        conf.int <- t(numeric(2))
        conf.int[1] <- hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, 0, tau)
        conf.int[2] <- hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, tau, 0)
        attr(conf.int, "conf.level") <- 1-alpha
      }
      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmtl
      res[[l,3]] <- conf.int[1]
      res[[l,4]] <- conf.int[2]
      res[[l,5]] <- teststats[l]
      res[[l,6]] <- quant^2
    }

    if(stepwise & (L > 1)){
      my_grid <- as.matrix(expand.grid(lapply(1:L, function(x) c(TRUE,FALSE))))
      # delete the FALSE ... FALSE row
      my_grid <- my_grid[rowSums(my_grid)>0,]
      pValues <- apply(my_grid, 1, function(ind){
        min(RMTL.asymptotic.test(time, status, group, hyp_mat = hyp_mat[ind],
                                 hyp_vec = hyp_vec[ind], M = M, tau = tau,
                                 stepwise = FALSE, alpha = alpha,
                                 Nres = Nres)$p.value)
      })

      # adjust p-values with closed testing procedure
      for(l in 1:L){
        pvalues[l] <- max(pValues[my_grid[,l]])
      }

      # delete the confidence interval and critical value
      res <- res[,-c(3,4,6)]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

    if(stepwise) method <- "Multiple asymptotic RMTL Wald-type tests with stepwise extension"
    if(!stepwise) method <- "Multiple asymptotic RMTL Wald-type tests"

    out <- list(method = method,
                test_stat = teststats,
                p.value = pvalues,
                res = res,
                alpha = alpha)
    class(out) <- "GFDrmst"

    return(out)

  }else{
    if(L == 1){

      pvalues <- stats::pchisq(teststats, df = qr(hyp_mat[[1]])$rank, lower.tail = FALSE)

      res <- matrix(list() , nrow = 1, ncol = 6)
      colnames(res) <- c("hyp_matrix", "estimator", "test_stat", "critical value", "adj_pvalue", "decision")
      cv <- stats::qchisq(1-alpha, df = qr(hyp_mat[[1]])$rank)


        # save the needed values for res
        res[[1,1]] <- hyp_mat[[1]]
        res[[1,2]] <- hyp_mat[[1]]%*%rmtl
        res[[1,3]] <- teststats[1]
        res[[1,4]] <- cv[1]


      res[,"adj_pvalue"] <- pvalues
      res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

      method <- "Asymptotic RMTL Wald-type test"

      out <- list(method = method,
                  test_stat = teststats,
                  p.value = pvalues,
                  res = res,
                  alpha = alpha)
      class(out) <- "GFDrmst"

      return(out)
    }else{
    random_numbers <- mvtnorm::rmvnorm(Nres, sigma = Sigma_hat)

    # determine values for approximating the limiting distribution
    random_values <- t(sapply(1:length(hyp_mat), function(mat_ind) apply(random_numbers, 1, function(z){
      vec <- hyp_mat[[mat_ind]]%*%z
      t(vec)%*%MASS::ginv(hyp_mat[[mat_ind]]%*% Sigma_hat %*%t(hyp_mat[[mat_ind]]))%*%vec
    } )))
    pvalues <- p.value(data_mat =  random_values, teststat=teststats)

    res <- matrix(list() , nrow = L, ncol = 6)
    colnames(res) <- c("hyp_matrix", "estimator", "test_stat", "critical value", "adj_pvalue", "decision")
    cv <- crit_values(random_values, alpha = alpha)

    for(l in 1:L){
      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmtl
      res[[l,3]] <- teststats[l]
      res[[l,4]] <- cv[l]

    }

    if(stepwise & (L > 1)){
      my_grid <- as.matrix(expand.grid(lapply(1:L, function(x) c(TRUE,FALSE))))
      # delete the FALSE ... FALSE row
      my_grid <- my_grid[rowSums(my_grid)>0,]
      pValues <- apply(my_grid, 1, function(ind){
        min(RMTL.asymptotic.test(time, status, group, hyp_mat = hyp_mat[ind],
                                 hyp_vec = hyp_vec[ind], tau = tau,
                                 stepwise = FALSE, alpha = alpha,
                                 Nres = Nres)$p.value)
      })

      # adjust p-values with closed testing procedure
      for(l in 1:L){
        pvalues[l] <- max(pValues[my_grid[,l]])
      }

      # delete the critical value
      res <- res[,-4]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

    if(stepwise) method <- "Multiple asymptotic RMTL Wald-type tests with stepwise extension"
    if(!stepwise) method <- "Multiple asymptotic RMTL Wald-type tests"

    out <- list(method = method,
                test_stat = teststats,
                p.value = pvalues,
                res = res,
                alpha = alpha)
    class(out) <- "GFDrmst"

    return(out)
  }}
}

# Multiple RMTL permutation tests
RMTL.permutation.test <- function(time = NULL, status = NULL, group = NULL,
                                 formula = NULL, event = NULL, data = NULL,
                                 hyp_mat, hyp_vec = NULL, M = NULL,
                                 tau, stepwise = FALSE, alpha = 0.05,
                                 Nres = 4999, seed = 1){
  set.seed(seed)

  if(is.null(time) | is.null(status) | is.null(group)){
    dat2 <- formula2input(formula, event, data)
    time   <- dat2$time
    status <- dat2$event
    group  <- dat2$group
  }

  k <- length(unique(group))
  group <- factor(group, labels = 1:k)
  status <- as.numeric(status)
  if(is.null(M)) M <- max(status)
  kM <- k*M
  rmtl <- numeric(kM)
  Sigma_hat <- matrix(0, ncol=kM,nrow=kM)

  # some possible values for hyp_mat
  if(hyp_mat %in% c("center", "Dunnett", "Tukey")){
    Y <- kronecker(GFDmcv::contr_mat(k, type = hyp_mat),diag(M))
    X <- data.frame(t(Y))
    colnames(X) <- rownames(Y)
    hyp_mat <- as.list(X)
    hyp_mat <- lapply(hyp_mat, t)
  }else{
    if(hyp_mat== "2by2"){
      if(k != 4) stop("More/less than 4 groups is not possible for 2-by-2 design.")
      nadat <- c("", "factor A", "factor B")
      if(!is.null(formula)){
        if(length(num_factor_levels(formula, event, data)) != 2 ||
           !all.equal(as.integer(c(2,2)), num_factor_levels(formula, event, data))){
          stop("Input data set does not have two factors with two levels each.")
        }
        formula2 <- paste0(formula, "*", event)
        dat <- model.frame(formula2, data)
        nadat <- names(dat)
      }
      HA <- cbind(diag(M), diag(M), -diag(M), -diag(M))
      HB <- cbind(diag(M), -diag(M), diag(M), -diag(M))
      HAB <- cbind(diag(M), -diag(M), -diag(M), diag(M))
      hyp_mat <- list(HA, HB, HAB)
      names(hyp_mat) <- c(paste("main effect of",nadat[2]),
                            paste("main effect of",nadat[3]), "interaction effect")
      }else{
      if(hyp_mat == "2by2 cause-wisely"){
        if(k != 4) stop("More/less than 4 groups is not possible for 2-by-2 design.")
        nadat <- c("", "factor A", "factor B")
        if(!is.null(formula)){
          if(length(num_factor_levels(formula, event, data)) != 2 ||
             !all.equal(as.integer(c(2,2)), num_factor_levels(formula, event, data))){
            stop("Input data set does not have two factors with two levels each.")
          }
          formula2 <- paste0(formula, "*", event)
          dat <- model.frame(formula2, data)
          nadat <- names(dat)
        }
        HA <- cbind(diag(M), diag(M), -diag(M), -diag(M))
        HB <- cbind(diag(M), -diag(M), diag(M), -diag(M))
        HAB <- cbind(diag(M), -diag(M), -diag(M), diag(M))
        Y <- rbind(HA, HB, HAB)
        X <- data.frame(t(Y))
        hyp_mat <- as.list(X)
        hyp_mat <- lapply(hyp_mat, t)
        names(hyp_mat) <- c(paste("main effect of",nadat[2],"on cause",1:M),
                            paste("main effect of",nadat[3],"on cause",1:M),
                            paste("interaction effect on cause",1:M))
      }else{
    if(is.null(hyp_vec)){
      warning("hyp_vec is chosen as zero vector since it is unspecified!")
    }}}
  }
  if(is.matrix(hyp_mat)) hyp_mat <- list(hyp_mat)
  L <- length(hyp_mat)
  if(is.null(hyp_vec)){
    hyp_vec <- lapply(hyp_mat, function(x) rep(0, nrow(x)))
    #warning("hyp_vec is chosen as zero vector since it is unspecified!")
  }
  if(!is.list(hyp_vec)) hyp_vec <- list(hyp_vec)

  # check whether all matrices are contrast matrices
  if(any(sapply(hyp_mat,
                function(x) is.character(all.equal(x %*% kronecker(rep(1,k),
                                                                   diag(M)),
                                                   matrix(0,nrow=nrow(x),ncol=M)))))){
    stop("Permutation approach is only valid for matrices with contrast property in terms of the different groups.")
  }

  # check whether all hypotheses have a solution in (0,tau]
  solution <- logical(length(hyp_mat))
  for(i in 1:length(hyp_mat)){
    solution[i] <- existence(hyp_mat[[i]], hyp_vec[[i]], tau, k)
  }
  if(any(!solution)){
    stop(paste0("Hypothesis ", which(!solution),
                " does not have a possible solution. "))
  }

  #time <- factor(time, exclude = c(NA, NaN))
  values <- data.frame(X=time, D=status, group=group)

  for(i in 1:k) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    Sigma_hat[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
  }
  teststats <- sapply(1:L, function(c_ind){
    c_mat <- hyp_mat[[c_ind]]
    vec <- c_mat %*% rmtl - hyp_vec[[c_ind]]
    out <- t(vec) %*% MASS::ginv(c_mat %*% Sigma_hat %*% t(c_mat)) %*% vec
    if(out == 0) out <- ifelse(all(vec == 0), 0, Inf)
    return(out)
  })
  attr(teststats, "names") <- paste0("W_n(H_", 1:L, ", c_", 1:L, ")")

  #calculate mu_hat for the pooled sample
  temp <- RMTL(values$time, values$status, tau, M = M, var = FALSE)
  mu_hat <- temp[["rmtl"]]

  # generate permuted samples
  perm_samples <- replicate(Nres,{
    perm_values <- values[,c("X", "D")]
    perm_values$group <- sample(group)
    perm_values
  }, simplify = FALSE)

  erg <- sapply(perm_samples, function(perm_values){
    perms <- test_stat_perm_bonf(perm_values=perm_values, tau=tau,
                                 mu_hat=mu_hat, n_group = k, M = M, kM = kM)
    # calculate the two pooled BS test statistics for all matrices
      erg_int <- sapply(hyp_mat, function(c_mat){
        vec <- c_mat %*% perms$vec1
        return( t(vec) %*% MASS::ginv(c_mat %*% perms$var %*% t(c_mat)) %*% vec )
      })
    erg_int
  })
  dim(erg) <- c(L,Nres)

  ### ab hier weiter #...
  # unadjusted p values
  pvalues <- rowMeans(erg > teststats)
  # Bonferroni adjustment
  if(!stepwise) pvalues <- p.adjust(pvalues, method = "bonferroni")
  # Holm adjustment
  if(stepwise) pvalues <- p.adjust(pvalues, method = "holm")

  # if all matrices are vectors
  if(all(sapply(hyp_mat, function(x) nrow(x)) == 1)){
    res <- matrix(list() , nrow = L, ncol = 8)
    colnames(res) <- c("hyp_matrix", "estimator", "lwr_conf", "upr_conf",
                       "test_stat", "critical value", "adj_pvalue", "decision")

    quant <- numeric(L)
    for(l in 1:L){
      quant[l] <- crit_values(t(erg[l,]), alpha = alpha/L)
    }

    for(l in 1:L){
      # determine simultaneous 1-alpha confidence intervals
      conf.int <- c(hyp_mat[[l]]%*%rmtl) + sqrt(c(hyp_mat[[l]]%*% Sigma_hat %*% t(hyp_mat[[l]])) * quant[l]) * c(-1,1)
      if(c(hyp_mat[[l]]%*% Sigma_hat %*% t(hyp_mat[[l]])) == 0) conf.int <- c(-Inf, Inf)
      # determine the smallest and largest possible value
      conf.int[1] <- max(conf.int[1], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, 0, tau))
      conf.int[2] <- min(conf.int[2], hyp_mat[[l]] %*% ifelse(t(hyp_mat[[l]]) >= 0, tau, 0))
      attr(conf.int, "conf.level") <- 1-alpha

      # save the needed values for res
      res[l,3:4] <- conf.int
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmtl
      res[[l,5]] <- teststats[l]
      res[[l,6]] <- quant[l]
    }
    if(stepwise & (L > 1)){# delete the confidence interval and critical value
      res <- res[,-c(3,4,6)]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")


  }else{
    res <- matrix(list() , nrow = L, ncol = 6)
    colnames(res) <- c("hyp_matrix", "estimator", "test_stat", "critical value", "adj_pvalue", "decision")

    quant <- numeric(L)
    for(l in 1:L){
      quant[l] <- crit_values(t(erg[l,]), alpha = alpha/L)
    }

    for(l in 1:L){

      # save the needed values for res
      res[[l,1]] <- hyp_mat[[l]]
      res[[l,2]] <- hyp_mat[[l]]%*%rmtl
    }

    res[,3] <- teststats
    res[,4] <- quant

    if(stepwise & (L > 1)){
      # delete the critical value
      res <- res[,-4]
    }

    res[,"adj_pvalue"] <- pvalues
    res[,"decision"] <- ifelse(pvalues <= alpha, "H1", "not significant")

  }

  if(stepwise) method <- "Permutation RMTL Wald-type tests with Holm correction"
  if(!stepwise) method <- "Permutation RMTL Wald-type tests with Bonferroni correction"

  out <- list(method = method,
              test_stat = teststats,
              p.value = pvalues,
              res = res,
              alpha = alpha)
  class(out) <- "GFDrmst"
  return(out)
}
