# one function for all methods
RMTL.test <- function(
    time = NULL,
    status = NULL,
    group = NULL,
    formula = NULL,
    event = NULL,
    data = NULL,
    hyp_mat,
    hyp_vec = NULL,
    M = NULL,
    tau,
    method = c("permutation", "asymptotic"),
    stepwise = FALSE,
    alpha = 0.05,
    Nres = 4999,
    seed = 1
){
  method <- match.arg(method)
  switch(method,
         permutation = RMTL.permutation.test(time, status, group, formula, event,
                                             data, hyp_mat, hyp_vec, M, tau, stepwise,
                                             alpha, Nres, seed),
         asymptotic  = RMTL.asymptotic.test(time, status, group, formula, event,
                                            data, hyp_mat, hyp_vec, M, tau, stepwise,
                                            alpha, Nres, seed))
}
