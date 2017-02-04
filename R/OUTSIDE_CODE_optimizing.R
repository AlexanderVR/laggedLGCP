#optimizing_new_method <- function (object, max_block_size = 50, draws=0, ...) {
    # Use rstan's optimizing function, but use faster (approximate)
    # method to take samples from N(0, -H^{-1}), where
    # H is the Hessian of log_p at the calculated (local) mode.
    # Assume that H is block-diagonal with blocks of size "max_block_size"

setGeneric(name = 'optimizing_new_method',
           def = function(object, ...) {
             standardGeneric("optimizing_new_method")
})


setMethod("optimizing_new_method", "stanmodel",
          function(object, data = list(), 
                   seed = sample.int(.Machine$integer.max, 1),
                   init = 'random', check_data = TRUE, sample_file = NULL, 
                   algorithm = c("LBFGS", "BFGS", "Newton"),
                   verbose = FALSE, hessian = FALSE, as_vector = TRUE, 
                   draws = 0, constrained = TRUE, max_block_size=40, use_sigma_points=FALSE, ...) {
            stan_fit_cpp_module <- object@mk_cppmodule(object)

            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- rstan:::parse_data(rstan:::get_cppcode(object))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(rstan:::new_empty_stanfit(object)))
              }
              for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
              data <- parsed_data
            } else if (is.character(data)) { # names of objects
              data <- try(rstan:::mklist(data))
              if (is(data, "try-error")) {
                message("failed to create the data; sampling not done")
                return(invisible(rstan:::new_empty_stanfit(object)))
              }
            }
            
            if (check_data) {
              data <- try(force(data))
              if (is(data, "try-error")) {
                message("failed to evaluate the data; sampling not done")
                return(invisible(NULL))
              }
              
              if (!missing(data) && length(data) > 0) {
                data <- try(rstan:::data_preprocess(data))
                if (is(data, "try-error")) {
                  message("failed to preprocess the data; optimization not done") 
                  return(invisible(list(stanmodel = object)))
                }
              } else data <- list()
            }
            cxxfun <- rstan:::grab_cxxfun(object@dso)
            sampler <- try(new(stan_fit_cpp_module, data, cxxfun))
            if (is(sampler, "try-error")) {
              message('failed to create the optimizer; optimization not done') 
              return(invisible(list(stanmodel = object)))
            } 
            m_pars <- sampler$param_names() 
            idx_wo_lp <- which(m_pars != "lp__")
            m_pars <- m_pars[idx_wo_lp]
            p_dims <- sampler$param_dims()[idx_wo_lp]
            if (is.numeric(init)) init <- as.character(init)
            if (is.function(init)) init <- init()
            if (!is.list(init) && !is.character(init)) {
              message("wrong specification of initial values")
              return(invisible(list(stanmodel = object)))
            } 
            seed <- rstan:::check_seed(seed, warn = 1)    
            if (is.null(seed))
              return(invisible(list(stanmodel = object)))
            args <- list(init = init, 
                         seed = seed, 
                         method = "optim", 
                         algorithm = match.arg(algorithm)) 
         
            if (!is.null(sample_file) && !is.na(sample_file)) 
              args$sample_file <- rstan:::writable_sample_file(sample_file) 
            dotlist <- list(...)
            rstan:::is_arg_recognizable(names(dotlist), 
                                c("iter",
                                  "save_iterations",
                                  "refresh",
                                  "init_alpha",
                                  "tol_obj",
                                  "tol_grad",
                                  "tol_param",
                                  "tol_rel_obj",
                                  "tol_rel_grad",
                                  "history_size"),
                                 pre_msg = "passing unknown arguments: ",
                                 call. = FALSE)
            if (!is.null(dotlist$method))  dotlist$method <- NULL
            if (!verbose && is.null(dotlist$refresh)) dotlist$refresh <- 0L
            optim <- sampler$call_sampler(c(args, dotlist))
            names(optim$par) <- rstan:::flatnames(m_pars, p_dims, col_major = TRUE)
            skeleton <- rstan:::create_skeleton(m_pars, p_dims)
            if (hessian || draws) {
              fn <- function(theta) {
                sampler$log_prob(theta, FALSE, FALSE)
              }
              gr <- function(theta) {
                sampler$grad_log_prob(theta, FALSE)
              }
              theta <- rstan:::rstan_relist(optim$par, skeleton)
              theta <- sampler$unconstrain_pars(theta)
              H <- optimHess(theta, fn, gr, control = list(fnscale = -1))
              colnames(H) <- rownames(H) <- sampler$unconstrained_param_names(FALSE, FALSE)
              if (hessian) optim$hessian <- H
              if (draws > 0) {
                K <- ncol(H)
                
                Nblocks <- max(1, ceiling(K / max_block_size))
                blocksize <- ceiling(K / Nblocks)

                # create indices of data
                pts <- 1 + c(0, blocksize*(1 : (Nblocks - 1)), K)
                # normal draws
                if (use_sigma_points) {
                  if (!(draws == (2 * K))) {
                    stop("incorrect number of sigma points", c(draws, K))
                  }
                  Z <- matrix(0., K, draws)
                  for (jj in 1:K)
                    Z[jj, 2*jj - 1] <- 1 #sqrt(K)
                    Z[jj, 2*jj] <- -1 #sqrt(K)
                }
                else {
                  Z <- matrix(rnorm(K * draws), K, draws)  
                }
                
                theta_tilde <- c() #matrix(0., K, draws)
                logdet_H_inv <- 0.
                for (j in 1:Nblocks) {
                  idx <- pts[j] : (pts[j+1] - 1)
                  Qb <- -H[idx, idx]
                  R <- try(chol(Qb))
                  if (inherits(R, "try-error"))
                    R <- diag(dim(Qb)[1])
                  
                  logdet_H_inv <- logdet_H_inv - sum(log(diag(R)))
                  
                  new_theta <- backsolve(R, Z[idx, ])
                  theta_tilde <- rbind(theta_tilde, new_theta)
                }

                theta_tilde <- t(theta_tilde + theta)
                # K <- ncol(R)
                # R_inv <- backsolve(R, diag(K))
                # Z <- matrix(rnorm(K * draws), K, draws)
                # theta_tilde <- t(theta + R_inv %*% Z)
                if (constrained) {
                  theta_tilde <- t(apply(theta_tilde, 1, FUN = function(theta) {
                    sampler$constrain_pars(theta)  
                  }))
                  colnames(theta_tilde) <- names(optim$par)
                }
                else {
                  log_p <- apply(theta_tilde, 1, FUN = function(theta) {
                    sampler$log_prob(theta, adjust_transform = TRUE, gradient = FALSE)
                  })
                  log_g <- colSums(dnorm(Z, log = TRUE)) - logdet_H_inv #sum(log(diag(R_inv)))
                  optim$log_p <- log_p
                  optim$log_g <- log_g
                  colnames(theta_tilde) <- colnames(H)
                  optim$log_prob <- sampler$log_prob
                  optim$grad_log_prob <- sampler$grad_log_prob
                }
                optim$theta_tilde <- theta_tilde
              }
            }
            if (!as_vector) optim$par <- rstan:::rstan_relist(optim$par, skeleton)
            return(optim)
}) 

#environment(optimizing_new_method) <- as.environment("package:rstan")
#     fit <- optimizing(object, hessian=TRUE, draws=0, ...)

#     
#     browser()

#     return(fit)
    # if (hessian || draws) {
    #     fn <- function(theta) {
    #         sampler$log_prob(theta, FALSE, FALSE)
    #     }
    #     gr <- function(theta) {
    #         sampler$grad_log_prob(theta, FALSE)
    #     }
    #     theta <- rstan_relist(optim$par, skeleton)
    #     theta <- sampler$unconstrain_pars(theta)
    #     H <- optimHess(theta, fn, gr, control = list(fnscale = -1))
    #     colnames(H) <- rownames(H) <- sampler$unconstrained_param_names(FALSE, 
    #         FALSE)
    #     if (hessian) 
    #         optim$hessian <- H
    #     if (draws > 0 && is.matrix(R <- try(chol(-H)))) {
    #         K <- ncol(R)
    #         R_inv <- backsolve(R, diag(K))
    #         Z <- matrix(rnorm(K * draws), K, draws)
    #         theta_tilde <- t(theta + R_inv %*% Z)
    #         if (constrained) {
    #           theta_tilde <- t(apply(theta_tilde, 1, FUN = function(theta) {
    #             sampler$constrain_pars(theta)
    #           }))
    #           colnames(theta_tilde) <- names(optim$par)
    #         }
    #         else {
    #           log_p <- apply(theta_tilde, 1, FUN = function(theta) {
    #             sampler$log_prob(theta, adjust_transform = TRUE, 
    #               gradient = FALSE)
    #           })
    #           log_g <- colSums(dnorm(Z, log = TRUE)) - sum(log(diag(R_inv)))
    #           optim$log_p <- log_p
    #           optim$log_g <- log_g
    #           colnames(theta_tilde) <- colnames(H)
    #           optim$log_prob <- sampler$log_prob
    #           optim$grad_log_prob <- sampler$grad_log_prob
    #         }
    #         optim$theta_tilde <- theta_tilde
    #     }
    # }
    # if (!as_vector) 
    #     optim$par <- rstan_relist(optim$par, skeleton)
    # return(optim)
