


set.seed(12345)
n <- 30
p <- 2
k = 3
X <- matrix(rnorm(n*p),n,p)
y <- sample(c(0,1,2),nrow(X),replace=TRUE)

beta <- matrix(c(1,2,3,4,3,2),p*k,1)
beta0 <- matrix(0,k*p,1) # initial starting vector
lambda <- 5

softmax <- function(x) {
  exp_x <- exp(x - max(x))
  return(exp_x / sum(exp_x))
}


f1<-function(beta) {
  
  arrays = cbind(X%*%beta[(1+0*p):(p+0*p)],X%*%beta[(1+1*p):(p+1*p)],X%*%beta[(1+2*p):(p+2*p)])
  return(-sum(log(t(apply(arrays,1,softmax)))*cbind(y==0,y==1,y==2))/n)
}

gradf1<-function(beta){
  X_T = t(X)
  arrays = cbind(X%*%beta[(1+0*p):(p+0*p)],X%*%beta[(1+1*p):(p+1*p)],X%*%beta[(1+2*p):(p+2*p)])
  p_arrays = (1-t(apply(arrays,1,softmax)))*cbind(y==0,y==1,y==2)
  # p_arrays = apply(p_arrays,1,sum)
  t= rbind(X_T%*%p_arrays,X_T%*%p_arrays,X_T%*%p_arrays)/n
  return(-t)
}

library(nnet)

f <- function(beta) {
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_one_hot <- class.ind(y)  # 需要nnet包
  loss <- -sum(y_one_hot * log(softmax_scores)) / n
  return(loss)
}

gradf <- function(beta) {
  beta_matrix <- matrix(beta, nrow = p, ncol = k)
  scores <- X %*% beta_matrix
  softmax_scores <- t(apply(scores, 1, softmax))
  y_one_hot <- class.ind(y)  # 需要nnet包
  grad_matrix <- t(X) %*% (softmax_scores - y_one_hot) / n
  return(matrix(grad_matrix, nrow = p * k, ncol = 1))
}

g <- function(beta) { lambda*norm(as.matrix(beta),'1') }
proxg <- function(beta, tau) { sign(beta)*(sapply(abs(beta) - tau*lambda,
  FUN=function(x) {max(x,0)})) }
x0 <- double(p*k) # initial starting iterate
tau1 <- 10
# 
# 
# set.seed(12345)
# n <- 100
# p <- 25
# X <- matrix(rnorm(n*p),n,p)
# beta <- matrix(rnorm(p),p,1)
# y <- X%*%beta + rnorm(n)
# beta0 <- matrix(0,p,1) # initial starting vector
# lambda <- 10
# 
# f <- function(beta){ 0.5*norm(X%*%beta - y, "F")^2 }
# gradf <- function(beta){ t(X)%*%(X%*%beta - y) }
# g <- function(beta) { lambda*norm(as.matrix(beta),'1') }
# proxg <- function(beta, tau) { sign(beta)*(sapply(abs(beta) - tau*lambda,
#   FUN=function(x) {max(x,0)})) }
# x0 <- double(p) # initial starting iterate
# tau1 <- 10




fasta <- function(f,gradf,g,proxg,x0,tau1,max_iters=1e2,w=10,backtrack=TRUE,recordIterates=FALSE,stepsizeShrink=0.5,eps_n=1e-15) {
  ## Allocate memory
  residual <- double(max_iters)          #  Residuals
  normalizedResid <- double(max_iters)   #  Normalized residuals
  taus       <- double(max_iters)        #  Stepsizes
  fVals      <- double(max_iters)        #  The value of 'f', the smooth objective term
  objective  <- double(max_iters+1)      #  The value of the objective function (f+g)
  totalBacktracks <- 0                   #  How many times was backtracking activated?
  backtrackCount  <- 0                   #  Backtracks on this iterations
  
  ## Intialize array values
  x1       <- x0
  d1       <- x1
  f1       <- f(d1)
  fVals[1] <- f1
  gradf1   <- gradf(d1)
  
  if (recordIterates) {
    iterates <- matrix(0,length(x0),max_iters+1)
    iterates[,1] <- x1
  } else {
    iterates <- NULL
  }
  
  ##  The handle non-monotonicity
  maxResidual       <- -Inf   #  Stores the maximum value of the residual that has been seen. Used to evaluate stopping conditions.
  minObjectiveValue <-  Inf   #  Stores the best objective value that has been seen.  Used to return best iterate, rather than last iterate
  
  objective[1] <- f1 + g(x0)
  
  ## Begin Loop
  for (i in 1:max_iters) {
    print(i)
    ##  Rename iterates relative to loop index.  "0" denotes index i, and "1" denotes index i+1
    x0 <- x1                       # x_i <--- x_{i+1}
    gradf0 <- matrix(gradf1)       # gradf0 is now $\nabla f(x_i)$
    tau0 <- tau1                   # \tau_i <--- \tau_{i+1}
    
    ##  FBS step: obtain x_{i+1} from x_i
    x1hat <- x0 - tau0 * c(gradf0) # Define \hat x_{i+1}
    x1 <- proxg(x1hat, tau0)       # Define x_{i+1}
    
    ##  Non-monotone backtracking line search
    Dx <- matrix(x1 - x0)
    d1 <- x1
    f1 <- f(d1)
    if (backtrack) {
      M <- max( fVals[max(i-w,1):max(i-1,1)] )  # Get largest of last 10 values of 'f'
      backtrackCount <- 0
      prop <- (f1 - 1e-12 > M + t(Dx)%*%gradf0 + 0.5*(norm(Dx,'f')**2)/tau0) && (backtrackCount < 20)
      #      print(paste("Pre AND: ", prop))
      #      print(paste("Sufficient Decrease: ",f1 - 1e-12 > M + t(Dx)%*%gradf0 + 0.5*(norm(Dx,'f')**2)/tau0) )
      #      print(paste("f1: ", f1, "M: ", M, "<Dx,gradf0>: ", t(Dx)%*%gradf0, "tau0: ", tau0))
      #      print(paste("backtrack: ",backtrackCount < 20))
      #  Note: 1e-12 is to quench rounding errors
      while ( prop ) {# The backtracking loop
        tau0 <- tau0*stepsizeShrink      # shrink stepsize
        x1hat <- x0 - tau0*c(gradf0)     # redo the FBS
        x1 <- proxg(x1hat, tau0)
        d1 <- x1
        f1 <- f(d1)
        Dx <- matrix(x1 - x0)
        backtrackCount <- backtrackCount + 1
        #        print(paste("Sufficient Decrease: ",f1 - 1e-12 > M + t(Dx)%*%gradf0 + 0.5*(norm(Dx,'f')**2)/tau0) )
        #        print(paste("backtrack: ",backtrackCount < 20))
        prop <- (f1 - 1e-12 > M + t(Dx)%*%gradf0 + 0.5*(norm(Dx,'f')**2)/tau0) && (backtrackCount < 20)
        #        print(paste("AND: ", prop))
        #        print(paste("backtrackCount: ",backtrackCount))
      }
      totalBacktracks <- totalBacktracks + backtrackCount
    }
    
    ## Record information
    taus[i] <- tau0  # stepsize
    residual[i]  <-  norm(Dx,'f')/tau0  # Estimate of the gradient, should be zero at solution
    maxResidual <- max(maxResidual, residual[i])
    normalizer <- max(norm(gradf0,'f'),norm(as.matrix(x1 - x1hat),'f')/tau0) + eps_n
    normalizedResid[i] <- residual[i]/normalizer   # Normalized residual:  size of discrepancy between the two derivative terms, divided by the size of the terms
    fVals[i] <- f1
    #  record function values
    objective[i+1] <- f1 + g(x1)
    newObjectiveValue <- objective[i+1]
    
    if (recordIterates) {  #  record iterate values
      iterates[,i+1] <- x1
    }
    
    if (newObjectiveValue < minObjectiveValue) { # Methods is non-monotone:  Make sure to record best solutiom
      bestObjectiveIterate <- x1
      minObjectiveValue <- min(minObjectiveValue, newObjectiveValue)
    }
    
    ## Compute stepsize needed for next iteration using BB/spectral method
    gradf1 <- gradf(d1)
    Dg <- matrix(gradf1 + (x1hat - x0)/tau0) # Delta_g, note that Delta_x was recorded above during backtracking
    dotprod <- t(Dx)%*%Dg
    #    print(paste("dotprod: ", t(Dx)%*%Dg))
    tau_s <- norm(Dx,'f')^2 / dotprod   #  First BB stepsize rule
    tau_m <- dotprod / norm(Dg,'f')^2   #  Alternate BB stepsize rule
    tau_m <- max(tau_m,0)
    if (abs(dotprod) < 1e-15) break
    if (2*tau_m > tau_s) {   #  Use "Adaptive"  combination of tau_s and tau_m
      tau1 <- tau_m
    } else {
      tau1 <- tau_s - .5*tau_m  #  Experiment with this param
    }
    if ( (tau1 <= 0) || is.infinite(tau1) || is.nan(tau1) ) {
      tau1 <- tau0*1.5
    }
  }
  
  if (recordIterates) {  #  record iterate values
    iterates <- iterates[,1:(i+1),drop=FALSE]
  }
  
  return(list(x=bestObjectiveIterate,objective=objective[1:(i+1)],fVals=fVals[1:i],
              totalBacktracks=totalBacktracks,
              residual=residual[1:i],taus=taus[1:i],iterates=iterates))
}


sol <- fasta(f,gradf,g,proxg,x0,tau1)
