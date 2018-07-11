geninv = function(G){
  # Returns the Moore-Penrose inverse of the argument 
  # Transpose if m < n
  # Translation of Courrieu 2018 Fast Computation of Moore-Penrose Inverse Matrices arXiv:0804.4809 [cs]
  
  m <- nrow(G)
  n <- ncol(G)
  transpose = FALSE
  if (nrow(G) < ncol(G)){
    transpose = T
    A = G%*%t(G)
    n = m
  }else{
    A = t(G)%*%G
  }

  # Full rank Cholesky factorization of A 
  dA <- diag(A); 
  tol <- min(dA[dA>0])*1e-9; 
  L <- matrix(0, nrow(A), ncol(A))
  r <- 0
  for (k in 1:n){
    r <- r+1
    L[k:n,r] <- A[k:n,k, drop = FALSE] - L[k:n,1:(r-1), drop = FALSE] %*% t(L[k,1:(r-1), drop = FALSE])
    # Note: for r=1, the substracted vector is zero
    if (L[k,r] > tol){
      L[k,r] <- sqrt(L[k,r])
      if (k < n){
        L[(k+1):n,r] <- L[(k+1):n,r, drop = F]/L[k,r]
      }
    }else{
      r = r - 1
    }
  }
  L <- L[,1:r]
  # Computation of the generalized inverse of G
  M = solve( t(L)%*%L )
  if (transpose) {
    Y <- t(G) %*% L %*% M %*% M %*% t(L)
  }else{
    Y <- L %*% M %*% M %*% t(L)%*% t(G)
  }
  return(Y)
}