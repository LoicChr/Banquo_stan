// Translation of Courrieu 2018 Fast Computation of Moore-Penrose Inverse Matrices arXiv:0804.4809 [cs]
// Written for square matrices currently
matrix geninv(matrix G){
    int dim[2];
    int m;
    int n;
    matrix[m,n] A;
    vector[m] dA;
    real tol;
    matrix[m,n] L;
    int r;
    matrix[m,r] LL;
    matrix[r,r] M;
    matrix[n,m] Y;

    dim = dims(G);
    m = dim[1];
    n = dim[2];

    A = (G') * G;

    # Full rank Cholesky factorization of A 
    dA = diagonal(A);
    tol = max(dA);
    for( i in 1:m ){
        if(dA[i] > 0 && dA[i] < tol){
            tol = dA[i];
        }
    }
    tol *= 1E-9;
    
    // I don't know if this is necessary
    for( i in 1:m ){
        for( j in 1:n ){
            L[i,j] = 0;
        }
    }

    r = 0;
    for( k in 1:n ){
        r += 1;
        L[k:n,r] = A[k:n,k] - L[k:n,1:(r-1)] * L[k,1:(r-1)]';
    
        # Note: for r=1, the substracted vector is zero
        if(L[k,r] > tol){
            L[k,r] = sqrt(L[k,r]);
            if (k < n){
                L[(k+1):n,r] = L[(k+1):n,r] / L[k,r];
            }
        }else{
             r = r - 1;
        }
    }

    LL = L[:,1:r];

    # Computation of the generalized inverse of G
    M = inverse(LL' * LL);
    Y = LL * M * M * LL' * G';

    return Y;
}