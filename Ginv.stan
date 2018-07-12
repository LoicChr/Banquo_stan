// Translation of Courrieu 2018 Fast Computation of Moore-Penrose Inverse Matrices arXiv:0804.4809 [cs]
// Written for square matrices currently
matrix geninv(matrix G, int n){
    vector[n] dA;
    matrix[n,n] A;
    matrix[n,n] L;
    matrix[n,n] LL;
    matrix[n,n] Y;
    real tol;
    int r;

    // classic starting move
    A = (G') * G;

    // Full rank Cholesky factorization of A 
    dA = diagonal(A);
    tol = max(dA);
    for( i in 1:n ){
        if(dA[i] > 0 && dA[i] < tol){
            tol = dA[i];
        }
    }
    tol *= 1E-9;
    
    // I don't know if this is necessary
    for( i in 1:n ){
        for( j in 1:n ){
            L[i,j] = 0;
        }
    }

    r = 0;
    for( k in 1:n ){
        r += 1;
        if(r == 1){
            L[k:n,r] = A[k:n,k];
        }
        // stan complains about having L on the LHS and RHS
        else{
            L[k:n,r] = A[k:n,k] - block(L,k,1,n-k+1,r-1) * sub_row(L,k,1,r-1)'; //L[k:n,1:(r-1)] * L[k,1:(r-1)]';
        }

        // Note: for r=1, the substracted vector is zero
        if(L[k,r] > tol){
            L[k,r] = sqrt(L[k,r]);
            if (k < n){
                for( j in (k+1):n ){
                    L[j,r] /= L[k,r];
                }
            }
        }else{
             r = r - 1;
        }
    }

    // LL is the reduced form and we add zeros to the final few columns just in case
    LL[:,1:r] = L[:,1:r];
    for( i in 1:n ){
        for( j in (r+1):n ){
            LL[i,j] = 0;
        }
    }

    // Computation of the generalized inverse of G
    // We regretably need to do the two separate inverses since I don't know how to dynamically set a matrix size
    Y = block(LL,1,1,n,r) * inverse(block(LL,1,1,n,r)' * block(LL,1,1,n,r)) * inverse(block(LL,1,1,n,r)' * block(LL,1,1,n,r)) * block(LL,1,1,n,r)' * G';

    return Y;
}
