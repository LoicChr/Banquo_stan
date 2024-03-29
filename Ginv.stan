// Translation of Courrieu 2018 Fast Computation of Moore-Penrose Inverse Matrices arXiv:0804.4809 [cs]
// Written for square matrices currently
matrix geninv(matrix G, int n, real epsilon){
    vector[n] dA;
    matrix[n,n] A;
    matrix[n,n] L;
    matrix[n,n] LL;
    matrix[n,n] Y;
    int r;
    real tol;

    // classic Nimzo-Indian start to everything
    A = (G') * G;

    // Full rank Cholesky factorization of A 
    dA = diagonal(A);
    tol = max(dA);
    for( i in 1:n ){
        if(dA[i] > 0 && dA[i] < tol){
            tol = dA[i];
        }
    }
    tol *= epsilon;
    
    // I don't know if this is necessary but didn't want L to not start with all zeroes
    // DEBUG check what happens if we comment it out
    for( i in 1:n ){
        for( j in 1:n ){
            L[i,j] = 0;
        }
    }

    // go column by column to get the Cholesky factorization
    r = 0;
    for( k in 1:n ){
        r += 1;
        // the first column is just a normalized version of the first column of A; nothing fancy required
        if(r == 1){
            L[k:n,r] = A[k:n,k];
        }
        // stan complains about having L on the LHS and RHS; I believe we could decouple these somehow
        else{
            L[k:n,r] = A[k:n,k] - block(L,k,1,n-k+1,r-1) * sub_row(L,k,1,r-1)'; //L[k:n,1:(r-1)] * L[k,1:(r-1)]';
        }

        // if the diagonal is larger than our tolerance we will keep it and must normalize
        if(L[k,r] > tol){
            L[k,r] = sqrt(L[k,r]);
            if (k < n){
                for( j in (k+1):n ){
                    L[j,r] /= L[k,r];
                }
            }
        }
        // otherwise the matrix is not full rank and hence we will not use these columns in the transformation later
        else{
             r = r - 1;
        }

        // DEBUG as a potential speed up, we could probably break the loop once we've reached the rank of the matrix
    }

    // LL is the reduced-rank form of L
    // we again add zeros to the final few columns just in case LL is not initialized as all zeroes
    // DEBUG as a potential speed up, we can check that this is not necessary
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
