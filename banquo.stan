functions{
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
    // if we want to define any functions we can do that here; this example might not work though
    row_vector traitspace_to_banquo(row_vector K, matrix sampled_alpha, int S){
        row_vector[S] N;
        row_vector[S] species;
        matrix[S,S] alpha_inv;
        matrix[S,S] sampled_alpha2;
        row_vector[S] K2;
        real minVal;
        alpha_inv = geninv(sampled_alpha, S, 10e-8); 
        N = (alpha_inv * K')';
        species = to_row_vector(rep_vector(1,S));
        minVal = min(N .* species);
        
        while (minVal < 0) { //Tol here?
         for (i in 1:S){
           if (N[i] == minVal){ 
             species[i] *= 0;
          }
         }
        sampled_alpha2 = quad_form_diag(sampled_alpha, species);
        alpha_inv = geninv(sampled_alpha2, S, 10e-8);
        K2 = K .* species;
        N = (alpha_inv * K2')';
        minVal = min(N .* species);
        
       }
        return N;
    }
    real hellinger(vector obs, vector pred){
      real hell;
      hell = dot_product(sqrt(obs)-sqrt(pred), sqrt(obs)-sqrt(pred));
      hell = 1-1*sqrt(hell)/sqrt(2);
      return hell;
    }
    
}
data{
    // input data passed to stan
    int<lower=1> N;
    int<lower=1> S;
    int<lower =1> Nsites;
    int sites[N];
    // int<lower=1> T;
    matrix[Nsites,S] observed;
    matrix[N,S] traitspace;
    matrix[S,1] traits; // Uncomment this
}
parameters{
    // these are the core parameters being inferred
    // noise around the predicted relative abundances
    real<lower=0> sigma_obs;

    // controls the intraspecific values
    real logmean_alphaii;
    real<lower=0> sigma_alphaii;
    real<lower=0> aii;
    
    // important output
  //  matrix[Nsites, S] banquo_agg;
  //  real<lower=0> alphaii[S];
    
    // controls the interspecific values
    // //real mean_alphaij;
   real<lower=0> alphaij_intercept;
     real alphaij_center;
     real<lower=0> alphaij_width;

    // reluctantly required parameters since they're used in multiple sections of the code

 }
 transformed parameters{
     //  int site_id;
      matrix[Nsites, S] banquo_agg;
      matrix[S,S] sampled_alpha;
          matrix[N,S] banquo;
          for( i in 1:S ){
        for( j in i:S ){
            if( i == j ){
                sampled_alpha[i,i] = aii;// pass as these are sampled by alphaii variable above
             }else{
                  // sampled_alpha[i,j] = alphaij_intercept*alphaij_width*sqrt(2*pi());
                   sampled_alpha[i,j] = (traits[i,1] - traits[j,1] - alphaij_center)/alphaij_width ;
                  sampled_alpha[i,j] = alphaij_intercept*exp(-sampled_alpha[i,j]*sampled_alpha[i,j]);
               //   sampled_alpha[i,j] *= exp(normal_lpdf(traits[i,1] - traits[j,1] | alphaij_center, alphaij_width));
                     //        sampled_alpha[j,i] = alphaij_intercept*sqrt(2*pi()); 
                   sampled_alpha[j,i] = (traits[j,1] - traits[i,1] - alphaij_center)/alphaij_width ;
                  sampled_alpha[j,i] = alphaij_intercept*exp(-sampled_alpha[j,i]*sampled_alpha[j,i]);
          //  sampled_alpha[j,i] =0;
          //  sampled_alpha[i,j] =0;
            }
        }
    }
    
     for( i in 1:N ){
      banquo[i] = traitspace_to_banquo(traitspace[i], sampled_alpha, S);
    }
        // // Aggregation of banquo output
    banquo_agg = rep_matrix(rep_row_vector(0,S),Nsites);
    for (i in 1:N){
     // site_id = ;
      banquo_agg[sites[i]] += banquo[i];
    }

    
    //Normalisation
    for (i in 1:Nsites){
      banquo_agg[i] = banquo_agg[i]/sum(banquo_agg[i]);
    }
//     // put here items that get calculated along the way that you would like to save and examine later as part of the posteriors
//     matrix[S,S] alpha;
// 
//     // perform transformations etc to create "effective" and inverse alpha matrices
//     for( i in 1:S ){
//         for( j in i:S ){
//             // intraspecific (diagonal) values
//             if( i == j ){
//                 alpha[i,j] = alphaii[i];
//             }
//             // interspecific (offdiagonal) values
//             else{
//                 // alpha[i,j] = 0;j
//                 // alpha[j,i] = 0;
//                 alpha[i,j] = sampled_alpha[i,j];
//                 alpha[j,i] = sampled_alpha[j,i];
//             }
//         }
//     }
// 
 }
model{
    // convenience parameters that get used below but that we don't want to save as part of the posteriors
    matrix[S,S] alpha_inv;





    // priors on core inferred parameters
    sigma_obs ~ cauchy(0, 2);
    logmean_alphaii ~ normal(0, 10);
    sigma_alphaii ~ cauchy(0, 2);
  //  mean_alphaij ~ normal(0, 10);
  //  sigma_alphaij ~ cauchy(0, 2);
    alphaij_center ~ normal(0, 10);
    alphaij_width ~ cauchy(0, 2);
    alphaij_intercept ~ cauchy(0, 2);
    
    aii ~ lognormal(logmean_alphaii, sigma_alphaii);
    // sample intraspecific effects as a random effect
    // I use a lognormal to keep them as positive deviations from 1; there might be a better way to approach this

    // generate the "sampled" interspecific interaction matrix


    
    // sample predicted values

    for( i in 1:Nsites ){
       observed[i] ~ normal(banquo_agg[i], sigma_obs);
       // target+= hellinger(to_vector(observed[i]), to_vector(banquo_agg[i]));
    }
}
// generated quantities{
//    // this section is useful for some metrics like WAIC which require tracking the likelihood of each observation separately
//     vector[N] log_lik;
//     real dev;
    
//     // model predictions
//     for ( i in 1:N ) {
//         // the linear component
//         blah blah blah

//         // calculate each observation's likelihood in order to use WAIC from rethinking
//         log_lik[i] = normal_lpdf( relfit[i] | mu[i], sigma );
//     }

//     // the deviance is -2 * the sum of all loglikelihood contributions across observations    
//     dev = 0;
//     dev = dev + (-2)*sum(log_lik);
// }
