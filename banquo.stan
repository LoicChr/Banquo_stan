functions{
    matrix geninv(matrix G, int n, real epsilon){
      vector[n] dA;
      matrix[n,n] A;
      matrix[n,n] L;
      matrix[n,n] LL;
      matrix[n,n] M;
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
      
      // go column by column to get the Cholesky factorization
      L = rep_matrix(0, n, n);
      LL = rep_matrix(0, n, n);
      r = 0;
      for( k in 1:n ){
          r += 1;
          // the first column is just a normalized version of the first column of A; nothing fancy required
          if(r == 1){
              L[k:n,r] = A[k:n,k];
          }
          // stan complains about having L on the LHS and RHS; I believe we could decouple these somehow
          else{
              L[k:n,r] = A[k:n,k] - block(LL,k,1,n-k+1,r-1) * sub_row(LL,k,1,r-1)'; //L[k:n,1:(r-1)] * L[k,1:(r-1)]';
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

          // to avoid using L on lhs and rhs above  
          LL = L;

          // DEBUG as a potential speed up, we could probably break the loop once we've reached the rank of the matrix
      }
  
      // LL is the reduced-rank form of L with some empty columns on the end
      LL = rep_matrix(0,n,n);
      LL[:,1:r] = L[:,1:r];

      // use the top left of the LLinv matrix to hold the inverse of LL' LL
      M = rep_matrix(0,n,n);
      M[1:r,1:r] = inverse(block(LL,1,1,n,r)' * block(LL,1,1,n,r));

      // final computation of the generalized inverse of G
      Y = block(LL,1,1,n,r) * block(M,1,1,r,r) * block(M,1,1,r,r) * block(LL,1,1,n,r)' * G';
  
      return Y;
    }

    // convert a vector of carrying capacities (Ks) to equilibrium abundances given an alpha matrix
    row_vector traitspace_to_banquo(row_vector K, matrix alpha, int S){
        row_vector[S] species;
        matrix[S,S] alpha_inv;
        row_vector[S] N;
        // row_vector[S] K_tmp;
        // matrix[S,S] alpha_tmp;
        real minN;

        // we start by assuming all species could coexist
        species = to_row_vector(rep_vector(1,S));
        
        // solve for the equilibrium abundances
        alpha_inv = geninv(alpha, S, 10e-8); 
        N = (alpha_inv * K')';

        // remove 'negative' species until all species have postive abundance
        minN = min(N .* species);
        while(minN < 0) { //Tol here?
          // "remove" the species with the lowest abundance
          for(i in 1:S){
            if(N[i] == minN){
              species[i] = 0;
              break;
            }
          }
        // // attempt to do the above in a more clever (and possibly faster) way
        // while(min(N .* species) < 0){  
          // // "remove" the species with the lowest abundance
          // species[sort_indices_asc(N)[1]] = 0;

          // recalculate the predicted abundances with the "reduced" community
          alpha_inv = geninv(quad_form_diag(alpha, species), S, 10e-8);
          N = (alpha_inv * (K .* species)')';

          minN = min(N .* species);
        }

        return N;
    }

    // compute the hellinger distance between two vectors
    real hellinger(vector obs, vector pred){
      real hell;
      hell = dot_product(sqrt(obs)-sqrt(pred), sqrt(obs)-sqrt(pred));
      hell = sqrt(hell / 2.0);
      return hell;
    }
    
}
data{
    // input data passed to stan
    int<lower=1> N;
    int<lower=1> S;
    int<lower=1> Nsites;
    int<lower=1> T;
    int sites[N];
    matrix[Nsites,S] observed;
    matrix[N,S] traitspace;
    matrix[S,T] traits; // Uncomment this
}
parameters{
    // control the intraspecific values
    real logmean_alphaii;
    real<lower=0> sigma_alphaii;
    real<lower=0> aii[S];
    
    // control the interspecific values
    real<lower=0> alphaij_intercept;
    real alphaij_center;
    real<lower=0> alphaij_width;
}
transformed parameters{
    matrix[S,S] sampled_alpha;
    matrix[N,S] banquo;
    matrix[Nsites, S] banquo_agg;

    // construct the alpha matrix given the sampled parameters
    for( i in 1:S ){
      for( j in i:S ){
        if( i == j ){
          // intraspecific alphas
          sampled_alpha[i,i] = aii[i];
        }else{
          // // for testing purposes only          
          // sampled_alpha[i,j] = 0;
          // sampled_alpha[j,i] = 0;

          // interspecific alpha based on distance from i to j
          sampled_alpha[i,j] = alphaij_intercept * exp(-pow(traits[i,1] - traits[j,1] - alphaij_center, 2)/alphaij_width);

          // interspecific alpha based on distance from j to i
          sampled_alpha[j,i] = alphaij_intercept * exp(-pow(traits[j,1] - traits[i,1] - alphaij_center, 2)/alphaij_width);
        }
      }
    }

    // calculate banquo predictions across all observations
    banquo = (geninv(sampled_alpha, S, 10e-8) * traitspace')';

    // for aggregation of banquo output within sites
    banquo_agg = rep_matrix(0, Nsites, S);
    
    // "correct" observations for which banquo predicts negative abundances (competitive exclusion)
    for( i in 1:N ){
      if (min(banquo[i]) < 0){
        banquo[i] = traitspace_to_banquo(traitspace[i], sampled_alpha, S);
      }
      banquo_agg[sites[i]] += banquo[i];
    }
    
    // normalize site-specific banquo results so that they represent relative abundances
    for (i in 1:Nsites){
      banquo_agg[i] = banquo_agg[i]/sum(banquo_agg[i]);
    }
}
model{
    // priors on core inferred parameters
    // DEBUG: we should consider regularizing the priors to get better behavior
    logmean_alphaii ~ normal(0, 1);
    sigma_alphaii ~ cauchy(0, 2);
    alphaij_intercept ~ cauchy(0, 2);
    alphaij_center ~ normal(0, 10);
    alphaij_width ~ cauchy(0, 2);

    // intraspecific interaction coefficients as lognormal "random effects"
    // aii ~ lognormal(logmean_alphaii, sigma_alphaii);
    for( i in 1:S ){
        aii[i] ~ lognormal(logmean_alphaii, sigma_alphaii);
    }

    // sample predicted values using banquo_agg computed in transformed_parameters block
    for( i in 1:Nsites ){
      // observed[i] ~ normal(banquo_agg[i], sigma_obs);
      target += 1 - hellinger(to_vector(observed[i]), to_vector(banquo_agg[i]));
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
