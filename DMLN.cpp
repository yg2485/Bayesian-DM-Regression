#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}



// [[Rcpp::export]]
List DMLN(NumericMatrix X,
          NumericMatrix Y,
          int T = 10000,
          double tau = 1, // variance for the proposal of beta_0
          double tau_beta_1 = 1, // variance for the proposal of beta in add_delete
          double tau_beta_2 = 0.5, // variance for the proposal of beta by itself
          double a_p = 0.4,
          double b_p = 1.6,
          double a_beta = 2,
          double b_beta = 10,
          int d = 3
){
  set_seed(d);  
  
  // Sample statistics
  int N = Y.nrow() ;
  int P = Y.ncol() ;
  int R = X.ncol() ;
  
  // acceptance rate
  NumericMatrix acc_1(R,P); // for delta and beta
  NumericVector acc_2(P);
  
  // Initialize
  double p_delta = a_p / (a_p + b_p);
  
  NumericMatrix Delta_temp(R, P);
  for(int i = 0; i < R; i++){
      for(int j = 0; j < P; j++){
        Delta_temp(i, j) = rbinom(1, 1, p_delta)[0];
    }
  }
  
  NumericMatrix Beta_temp(R, P);
  for(int i = 0; i < R; i++){
    for(int j = 0; j < P; j++){
      if (Delta_temp(i, j) == 1)
        Beta_temp(i, j) = R::rnorm(0, tau);
    }
  }
  
  NumericVector Beta_0_temp = rnorm(P, 0, 1);
  
  NumericMatrix X_design(N, R + 1);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < R + 1; j++){
      if(j == 0){
        X_design(i, j) = 1;
      }else{
        X_design(i, j) = X(i, j-1);
      }
    }
  }
  
  NumericMatrix Beta(R + 1, P);
  for(int i = 0; i < R + 1; i++){
    for(int j = 0; j < P; j++){
      if(i == 0){
        Beta(i, j) = Beta_0_temp[j];
      }else{
        Beta(i, j) = Beta_temp(i-1, j);
      }
    }
  }
    
  NumericMatrix log_alpha(N, P), Alpha(N, P);
  for(int i = 0; i < N; i++){
    for(int j = 0; j < P; j++){
      for(int k = 0; k < R + 1; k++){
        log_alpha(i,j) += X_design(i,k) * Beta(k,j);
      }
      Alpha(i,j) = exp(log_alpha(i,j));
    }
  }
  
  // Create space to store the result
  NumericMatrix Delta_store(T*R*P, 4);
  NumericMatrix Beta_store(T*R*P, 4);
  NumericMatrix Beta_0_store(T, P);
  int count = 0;
  
  
  ///////////////////// Start the MCMC Samplying /////////////////////
  for (int t = 0; t < T; t++){
    if(t*100/T == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    
    
    //Jointly update Beta and Delta
    for (int j = 0; j < P; j++){
      for (int r = 0; r < R; r++){
        Delta_temp(r,j) = 1 - Delta_temp(r,j);
        if (Delta_temp(r,j) == 1){
          
          //add
          double beta_star = rnorm(1, 0, tau_beta_1)[0];
          NumericMatrix log_alpha_star = clone(log_alpha);
          NumericMatrix Alpha_star(N,P);
          
          for (int i = 0; i < N; i++){
            log_alpha_star(i,j) = log_alpha_star(i,j) + beta_star*X(i,r);
          }
          for(int i = 0; i < N; i++){
            for(int k = 0; k < P; k++){
              Alpha_star(i,k) = exp(log_alpha_star(i,k));
            }
          }
            
          double r_mh = 0;
          for (int l = 0; l < N; l++){
            NumericVector alpha_old = Alpha(l,_), alpha_star = Alpha_star(l,_), y=Y(l,_);
            
            r_mh = r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y));
            r_mh = r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j]);
          }
            double sum_delta = sum(Delta_temp(_,j));
            // for (int h = 0; h < R; h++){
            //   sum_delta += Delta_temp(h, j);
            // }
            double lb_star = R::lbeta(a_p + sum_delta, R-sum_delta + b_p) - R::lbeta(a_p, b_p);
            double lb_old = R::lbeta(a_p + sum_delta - 1, R - sum_delta + 1 + b_p) - R::lbeta(a_p, b_p);
            r_mh = r_mh + lb_star - log(1-exp(lb_old));
            r_mh = r_mh + R::dnorm(beta_star, 0, 10, TRUE);
            r_mh = r_mh - R::dnorm(beta_star,0, 1, TRUE);
            
            if(r_mh > log(unif_rand())){
              Beta_temp(r, j) = beta_star;
              log_alpha = clone(log_alpha_star);
              Alpha = clone(Alpha_star);
              acc_1(r,j) += 1;
            }else{
              Delta_temp(r, j) = 0;
            }
        }else{
          //delete
          double beta_old = Beta_temp(r, j);
          NumericMatrix log_alpha_star = clone(log_alpha);
          NumericMatrix Alpha_star(N,P);
          
          for (int i = 0; i < N; i++){
            log_alpha_star(i,j) = log_alpha_star(i,j) - beta_old*X(i,r);
          }
          for(int i = 0; i < N; i++){
            for(int k = 0; k < P; k++){
              Alpha_star(i,k) = exp(log_alpha_star(i,k));
            }
          }
          double r_mh = 0;
          for (int l = 0; l < N; l++){
            NumericVector alpha_old = Alpha(l,_), alpha_star = Alpha_star(l,_), y=Y(l,_);
            for (int k = 0; k < P; k++){
              alpha_old[k] = Alpha(l,k);
              alpha_star[k] = Alpha_star(l,k);
              y[k] = Y(l,k);
            }
            r_mh = r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y));
            r_mh = r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j]);
          }
          double sum_delta = sum(Delta_temp(_,j));
          // for (int h = 0; h < R; h++){
          //   sum_delta = sum_delta + Delta_temp(h, j);
          // }
          double lb_star = R::lbeta(a_p + sum_delta, R-sum_delta + b_p) - R::lbeta(a_p, b_p);
          double lb_old = R::lbeta(a_p + sum_delta + 1, R - sum_delta - 1 + b_p) - R::lbeta(a_p, b_p);
          r_mh = r_mh + log(1-exp(lb_star)) - lb_old;
          r_mh = r_mh - R::dnorm(beta_old, 0, 10, TRUE);
          r_mh = r_mh + R::dnorm(beta_old,0, 1, TRUE);
          if(r_mh > log(unif_rand())){
            Beta_temp(r, j) = 0;
            log_alpha = clone(log_alpha_star);
            Alpha = clone(Alpha_star);
            acc_1(r,j) += 1;
          }else{
            Delta_temp(r, j) = 1;
            }
          }
        }
      }
    // end of updating Beta and Delta
    
    // Update Beta
    for (int j = 0; j < P; j++){
      for (int r = 0; r < R; r++){
        if (Delta_temp(r, j) == 1){
          double beta_old = Beta_temp(r,j);
          double beta_star = R::rnorm(beta_old, tau_beta_2);
          NumericMatrix log_alpha_star = clone(log_alpha);
          NumericMatrix Alpha_star(N,P);
          
          for (int i = 0; i < N; i++){
            log_alpha_star(i,j) = log_alpha_star(i,j) - beta_old*X(i,r) + beta_star*X(i,r);
          }
          for(int i = 0; i < N; i++){
            for(int k = 0; k < P; k++){
              Alpha_star(i,k) = exp(log_alpha_star(i,k));
            }
          }
          double r_mh = 0;
          for (int l = 0; l < N; l++){
            NumericVector alpha_old = Alpha(l,_), alpha_star = Alpha_star(l,_), y=Y(l,_);

            r_mh = r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y));
            r_mh = r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j]);
          }
          r_mh = r_mh + (a_beta+0.5)*(log(b_beta + beta_old*beta_old/2) - log(b_beta + beta_star*beta_star/2));
          if(r_mh > log(unif_rand())){
            Beta_temp(r,j) = beta_star;
            log_alpha = clone(log_alpha_star);
            Alpha = clone(Alpha_star);
          }
        }
      }
    }// end of updating Beta
    
    
    // update Beta_0
    for (int j = 0; j < P; j++){
      double beta_0_old = Beta_0_temp[j];
      double beta_0_star = R::rnorm(beta_0_old, tau);
      NumericMatrix log_alpha_star = clone(log_alpha);
      NumericMatrix Alpha_star(N,P);
      
      for (int i = 0; i < N; i++){
        log_alpha_star(i,j) = log_alpha_star(i,j) - beta_0_old + beta_0_star;
      }
      for(int i = 0; i < N; i++){
        for(int k = 0; k < P; k++){
          Alpha_star(i,k) = exp(log_alpha_star(i,k));
        }
      }
      
      double r_mh = 0;
      for (int l = 0; l < N; l++){
        NumericVector alpha_old = Alpha(l,_), alpha_star = Alpha_star(l,_), y=Y(l,_);
        r_mh = r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y));
        r_mh = r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j]);
      }
      r_mh = r_mh + (beta_0_old*beta_0_old - beta_0_star*beta_0_star)/ 20;
      if(r_mh > log(unif_rand())){
        Beta_0_temp[j] = beta_0_star;
        log_alpha = clone(log_alpha_star);
        Alpha = clone(Alpha_star);
        acc_2(j) += 1;
      }
    }//end of updating Beta_0
    
    for(int r=0; r<R; r++){
      for(int p=0; p<P; p++){
        Delta_store(t*R*P + (r*P + p),_) = NumericVector::create(t+1, r+1, p+1, Delta_temp(r,p));
        Beta_store(t*R*P + (r*P + p),_) = NumericVector::create(t+1, r+1, p+1, Beta_temp(r,p));
      }
    }
    Beta_0_store(t,_) = Beta_0_temp;
    
  }// end of MCMC
  List result;
  result["Delta_store"] = Delta_store;
  result["Beta_store"] = Beta_store;
  result["Beta_0_store"] = Beta_0_store;
  result["accept_1"] = acc_1/T;
  result["accept_2"] = acc_2/T;
  // result["T"] = T;
  
  return result;
  
}// end of the function




