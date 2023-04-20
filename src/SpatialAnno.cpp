// This script considers the case with the hetergenous error with/without spatial information.
// lightly edited version of question, working fine for me
// #define ARMA_64BIT_WORD 1

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "utilSimulDRcluster.h"

//#define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);


//' @title
//' getneighborhood_fast
//' @description
//' an efficient function to find the neighborhood based on the matrix of position and a pre-defined cutoff
//'
//' @param x is a n-by-2 matrix of position.
//' @param cutoff is a threashold of Euclidean distance to decide whether a spot is an neighborhood of another spot. For example, if the Euclidean distance between spot A and B is less than cutoff, then A is taken as the neighbourhood of B. 
//' @return A sparse matrix containing the neighbourhood
//'
//' @export
// [[Rcpp::export]]
arma::sp_umat getneighborhood_fast(const arma::mat x, double cutoff)	{
	int N = x.n_rows;
    arma::sp_umat D(N, N);
    double dis;
    uvec idx, idx2;
	for (int j = 0; j < N-1; ++j)
	{    
        idx = find(abs(x(j,0) - x.col(0))<cutoff); 
        idx2 = find(idx>j);
        int p = idx2.n_elem;
		for (int i = 0; i < p; ++i)
		{
            dis = norm(x.row(idx(idx2(i))) - x.row(j), 2);
            if (dis < cutoff){
                D(idx(idx2(i)),j) = 1;
                D(j,idx(idx2(i))) = 1;
            }
		}
	}
	return D;
}
 

///////////////////////////////////////////////////////////////

void multi_det_SkCpp2(const arma::mat& X, const arma::vec& Lam_vec0,
                      const arma::mat& W0, const arma::mat& Ck, 
                      const arma::rowvec Muk, const arma::mat& Sigmak,
                      double& logdSk, arma::vec& mSk){
  //int p = X.n_cols;
  int n = X.n_rows;
  
  
  mat WC12,  tmp2;
  vec tmp1, s, tmp3;
  mat U, V, X_tk;
  
  // method B: use SVD to compute |Srk.i()|
  svd(U, s, V, Sigmak);
  WC12 = W0 * (U * diagmat(sqrt(s)));
  WC12 = sp_mat(diagmat(1.0/sqrt(Lam_vec0))) * WC12;  // change to sparse matrix multiplication.
  vec d = svd(WC12);
  logdSk = -accu(log(1 +  d%d)) - accu(log(Lam_vec0));
  // method A: directly compuate log|Srki|
  // mat Srki = W0*Sigmak*W0.t()+ sp_mat(diagmat(Lam_vec0));
  // s = eig_sym(Srki);
  // //svd(U,s, V, Srki);
  // logdSk = -accu(log(s));
  // WC12 = U* diagmat(sqrt(1.0/ s));
  // X_tk = (X - repmat(Muk* W0.t(), n, 1));  // change to sparse matrix multiplication.
  // tmp2 = X_tk * WC12;
  // mSk  = sum(tmp2 % tmp2, 1);
  svd(U, s, V, Ck.i());
  WC12 = W0 * (U * diagmat(sqrt(s)));
  WC12 = sp_mat(diagmat(1.0/sqrt(Lam_vec0))) * WC12;  
  X_tk = (X - repmat(Muk* W0.t(), n, 1)) * sp_mat(diagmat(1/sqrt(Lam_vec0))) ;  // change to sparse matrix multiplication.
  tmp1 = sum(X_tk % X_tk, 1);
  tmp2 = X_tk * WC12;
  tmp3 = sum(tmp2 % tmp2, 1);
  mSk = tmp1 - tmp3;
  
}

void multi_det_SkCpp(const arma::mat& X, const arma::vec& Lam_vec0, const arma::mat& W0, const arma::mat& Ck, 
                     const arma::rowvec Muk, 
                     double& logdSk, arma::vec& mSk){
  //int p = X.n_cols;
  int n = X.n_rows;
  // int p = X.n_cols;
  // // mSk = zeros(n);
  // S2k = zeros(p);
  // dSk = 0;
  
  mat WC12,  tmp2;
  vec tmp1, s, tmp3;
  mat U, V, X_tk;
  
  svd(U, s, V, Ck.i());
  
  WC12 = W0 * (U * diagmat(sqrt(s)));
  WC12 = sp_mat(diagmat(1/sqrt(Lam_vec0))) * WC12; // sparse matrix multiplication to speed up
  vec d = svd(WC12);
  //dSk = arma::as_scalar(prod(1- d % d)) / prod(Lam_vec0);
  logdSk = accu(log(1 - d%d)) - accu(log(Lam_vec0));
  X_tk = (X - repmat(Muk* W0.t(), n, 1)) * sp_mat(diagmat(1/sqrt(Lam_vec0))) ;
  tmp1 = sum(X_tk % X_tk, 1);
  tmp2 = X_tk * WC12;
  tmp3 = sum(tmp2 % tmp2, 1);
  mSk = tmp1 - tmp3;
}


vec decomp(const arma::mat& Cki, const arma::mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s)));
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  

// [[Rcpp::export]]  
arma::mat update_W_l(const arma::mat& x_l, const arma::mat& y_l, const arma::mat& Ez_l, const arma::mat& Ci_l){
  int n = y_l.n_rows;

  mat Ezzt = Ez_l.t() * Ez_l + n * Ci_l;
  mat Ezzt_i = inv_sympd(Ezzt);
  return x_l.t() * Ez_l * Ezzt_i;
}


// [[Rcpp::export]]  
arma::mat update_W_u(const arma::mat& x_u, const arma::mat& R_u, const arma::cube& Ez_u, const arma::mat& C_ui,  const arma::vec& N_u){
  int k, K= R_u.n_cols, q= Ez_u.n_cols, n= R_u.n_rows;
  mat tmpMat(n,q, fill::zeros), Ezzt(q,q, fill::zeros), tmpMat2;
  // vec N = sum(R.t(), 1);
  for(k=0; k<K; ++k){
    tmpMat2= repmat(R_u.col(k), 1, q) % Ez_u.slice(k);
    tmpMat += tmpMat2;
    Ezzt += tmpMat2.t() * Ez_u.slice(k) + N_u(k) * C_ui;
  }
  return x_u.t() * tmpMat * inv_sympd(Ezzt);
}




//' 
//'  
//' 
//' 
//' 
//' 
//' 
// [[Rcpp::export]]  
arma::mat update_W0(const arma::mat& X, const arma::mat& R, const arma::cube& Ez, const arma::cube& Ci,  const arma::vec& N){
  int k, K= R.n_cols, q= Ez.n_cols, n= R.n_rows;
  mat tmpMat(n,q, fill::zeros), Ezzt(q,q, fill::zeros), tmpMat2;
  // vec N = sum(R.t(), 1);
  for(k=0; k<K; ++k){
    tmpMat2= repmat(R.col(k), 1, q) % Ez.slice(k);
    tmpMat += tmpMat2;
    Ezzt+= tmpMat2.t() * Ez.slice(k) + N(k) * Ci.slice(k);
  }
  return X.t() * tmpMat * Ezzt.i();
}


 
// initallize sigma
// [[Rcpp::export]]
arma::mat init_Sgm(const arma::mat& y_l, const arma::mat& Ez_l, const arma::mat& Ci_l, const arma::mat& Mu_l,  
                    const bool& diagSigma){
  int k, K= y_l.n_cols, q= Mu_l.n_cols, n= y_l.n_rows;
  mat sigma(q, q, fill::zeros);
  for(k = 0; k<K; ++k){
    mat res = Ez_l - repmat(Mu_l.row(k), n, 1);
    for (int i = 0; i < n; ++i)
    {
      sigma += trans(res.row(i)) * y_l(i, k) * res.row(i);
    }
  }
  mat Sgm = sigma/n + Ci_l;
  if(diagSigma){
      Sgm = diagmat(Sgm);
  } 
  return(Sgm);
}



 
// [[Rcpp::export]]  
arma::mat update_Sgm_u(const arma::mat& R_u, const arma::cube& Ez_u, const arma::mat& C_ui, const arma::mat& Mu_u,  
                   const arma::vec&  N_u, const bool& diagSigma){
  int k, K= R_u.n_cols, q= Mu_u.n_cols, n= R_u.n_rows;
  mat sigma(q, q, fill::zeros);
  for(k = 0; k<K; ++k){
    mat res = Ez_u.slice(k) - repmat(Mu_u.row(k), n, 1);
    for (int i = 0; i < n; ++i)
    {
      sigma += trans(res.row(i)) *R_u(i, k) * res.row(i);
    }
    
  }
  mat Sgm = sigma/accu(N_u) + C_ui; 
  if(diagSigma){
    Sgm = diagmat(Sgm);
  } 
  return(Sgm);
}

  

// update Lambda in labled data
// [[Rcpp::export]]
arma::vec update_Lam_l(const arma::mat& x_l, const arma::mat& W_l, const arma::mat& Ez_l, const arma::mat & C_li,const bool& homo){
  int p = x_l.n_cols;
  vec Lam;
  mat tmpXk = (x_l - Ez_l * W_l.t() );
  mat tmpXk_2 = tmpXk % tmpXk;
  rowvec Lsum1 = sum(tmpXk_2)/(x_l.n_rows*1.0);
  vec Lsum2 = decomp(C_li, W_l);
  Lam = Lsum1.t() + Lsum2;
  if (homo){
    Lam = mean(Lam) * ones(p,1);
  }  
  
  return Lam; 
}

// update Lambda
// [[Rcpp::export]]
arma::vec update_Lam_u(const arma::mat& R_u, const arma::mat& X_u, const arma::mat& W_u, const arma::cube& Ez_u, const arma::mat& C_ui,const bool& homo){
  int k, K= R_u.n_cols, p = X_u.n_cols;
  rowvec N_u = sum(R_u);
  vec Lsum(p, fill::zeros), Lam;
  mat tmpXk;
  vec Lsum1 = decomp(C_ui, W_u);
  for(k=0; k<K; ++k){
    tmpXk = (X_u - Ez_u.slice(k) * W_u.t() );
    Lsum += trans(sum(tmpXk % tmpXk % repmat(R_u.col(k), 1, p)));
    Lsum += N_u(k) * Lsum1;
  }
  if(homo){
    Lam = mean(Lsum)* ones(p,1) / (X_u.n_rows*1.0);
  }else{
    Lam = Lsum/(X_u.n_rows*1.0);
  }
  return Lam; 
}



// update Lambda
// [[Rcpp::export]]
arma::vec update_Lam(const arma::mat& R, const arma::mat& X, const arma::mat& W, const arma::cube& Ez, const arma::cube& Ci,const bool& homo){
  int k, K= R.n_cols, p = X.n_cols;
  rowvec N = sum(R);
  vec Lsum(p,fill::zeros), Lam;
  mat tmpXk;
  for(k=0; k<K; ++k){
    tmpXk = (X - Ez.slice(k) * W.t() );
    Lsum += trans(sum(tmpXk % tmpXk % repmat(R.col(k), 1, p)));
    Lsum += N(k) * decomp(Ci.slice(k), W);
  }
  if(homo){
    Lam = mean(Lsum)* ones(p,1) / (X.n_rows*1.0);
  }else{
    Lam = Lsum/(X.n_rows*1.0);
  }
  return Lam; 
}



// // update sigma20: homo variance
// vec update_Lam2(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci){
//   int k, K= R.n_cols, p = X.n_cols;
//   mat tmpXk;
//   double term1=0, term2=0;
//   for(k=0; k <K; ++k){
//     tmpXk = (X - Ez.slice(k) * W.t() );
//     term1 += accu(sum(tmpXk % tmpXk, 1) % R.col(k));
//     term2  += accu(decomp(Ci.slice(k), W))* accu(R.col(k));
//   }
//   double sigma20 = (term1+ term2)/ X.n_elem;
//   return(ones(p)* sigma20);
// }

//Calculate Q function
double Q_fun(const arma::mat& X, const arma::mat& R,  const arma::cube& Ez, const arma::cube& Ci, 
             const arma::mat& W0, const arma::mat& Mu0, const arma::cube& Sigma0, const arma::vec& Pi0, 
             const arma::vec& Lam_vec0){
  
  double Q = 0, tmp_scalar =0;
  int  K = Pi0.n_elem, n = X.n_rows;
  int q = Mu0.n_cols;
  mat tmpSig(q, q, fill::zeros), tmpMat;
  colvec Ezik;
  mat Ezzt(q,q, fill::zeros);
  for(int k=0; k<K; k++){
    tmpSig = Sigma0.slice(k);
    tmpMat = Ez.slice(k);
    for(int i = 0; i<n; i++){
      Ezik =  trans(tmpMat.row(i)); 
      Ezzt = Ci.slice(k) + Ezik * Ezik.t();
      
      tmp_scalar =  0.5* accu(log(Lam_vec0)) + 0.5* accu( (X.row(i) % X.row(i)) / Lam_vec0.t()) +
        0.5 * arma::as_scalar(trace(W0.t()* diagmat(1.0/Lam_vec0)*W0* Ezzt)- 2* X.row(i)* diagmat(1.0/Lam_vec0)*W0*Ezik);
      Q +=  - R(i,k) * tmp_scalar + R(i,k)*(log(Pi0(k)) -  0.5* log(det(tmpSig))- 0.5* trace(tmpSig.i()*
        Ezzt)+ arma::as_scalar(Mu0.row(k) * tmpSig.i()*(Ezik- 0.5*trans(Mu0.row(k)))));
    }
  }
  return Q;
}

///////////////////////////////////////////////////////////////


//[[Rcpp::export]]
arma::mat getPairDist(const arma::mat x)	{
  int N = x.n_rows;
  arma::mat D(N, N);
  for (int j = 0; j < N; ++j)
  {
    for (int i = j; i < N; ++i)
    {
      D(i, j)	= norm(x.row(i) - x.row(j), 2);
      D(j, i) = D(i, j);
    }
  }
  
  return D;
}



arma::sp_mat get_spNbs(arma::ivec y, const arma::sp_mat& Adj) {   // ivec是索引型向量
  // row is for pixel.
  //output a sparse matrix, i-th row contains labels of neighbor_i. 
  // Make const iterator
  arma::sp_mat::const_iterator start = Adj.begin(); //构造一个sp_mat的常数迭代器,常数迭代器只可读，不可写，实现对矩阵按照列对每个非零元素进行访问。
  //arma::sp_mat::const_iterator end   = Adj.end();
  
  // Calculate number of nonzero points
  //int n = std::distance(start, end);
  int n = Adj.n_nonzero; // 计算Adj的所有非零元的个数
  //cout << "n=" << n << endl;
  //cout << "n=" << Adj.n_nonzero << endl;
  
  sp_mat spNbs(y.n_elem, y.n_elem);    // neiborhood state matrix, matched with Adj.
  
  
  arma::sp_mat::const_iterator it = start; // Note spNbs is not a symmetric matrix, the nonzero in i-th row is the class label of sample i.
  for(int i = 0; i < n; ++i)
  {
    //temp(0) = it.row();
    //temp(1) = it.col();
    spNbs(it.row(), it.col()) = y(it.col()); // it只自加非零元个数次，得到每个i对应的邻居的状态
    ++it; // increment
  }
  
  return spNbs.t(); // return the class label of neighbor matrix, i-th column is the neighbor label of sample i
}

// [[Rcpp::export]]
arma::mat calYenergy2D_sp(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  // Calculate the value of energy function of y, which is equal to negative logliklihood up to a constant
  int n = y.n_rows;
  arma::sp_mat spNbs_t = get_spNbs(y, Adj); // transform spNbs to iterate by column.
  arma::mat Uy(n, K);
  double n_sameS;
  int i, k, nn;
  for (k = 0; k < K; k++)
  {
    for (i = 0; i < n; i++)
    {
      arma::sp_mat col(spNbs_t.col(i)); // the class label of neighbors of i-th sample.
      n_sameS = 0;
      
      nn = col.n_nonzero; // the number of neighbors of i-th sample
      for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
        n_sameS += ((*j) == (k+1));
        
      }
      Uy(i, k) = alpha(k) + beta * (nn - n_sameS)/2;
      
      
    }
  }
  
  arma::mat C_mat = normalise(exp(-Uy), 1, 1); // pseudo likelihood of Y.
  Uy = -log(C_mat); // normalized Uy, this is the energy of y.
  return Uy;
  
}
// Suppose we know true y, can we estimate beta correctly by maximizing pseudo loglikelihood  
// //[[Rcpp::export]] 
// double obj_beta2(const arma::ivec& y, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
//   
//   mat Uy = calYenergy2D_sp(y, Adj, K, alpha, beta);
//   arma::mat C_mat = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE to get the likelihood
//   double loglike = 0;
//   int n = y.n_elem;
//   int i = 0;
//   for(; i < n; ++i)
//   {
//     
//     loglike += log(C_mat(i, y(i)-1));
//   }
//   return loglike;
// }

//[[Rcpp::export]]  
double obj_beta(const arma::ivec& y, const arma::mat& R, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
  
  mat Uy = calYenergy2D_sp(y, Adj, K, alpha, beta); // Uy was normalized, so there is no need to normalized Uy. 
  //arma::mat C_mat = normalise(exp(-Uy), 1, 1); // set all rowSums to be ONE to get the likelihood
  //return accu(R % log(C_mat)); 
  return -accu(R % Uy);
}

           
            

////////////////////////////////////////////////////////////////////////////    

// update variance of marker genes
// [[Rcpp::export]]
arma::vec update_sigma(arma::mat X_m, arma::mat rho, arma::mat R, arma::vec alpha, arma::mat bet){
  int k, K = R.n_cols, n = X_m.n_rows, m = X_m.n_cols;
  vec V(m, fill::zeros), sigma;
  mat tmpXk;
  rowvec rho_bet;
  mat A = repmat(trans(alpha), n, 1);
  for(k=0; k<K; ++k){
    rho_bet = trans(rho.col(k) % bet.col(k));
    tmpXk = (X_m - A - repmat(rho_bet, n, 1));
    V += trans(sum(tmpXk % tmpXk % repmat(R.col(k), 1, m)));
  }
  sigma = V/n;   
  return sigma; 
}

// [[Rcpp::export]]
arma::vec update_Rsigma(arma::mat X_m, arma::mat rho, arma::mat R, arma::vec alpha, arma::mat bet, arma::cube eta){
  int k, K = R.n_cols, n = X_m.n_rows, m = X_m.n_cols;
  vec V(m, fill::zeros), sigma;
  mat tmpXk;
  rowvec rho_bet;
  mat A = repmat(trans(alpha), n, 1);
  for(k=0; k<K; ++k){
    rho_bet = trans(rho.col(k) % bet.col(k));
    tmpXk = (X_m - A - repmat(rho_bet, n, 1) % eta.slice(k));
    V += trans(sum(tmpXk % tmpXk % repmat(R.col(k), 1, m)));
  }
  sigma = (V+1e-6)/n;   
  return sigma; 
}



// update base mean of marker genes
// [[Rcpp::export]]
arma::vec update_alpha(arma::mat X_m, arma::mat rho, arma::mat R, arma::mat bet){
  int k, K = R.n_cols, n = X_m.n_rows, m = X_m.n_cols;
  vec a(m, fill::zeros);
  mat tmpXk;
  rowvec rho_bet;
  
  for(k=0; k<K; ++k){
    rho_bet = trans(rho.col(k) % bet.col(k));
    tmpXk = (X_m - repmat(rho_bet, n, 1));
    a += trans(sum(tmpXk % repmat(R.col(k), 1, m)));
  }
  
  return a/n; 
}

// [[Rcpp::export]]
arma::vec update_Ralpha(arma::mat X_m, arma::mat rho, arma::mat R, arma::mat bet, arma::cube eta){
  int k, K = R.n_cols, n = X_m.n_rows, m = X_m.n_cols;
  vec a(m, fill::zeros);
  mat tmpXk;
  rowvec rho_bet;
  
  for(k=0; k<K; ++k){
    rho_bet = trans(rho.col(k) % bet.col(k));
    tmpXk = (X_m - repmat(rho_bet, n, 1) % eta.slice(k));
    a += trans(sum(tmpXk % repmat(R.col(k), 1, m)));
  }
  
  return a/n; 
}

 
// update specific mean of marker genes
// [[Rcpp::export]]
arma::rowvec update_beta_j(arma::mat X_m, arma::mat rho, arma::mat R, arma::vec alpha, arma::vec N, int j, double lfc){
  int k, K = R.n_cols;

  rowvec beta_j(K, fill::zeros);
  for (k=0; k<K; ++k){
      if(rho(j, k) != 0) {
        beta_j(k) =  sum(R.col(k) % (X_m.col(j) - alpha(j))) / N(k);
        if(beta_j(k) < lfc)  beta_j(k) = lfc;
      }
  }
  
  return beta_j; 
}



// update specific mean of marker genes
// [[Rcpp::export]]
arma::rowvec update_Rbeta_j(arma::mat X_m, arma::mat rho, arma::mat R, arma::vec alpha, arma::cube eta, int j, double lfc){
  int k, K = R.n_cols;

  rowvec beta_j(K, fill::zeros);

  for (k=0; k<K; ++k){
    double N_eta = sum(R.col(k) % eta.slice(k).col(j));
    if((rho(j, k) * N_eta) != 0) { 
      beta_j(k) =  sum(R.col(k) % (X_m.col(j) - alpha(j)) % eta.slice(k).col(j)) / N_eta;
      if(beta_j(k) < lfc) beta_j(k) = lfc;
    }
  }
  
  return beta_j; 
}


// [[Rcpp::export]]
void Mstep_marker(const arma::mat& X_m, const arma::mat& rho, const arma::mat& R, 
                arma::vec& alpha, arma::mat& bet, arma::vec& sigma, const arma::vec& N, const double& lfc){
    int j, m = X_m.n_cols;
    int K = rho.n_cols;
    for (j = 0; j < m; ++j) {
        for (int it = 0; it < 20; it++) {
          mat mu_old = repmat(alpha, 1, K) + bet;
          alpha = update_alpha(X_m, rho, R, bet);
          bet.row(j) = update_beta_j(X_m, rho, R, alpha, N, j, lfc);
          mat mu = repmat(alpha, 1, K) + bet;
          if (accu(abs(mu - mu_old)) < 1e-4) {
            //cout << it << endl;
            break;}
        }
    }

    sigma = update_sigma(X_m, rho, R, alpha, bet);
}


// update specific mean of marker genes
// [[Rcpp::export]]
arma::vec update_etajk(arma::mat X_m, arma::vec alpha, arma::mat bet, int j, int k){
  int i, n = X_m.n_rows;

  vec eta_jk(n, fill::zeros);

  for (i=0; i<n; ++i){
    eta_jk(i) = 1 * (X_m(i, j) > (alpha(j) + 0.5*bet(j, k)));
  }
  
  return eta_jk; 
}


// [[Rcpp::export]]
arma::cube add_mu(const arma::mat& rho, const arma::vec& alpha,  const arma::mat& bet, const arma::cube& eta) {
    int n = eta.n_rows;
    int m = bet.n_rows;
    int K = bet.n_cols;
    cube mu(eta);
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < m; ++j)
      {
        for (int k = 0; k < K; ++k)
        {
          mu(i, j, k) = alpha(j) + rho(j, k) * bet(j, k) * eta(i, j, k);
        }
      }
    }
    return mu;
}

// [[Rcpp::export]]
void Mstep_Rmarker(const arma::mat& X_m, const arma::mat& rho, const double& lfc, const arma::mat& R, 
                arma::vec& alpha, arma::mat& bet, arma::cube& eta, arma::cube& mu, arma::vec& sigma){
    int j, m = X_m.n_cols;
    int k, K = rho.n_cols;
    cube mu_old(mu); 
    for (j = 0; j < m; ++j) {
        for (int it = 0; it < 20; it++) {
          mu_old = add_mu(rho, alpha, bet, eta); 
          alpha = update_Ralpha(X_m, rho, R, bet, eta);
          bet.row(j) = update_Rbeta_j(X_m, rho, R, alpha, eta, j, lfc);
          for (k = 0; k < K; ++k) {
            if(rho(j, k) != 0) eta.slice(k).col(j) = update_etajk(X_m, alpha, bet, j, k);
          }
          cube mu = add_mu(rho, alpha, bet, eta); 
          if (accu(abs(mu - mu_old)) < 1e-4) {
            //cout << it << endl;
            break;}
        }
    }
    mu = add_mu(rho, alpha, bet, eta); 
    sigma = update_Rsigma(X_m, rho, R, alpha, bet, eta);
}



// [[Rcpp::export]]
arma::vec log_cp_marker(arma::mat X_m, arma::vec mu_k_m, arma::vec sigma){
  int n = X_m.n_rows;
  int m = X_m.n_cols;

  vec mSk(n, fill::zeros);
  for (int i = 0; i < n; i++){
      for(int j = 0; j < m; j++)  {
        mSk(i) += (X_m(i, j) - mu_k_m(j)) * (X_m(i, j) - mu_k_m(j)) / sigma(j);
      }
  }
  return(mSk);
}

// for mathod which is robust to markers
// [[Rcpp::export]]
arma::vec log_cp_Rmarker(arma::mat X_m, arma::mat mu_k_m, arma::vec sigma){
  int n = X_m.n_rows;
  int m = X_m.n_cols;

  vec mSk(n, fill::zeros);
  for (int i = 0; i < n; i++){
      for(int j = 0; j < m; j++)  {
        mSk(i) += (X_m(i, j) - mu_k_m(i, j)) * (X_m(i, j) - mu_k_m(i, j)) / sigma(j);
      }
  }
  return(mSk);
}
 
// clustering with marker genes
// [[Rcpp::export]]
Rcpp::List em_marker(const arma::mat& X_m, const arma::mat& rho,  const double& lfc,       // input
                     const int& maxIter, const double& eps, const bool& verbose){
    // basic info
    int n = X_m.n_rows;
    int m = X_m.n_cols;
    int K = rho.n_cols;

    
    mat R(n, K, fill::randu);  
    vec Rsum = sum(R, 1);
    R = R / repmat(Rsum, 1, K);

    vec N = arma::sum(R.t(), 1);
    vec Pi = N/ accu(N);
    mat A(R);

    // valid initial values
    mat bet(m, K, fill::zeros);
    vec alpha = update_alpha(X_m, rho, R, bet);
    vec sigma(K, fill::ones);    
    Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
    mat mu = repmat(alpha, 1, K) + bet;

    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec diff(maxIter);
    diff(0) = INT_MAX;

    vec maxA(n, fill::zeros);
    vec loglik_more_vec = maxA;

    List Mstep;
    int k, iter;
    //mat mu_old = mu;
    //mat M(n, K, fill::randu);
    //ucolvec type = index_max(M, 1);
    for (iter = 1; iter < maxIter; ++iter)
    {

      //ucolvec type_old = type;
     
      // E-step for marker genes
      double logdSigma = accu(log(sigma));
      for (k = 0; k < K; k++) {
        vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
        A.col(k) = log(Pi(k)) - 0.5 * logdSigma - 0.5 * mSk;
      }
      //cout << "checkpoint" << endl;
      maxA = max(A, 1);
      A = (A - repmat(maxA, 1, K));
      loglik_more_vec = sum(exp(A), 1);
      loglik(iter) = sum(log(loglik_more_vec) + maxA) - n * m /2.0 * log(2* M_PI); 
      R = exp(A) / repmat(loglik_more_vec, 1, K);
 
      //type = index_max(R, 1);
      // M-step for marker genes
      N = arma::sum(R.t(), 1);
      Pi = N/ accu(N); 
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;


      //cout << "mu" <<  mu(0,0) << endl;
      //cout << "mu_old" <<  mu_old(0,0) << endl;
      //mat D = (mu - mu_old);
      //diff(iter) = accu(D % D);
 
      //mu_old = mu;
      if(abs(loglik(iter)  - loglik(iter-1))   < -1e-7){
         perror("The likelihood failed to increase!");
         break;
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
      

    }
    List resList = List::create(
      Rcpp::Named("iter") = iter,
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1),
      //Rcpp::Named("diff") = diff.subvec(0, iter-1),
      Rcpp::Named("type") = index_max(R, 1) + 1,
      Rcpp::Named("R") = R,
      Rcpp::Named("Pi") = Pi,
      Rcpp::Named("alpha") = alpha,
      Rcpp::Named("bet") = bet,
      Rcpp::Named("sigma") = sigma,
      Rcpp::Named("mu") = mu
   );
   return(resList);
} 



// clustering with marker genes
// [[Rcpp::export]]
Rcpp::List em_Rmarker(const arma::mat& X_m, const arma::mat& rho, const double& lfc,       // input
                     const int& maxIter, const double& eps, const bool& verbose){
    // basic info
    int n = X_m.n_rows;
    int m = X_m.n_cols;
    int K = rho.n_cols;

    
    mat R(n, K, fill::randu);  
    vec Rsum = sum(R, 1);
    R = R / repmat(Rsum, 1, K);

    vec N = arma::sum(R.t(), 1);
    vec Pi = N/ accu(N);
    mat A(R);

    // valid initial values
    mat bet(m, K, fill::zeros);
    vec alpha = update_alpha(X_m, rho, R, bet);
    cube eta(n, m, K, fill::zeros);
    vec sigma(K, fill::ones);    
    Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
    cube mu = add_mu(rho, alpha, bet, eta); 

    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec diff(maxIter);
    diff(0) = INT_MAX;

    vec maxA(n, fill::zeros);
    vec loglik_more_vec = maxA;

    List Mstep;
    int k, iter;
    //mat mu_old = mu;
    //mat M(n, K, fill::randu);
    //ucolvec type = index_max(M, 1);
    for (iter = 1; iter < maxIter; ++iter)
    {
      // E-step for marker genes
      double logdSigma = accu(log(sigma));
      for (k = 0; k < K; k++) {
        vec mSk = log_cp_Rmarker(X_m, mu.slice(k), sigma);
        A.col(k) = log(Pi(k)) - 0.5 * logdSigma - 0.5 * mSk;
      }
      //cout << "checkpoint" << endl;
      maxA = max(A, 1);
      A = (A - repmat(maxA, 1, K));
      loglik_more_vec = sum(exp(A), 1);
      loglik(iter) = sum(log(loglik_more_vec) + maxA) - n * m /2.0 * log(2* M_PI); 
      R = exp(A) / repmat(loglik_more_vec, 1, K);
 
      // M-step for marker genes
      Mstep_Rmarker(X_m, rho, lfc, R, alpha, bet, eta, mu, sigma);
      

      //cout << "mu" <<  mu(0,0) << endl;
      if(abs(loglik(iter)  - loglik(iter-1))   < -1e-7){
         perror("The likelihood failed to increase!");
         break;
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
      
    }
    N = arma::sum(R.t(), 1);
    Pi = N/ accu(N); 

    List resList = List::create(
      Rcpp::Named("iter") = iter,
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1),
      //Rcpp::Named("diff") = diff.subvec(0, iter-1),
      Rcpp::Named("type") = index_max(R, 1) + 1,
      Rcpp::Named("R") = R,
      Rcpp::Named("Pi") = Pi,
      Rcpp::Named("alpha") = alpha,
      Rcpp::Named("bet") = bet,
      Rcpp::Named("eta") = eta,
      Rcpp::Named("sigma") = sigma,
      Rcpp::Named("mu") = mu
   );
   return(resList);
} 

 
  
 
// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List em(const arma::mat& X_m, const arma::mat& X_u, const arma::mat& rho, const double& lfc,
    const arma::mat& R_int, const arma::vec& Pi_int, const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int m = X_m.n_cols;
    int K = rho.n_cols;
    int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    mat R(R_int);  
    vec N = arma::sum(R.t(), 1);
    vec Pi = N/ accu(N);
    mat A(R);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    double  logdS_u;
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    cube Ez_u(n, q, K, fill::zeros);
    
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      // compute loglikelihood
      for (k=0; k<K; k++){
        vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
        multi_det_SkCpp2(X_u, Lam_u, W_u, C_u, Mu_u.row(k), Sgm_u, logdS_u, mS_u);
        // cout<<"dSk="<<exp(logdSk_u)<<"mSk=" <<mS_u(0)<<endl;
        A.col(k) = log(Pi(k)) - 0.5 * logdSigma - 0.5 * mSk + 0.5*logdS_u  -0.5 * mS_u;
      
        Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
                        repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
      } 


      maxA = max(A, 1);
      A = (A - repmat(maxA, 1, K));
      loglik_vec = sum(exp(A),1);
      loglik(iter) = sum(log(loglik_vec) + maxA) - n * (p + m) /2.0 * log2pi; 
      R = exp(A) / repmat(loglik_vec, 1, K);
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Pi = N/ accu(N); 
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;

      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 

      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("type") = index_max(R, 1) + 1,
      Rcpp::Named("Pi") = Pi,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ez_u,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 


// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List emR(const arma::mat& X_m, const arma::mat& X_u, const arma::mat& rho, const double& lfc,
    const arma::mat& R_int, const arma::vec& Pi_int, const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int m = X_m.n_cols;
    int K = rho.n_cols;
    int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    mat R(R_int);  
    vec N = arma::sum(R.t(), 1);
    vec Pi = N/ accu(N);
    mat A(R);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int); 
    cube eta(n, m, K, fill::zeros);   
    cube mu = add_mu(rho, alpha, bet, eta); 
    //mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    double  logdS_u;
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    cube Ez_u(n, q, K, fill::zeros);
    
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      // compute loglikelihood
      for (k=0; k<K; k++){
        vec mSk = log_cp_Rmarker(X_m, mu.slice(k), sigma);
        multi_det_SkCpp2(X_u, Lam_u, W_u, C_u, Mu_u.row(k), Sgm_u, logdS_u, mS_u);
        // cout<<"dSk="<<exp(logdSk_u)<<"mSk=" <<mS_u(0)<<endl;
        A.col(k) = log(Pi(k)) - 0.5 * logdSigma - 0.5 * mSk + 0.5*logdS_u  -0.5 * mS_u;
      
        Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
                        repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
      } 


      maxA = max(A, 1);
      A = (A - repmat(maxA, 1, K));
      loglik_vec = sum(exp(A),1);
      loglik(iter) = sum(log(loglik_vec) + maxA) - n * (p + m) /2.0 * log2pi; 
      R = exp(A) / repmat(loglik_vec, 1, K);
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Pi = N/ accu(N); 
      Mstep_Rmarker(X_m, rho, lfc, R, alpha, bet, eta, mu, sigma);

      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 

      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("type") = index_max(R, 1) + 1,
      Rcpp::Named("Pi") = Pi,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("eta_m") = eta,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ez_u,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 


// [[Rcpp::export]]
Rcpp::List runICM_sp (const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::mat& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {

    vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 

// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmem(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    //int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      // compute loglikelihood
      ICM_fit = runICM_sp(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;

      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 

      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 



// [[Rcpp::export]]
Rcpp::List runICM_sp_anno (const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::mat& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {

    vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
      /*
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
      */
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  /*
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  */
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 


// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmem_anno(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    //int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      // compute loglikelihood
      ICM_fit = runICM_sp_anno(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;

        
      /*
      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 
      */
        
        
        
      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 






// [[Rcpp::export]]
Rcpp::List runICM_sp_factor(const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::mat& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {

    vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  /*
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  */
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 



// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmem_factor(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    //int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      
      // compute loglikelihood
      ICM_fit = runICM_sp_factor(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
       
        
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      /*
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;
      */
        
      
      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 
      
        
        
        
      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 




// [[Rcpp::export]]
Rcpp::List runICM_sp_spatialfactor(const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::mat& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {

    vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 



// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmem_spatialfactor(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    //int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      
      // compute loglikelihood
      ICM_fit = runICM_sp_spatialfactor(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
       
        
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      /*
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;
      */
        
      
      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 
      
        
        
        
      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 






// [[Rcpp::export]]
Rcpp::List runICM_sp_spatialanno (const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::mat& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {

    vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
      /*
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
      */
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 


// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmem_spatialanno(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    //int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      // compute loglikelihood
      ICM_fit = runICM_sp_spatialanno(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;

        
      /*
      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 
      */
        
        
        
      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 





// [[Rcpp::export]]
Rcpp::List runICM_sp_annofactor (const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::mat& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {

    vec mSk = log_cp_marker(X_m, mu.col(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
    /*
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  */
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 

// without considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmem_annofactor(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    //int p = X_u.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    mat mu = repmat(alpha, 1, K) + bet;

    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      // compute loglikelihood
      ICM_fit = runICM_sp_annofactor(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Mstep_marker(X_m, rho, R, alpha, bet, sigma, N, lfc);
      mu = repmat(alpha, 1, K) + bet;

      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 

      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 






// [[Rcpp::export]]
Rcpp::List runICM_Rsp (const arma::mat& X_m, const arma::mat& X_u,  arma::ivec& y, const arma::sp_mat& Adj, const arma::vec& Pi, 
                      const arma::cube& mu, const arma::vec& sigma, double logdSigma,
                      const arma::mat& W_u, const arma::vec& Lam_u, const arma::mat& Mu_u,
                      const arma::mat& Sgm_u, const arma::mat& Sgm_ui, const arma::mat& C_u, const arma::mat& C_ui,
                      const arma::vec& beta_grid, double beta, int maxIter_ICM) {
  // Target: estimate Y, evaluate R, Ez, and update beta by using grid search.
  
  // basic info.
  int n = X_u.n_rows, m = X_m.n_cols, p=X_u.n_cols, K = Mu_u.n_rows, q= Mu_u.n_cols;
  int iter, k;
  
  // two cached objects used for parameters update.
  cube Ez_u(n,q,K, fill::zeros);
  double  logdS_u;
  vec mS_u(n);
  
  // evaluate energy of x, Ux
  arma::mat Ux(n, K);
  for (k = 0; k < K; k++) {
    vec mSk = log_cp_Rmarker(X_m, mu.slice(k), sigma);
    multi_det_SkCpp2(X_u, Lam_u,W_u, C_u, Mu_u.row(k), Sgm_u, // Use SVD to speed up.
                    logdS_u, mS_u);
    // cout<<"dSk="<<exp(logdSk)<<"mSk=" <<mSk(0)<<endl;
    Ux.col(k) =  0.5 * logdSigma + 0.5 * mSk - 0.5*logdS_u  + 0.5 * mS_u; // calculate energy by column.
    
    Ez_u.slice(k) = (X_u * (repmat(1.0/ Lam_u, 1, q) % W_u) + 
      repmat(Mu_u.row(k) * Sgm_ui, n, 1)) * C_ui;
  }
  
  // Estimate Y by ICM
  arma::vec Energy(maxIter_ICM);
  Energy(0) = INFINITY;
  arma::mat Uy(n, K);
  arma::mat U(n, K);
  //--------------------------------------------------------------------------------  
  // ICM algrithm to estimate Y
  //--------------------------------------------------------------------------------
  // int Iteration = 1;
  for (iter = 1; iter < maxIter_ICM; iter ++ ) {
    
    Uy = calYenergy2D_sp(y, Adj, K, Pi, beta);
    
    U = Uy + Ux; // log likelihood of (x, y).
    arma::vec Umin = min(U, 1);
    arma::uvec y_u = index_min(U, 1);
    y = conv_to< ivec >::from(y_u) + 1;
    //y = index_min(U, 1);

    Energy(iter) = sum(Umin);
    if (Energy(iter) - Energy(iter - 1) > 1e-5) {
      //cout << "diff Energy = " << Energy(iter) - Energy(iter - 1)  << endl;
      break;
    }
    
    if (Energy(iter-1) - Energy(iter) < 1e-5)
    {
      //cout << "ICM Converged at Iteration = " << iter  << endl;
      break;
    }
  }
  
  // if (iter == maxIter_ICM) {
  //   Iteration = iter - 1;
  // } else {
  //   Iteration = iter;
  // }
  
  // calculate R and pseudo observed loglikelihood
  vec maxA1 = max(-U, 1);
  U = (-U - repmat(maxA1, 1, K));
  vec loglik_more_vec = sum(exp(U),1);
  double loglik = sum(log(loglik_more_vec) + maxA1) - n * (p + m) /2.0 * log(2* M_PI); 
  arma::mat R = exp(U) / repmat(loglik_more_vec, 1, K);
  
  
  // vec energy = Energy.subvec(1, Iteration);
  
  // update beta: grid search.
  int ng_beta = beta_grid.n_elem;
  vec objBetaVec(ng_beta);
  for(k=0; k < ng_beta; ++k){
    objBetaVec(k) = obj_beta(y, R, Adj, K, Pi, beta_grid(k));
  }
  //cout <<  objBetaVec  << endl;
  List output = List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("R") = R,
    Rcpp::Named("Ez_u") = Ez_u,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("hbeta") = beta_grid(index_max(objBetaVec))); // get the maximum beta.
  
  return output; 
  
}  
 

// considering the spatial information
// [[Rcpp::export]]
Rcpp:: List icmemR(const arma::mat& X_m, const arma::mat& X_u, 
    const arma::sp_mat& Adj, const arma::mat& rho, const double& lfc,
    const arma::ivec& y_int, const arma::vec& Pi, const double& xi_int, const arma::vec& xi_grid,
    const arma::vec& alpha_int, const arma::mat& bet_int, const arma::mat& sigma_int,
    const arma::mat& Mu_u_int, const arma::mat& W_u_int, const arma::mat& Sgm_u_int, const  arma::vec& Lam_u_int,
    const int& maxIter, const int& maxIter_ICM, const double& eps, const bool& verbose, 
    const bool& homo = false, const bool& diagSigmak = false){
    // basic info
    int n = X_m.n_rows;
    int K = rho.n_cols;
    int m = X_m.n_cols;
    int q = Mu_u_int.n_cols;


    // Initialize the iterative parameters of marker genes
    ivec y(y_int);
    vec N(K, fill::zeros);
    //mat A(R_int);
    double xi(xi_int);

    // valid initial values
    mat bet(bet_int);
    vec alpha(alpha_int);
    vec sigma(sigma_int);    
    //mat mu = repmat(alpha, 1, K) + bet;
    cube eta(n, m, K, fill::zeros); 
    cube mu = add_mu(rho, alpha, bet, eta); 
    // Initialize the iterative parameters of non-marker genes
    mat Mu_u(Mu_u_int), W_u(W_u_int);
    vec Lam_u(Lam_u_int);
    mat Sgm_u(Sgm_u_int);
   
    // If p is sufficient large, loglik can not be computed.
    // But this can be solved by some programming tricks.
    vec loglik(maxIter);
    loglik(0) = INT_MIN;
    vec maxA(n,fill::zeros);
    vec loglik_vec = maxA;
  
    // Define the variables that will be used in algorithm
    // variables usded in updating Pi0
    mat C_u(q, q, fill::zeros), C_ui(q,q, fill::zeros);
    vec mS_u(n);
    int k, iter;
  
    // variables usded in updating Mu0
    //cube Ez_u(n, q, K, fill::zeros);
    //mat R(R_int);  
    // cout<<"start EM algorithm in CPP::"<<endl;
    // begin algorithm
    List ICM_fit;
    for (iter = 1; iter < maxIter; iter++){
      
      // maker genes
      double logdSigma = accu(log(sigma));
      // non-marker genes
      mat Sgm_ui = inv_sympd(Sgm_u);
      C_u = W_u.t() * sp_mat(diagmat(1.0/ Lam_u)) * W_u +  Sgm_ui;
      C_ui = inv_sympd(C_u);  

      
      // compute loglikelihood
      ICM_fit = runICM_Rsp(X_m, X_u, y, Adj, Pi, mu, sigma, logdSigma, W_u, Lam_u, Mu_u, Sgm_u, Sgm_ui,
                          C_u, C_ui, xi_grid, xi, maxIter_ICM);
      loglik(iter) = ICM_fit["loglik"];
      mat R = ICM_fit["R"];
      cube Ez_u = ICM_fit["Ez_u"];
      xi = ICM_fit["hbeta"];
      // cout<<"Finish R computing in CPP::"<<endl;


      // M-step
      N = arma::sum(R.t(), 1);
      Mstep_Rmarker(X_m, rho, lfc, R, alpha, bet, eta, mu, sigma);

      // update Mu_u 
      for (k = 0; k<K; ++k){
        Mu_u.row(k) = trans(R.col(k)) * Ez_u.slice(k) / N(k);
      }
  
      // update Sgm_u
      Sgm_u = update_Sgm_u(R, Ez_u, C_ui, Mu_u, N, diagSigmak);
  
      // update W_u
      W_u = update_W_u(X_u, R, Ez_u, C_ui, N);
      // update  Lambda
      Lam_u = update_Lam_u(R, X_u, W_u, Ez_u, C_ui, homo); 

      /////////////////////////////////////////////////////////////////////
      // calculate loglikelihood
      // output return value
      // calculate loglikelihood
      if(loglik(iter)  - loglik(iter-1)   < -1e-7){
         //perror("The likelihood failed to increase!");
         break;
      }
    
      // output algorithm info.
      if(verbose){
         // cout<<"iter = "<< iter +1 <<", loglik="<<loglik(iter)<<", dobj ="<<(loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1))<<endl;
         // cout<<"iter = "<< iter+1<<", Qval="<<Q<<", dQ ="<<Q  - tmp_Q <<endl;
         Rprintf("iter = %d, loglik= %4f, dloglik=%4f \n", 
               iter +1, loglik(iter), (loglik(iter)  - loglik(iter-1))/ abs(loglik(iter-1)));
      }
      if(abs((loglik(iter)  - loglik(iter-1))/ loglik(iter-1)) < eps) break;
   }
   
  
  mat R = ICM_fit["R"]; // R's existence is temporary, so we require to redefine it.
  cube Ez_u = ICM_fit["Ez_u"]; // Ez is also the problem.
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    
    Ezz +=  Ez_u.slice(k) % repmat(R.col(k), 1, q);
  }
  

   List resList = List::create(
      Rcpp::Named("R") = R,
      Rcpp::Named("xi") = xi,
      Rcpp::Named("type") = y,
      Rcpp::Named("alpha_m") = alpha,
      Rcpp::Named("bet_m") = bet,
      Rcpp::Named("eta_m") = eta,
      Rcpp::Named("mu_m") = mu,
      Rcpp::Named("sigma_m") = sigma,
      Rcpp::Named("Ez_u") = Ezz,
      Rcpp::Named("Mu_u") = Mu_u,
      Rcpp::Named("Sgm_u") = Sgm_u,
      Rcpp::Named("W_u") = W_u,
      Rcpp::Named("Lam_u") = Lam_u,
      //Rcpp::Named("loglik") = loglik(iter-1),
      //Rcpp::Named("dLogLik") = loglik(iter-1)  - loglik(iter-2),
      Rcpp::Named("loglik") = loglik.subvec(0, iter-1)
   );
   return(resList);
} 

