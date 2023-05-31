
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;







// [[Rcpp::export(getH_ik_Wald_RCPP)]] //getH_ik_Wald
arma::vec H_ik_Wald(arma::mat& x, int i, int k, int p, arma::vec a)
{
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec ai = a.rows( span((i-1)*(d*p+1)+1 -1 , (i-1)*(d*p+1)+ d*p +1 -1) ); //intercept - index

  arma::mat X = x.rows( span(k-p-1,k-1-1)) ;
  //arma::mat XX = x.rows( span(k-p,k-1));
  //XX.attr("dim") = Dimension(1, d*p);
  arma::vec V = vectorise(flipud(X),1).t();
  arma::vec O = ones(1);
  arma::vec VV = join_cols(O,V);
  //V.rows(span(1  ,d*p +1) ) = XX;
  //Numericarma::vector Xx = V.push_front(1); //add intercept
  //arma::vec xr = x.row(k);
  //for(int ii=1; ii< (p+1); i++){
  // xr = x.row(k-ii);
  // V.rows( span((ii-1)*(d+1)+1  , (ii-1)*(d+1)+ d +1 ) ) =  xr.t();
  //}

  double y = x(k-1,i-1); double aV = dot(ai, VV);
  double  e = y - aV;  //residual
  //vec Vt = transpose(V)
  //double vv =dot(VV,VV);
  arma::vec  H_ik_Wald = - y*VV +  VV*VV.t()*ai - e*VV ;
  //vec Y(size(H_ik_Wald));
  //Y.fill(vv);
  return H_ik_Wald;

}


// [[Rcpp::export(makeH_k_Wald_RCPP)]] //makeH_k_Wald
arma::vec H_k_Wald(arma::mat& x, int k, int p, arma::vec a)
{
  int d = x.n_cols;
  arma::vec H = zeros(d+  d*d*p); //accounts for intercept

  for (int ii=1; ii< d+1; ii++) {
  //H.subvec((1-1)*(d*p+1)+1-1, (1-1)*(d*p+1)+ d*p +1-1 )= H_ik_Wald(x,1,k,p,a) ;
      H.subvec((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1 )= H_ik_Wald(x,ii,k,p,a) ;
      //H = join_cols(H1, H_ik_Wald(x,ii,k,p,a));
  };
  return H;
}


// [[Rcpp::export(makeH_l_u_RCPP)]] //makeH_l_u
arma::mat H_l_u(arma::mat& x,  int p, int l,int u, arma::vec a)
{
  int n = x.n_rows;
  int d = x.n_cols;
  int nr = d+ d*d*p;
  int nc = u-l+1;
  arma::mat H; H.zeros(nr,nc); //matrix of H values #accounts for intercept
  for (int t=0; t <(u-l+1); t++ ) {
          H.col(t) = H_k_Wald(x, l+t, p, a) ;//-1-1
  };
  return H;

}

struct Matrix_types {
  arma::mat m;
  arma::sp_mat M;
};
// [[Rcpp::export]] // INTERLEAVE rows of input list
arma::mat write_rows2(Rcpp::List data,  int nrows, int ncols) {//Rcpp::CharacterVector clss,

  const int len = data.length();
  std::vector<Matrix_types> matr(len);
  std::vector<bool> is_dense(len);
  arma::mat result(nrows*len, ncols);

  // populate the structs
  for (int j = 0; j < len; j++) {
    //is_dense[j] = (clss[j] == "matrix");
    //if (is_dense[j]) {
      matr[j].m = Rcpp::as<arma::mat>(data[j]);
    //}
    //else {
    //  matr[j].M = Rcpp::as<arma::sp_mat>(data[j]);//SPARSE M
    //}
  }

  // populate the result
  for (int i = 0, k = 0; i < nrows; i++) {
    for (int j = 0; j < len; j++, k++) {
      //if (is_dense[j]) {
        result.row(k) = matr[j].m.row(i);
      //}
      //else {
      //  arma::rowvec r(matr[j].M.row(i));
      //  result.row(k) = r;
      //}
    }
  }
  return result;
}

// [[Rcpp::export(get_a_lu_i_RCPP)]] //get_a_lu_i
arma::vec a_lu_i(arma::mat& x, int i, int p, int l, int u)
{
    int d = x.n_cols;
    arma::vec y =   x( span(l-1,u-1),i-1 )  ;
    arma::mat Y = diagmat(y);
    arma::mat O = ones(y.n_elem,1) ;
    arma::mat X;
    arma::rowvec y_soln;
    if(p==1)  {
       X = join_rows(O, x.rows( l-1-1,u-1-1) ); //regressors ##p = 1 case only
       y_soln = sum(Y * X);
    }
     if(p==2){
       List L = List::create(x.rows( l-1-1,u-1-1).t(), x.rows( l-2-1,u-2-1).t() ); //construct data list
       X = join_rows(O, write_rows2(L,d,y.n_elem ).t() ) ;
       y_soln = sum(Y * X);
     }
     if(p==3){
       List L = List::create(x.rows( l-1-1,u-1-1).t(), x.rows( l-2-1,u-2-1).t(),x.rows( l-3-1,u-3-1).t() ); //construct data list
       X = join_rows(O, write_rows2(L,d,y.n_elem ).t() ) ;
       y_soln = sum(Y * X);
     }
     if(p==4){
       List L = List::create(x.rows( l-1-1,u-1-1).t(), x.rows( l-2-1,u-2-1).t(),x.rows( l-3-1,u-3-1).t(),x.rows( l-4-1,u-4-1).t() ); //construct data list
       X = join_rows(O, write_rows2(L,d,y.n_elem ).t() ) ;
       y_soln = sum(Y * X);
     }
     if(p==5){
       List L = List::create(x.rows( l-1-1,u-1-1).t(), x.rows( l-2-1,u-2-1).t(),x.rows( l-3-1,u-3-1).t(),x.rows( l-4-1,u-4-1).t(),x.rows( l-5-1,u-5-1).t() ); //construct data list
       X = join_rows(O, write_rows2(L,d,y.n_elem ).t() ) ;
       y_soln = sum(Y * X);
     }
     if(p==6){
       List L = List::create(x.rows( l-1-1,u-1-1).t(), x.rows( l-2-1,u-2-1).t(),x.rows( l-3-1,u-3-1).t(),x.rows( l-4-1,u-4-1).t(),x.rows( l-5-1,u-5-1).t(),x.rows( l-6-1,u-6-1).t() ); //construct data list
       X = join_rows(O, write_rows2(L,d,y.n_elem ).t() ) ;
       y_soln = sum(Y * X);
     }
    arma::mat X_soln = (X.t() * X ).i();
    arma::vec a_out = (y_soln * X_soln).t();
    return a_out;
}

// [[Rcpp::export(make_a_lu_RCPP)]] //make_a_lu
arma::vec make_a_lu(arma::mat& x, int p, int l, int u)
{
  int d = x.n_cols;
  arma::vec a_lu = zeros(d + d*d*p);
  for (int ii=1; ii< d+1; ii++) {
    a_lu.subvec( (ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1 ) = a_lu_i(x,ii,p,l,u)  ;
  };
  return a_lu;
}


// // [[Rcpp::export]]
arma::sp_mat blockDiag_Wald( arma::field<arma::mat>& xlist ) {

  //xlist: list of matrices

  unsigned int n = xlist.n_rows ;
  int dimen = 0 ;
  arma::ivec dimvec(n) ;

  for(unsigned int i=0; i<n; i++) {
    dimvec(i) = xlist(i,0).n_rows ;
    dimen += dimvec(i) ;
  }

  sp_mat X(dimen,dimen);
  int idx=0;

  for(unsigned int i=0; i<n; i++) {
    X.submat( idx, idx, idx + dimvec(i) - 1, idx + dimvec(i) - 1 ) = xlist(i,0) ;
    idx = idx + dimvec(i) ;
  }
  return(X);
}


// [[Rcpp::export(get_V_nk_RCPP)]] //get_V_nk
arma::sp_mat V_nk(arma::mat x, int p, int l, int u)
 {
  int d = x.n_cols;
  arma::mat xk;
  arma::mat O = ones(u-l+1,1) ;
  if(p==1) xk = join_rows(O, x.rows( l-1-1,u-1-1) );
  if(p==2) xk = join_rows(O, x.rows( l-1-1,u-1-1) , x.rows( l-2-1,u-2-1));
  if(p==3) xk = join_rows(O, x.rows( l-1-1,u-1-1) , x.rows( l-2-1,u-2-1), x.rows( l-3-1,u-3-1));
  if(p==4) {
    xk = join_rows(x.rows( l-1-1,u-1-1) , x.rows( l-2-1,u-2-1), x.rows( l-3-1,u-3-1), x.rows( l-4-1,u-4-1));
    xk.insert_cols(0,O);
  };
  if(p==5) {
    xk = join_rows( x.rows( l-2-1,u-2-1), x.rows( l-3-1,u-3-1), x.rows( l-4-1,u-4-1), x.rows( l-5-1,u-5-1));
    xk.insert_cols(0, join_rows(O,x.rows( l-1-1,u-1-1)  ));
  };
  if(p==6) {
    xk = join_rows(  x.rows( l-3-1,u-3-1), x.rows( l-4-1,u-4-1), x.rows( l-5-1,u-5-1), x.rows( l-6-1,u-6-1));
    xk.insert_cols(0, join_rows(O,x.rows( l-1-1,u-1-1),x.rows( l-2-1,u-2-1)  ));
  };
  arma::mat C =  xk.t() * xk /(u-l);
  arma::field<arma::mat> Vlist(d);
  for (int ii=0; ii< d; ii++) {
    Vlist(ii) = C; //coerce into block diagonal form
  };
  arma::sp_mat out = blockDiag_Wald(Vlist);
  return out;
 }

// [[Rcpp::export(getsigma_i_kLOCAL1_RCPP)]] //getsigma_i_kLOCAL1
arma::vec sigma_i_k(arma::mat x, int i,int k,int G,int p,arma::vec a_upper) //,vec a_lower)
{
  arma::mat x_upper; arma::mat x_lower;
  arma::mat O = ones(G,1) ;
  if(p==1) x_upper = join_rows(O, x.rows( k-1,k+G-1-1) ); //upper sample
  if(p==2) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ) );
  if(p==3) x_upper = join_rows(O, x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ) );
  if(p==4){
    x_upper = join_rows(x.rows( k-1,k+G-1-1), x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ),x.rows(k-3-1, k+G-4-1 ) );
    x_upper.insert_cols(0,O);
  } ;
  if(p==5){
    x_upper = join_rows( x.rows(k-1-1, k+G-2-1 ), x.rows(k-2-1, k+G-3-1 ),x.rows(k-3-1, k+G-4-1 ),x.rows(k-4-1, k+G-5-1 ) );
    x_upper.insert_cols(0, join_rows(O,x.rows( k-1,k+G-1-1)) );
  } ;
  if(p==6){
    x_upper = join_rows(  x.rows(k-2-1, k+G-3-1 ),x.rows(k-3-1, k+G-4-1 ),x.rows(k-4-1, k+G-5-1 ),x.rows(k-5-1, k+G-6-1 ) );
    x_upper.insert_cols(0, join_rows(O,x.rows( k-1,k+G-1-1),x.rows(k-1-1, k+G-2-1 )   ));
  } ;
  arma::rowvec res_upper =  x( span(k+1-1,k+G-1), i-1).t() - a_upper.t() * x_upper.t(); //upper residuals

  // if(p==1)x_lower = join_rows(O, x.rows( k-G-1,k-1-1) );
  // if(p==2)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1));
  // if(p==3)x_lower  = join_rows(O, x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1), x.rows(k-G-3, k-3-1));
  // if(p==4){
  //   x_lower = join_rows(x.rows( k-G-1,k-1-1), x.rows(k-G-2, k-2-1), x.rows(k-G-3, k-3-1),x.rows(k-G-4, k-4-1) );
  //   x_lower.insert_cols(0,1);
  // };
  // rowvec res_lower =  x( span(k-G+1-1,k-1), i-1).t() - a_lower.t() * x_lower.t(); //lower residuals
 // double sigma_i =  sum(square(res_upper))/G;  //+ sum(square(res_lower)) ) /(2*G);
  return res_upper.t();

}

// [[Rcpp::export(getsigma_d_kLOCAL1_RCPP)]] //getsigma_d_kLOCAL1
arma::mat sigma_d_k(arma::mat x, int k,int G,int p,arma::vec a_upper, arma::vec a_lower)
{
  int d = x.n_cols;
  arma::mat sigma_u = zeros(G,d); arma::mat sigma_l = zeros(G,d);
   for (int ii=1; ii< d+1; ii++) {
       sigma_u.col(ii-1) = sigma_i_k(x,ii,k,G,p, a_upper( span((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1)) );
       sigma_l.col(ii-1) = sigma_i_k(x,ii,k-G,G,p, a_lower( span((ii-1)*(d*p+1)+1-1, (ii-1)*(d*p+1)+ d*p +1-1)) );
   };
    return (sigma_u.t() *sigma_u + sigma_l.t()* sigma_l)/(G);// cov(sigma_u)+cov(sigma_l);
}




// SIGMA ESTIMATORS ------------------------

// [[Rcpp::export(get_DiagH_Wald_RCPP)]] //get_DiagH_Wald
arma::sp_mat diagH_Wald(arma::mat x, int G,int p,arma::mat H_l, arma::mat H_u)//sp_mat
{
  int n = x.n_rows;
  int d = x.n_cols;
  arma::field<arma::mat> H_list(d);
  arma::mat H_u_i; arma::mat H_l_i; arma::mat H_out;
  arma::vec Hbar_u; arma::vec Hbar_l;
  arma::vec evals; arma::mat evecs;
  for (int ii=1; ii< d+1; ii++) {
        H_u_i = H_u.rows( ((ii-1)*(d*p+1)+1-1), ((ii-1)*(d*p+1)+ d*p +1-1) );
        H_l_i = H_l.rows( ((ii-1)*(d*p+1)+1-1), ((ii-1)*(d*p+1)+ d*p +1-1) );
        Hbar_u = mean(H_u_i,1) ;
        H_u_i.each_col() -= Hbar_u ;//H_u_centred
        Hbar_l = mean(H_l_i,1) ;
        H_l_i.each_col() -= Hbar_l; //H_l_centred
        H_out = H_l_i * H_l_i.t() + H_u_i * H_u_i.t();
        //eigen decomposition
        eig_sym(evals,evecs,H_out);
        H_out = evecs * diagmat(pow(evals,-0.5))*  evecs.t();
        H_list(ii-1) = H_out;
  }
  arma::sp_mat Sig_ =  sqrt(2*G) * blockDiag_Wald(H_list); //coerce into block diagonal form
  return Sig_;
}



// [[Rcpp::export(get_FullH_Wald_RCPP)]] //get_FullH_Wald
arma::mat FullH_Wald(arma::mat x, int G,arma::mat H_l, arma::mat H_u)
{
  int n = x.n_rows;
  int d = x.n_cols;
  //field<mat> H_list(d);
  arma::mat H_out;
  arma::vec Hbar_u; arma::vec Hbar_l;
  arma::vec evals; arma::mat evecs;

    Hbar_u = mean(H_u,1) ;
    H_u.each_col() -= Hbar_u ;//H_u_centred
    Hbar_l = mean(H_l,1) ;
    H_l.each_col() -= Hbar_l; //H_l_centred
    H_out = H_l * H_l.t() + H_u * H_u.t();
    //eigen decomposition
    eig_sym(evals,evecs,H_out);
    H_out = evecs * diagmat(pow(evals/(2*G),-0.5))*  evecs.t();


    //DEAL WITH OVERFLOW
    // int n = x.n_rows;
    // int d = x.n_cols;
    // //field<mat> H_list(d);
    // mat H_; Mat<long double> H_out;
    // vec Hbar_u; vec Hbar_l;
    // vec evals; mat evecs; std::vector<long double> evld;
    //
    // Hbar_u = mean(H_u,1) ;
    // H_u.each_col() -= Hbar_u ;//H_u_centred
    // Hbar_l = mean(H_l,1) ;
    // H_l.each_col() -= Hbar_l; //H_l_centred
    // H_ = H_l * H_l.t() + H_u * H_u.t();
    // //eigen decomposition
    // eig_sym(evals,evecs,H_);
    // evld = conv_to<std::vector<long double>>::from(evals);
    // //Col<long double> pp = conv_to<Col<long double>>::from(powl(evld/(2*G),-0.5) );
    // Mat<long double> D = diagmat()evld;
    // H_out = evecs * (D)*  evecs.t();
  return H_out; // RETURNS NaNs
}


// [[Rcpp::export(get_DiagC_Wald_RCPP)]] //get_DiagC_Wald
arma::mat DiagC_Wald(arma::mat x, int p, arma::mat sigma_d, int k, int G)
{
  int n = x.n_rows;
  int d = x.n_cols;
  arma::mat xk;arma::vec evals; arma::mat evecs;arma::vec evalS; arma::mat evecS;
  arma::mat O = ones(2*G,1) ;
  if(p==1) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) );
  if(p==2) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1));
  if(p==3) xk = join_rows(O, x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1));
  if(p==4) {
    xk = join_rows(x.rows( k-G-1,k+G-1-1) , x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1));
    xk.insert_cols(0,O);
  };
  if(p==5) {
    xk = join_rows(x.rows( k-G-1-1,k+G-2-1),x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1),x.rows( k-G-4-1,k+G-5-1));
    xk.insert_cols(0, join_rows(O,x.rows( k-G-1,k+G-1-1)));
  };
  if(p==6) {
    xk = join_rows(x.rows( k-G-2-1,k+G-3-1), x.rows( k-G-3-1,k+G-4-1),x.rows( k-G-4-1,k+G-5-1),x.rows( k-G-5-1,k+G-6-1));
    xk.insert_cols(0, join_rows(O,x.rows( k-G-1,k+G-1-1), x.rows( k-G-1-1,k+G-2-1) ));
  };

  arma::mat C =  xk.t() * xk /(2*G);
  //eigen decomposition
  eig_sym(evals,evecs,C);
  arma::mat C_ = evecs * diagmat(1/sqrt(evals) )*  evecs.t(); //sqrt
  eig_sym(evalS,evecS,sigma_d);
  arma::mat S_ = evecS * diagmat(pow(2*evalS, -0.5) )*  evecS.t(); //inverse sqrt
  // field<mat> Clist(d);
  // for (int ii=0; ii< d; ii++) {
  //   Clist(ii) =   C_ /sqrt(2*sigma_d(ii));
  // };
  arma::mat out = kron(S_, C_);//blockDiag_Wald(Clist); //coerce into block diagonal form

  return out;
}


// [[Rcpp::export(get_Wkn_RCPP)]] //get_Wkn
double Wkn(arma::mat x, int p, int k, int G, String estim = "C")
{
  int n = x.n_rows;
  int d = x.n_cols;

  arma::vec a_upper = make_a_lu(x, p, k+1, k+G);
  arma::vec a_lower = make_a_lu(x, p, k-G+1, k);
  arma::vec W_mat;
  //double v;
//Sigma estimator options------
  arma::sp_mat V = V_nk(x, p, k-G+1, k);
    if(estim == "C"){
      arma::mat sigma_d = sigma_d_k(x,k,G,p,a_upper, a_lower);
      arma::mat Sig_ = DiagC_Wald(x,p,sigma_d,k,G) ;
      W_mat = Sig_ * V * (a_upper-a_lower) ;
    } else{
       arma::mat H_l = H_l_u(x, p, k-G+1, k , a_lower);
       arma::mat H_u = H_l_u(x, p, k+1, k+G , a_upper);//v= H_u(0,0);
       if(estim == "H") {
         arma::sp_mat Sig_ = diagH_Wald(x,G,p,H_l,H_u); //DiagH estimator for Sigma
        W_mat = Sig_ * V * (a_upper-a_lower) ;
       }
       // if(estim == "FullH") {
       //   arma::mat Sig_ = FullH_Wald(x,G,H_l,H_u);  //FullH estimator
       //   W_mat = Sig_ * V * (a_upper-a_lower) ;
       // }
    }
//------------------------------
  double W = sqrt(G/2) * norm(W_mat, "fro");
  return W;
}


double floop(int k) {
  arma::mat x; int p; int G; String estim;
  return Wkn(x,p,k,G,estim);
} //wrapper

// [[Rcpp::export(get_W_RCPP)]] //get_W
arma::vec W(arma::mat x, int p, int G, String estim = "C") //int ncores=1)
{

  int n = x.n_rows;
  arma::vec out = zeros( n);
  //RcppParallel::RVector<double> wo(out);
  //RcppParallel::RMatrix<double> wx(x);
  //
  // #if defined(_OPENMP)
  // #pragma omp parallel for num_threads(ncores)
  // #endif
  // for(size_t ii = 0; ii < n; ii++)
  // {
  //   wo[ii] = boost::math::erf(wx[ii]);
  // }

  for (int k=G+ p+1; k< n-G-p; k++) {
      //wo(k) = Wkn(wx,p,k,G,estim); }
      out(k) = Wkn(x,p,k+p-1,G,estim); }
  // Col<int> K = conv_to<Col<int>>::from(linspace(G+p, n-G-2) );
  // std::transform(K.begin(),K.end(), out(G+p), floop );
  return out; //wo
}



//SIMULATION -------------------------------------


// [[Rcpp::export(get_cps_RCPP)]] //get_cps
arma::vec cps_Wald(arma::vec Wn, double D_n, int G, double nu = 0.25)
{
    int n = Wn.size(); arma::vec out;
    //rshift <- c(Tn[-1],0); lshift <- c(0,Tn[-n]);
    arma::vec rshift = shift(Wn,-1); arma::vec lshift = shift(Wn,1);
    arma::uvec over = find(Wn >D_n); //indices are greater than D_n?
    arma::uvec lunder = find(lshift < D_n);
    arma::uvec v = intersect(over, lunder); //lowers

    arma::uvec runder = find(rshift < D_n);
    arma::uvec w = intersect(over, runder); //uppers
    arma::uvec nu_remove = find(w-v >= nu*G); //(epsilon) test for distance between
    arma::uvec v_nu = v.elem(nu_remove); arma::uvec w_nu = w.elem(nu_remove); //w_nu.insert_rows(w_nu.size(),1); w_nu(w_nu.size()-1)=n;
    int q = nu_remove.size(); //number of CPs
    if(q>0){
       out = zeros(q);
       for (int ii=0; ii< q; ii++) {
             out(ii) = v_nu(ii) + Wn( span(v_nu(ii)-1,w_nu(ii)-1)).index_max() ;
       };
       };
    return out;
}

// [[Rcpp::export(test_Wald_RCPP)]] //test_Wald
List test_Wald(arma::mat x, int p, int G, double alpha =0.05, String estim = "C"){
  int n = x.n_rows;
  int d = x.n_cols;
  arma::vec cp; //double nu=1/4;
//Test setup----------------------------
      double c_alpha = -log(log( pow((1-alpha),-0.5)) ); //critical value
      double a = sqrt(2*log(n/G)); //test transform multipliers
      double g = lgamma(d*(d*p+1)/2); double l23 = log(2)-log(3);
      double b = 2*log(n/G) + (d*(d*p+1)/2) * log(log(n/G)) -g - l23   ;
      // D_n =  ;//threshold
      double D_n = max((b + c_alpha)/a, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) ); //##ASYMPTOTIC correction
      int Reject = 0; //
//Run test-----------------------------
      arma::vec Wn = W(x,p,G,estim); //evaluate statistic at each time k
      double test_stat = Wn.max();
      if(test_stat > D_n){ //compare test stat with threshold
        Reject = 1; //test true
         cp = cps_Wald(Wn,D_n,G);
         if( cp.size()==0 ) Reject = 0 ;//doesn't pass eps-test
      } ;
//Plot------------------------------------
//       plot(Wn) # plot test statistic
//         abline(h = D_n, col = "blue") #add threshold
//         if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
//           pl <- recordPlot()
// #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
//Output------------------------------------
      List out = List::create(Named("Reject") = Reject,  _["ChangePoints"] = cp);//, _["D_n"]=D_n, _["Wn"] = Wn,);
      return out ;
}


// [[Rcpp::export(sim_data_RCPP)]]
arma::mat sim_data(List pars, int n=1000, int d=5, double sd = 0.2){
  arma::mat errors = randn(n,d)*sd;
  arma::mat simdata = errors;
  int p = pars.length();
  //simdata.rows(0,p-1) = errors.rows(0,p-1);
  for(int r =p; r < n; r++){
    //simdata.row(r)   errors[row,]
    for (int ii=0; ii<p;  ii++) {
      arma::mat mp = pars[ii];
      arma::vec pp = mp * simdata.row(r-ii).t();
      simdata.row(r) = simdata.row(r) +  pp.t() ;
    };
  };
  return simdata;
}



// // [[Rcpp::export(var_simulate_RCPP)]]
// NumericVector var_sim(List pars, int reps =100, int p=2, int G=200, double alpha =0.05, String estim = "DiagC",int ncores =1){
//   vec cp={500,1000,1500};
//   NumericVector out(reps) ;
//   RcppParallel::RVector<double> wo(out);
//   //RcppParallel::RVector<double> wx(x);
//
//   #if defined(_OPENMP)
//   #pragma omp parallel for num_threads(ncores)
//   #endif
//    for(int repl = 0; repl < reps; repl++ ){
//    List p1 = List::create(pars[0], pars[1]);List p2 = List::create(pars[2], pars[3]);List p3 = List::create(pars[4], pars[5]);
//    mat r1 = sim_data(p1, cp(0), 5);mat r2 = sim_data(p2, cp(1)-cp(0), 5);mat r3 = sim_data(p3, cp(2)-cp(1), 5); //in-regime data
//    mat r = join_cols(r1,r2,r3); //full data
//    List t = test_Wald(r, p, G, alpha,estim);
//    wo[repl] = t[0];
//    };
//   //Output------------------------------------
//   // List out = List::create(Named("Reject") = Reject, _["Wn"] = Wn, _["ChangePoints"] = cp, _["D_n"]=D_n);
//   return out ;
// }




// [[Rcpp::export(MFA_Wald)]]
List MFA_Wald(arma::mat x, int p, arma::vec Gset, String estim = "C", double alpha = 0.05, double nu = 0.25){
  Gset = sort(Gset); //into ascending order
  int Glen = Gset.size();
  NumericVector cps;
  bool Reject = FALSE;
  NumericVector v(Glen); List tests(Glen);
  List t;
  for(int ii =0; ii < Glen; ii++){
    t =  test_Wald(x, p, Gset[ii], alpha, estim);
    tests[ii] = t;
    if(t["Reject"]){ Reject = TRUE;}
  }
  if(Reject){
   cps = as<List>(tests[0])["ChangePoints"];
   if(Glen > 1){
     for(int ii=1;  ii<Glen; ii++){
       NumericVector K = as<List>(tests[ii])["ChangePoints"];
       for(int jj=0; jj < K.size(); jj++) {
         if(min(abs(K[jj] - cps)) > nu*Gset[ii]) {cps.push_back(K[jj]);}
       }
     }
   }
  } //list(Reject, cps, q= length(cps))
  return(List::create(Named("Reject") = Reject, _["ChangePoints"]=cps, _["q"] = cps.size() ));
}
