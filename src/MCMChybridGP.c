#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

SEXP calcSigma(SEXP X, SEXP dsq_eta)
{
  int n, d;
  int *dims;
  double *Sigma_p, *X_p, *dsq_eta_p;
  SEXP Sigma;
    
  dims = INTEGER(coerceVector(getAttrib(X, R_DimSymbol), INTSXP));
  n = dims[0];
  d = dims[1];
    
  PROTECT(X = coerceVector(X, REALSXP));
  PROTECT(dsq_eta = coerceVector(dsq_eta, REALSXP));
  PROTECT(Sigma = allocMatrix(REALSXP, n, n));
  X_p = REAL(X);
  dsq_eta_p = REAL(dsq_eta);
  Sigma_p = REAL(Sigma);
    
  int i, j, k;
  double prod, h;
    
  for(i=0; i<n; i++) for(j=0; j<n; j++) {
    prod = 1.0;
    for(k=0; k<d; k++) {
      h = X_p[i + k*n] - X_p[j + k*n];
      h = fabs(h);
      prod *= (1 + dsq_eta_p[k+1]*h) * exp(-dsq_eta_p[k+1]*h);
      }
    Sigma_p[i + j*n] = dsq_eta_p[0] * prod;
    }
  UNPROTECT(3);
  return(Sigma);
}

SEXP Leapfrog(SEXP X_mat, SEXP y_vec, SEXP x_vec, SEXP p_vec, SEXP L_int,
    SEXP delta_num, SEXP Sigma_inv_mat, SEXP dsq_eta_vec, SEXP T_num, SEXP max_sig,
    SEXP lb_vec, SEXP ub_vec)
{
  int n, d;
  int *Xdims;
  double *X, *y, *x, *p, *lb, *ub;
  int L;
  double delta, maxsig, T, *Sigma_inv, *dsq_eta;
  Xdims = INTEGER(coerceVector(getAttrib(X_mat, R_DimSymbol), INTSXP));
  n = Xdims[0]; d = Xdims[1];
  SEXP c_x, cx_vec, ans_1, ans_2;
  double *cx, *cxvec, *ans1, *ans2;
  double prod, h, h1, val, sum;
  int i,j, i1, ii, jj , kk, k;
  int OUTSIDE;
  
  SEXP X_i, x_p_E;
  SEXP y_star;
  
  PROTECT(X_mat = coerceVector(X_mat, REALSXP));
  PROTECT(y_vec = coerceVector(y_vec, REALSXP));
  PROTECT(x_vec = coerceVector(x_vec, REALSXP));
  PROTECT(p_vec = coerceVector(p_vec, REALSXP));
  PROTECT(L_int = coerceVector(L_int, INTSXP));
  PROTECT(delta_num = coerceVector(delta_num, REALSXP));
  PROTECT(Sigma_inv_mat = coerceVector(Sigma_inv_mat, REALSXP));
  PROTECT(dsq_eta_vec = coerceVector(dsq_eta_vec, REALSXP));
  PROTECT(T_num = coerceVector(T_num, REALSXP));
  PROTECT(max_sig = coerceVector(max_sig, REALSXP));
  PROTECT(lb_vec = coerceVector(lb_vec, REALSXP));
  PROTECT(ub_vec = coerceVector(ub_vec, REALSXP));
  PROTECT(c_x = allocVector(REALSXP, n));
  PROTECT(cx_vec = allocVector(REALSXP, n));
  PROTECT(ans_1 = allocVector(REALSXP, n));
  PROTECT(ans_2 = allocVector(REALSXP, 1));
  PROTECT(X_i = allocVector(REALSXP, d));
  PROTECT(x_p_E = allocMatrix(REALSXP, d, 3));
  PROTECT(y_star = allocMatrix(REALSXP, n, 1));
  
  X = REAL(X_mat);
  y = REAL(y_vec);
  x = REAL(x_vec);
  p = REAL(p_vec);
  L = INTEGER(L_int)[0];
  delta = REAL(delta_num)[0];
  Sigma_inv = REAL(Sigma_inv_mat);
  dsq_eta = REAL(dsq_eta_vec);
  T = REAL(T_num)[0];
  maxsig = REAL(max_sig)[0];
  lb = REAL(lb_vec);
  ub = REAL(ub_vec);
  cx = REAL(c_x);
  cxvec = REAL(cx_vec);
  ans1 = REAL(ans_1);
  ans2 = REAL(ans_2);
  
  double *Xi, *xpE;
  Xi = REAL(X_i);
  xpE = REAL(x_p_E);
  
  for(i=0; i<d; i++) {
    xpE[i +d* 0] = x[i];
    xpE[i +d* 1] = p[i];
    }
  
  double meanY = 0.0;
  double *ystar;
  ystar = REAL(y_star);
  for(i=0; i<n; i++) meanY += y[i];
  meanY /= n;
  for(i=0; i<n; i++) ystar[i] = y[i] - meanY;
  
  
  SEXP x_2, p_2, dZ_old, dZ_new;
  double *x2, *p2, *dZold, *dZnew;
  int Lexp;
  PROTECT(x_2 = allocVector(REALSXP, d));
  PROTECT(p_2 = allocVector(REALSXP, d));
  PROTECT(dZ_old = allocVector(REALSXP, d));
  PROTECT(dZ_new = allocVector(REALSXP, d));
  x2 = REAL(x_2);
  p2 = REAL(p_2);
  dZold = REAL(dZ_old);
  dZnew = REAL(dZ_new);
  
  for(i=0; i<d; i++) {
    x2[i] = x[i];
    p2[i] = p[i];
    }
  double sigmaf_x = 0.0;
  
  
//  dZ(x=x2, val=dZold);
  for(ii=0; ii<n; ii++) {
    for(jj=0; jj<1; jj++) {
      sum = 0;
      for(kk=0; kk<n; kk++)
        sum = sum + Sigma_inv[ii + n*kk] * ystar[kk + n*jj];
      ans1[ii+n*jj] = sum;
    }
  }
  for(k=0; k<d; k++) {
    for(i=0; i<n; i++) {
      for(j=0; j<d; j++) Xi[j] = X[i + n*j];
      prod=1.0;
      h = x2[k] - Xi[k];
      h = fabs(h);
      for(i1=0; i1<d; i1++) {
        h1 = x2[i1] - Xi[i1];
        h1 = fabs(h1);
        prod *= (1.0 + dsq_eta[1+i1]*h1) * exp(-dsq_eta[1+i1]*h1);
      }
      cxvec[i] = dsq_eta[0] * prod *
          dsq_eta[1+k] * dsq_eta[1+k] *
          (x2[k]-Xi[k]) / (1.0 + dsq_eta[1+k]*h);
    }

    for(ii=0; ii<1; ii++) {
      for(jj=0; jj<1; jj++) {
        sum = 0;
        for(kk=0; kk<n; kk++)
          sum = sum + cxvec[ii + 1*kk] * ans1[kk + n*jj];
        ans2[ii+1*jj] = sum;
      }
    }
    dZold[k] = ans2[0]/T;
  }

  for(Lexp=0; Lexp<L; Lexp++) {
    for(i=0; i<d; i++) x2[i] += delta*p2[i] - delta*delta/2.0*dZold[i];
//    dZ(x=x2, val=dZnew);
    double prod, h, h1, sum;
    int ii, jj , kk, k, i,j, i1;
    for(ii=0; ii<n; ii++) {
      for(jj=0; jj<1; jj++) {
        sum = 0;
        for(kk=0; kk<n; kk++)
          sum = sum + Sigma_inv[ii + n*kk] * ystar[kk + n*jj];
        ans1[ii+n*jj] = sum;
      }
    }
    for(k=0; k<d; k++) {
      for(i=0; i<n; i++) {
        for(j=0; j<d; j++) Xi[j] = X[i + n*j];
        prod=1.0;
        h = x2[k] - Xi[k];
        h = fabs(h);
        for(i1=0; i1<d; i1++) {
          h1 = x2[i1] - Xi[i1];
          h1 = fabs(h1);
          prod *= (1.0 + dsq_eta[1+i1]*h1) * exp(-dsq_eta[1+i1]*h1);
        }
        cxvec[i] = dsq_eta[0] * prod *
            dsq_eta[1+k] * dsq_eta[1+k] *
            (x2[k]-Xi[k]) / (1.0 + dsq_eta[1+k]*h);
      }
  
      for(ii=0; ii<1; ii++) {
        for(jj=0; jj<1; jj++) {
          sum = 0;
          for(kk=0; kk<n; kk++)
            sum = sum + cxvec[ii + 1*kk] * ans1[kk + n*jj];
          ans2[ii+1*jj] = sum;
        }
      }
      dZnew[k] = ans2[0]/T;
    }

    for(i=0; i<d; i++) p2[i] -= delta/2.0 * (dZold[i] + dZnew[i]);
    for(i=0; i<d; i++) dZold[i] = dZnew[i];
//  sigmaf_x = sigmaf(x=x2);
    for(i=0; i<n; i++) {
      for(j=0; j<d; j++) Xi[j] = X[i + n*j];
      prod = 1.0;
      for(i1=0; i1<d; i1++) {
        h = x2[i1] - Xi[i1];
        h = fabs(h);
        prod *= (1.0 + dsq_eta[1+i1]*h) * exp(-dsq_eta[1+i1]*h);
      }
      cxvec[i] = dsq_eta[0] * prod;
    }

    for(ii=0; ii<n; ii++) {
      for(jj=0; jj<1; jj++) {
        sum = 0;
        for(kk=0; kk<n; kk++)
          sum = sum + Sigma_inv[ii + n*kk] * cxvec[kk + n*jj];
        ans1[ii+n*jj] = sum;
      }
    }

    for(i=0; i<n; i++) {
      for(j=0; j<d; j++) Xi[j] = X[i + n*j];
      prod = 1.0;
      for(i1=0; i1<d; i1++) {
        h = x2[i1] - Xi[i1];
        h = fabs(h);
        prod *= (1.0 + dsq_eta[1+i1]*h) * exp(-dsq_eta[1+i1]*h); }
      cx[i] = dsq_eta[0] * prod;
    }

    for(ii=0; ii<1; ii++) {
      for(jj=0; jj<1; jj++) {
        sum = 0;
        for(kk=0; kk<n; kk++)
          sum = sum + cx[ii + 1*kk] * ans1[kk + n*jj];
        ans2[ii+1*jj] = sum;
      }
    }

    val = dsq_eta[0] - ans2[0];
    sigmaf_x = sqrt(val);

    if(sigmaf_x > maxsig) {
      sigmaf_x = -sigmaf_x;
      break;
    }
  }

  for(i=0; i<d; i++) do {
    OUTSIDE = 0;
    if(x2[i] < lb[i]) x2[i] = 2.0*lb[i] - x2[i];
    if(x2[i] > ub[i]) x2[i] = 2.0*ub[i] - x2[i];
    if(x2[i] < lb[i]) OUTSIDE = 1;
    if(x2[i] > ub[i]) OUTSIDE = 1;
  }while(OUTSIDE);
  
  for(i=0; i<d; i++) {
    xpE[i +d* 0] = x2[i];
    xpE[i +d* 1] = p2[i];
    xpE[i +d* 2] = sigmaf_x;
  }
  
  UNPROTECT(23);
  return(x_p_E);
  }


