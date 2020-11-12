#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

void GPcovar(double* X, int* n_ptr, int* d_ptr, double* dsq_eta,
             double* Sigma)
{
  int n = n_ptr[0], d = d_ptr[0];
  int i, j, k;
  double prod, h;
  for(i=0; i<n; i++) for(j=0; j<n; j++) {
    prod = 1.0;
    for(k=0; k<d; k++) {
      h = X[i +k* n] - X[j +k* n];
      h = fabs(h);
      prod *= (1 + dsq_eta[k+1]*h) * exp(-dsq_eta[k+1]*h);
    }
    Sigma[i +j* n] = dsq_eta[0] * prod;
  }
}

void Leap(double* X, int* n_ptr, int* d_ptr, double* y, double* x,
          double* p, int* L_ptr, double* delta_ptr, double* Sigma_inv,
          double* dsq_eta, double* T_ptr, double* maxsig_ptr,
          double* lb, double* ub, double* cx, double* ans1,
          double* dZold, double* dZnew, double* xpE)
{
  int n = n_ptr[0], d = d_ptr[0], L = L_ptr[0];
  double delta = delta_ptr[0], maxsig = maxsig_ptr[0], T = T_ptr[0];
  double prod = 1.0, h = 0.0, h1 = 0.0, val = 0.0, sum = 0.0;
  int i = 0, j = 0, ii = 0, jj = 0;
  int OUTSIDE;
  double meanY = 0.0;
  int Lexp;
  double sigmaf_x = 0.0;
  
  for(i=0; i<n; i++) meanY += y[i];
  meanY /= n;

//  dZ(x=x2, val=dZold);
  for(i=0; i<n; i++) {
      sum = 0;
      for(j=0; j<n; j++)
        sum = sum + Sigma_inv[i + n*j] * (y[j]-meanY);
      ans1[i] = sum;
  }
  for(j=0; j<d; j++) {
    for(i=0; i<n; i++) {
      prod=1.0;
      h = x[j] - X[i +n* j];
      h = fabs(h);
      for(ii=0; ii<d; ii++) {
        h1 = x[ii] - X[i + n*ii];
        h1 = fabs(h1);
        prod *= (1.0 + dsq_eta[1+ii]*h1) * exp(-dsq_eta[1+ii]*h1);
      }
      cx[i] = dsq_eta[0] * prod *
          dsq_eta[1+j] * dsq_eta[1+j] *
          (x[j]-X[i +n* j]) / (1.0 + dsq_eta[1+j]*h);
    }

    sum = 0;
    for(jj=0; jj<n; jj++)
      sum = sum + cx[1*jj] * ans1[jj];
    dZold[j] = sum/T;
  }

  for(Lexp=0; Lexp<L; Lexp++) {
    for(i=0; i<d; i++) x[i] += delta*p[i] - delta*delta/2.0*dZold[i];
//    dZ(x=x2, val=dZnew);
    for(i=0; i<n; i++) {
        sum = 0;
        for(j=0; j<n; j++)
          sum = sum + Sigma_inv[i + n*j] * (y[j]-meanY);
        ans1[i] = sum;
    }
    for(j=0; j<d; j++) {
      for(i=0; i<n; i++) {
        prod=1.0;
        h = x[j] - X[i +n* j];
        h = fabs(h);
        for(ii=0; ii<d; ii++) {
          h1 = x[ii] - X[i + n*ii];
          h1 = fabs(h1);
          prod *= (1.0 + dsq_eta[1+ii]*h1) * exp(-dsq_eta[1+ii]*h1);
        }
        cx[i] = dsq_eta[0] * prod *
            dsq_eta[1+j] * dsq_eta[1+j] *
            (x[j]-X[i +n* j]) / (1.0 + dsq_eta[1+j]*h);
      }
  
      sum = 0;
      for(jj=0; jj<n; jj++) sum = sum + cx[jj] * ans1[jj];
      dZnew[j] = sum/T;
    }

    for(i=0; i<d; i++) p[i] -= delta/2.0 * (dZold[i] + dZnew[i]);
    for(i=0; i<d; i++) dZold[i] = dZnew[i];
//  sigmaf_x = sigmaf(x=x2);
    for(i=0; i<n; i++) {
      prod = 1.0;
      for(ii=0; ii<d; ii++) {
        h = x[ii] - X[i + n*ii];
        h = fabs(h);
        prod *= (1.0 + dsq_eta[1+ii]*h) * exp(-dsq_eta[1+ii]*h);
      }
      cx[i] = dsq_eta[0] * prod;
    }

    for(i=0; i<n; i++) {
        sum = 0;
        for(j=0; j<n; j++)
          sum = sum + Sigma_inv[i + n*j] * cx[j];
        ans1[i] = sum;
    }

    for(i=0; i<n; i++) {
      prod = 1.0;
      for(ii=0; ii<d; ii++) {
        h = x[ii] - X[i + n*ii];
        h = fabs(h);
        prod *= (1.0 + dsq_eta[1+ii]*h) * exp(-dsq_eta[1+ii]*h); }
      cx[i] = dsq_eta[0] * prod;
    }

    sum = 0;
    for(j=0; j<n; j++) sum = sum + cx[j] * ans1[j];
    val = dsq_eta[0] - sum;
    sigmaf_x = sqrt(val);

    if(sigmaf_x > maxsig) {
      sigmaf_x = -sigmaf_x;
      break;
    }
  }

  for(i=0; i<d; i++) do {
    OUTSIDE = 0;
    if(x[i] < lb[i]) x[i] = 2.0*lb[i] - x[i];
    if(x[i] > ub[i]) x[i] = 2.0*ub[i] - x[i];
    if(x[i] < lb[i]) OUTSIDE = 1;
    if(x[i] > ub[i]) OUTSIDE = 1;
  }while(OUTSIDE);
  
  for(i=0; i<d; i++) {
    xpE[i +d* 0] = x[i];
    xpE[i +d* 1] = p[i];
    xpE[i +d* 2] = sigmaf_x;
  }
}

