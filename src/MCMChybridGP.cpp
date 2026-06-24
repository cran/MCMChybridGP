#include <Rcpp.h>
using namespace Rcpp;

// ======================================================
// GPcovar
// ======================================================

// [[Rcpp::export]]
NumericMatrix GPcovar(const NumericMatrix& X,
                      const NumericVector& dsq_eta)
{
    const int n = X.nrow();
    const int d = X.ncol();

    if (dsq_eta.size() != d + 1)
        stop("length(params) must equal ncol(X)+1");

    NumericMatrix Sigma(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            double prod = 1.0;

            for (int k = 0; k < d; k++) {
                double h = std::fabs(X(i, k) - X(j, k));
                double eta = dsq_eta[k + 1];

                prod *= (1.0 + eta * h) * std::exp(-eta * h);
            }

            Sigma(i, j) = dsq_eta[0] * prod;
        }
    }

    return Sigma;
}


// ======================================================
// Leap
// ======================================================

// [[Rcpp::export]]
NumericMatrix Leap(const NumericMatrix& X,
                   const NumericVector& y,
                   NumericVector x,
                   NumericVector p,
                   int L,
                   double delta,
                   const NumericMatrix& Sigma_inv,
                   const NumericVector& dsq_eta,
                   double T,
                   double maxsig,
                   const NumericVector& lb,
                   const NumericVector& ub)
{
    const int n = X.nrow();
    const int d = X.ncol();

    NumericVector cx(n);
    NumericVector ans1(n);
    NumericVector dZold(d);
    NumericVector dZnew(d);

    double meanY = 0.0;
    double prod = 1.0;
    double h = 0.0;
    double h1 = 0.0;
    double sum = 0.0;
    double val = 0.0;
    double sigmaf_x = 0.0;

    int i, j, ii, jj;

    // ---------------------------------------
    // mean(y)
    // ---------------------------------------

    for(i = 0; i < n; i++)
        meanY += y[i];

    meanY /= (double)n;

    // ---------------------------------------
    // ans1 = Sigma_inv %*% (y - meanY)
    // ---------------------------------------

    for(i = 0; i < n; i++) {

        sum = 0.0;

        for(j = 0; j < n; j++)
            sum += Sigma_inv(i,j) * (y[j] - meanY);

        ans1[i] = sum;
    }

    // ---------------------------------------
    // dZold
    // ---------------------------------------

    for(j = 0; j < d; j++) {

        for(i = 0; i < n; i++) {

            prod = 1.0;

            h = std::fabs(x[j] - X(i,j));

            for(ii = 0; ii < d; ii++) {

                h1 = std::fabs(x[ii] - X(i,ii));

                prod *= (1.0 + dsq_eta[ii+1]*h1) *
                    std::exp(-dsq_eta[ii+1]*h1);
            }

            cx[i] =
                dsq_eta[0] *
                prod *
                dsq_eta[j+1] *
                dsq_eta[j+1] *
                (x[j] - X(i,j)) /
                (1.0 + dsq_eta[j+1]*h);
        }

        sum = 0.0;

        for(jj = 0; jj < n; jj++)
            sum += cx[jj] * ans1[jj];

        dZold[j] = sum / T;
    }

    // ---------------------------------------
    // Leapfrog iterations
    // ---------------------------------------

    for(int Lexp = 0; Lexp < L; Lexp++) {

        for(i = 0; i < d; i++)
            x[i] += delta * p[i]
                  - 0.5 * delta * delta * dZold[i];

        // dZnew

        for(i = 0; i < n; i++) {

            sum = 0.0;

            for(j = 0; j < n; j++)
                sum += Sigma_inv(i,j) * (y[j] - meanY);

            ans1[i] = sum;
        }

        for(j = 0; j < d; j++) {

            for(i = 0; i < n; i++) {

                prod = 1.0;

                h = std::fabs(x[j] - X(i,j));

                for(ii = 0; ii < d; ii++) {

                    h1 = std::fabs(x[ii] - X(i,ii));

                    prod *=
                        (1.0 + dsq_eta[ii+1]*h1) *
                        std::exp(-dsq_eta[ii+1]*h1);
                }

                cx[i] =
                    dsq_eta[0] *
                    prod *
                    dsq_eta[j+1] *
                    dsq_eta[j+1] *
                    (x[j] - X(i,j)) /
                    (1.0 + dsq_eta[j+1]*h);
            }

            sum = 0.0;

            for(jj = 0; jj < n; jj++)
                sum += cx[jj] * ans1[jj];

            dZnew[j] = sum / T;
        }

        for(i = 0; i < d; i++)
            p[i] -= 0.5 * delta *
                (dZold[i] + dZnew[i]);

        for(i = 0; i < d; i++)
            dZold[i] = dZnew[i];

        // sigmaf(x)

        for(i = 0; i < n; i++) {

            prod = 1.0;

            for(ii = 0; ii < d; ii++) {

                h = std::fabs(x[ii] - X(i,ii));

                prod *=
                    (1.0 + dsq_eta[ii+1]*h) *
                    std::exp(-dsq_eta[ii+1]*h);
            }

            cx[i] = dsq_eta[0] * prod;
        }

        for(i = 0; i < n; i++) {

            sum = 0.0;

            for(j = 0; j < n; j++)
                sum += Sigma_inv(i,j) * cx[j];

            ans1[i] = sum;
        }

        for(i = 0; i < n; i++) {

            prod = 1.0;

            for(ii = 0; ii < d; ii++) {

                h = std::fabs(x[ii] - X(i,ii));

                prod *=
                    (1.0 + dsq_eta[ii+1]*h) *
                    std::exp(-dsq_eta[ii+1]*h);
            }

            cx[i] = dsq_eta[0] * prod;
        }

        sum = 0.0;

        for(j = 0; j < n; j++)
            sum += cx[j] * ans1[j];

        val = dsq_eta[0] - sum;

        sigmaf_x = std::sqrt(val);

        if(sigmaf_x > maxsig) {
            sigmaf_x = -sigmaf_x;
            break;
        }
    }

    // ---------------------------------------
    // Reflective bounds
    // ---------------------------------------

    for(i = 0; i < d; i++) {

        bool outside;

        do {

            outside = false;

            if(x[i] < lb[i])
                x[i] = 2.0 * lb[i] - x[i];

            if(x[i] > ub[i])
                x[i] = 2.0 * ub[i] - x[i];

            if(x[i] < lb[i]) outside = true;
            if(x[i] > ub[i]) outside = true;

        } while(outside);
    }

    // ---------------------------------------
    // Return d x 3 matrix
    // ---------------------------------------

    NumericMatrix xpE(d,3);

    for(i = 0; i < d; i++) {
        xpE(i,0) = x[i];
        xpE(i,1) = p[i];
        xpE(i,2) = sigmaf_x;
    }

    return xpE;
}
