

/*
%----------------------------------------------------------------------------
%
% Purpose:
%   Numerical integration methods for ordinaray differential equations
%
%   This module provides implemenation of the variable order variable
%   stepsize multistep method of Shampine & Gordon.
%
% Last modified:   2015/08/25   M. Mahooti
%
% Reference:
%
%   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
%   Freeman and Comp., San Francisco (1975).
%
%----------------------------------------------------------------------------*/


#include "DEInteg.h"
#include "sign_.h"

// Inicialización de constantes estáticas
const std::vector<double> DEInteg::two = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0,
                                          256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};

const std::vector<double> DEInteg::gstr = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188,
                                           0.0143, 0.0114, 0.00936, 0.00789, 0.00679,
                                           0.00592, 0.00524, 0.00468};

// Constructor
DEInteg::DEInteg() : PermitTOUT(true), told(0), OldPermit(false), delsgn(0), start(false),
                     x(0), yy(Matrix(1,1)), h(0), hold(0), hnew(0), k(0), kold(0),
                     phase1(false), nornd(false), stiff(false), nostep(0), kle4(0),
                     ifail(0), crash(false), wt(Matrix(1,1)), p(Matrix(1,1)),
                     yp(Matrix(1,1)), phi(Matrix(1,17)), g(14, 0.0), sig(14, 0.0),
                     rho(14, 0.0), w(13, 0.0), alpha(13, 0.0), beta(13, 0.0),
                     v(13, 0.0), psi_(13, 0.0), State_(DE_INIT) {}

// Método principal de integración
Matrix DEInteg::integrate(std::function<Matrix(double, const Matrix&)> func,double t, double tout, double relerr, double abserr,Matrix& y) {

    int neq = y.getFilas();
    if (y.getColumnas()!=1) throw std::invalid_argument("state must be Nx1");
    int kp1,kp2,km1, km2,knew,ns;
    double erk,erkm1,absh;
    int n_eqn = y.getFilas();

    double twou = 2.0 * std::numeric_limits<double>::epsilon();
    double fouru = 4.0 * std::numeric_limits<double>::epsilon();

    // Redimensionar matrices según n_eqn
    yy = Matrix(n_eqn, 1);
    wt = Matrix(n_eqn, 1);
    p = Matrix(n_eqn, 1);
    yp = Matrix(n_eqn, 1);
    phi = Matrix(n_eqn, 17);

    // Return if output time equals input time
    if (t == tout) {
        return y;
    }

    // Test for improper parameters
    double epsilon = std::max(relerr, abserr);

    if ((relerr < 0.0) || (abserr < 0.0) || (epsilon <= 0.0) ||
        (State_ > DE_INVPARAM) || ((State_ != DE_INIT) && (t != told))) {
        State_ = DE_INVPARAM;
        return y;
    }

    // Set interval of integration and counter for number of steps
    double del = tout - t;
    double absdel = std::abs(del);

    double tend = t + 100.0 * del;
    if (!PermitTOUT) {
        tend = tout;
    }

    nostep = 0;
    kle4 = 0;
    stiff = false;
    double releps = relerr / epsilon;
    double abseps = abserr / epsilon;

    if ((State_ == DE_INIT) || (!OldPermit) || (delsgn * del <= 0.0)) {
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0, del);
        h = sign_(std::max( fouru * std::fabs(x),std::fabs(tout - x) ),tout - x);
    }



    while (true) { // Start step loop
        // If already past output point, interpolate solution and return
        if (std::abs(x - t) >= absdel) {
            Matrix yout(n_eqn, 1);
            Matrix ypout(n_eqn, 1);
            g[1] = 1.0;
            rho[1] = 1.0;
            double hi = tout - x;
            int ki = kold + 1;

            // Initialize w[*] for computing g[*]
            for (int i = 0; i < ki; i++) {
                double temp1 = i + 1;
                w[i] = 1.0 / temp1;
            }

            // Compute g[*]
            double term = 0.0;
            for (int j = 1; j < ki; j++) {
                double psijm1 = psi_[j];
                double gamma = (hi + term) / psijm1;
                double eta = hi / psijm1;
                for (int i = 0; i < ki - j; i++) {
                    w[i] = gamma * w[i] - eta * w[i + 1];
                }
                g[j + 1] = w[0];
                rho[j + 1] = gamma * rho[j];
                term = psijm1;
            }

            // Interpolate for the solution yout and for the derivative ypout
            for (int j = 0; j < ki; j++) {
                int i = ki - j;
                yout = yout + g[i] * phi.getColumn(i);
                ypout = ypout + rho[i] * phi.getColumn(i);
            }
            yout = y + hi * yout;
            State_ = DE_DONE;
            told = tout;
            OldPermit = PermitTOUT;
            return yout;
        }

        // If cannot go past output point and sufficiently close, extrapolate and return
        if (!PermitTOUT && (std::abs(tout - x) < fouru * std::abs(x))) {
            h = tout - x;
            yp = func(x, yy);
            Matrix result = yy + h * yp;
            State_ = DE_DONE;
            told = tout;
            OldPermit = PermitTOUT;
            return result;
        }

        // Limit step size, set weight vector and take a step
        h = sign_(std::min(std::abs(h), std::abs(tend - x)), h);
        for (int l = 1; l <= n_eqn; ++l) {
            wt(l, 1) = releps * std::abs(yy(l, 1)) + abseps;
        }

        // Begin block 0
        if (std::abs(h) < fouru * std::abs(x)) {
            h = sign_(fouru * std::abs(x), h);
            crash = true;
            State_ = DE_BADACC;
            relerr = epsilon * releps;
            abserr = epsilon * abseps;
            told = x;
            OldPermit = true;
            return yy;
        }

        double p5eps = 0.5 * epsilon;
        crash = false;
        g[1] = 1.0;
        g[2] = 0.5;
        sig[1] = 1.0;

        double round = 0.0;
        for (int l = 1; l <= n_eqn; ++l) {
            double yi = y(l, 1);
            double wti = wt(l, 1);
            round += (yi * yi) / (wti * wti);
        }
        round = twou * std::sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0 * round * (1.0 + fouru);
            crash = true;
            State_ = DE_BADACC;
            relerr = epsilon * releps;
            abserr = epsilon * abseps;
            told = x;
            OldPermit = true;
            return yy;
        }

        if (start) {
            // Initialize. Compute appropriate step size for first step.
            yp = func(x, y);
            double sum = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum += (yp(l, 1) * yp(l, 1))
                       / (wt(l, 1) * wt(l, 1));
            }

            sum = std::sqrt(sum);
            absh = std::abs(h);
            if (epsilon < 16.0 * sum * h * h) {
                absh = 0.25 * std::sqrt(epsilon / sum);
            }
            h = sign_(std::max(absh, fouru * std::abs(x)), h);
            hold = 0.0;
            hnew = 0.0;
            k = 0;
            kold = 0;
            start = false;
            phase1 = true;
            nornd = true;
            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l, 16) = 0.0;
                }

            }
        }
        // End block 0

        bool success = false;
        while (!success) {
            // Begin block 1
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            ns = 0;
            if (h != hold) {
                ns = 0;
            }
            if (ns <= kold) {
                ns = ns + 1;
            }
            int nsp1 = ns + 1;

            if (k >= ns) {
                // Compute components of alpha, beta, psi, sig which are changed
                beta[ns] = 1.0;
                double realns = ns;
                alpha[ns] = 1.0 / realns;
                double temp1 = h * realns;
                sig[nsp1] = 1.0;

                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; i++) {
                        int im1 = i - 1;
                        double temp2 = psi_[im1];
                        psi_[im1] = temp1;
                        beta[i] = beta[im1] * psi_[im1] / temp2;
                        temp1 = temp2 + h;
                        alpha[i] = h / temp1;
                        double reali = i;
                        sig[i + 1] = reali * alpha[i] * sig[i];
                    }
                }
                psi_[k] = temp1;

                // Compute coefficients g[*]; initialize v[*] and set w[*]
                if (ns > 0) {
                    // If order was raised, update diagonal part of v[*]
                    if (k > kold) {
                        double temp4 = k * kp1;
                        v[k] = 1.0 / temp4;
                        int nsm2 = ns - 2;
                        for (int j = 0; j < nsm2; j++) {
                            int i = k - j - 1;
                            v[i] = v[i] - alpha[j + 1] * v[i + 1];
                        }
                    }

                    // Update V[*] and set W[*]
                    int limit1 = kp1 - ns;
                    double temp5 = alpha[ns];
                    for (int iq = 0; iq < limit1; iq++) {
                        v[iq] = v[iq] - temp5 * v[iq + 1];
                        w[iq] = v[iq];
                    }
                    g[nsp1] = w[0];
                } else {
                    for (int iq = 0; iq < k; iq++) {
                        double temp3 = (iq + 1) * (iq + 2);
                        v[iq] = 1.0 / temp3;
                        w[iq] = v[iq];
                    }
                }

                // Compute the g[*] in the work vector w[*]
                int nsp2 = ns + 1;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; i++) {
                        int limit2 = kp2 - i;
                        double temp6 = alpha[i - 1];
                        for (int iq = 0; iq < limit2; iq++) {
                            w[iq] = w[iq] - temp6 * w[iq + 1];
                        }
                        g[i] = w[0];
                    }
                }
            } // if K>=NS
            // End block 1

            // Begin block 2
            // Predict a solution, evaluate derivatives using predicted solution
            if (k >= nsp1) {
                for (int i = nsp1; i <= k; i++) {
                    double temp1 = beta[i];
                    for (int l = 1; l <= n_eqn; ++l) {
                        phi(l, i+1) = temp1 * phi(l, i+1);
                    }
                }
            }

            // Predict solution and differences
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, kp1+1) = phi(l, k+1);
                phi(l, k+1)   = 0.0;
                p(l, 1)       = 0.0;
            }


            for (int j = 0; j < k; ++j) {
                int i = kp1 - j - 1;
                double temp2 = g[i];
                for (int l = 1; l <= n_eqn; ++l) {
                    p(l, 1)    = p(l, 1) + temp2 * phi(l, i+1);
                    phi(l, i+1)= phi(l, i+1) + phi(l, i+2);
                }
            }


            if (nornd) {
                p = y + h * p;
            } else {
                for (int l = 1; l <= n_eqn; ++l) {
                    double tau   = h * p(l, 1) - phi(l, 16);
                    p(l, 1)      = y(l, 1) + tau;
                    phi(l, 17)   = (p(l, 1) - y(l, 1)) - tau;
                }

            }

            double xold = x;
            x = x + h;
            double absh = std::abs(h);
            yp = func(x, p);

            // Estimate errors at orders k, k-1, k-2
            double erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;

            for (int l = 1; l <= n_eqn; ++l) {
                double temp3 = 1.0 / wt(l, 1);
                double temp4 = yp(l, 1) - phi(l, 2);
                if (km2 >= 0) {
                    erkm2 += ((phi(l, km1+1) + temp4) * temp3) * ((phi(l, km1+1) + temp4) * temp3);
                }
                if (km1 >= 0) {
                    erkm1 += ((phi(l, k+1)   + temp4) * temp3) * ((phi(l, k+1)   + temp4) * temp3);
                }
                erk   += (temp4 * temp3) * (temp4 * temp3);
            }


            if (km2 >= 0) {
                erkm2 = absh * sig[km1] * gstr[km2] * std::sqrt(erkm2);
            }
            if (km1 >= 0) {
                erkm1 = absh * sig[k] * gstr[km1] * std::sqrt(erkm1);
            }

            double temp5 = absh * std::sqrt(erk);
            double err = temp5 * (g[k] - g[kp1]);
            erk = temp5 * sig[kp1] * gstr[k];
            knew = k;

            // Test if order should be lowered
            if (km2 >= 0) {
                if (std::max(erkm1, erkm2) <= erk) {
                    knew = km1;
                }
            }
            if (km2 == -1) {
                if (erkm1 <= 0.5 * erk) {
                    knew = km1;
                }
            }
            // End block 2

            success = (err <= epsilon);

            if (!success) {
                // Begin block 3
                // The step is unsuccessful. Restore x, phi, psi.
                phase1 = false;
                x = xold;
                for (int i = 0; i < k; ++i) {
                    double temp1 = 1.0 / beta[i + 1];  // stays zero-based for vector
                    for (int l = 1; l <= n_eqn; ++l) {
                        phi(l, i+2) = temp1 * (phi(l, i+2) - phi(l, i+3));
                    }
                }


                if (k >= 1) {
                    for (int i = 1; i <= k; i++) {
                        psi_[i - 1] = psi_[i] - h;
                    }
                }

                // On third failure, set order to one.
                // Thereafter, use optimal step size
                ifail = ifail + 1;
                double temp2 = 0.5;
                if (ifail > 2) {
                    if (p5eps < 0.25 * erk) {
                        temp2 = std::sqrt(p5eps / erk);
                    }
                }
                if (ifail >= 2) {
                    knew = 0;
                }
                h = temp2 * h;
                k = knew;
                if (std::abs(h) < fouru * std::abs(x)) {
                    crash = true;
                    h = sign_(fouru * std::abs(x), h);
                    epsilon = epsilon * 2.0;
                    State_ = DE_BADACC;
                    relerr = epsilon * releps;
                    abserr = epsilon * abseps;
                    told = x;
                    OldPermit = true;
                    return yy;
                }
                // End block 3
            }
        } // while (!success)

        // Begin block 4
        // The step is successful. Correct the predicted solution
        kold = k;
        hold = h;

        // Correct and evaluate
        double temp1 = h * g[kp1];
        if (nornd) {
            for (int l = 1; l <= n_eqn; ++l) {
                y(l, 1) = p(l, 1) + temp1 * (yp(l, 1) - phi(l, 2));
            }

        } else {
            for (int l = 1; l <= n_eqn; ++l) {
                double rho_val = temp1 * (yp(l, 1) - phi(l, 2)) - phi(l, 17);
                y(l, 1)        = p(l, 1) + rho_val;
                phi(l, 16)     = (y(l, 1) - p(l, 1)) - rho_val;
            }

        }
        yp = func(x, y);

        // Update differences for next step
        for (int l = 1; l <= n_eqn; ++l) {
            phi(l, kp1+1) = yp(l, 1) - phi(l, 2);
            phi(l, kp2+1) = phi(l, kp1+1) - phi(l, kp2+1);
        }

        for (int i = 0; i < k; ++i) {
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, i+2) = phi(l, i+2) + phi(l, kp1+1);
            }
        }


        // Estimate error at order k+1 unless:
        // - in first phase when always raise order
        // - already decided to lower order
        // - step size not constant so estimate unreliable
        double erkp1 = 0.0;
        if ((knew == km1) || (k == 11)) {
            phase1 = false;
        }

        if (phase1) {
            k = kp1;
            erk = erkp1;
        } else {
            if (knew == km1) {
                // lower order
                k = km1;
                erk = erkm1;
            } else {
                if (kp1 <= ns) {
                    for (int l = 1; l <= n_eqn; ++l) {
                        double tmp = phi(l, kp2+1) / wt(l, 1);
                        erkp1 += tmp * tmp;
                    }

                    erkp1 = absh * gstr[kp1] * std::sqrt(erkp1);

                    // Using estimated error at order k+1, determine appropriate order
                    if (k > 0) {
                        if (erkm1 <= std::min(erk, erkp1)) {
                            // lower order
                            k = km1;
                            erk = erkm1;
                        } else {
                            if ((erkp1 < erk) && (k != 11)) {
                                // raise order
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    } else if (erkp1 < 0.5 * erk) {
                        // raise order
                        k = kp1;
                        erk = erkp1;
                    }
                }
            }
        }

        // With new order determine appropriate step size for next step
        if (phase1 || (p5eps >= erk * two[k + 1])) {
            hnew = 2.0 * h;
        } else {
            if (p5eps < erk) {
                double temp2 = k + 1;
                double r = std::pow(p5eps / erk, 1.0 / temp2);
                hnew = absh * std::max(0.5, std::min(0.9, r));
                hnew = sign_(std::max(hnew, fouru * std::abs(x)), h);
            } else {
                hnew = h;
            }
        }
        h = hnew;
        // End block 4

        // Test for too small tolerances
        if (crash) {
            State_ = DE_BADACC;
            relerr = epsilon * releps;
            abserr = epsilon * abseps;
            told = x;
            OldPermit = true;
            return yy;
        }

        nostep = nostep + 1;

        // Count number of consecutive steps taken with order <= 4
        kle4 = kle4 + 1;
        if (kold > 3) {
            kle4 = 0;
        }
        if (kle4 >= 49) {
            stiff = true;
        }
    } // End step loop
}