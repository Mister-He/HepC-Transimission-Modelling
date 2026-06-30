// =============================================================================
// sim.cpp — HCV Compartmental Model among PWID in Singapore
// =============================================================================
// Rcpp interface for a continuous-time ODE system.
// Compile via:  Rcpp::sourceCpp("sim.cpp")
//
// Compartment index convention
// ─────────────────────────────
//   Stratum  s ∈ {D=0, J=1, F=2, X=3}   (D=never, J=current, F=ever, X=former)
//   Stage    k ∈ {1,2,3,4}               (NC, CC, DC, HCC)
//   State    h ∈ {u=0, a=1, c=2, t=3}   (susceptible/post-SVR, acute, chronic, treated)
//   Age grp  i ∈ {0,...,9}              (10 age groups)
//
// Flat index: idx(s,k,h,i) = s*4*4*10 + (k-1)*4*10 + h*10 + i
// Total compartments: 4*4*4*10 = 640
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <stdexcept>
using namespace Rcpp;

// ─── helper: flat compartment index ─────────────────────────────────────────
// s in {0,1,2,3}, k in {1,2,3,4}, h in {0,1,2,3}, i in {0,...,9}
inline int idx(int s, int k, int h, int i) {
    return s * 4 * 4 * 10 + (k - 1) * 4 * 10 + h * 10 + i;
}

// =============================================================================
// PARAMETER STRUCT
// =============================================================================
struct Params {

    // --- Transmission & natural history ------------------------------------
    double q;       // transmission prob per needle-sharing event
    double kappa;   // prob spontaneous clearance after acute infection
    double iota1;   // mean duration of acute infection (years)
    double iota2;   // mean duration of DAA treatment   (years)

    // --- Genotype ----------------------------------------------------------
    double rho;     // proportion with genotype 3

    // --- SVR12 rates (effective, genotype-weighted) -------------------------
    // Stages 1-2 (NC and CC), RBV-free
    double alpha_NC;     // = rho*alpha3_NC + (1-rho)*alphao_NC
    // Stage 3 (DC) with RBV+ (used for GT3 DC patients)
    double alpha_DC_pos; // = rho*alpha3_DC_pos + (1-rho)*alphao_DC_pos
    // Stage 3 (DC) RBV-free
    double alpha_DC_neg; // = rho*alpha3_DC_neg + (1-rho)*alphao_DC_neg
    // Stage 4 (HCC) — assumed same as DC+ (conservative)
    double alpha_HCC;    // = alpha_DC_pos  (assumption; see Section X)

    // --- Treatment initiation rates (per scenario) -------------------------
    // tau[k-1] for stage k = 1,2,3,4
    std::vector<double> tau; // length 4

    // --- Baseline progression rates (other GT, per year) -------------------
    double p_NC_CC;   // p3^{NC→CC}
    double p_CC_DC;   // p3^{CC→DC}
    double p_CC_HCC;  // p3^{CC→HCC}
    double p_DC_HCC;  // p3^{DC→HCC}

    // --- GT3 relative risks ------------------------------------------------
    double r3_NC_CC;   // r3^{NC→CC}  = 1.36
    double r3_CC_DC;   // r3^{CC→DC}  = 1.36
    double r3_CC_HCC;  // r3^{CC→HCC} = 1.93
    double r3_DC_HCC;  // r3^{DC→HCC} = 1.93

    // --- Effective progression rates (computed from above) -----------------
    // p_tilde_c_{A→B} = rho*r3*p + (1-rho)*p
    double ptc_NC_CC;
    double ptc_CC_DC;
    double ptc_CC_HCC;
    double ptc_DC_HCC;

    // --- Post-SVR progression modifiers ------------------------------------
    double phi_CC_DC;   // phi^{CC→DC}  = 0.07
    double phi_CC_HCC;  // phi^{CC→HCC} = 0.23
    double phi_DC_HCC;  // phi^{DC→HCC} = 1.00 (assumed no reduction)

    // --- Effective post-SVR progression rates (computed) -------------------
    double ptu_CC_DC;   // = phi_CC_DC  * ptc_CC_DC
    double ptu_CC_HCC;  // = phi_CC_HCC * ptc_CC_HCC
    double ptu_DC_HCC;  // = phi_DC_HCC * ptc_DC_HCC

    // --- Background mortality (age-varying, length 10) ----------------------
    std::vector<double> mu;   // mu[i] for age group i

    // --- Standardized mortality ratio for ever-PWIDs -----------------------
    double omega;  // SMR for ever-PWIDs (Degenhardt et al. 2011)

    // --- Disease-specific excess mortality ---------------------------------
    double mu_DC;   // additional mortality in decompensated cirrhosis
    double mu_HCC;  // additional mortality in HCC

    // --- SVR mortality modifiers -------------------------------------------
    double psi_DC;  // relative mortality risk in DC after SVR  = 0.45
    double psi_HCC; // relative mortality risk in HCC after SVR = 0.37

    // --- Incarceration rates (age-varying, length 10 each) ------------------
    // lambda1[i] is pre-computed in R as lambda3[i] * c_true[i] before each call
    std::vector<double> lambda1; // effective first-arrest rate lambda_i^(1)
    std::vector<double> lambda2; // release rate        lambda_i^(2)
    std::vector<double> lambda3; // re-arrest rate      lambda_i^(3)
    double pi_recid;             // recidivism probability pi

    // --- Needle-sharing contact rate ---------------------------------------
    // GUESS: scalar homogeneous mixing; replace with 10×10 matrix if needed
    arma::mat C_contact; // 10×10 age-structured contact matrix, C_contact(i,j)

    // --- Population entry (age-varying, time-constant approximation) -------
    // GUESS: constant over time; ideally beta_i(t) from calibration
    std::vector<double> beta; // entry rate into D_{u,1,i}  (GUESS)
};

// =============================================================================
// HELPER: compute force of infection gamma_{i,j}(t)
// Assumes proportionate (homogeneous) mixing across strata and age groups.
// gamma_i(t) = q * c_contact * (total infectious / total active PWID)
// Infectious = D_a + F_a + D_c + F_c (acute, active strata only; J treated as no sharing)
// =============================================================================
double forceOfInfection(int i, const std::vector<double>& y, const Params& p) {

    double lambda_i = 0.0;

    for (int j = 0; j < 10; ++j) {

        // infectious (acute and chronic) in strata D and F for age group j
        double infectious_j = y[idx(0,1,1,j)] + y[idx(0,1,2,j)]
                            + y[idx(0,2,1,j)] + y[idx(0,2,2,j)]
                            + y[idx(0,3,1,j)] + y[idx(0,3,2,j)]
                            + y[idx(0,4,1,j)] + y[idx(0,4,2,j)]
                            + y[idx(2,1,1,j)] + y[idx(2,1,2,j)]
                            + y[idx(2,2,1,j)] + y[idx(2,2,2,j)]
                            + y[idx(2,3,1,j)] + y[idx(2,3,2,j)]
                            + y[idx(2,4,1,j)] + y[idx(2,4,2,j)];

        // active PWID in strata D and F for age group j
        double active_j = 0.0;
        for (int k = 1; k <= 4; ++k)
            for (int h = 0; h <= 3; ++h) {
                active_j += y[idx(0,k,h,j)];
                active_j += y[idx(2,k,h,j)];
            }

        if (active_j <= 0.0) continue;

        lambda_i += p.C_contact(i, j) * (infectious_j / active_j);
    }

    return p.q * lambda_i;
}

// =============================================================================
// ODE RIGHT-HAND SIDE
// dydt has 640 entries matching the flat index convention above.
// =============================================================================
std::vector<double> rhs(double t,
                        const std::vector<double>& y,
                        const Params& p) {

    std::vector<double> dydt(640, 0.0);

    // SVR12 rates by stage
    // Stage 1,2: alpha_NC; Stage 3: alpha_DC_pos; Stage 4: alpha_HCC
    auto svr = [&](int k) -> double {
        if (k <= 2) return p.alpha_NC;
        if (k == 3) return p.alpha_DC_pos;
        return p.alpha_HCC;
    };

    // Effective mortality at stage k for state h
    // h=3 (treated): apply SVR mortality modifier for stages 3,4
    auto mu_eff = [&](int k, int h, int i) -> double {
        double base = p.mu[i] * p.omega; 
        if (k == 3) {
            double adjusted = (h == 3) ? p.psi_DC * (p.mu_DC + base) : (p.mu_DC + base);
            return adjusted;
        }
        if (k == 4) {
            double adjusted = (h == 3) ? p.psi_HCC * (p.mu_HCC + base) : (p.mu_HCC + base);
            return adjusted;
        }
        return base;
    };

    // Post-SVR progression rates (applies to u compartments)
    auto ptu = [&](int from_k, int to_k) -> double {
        if (from_k == 2 && to_k == 3) return p.ptu_CC_DC;
        if (from_k == 2 && to_k == 4) return p.ptu_CC_HCC;
        if (from_k == 3 && to_k == 4) return p.ptu_DC_HCC;
        return 0.0;
    };

    // Chronic progression rates (applies to c compartments)
    auto ptc = [&](int from_k, int to_k) -> double {
        if (from_k == 1 && to_k == 2) return p.ptc_NC_CC;
        if (from_k == 2 && to_k == 3) return p.ptc_CC_DC;
        if (from_k == 2 && to_k == 4) return p.ptc_CC_HCC;
        if (from_k == 3 && to_k == 4) return p.ptc_DC_HCC;
        return 0.0;
    };

    // =========================================================================
    //  LOOP OVER AGE GROUPS
    // =========================================================================
    for (int i = 0; i < 10; ++i) {

        double gam  = forceOfInfection(i, y, p);  // gamma_{i,j}(t)
        double mu_i = p.mu[i] * p.omega;          // mortality for PWIDs in age group i
        double l1   = p.lambda1[i];               // pre-computed as lambda3[i]*c_true[i] in R
        double l2   = p.lambda2[i];
        double l3   = p.lambda3[i];
        double pi   = p.pi_recid;

        // =====================================================================
        //  STRATUM D — Never incarcerated
        // =====================================================================
        // --- Stage 1: Non-cirrhosis ------------------------------------------
        {
            double Du1 = y[idx(0,1,0,i)];
            double Da1 = y[idx(0,1,1,i)];
            double Dc1 = y[idx(0,1,2,i)];
            double Dt1 = y[idx(0,1,3,i)];

            double a = svr(1);

            // D_{u,1,i}
            dydt[idx(0,1,0,i)] +=
                - (gam + l1 + mu_i) * Du1
                + (p.kappa / p.iota1) * Da1
                + (a / p.iota2) * Dt1
                + p.beta[i];

            // D_{a,1,i}
            dydt[idx(0,1,1,i)] +=
                  gam * Du1
                - (1.0/p.iota1 + l1 + mu_i) * Da1;

            // D_{c,1,i}
            dydt[idx(0,1,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Da1
                + ((1.0-a) / p.iota2) * Dt1
                - (p.tau[0] + l1 + mu_i + ptc(1,2)) * Dc1;

            // D_{t,1,i}
            dydt[idx(0,1,3,i)] +=
            p.tau[0] * Dc1
            - (1.0/p.iota2 + l1 + mu_eff(1,3,i)) * Dt1;
        }

        // --- Stage 2: Compensated cirrhosis ----------------------------------
        {
            double Du2 = y[idx(0,2,0,i)];
            double Da2 = y[idx(0,2,1,i)];
            double Dc2 = y[idx(0,2,2,i)];
            double Dt2 = y[idx(0,2,3,i)];
            double Dc1 = y[idx(0,1,2,i)];

            double a = svr(2);

            dydt[idx(0,2,0,i)] +=
                - (gam + l1 + mu_i + ptu(2,3) + ptu(2,4)) * Du2
                + (p.kappa / p.iota1) * Da2
                + (a / p.iota2) * Dt2;

            dydt[idx(0,2,1,i)] +=
                  gam * Du2
                - (1.0/p.iota1 + l1 + mu_i) * Da2;

            dydt[idx(0,2,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Da2
                + ((1.0-a) / p.iota2) * Dt2
                + ptc(1,2) * Dc1
                - (p.tau[1] + l1 + mu_i + ptc(2,3) + ptc(2,4)) * Dc2;

            dydt[idx(0,2,3,i)] +=
                  p.tau[1] * Dc2
                - (1.0/p.iota2 + l1 + mu_eff(2,3,i)) * Dt2;
        }

        // --- Stage 3: Decompensated cirrhosis --------------------------------
        {
            double Du3 = y[idx(0,3,0,i)];
            double Da3 = y[idx(0,3,1,i)];
            double Dc3 = y[idx(0,3,2,i)];
            double Dt3 = y[idx(0,3,3,i)];
            double Du2 = y[idx(0,2,0,i)];
            double Dc2 = y[idx(0,2,2,i)];

            double a = svr(3);

            dydt[idx(0,3,0,i)] +=
                - (gam + l1 + mu_eff(3,0,i) + ptu(3,4)) * Du3
                + (p.kappa / p.iota1) * Da3
                + (a / p.iota2) * Dt3
                + ptu(2,3) * Du2;

            dydt[idx(0,3,1,i)] +=
                  gam * Du3
                - (1.0/p.iota1 + l1 + mu_eff(3,1,i)) * Da3;

            dydt[idx(0,3,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Da3
                + ((1.0-a) / p.iota2) * Dt3
                + ptc(2,3) * Dc2
                - (p.tau[2] + l1 + mu_eff(3,2,i) + ptc(3,4)) * Dc3;

            dydt[idx(0,3,3,i)] +=
                  p.tau[2] * Dc3
                - (1.0/p.iota2 + l1 + mu_eff(3,3,i)) * Dt3;
        }

        // --- Stage 4: HCC ----------------------------------------------------
        {
            double Du4 = y[idx(0,4,0,i)];
            double Da4 = y[idx(0,4,1,i)];
            double Dc4 = y[idx(0,4,2,i)];
            double Dt4 = y[idx(0,4,3,i)];
            double Du2 = y[idx(0,2,0,i)];
            double Du3 = y[idx(0,3,0,i)];
            double Dc2 = y[idx(0,2,2,i)];
            double Dc3 = y[idx(0,3,2,i)];

            double a = svr(4);

            dydt[idx(0,4,0,i)] +=
                - (gam + l1 + mu_eff(4,0,i)) * Du4
                + (p.kappa / p.iota1) * Da4
                + (a / p.iota2) * Dt4
                + ptu(2,4) * Du2
                + ptu(3,4) * Du3;

            dydt[idx(0,4,1,i)] +=
                  gam * Du4
                - (1.0/p.iota1 + l1 + mu_eff(4,1,i)) * Da4;

            dydt[idx(0,4,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Da4
                + ((1.0-a) / p.iota2) * Dt4
                + ptc(2,4) * Dc2
                + ptc(3,4) * Dc3
                - (p.tau[3] + l1 + mu_eff(4,2,i)) * Dc4;

            dydt[idx(0,4,3,i)] +=
                  p.tau[3] * Dc4
                - (1.0/p.iota2 + l1 + mu_eff(4,3,i)) * Dt4;
        }

        // =====================================================================
        //  STRATUM J — Currently incarcerated
        //  Receives from D (first arrest) and F (re-arrest)
        //  Releases to F (pi) and X (1-pi)
        // =====================================================================
        for (int k = 1; k <= 4; ++k) {
            for (int h = 0; h <= 3; ++h) {
                double Jkh = y[idx(1,k,h,i)];
                double Dkh = y[idx(0,k,h,i)];
                double Fkh = y[idx(2,k,h,i)];

                // inflow from D (first arrest) and F (re-arrest)
                double inflow = l1 * Dkh + l3 * Fkh;

                // net outflow rate (release only; no force of infection in J)
                double out = l2 + mu_eff(k,h,i);

                dydt[idx(1,k,h,i)] += inflow - out * Jkh;

                // internal disease dynamics within J (same structure as D,
                // but no gamma infection and no lambda1 arrest term)
                // These are added below per-state
            }
        }

        // J internal disease transitions (mirrors D block minus gamma/lambda1)
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);

            double Ju = y[idx(1,k,0,i)];
            double Ja = y[idx(1,k,1,i)];
            double Jc = y[idx(1,k,2,i)];
            double Jt = y[idx(1,k,3,i)];

            // u: SVR return inflow, no infection in prison
            dydt[idx(1,k,0,i)] +=
                  (p.kappa / p.iota1) * Ja
                + (a / p.iota2) * Jt;

            // a: no new infection in J
            dydt[idx(1,k,1,i)] +=
                - (1.0/p.iota1) * Ja;  // recovery/chronification handled below

            // c: treatment initiation and disease entry
            dydt[idx(1,k,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Ja
                + ((1.0-a) / p.iota2) * Jt
                - p.tau[k-1] * Jc;

            // t: treatment
            dydt[idx(1,k,3,i)] +=
                  p.tau[k-1] * Jc
                - (1.0/p.iota2) * Jt;

            // progression between stages within J
            if (k < 4) {
                // u→ progression (post-SVR)
                if (k == 2) {
                    dydt[idx(1,2,0,i)] -= (ptu(2,3) + ptu(2,4)) * Ju;
                    dydt[idx(1,3,0,i)] += ptu(2,3) * Ju;
                    dydt[idx(1,4,0,i)] += ptu(2,4) * Ju;
                }
                if (k == 3) {
                    dydt[idx(1,3,0,i)] -= ptu(3,4) * Ju;
                    dydt[idx(1,4,0,i)] += ptu(3,4) * Ju;
                }
                // c→ progression (chronic)
                if (k == 1) {
                    dydt[idx(1,1,2,i)] -= ptc(1,2) * Jc;
                    dydt[idx(1,2,2,i)] += ptc(1,2) * Jc;
                }
                if (k == 2) {
                    dydt[idx(1,2,2,i)] -= (ptc(2,3) + ptc(2,4)) * Jc;
                    dydt[idx(1,3,2,i)] += ptc(2,3) * Jc;
                    dydt[idx(1,4,2,i)] += ptc(2,4) * Jc;
                }
                if (k == 3) {
                    dydt[idx(1,3,2,i)] -= ptc(3,4) * Jc;
                    dydt[idx(1,4,2,i)] += ptc(3,4) * Jc;
                }
            }
        }

        // =====================================================================
        //  STRATUM F — Ever-incarcerated (active PWID, at risk of infection)
        // =====================================================================
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);

            double Fu = y[idx(2,k,0,i)];
            double Fa = y[idx(2,k,1,i)];
            double Fc = y[idx(2,k,2,i)];
            double Ft = y[idx(2,k,3,i)];
            double Ju = y[idx(1,k,0,i)];
            double Ja = y[idx(1,k,1,i)];
            double Jc = y[idx(1,k,2,i)];
            double Jt = y[idx(1,k,3,i)];

            // F_u: released from J with recidivism, minus infection and re-arrest
            dydt[idx(2,k,0,i)] +=
                  pi * l2 * Ju
                + (p.kappa / p.iota1) * Fa
                + (a / p.iota2) * Ft
                - (gam + l3 + mu_eff(k,0,i)) * Fu;

            dydt[idx(2,k,1,i)] +=
                  pi * l2 * Ja
                + gam * Fu
                - (1.0/p.iota1 + l3 + mu_eff(k,1,i)) * Fa;

            dydt[idx(2,k,2,i)] +=
                  pi * l2 * Jc
                + ((1.0-p.kappa) / p.iota1) * Fa
                + ((1.0-a) / p.iota2) * Ft
                - (p.tau[k-1] + l3 + mu_eff(k,2,i)) * Fc;

            dydt[idx(2,k,3,i)] +=
                  pi * l2 * Jt
                + p.tau[k-1] * Fc
                - (1.0/p.iota2 + l3 + mu_eff(k,3,i)) * Ft;
        }

        // F progression between stages
        for (int k = 1; k <= 3; ++k) {
            double Fu = y[idx(2,k,0,i)];
            double Fc = y[idx(2,k,2,i)];

            if (k == 1) {
                dydt[idx(2,1,2,i)] -= ptc(1,2) * Fc;
                dydt[idx(2,2,2,i)] += ptc(1,2) * Fc;
            }
            if (k == 2) {
                dydt[idx(2,2,0,i)] -= (ptu(2,3) + ptu(2,4)) * Fu;
                dydt[idx(2,3,0,i)] += ptu(2,3) * Fu;
                dydt[idx(2,4,0,i)] += ptu(2,4) * Fu;

                dydt[idx(2,2,2,i)] -= (ptc(2,3) + ptc(2,4)) * Fc;
                dydt[idx(2,3,2,i)] += ptc(2,3) * Fc;
                dydt[idx(2,4,2,i)] += ptc(2,4) * Fc;
            }
            if (k == 3) {
                dydt[idx(2,3,0,i)] -= ptu(3,4) * Fu;
                dydt[idx(2,4,0,i)] += ptu(3,4) * Fu;

                dydt[idx(2,3,2,i)] -= ptc(3,4) * Fc;
                dydt[idx(2,4,2,i)] += ptc(3,4) * Fc;
            }
        }

        // =====================================================================
        //  STRATUM X — Former PWID (no infection risk, no re-arrest)
        // =====================================================================
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);

            double Xu = y[idx(3,k,0,i)];
            double Xa = y[idx(3,k,1,i)];
            double Xc = y[idx(3,k,2,i)];
            double Xt = y[idx(3,k,3,i)];
            double Ju = y[idx(1,k,0,i)];
            double Ja = y[idx(1,k,1,i)];
            double Jc = y[idx(1,k,2,i)];
            double Jt = y[idx(1,k,3,i)];

            // inflow from J on release (1-pi fraction)
            dydt[idx(3,k,0,i)] +=
                  (1.0-pi) * l2 * Ju
                + (p.kappa / p.iota1) * Xa
                + (a / p.iota2) * Xt
                - mu_eff(k,0,i) * Xu;

            dydt[idx(3,k,1,i)] +=
                  (1.0-pi) * l2 * Ja
                - (1.0/p.iota1 + mu_eff(k,1,i)) * Xa;

            dydt[idx(3,k,2,i)] +=
                  (1.0-pi) * l2 * Jc
                + ((1.0-p.kappa) / p.iota1) * Xa
                + ((1.0-a) / p.iota2) * Xt
                - (p.tau[k-1] + mu_eff(k,2,i)) * Xc;

            dydt[idx(3,k,3,i)] +=
                  (1.0-pi) * l2 * Jt
                + p.tau[k-1] * Xc
                - (1.0/p.iota2 + mu_eff(k,3,i)) * Xt;
        }

        // X progression between stages
        for (int k = 1; k <= 3; ++k) {
            double Xu = y[idx(3,k,0,i)];
            double Xc = y[idx(3,k,2,i)];

            if (k == 1) {
                dydt[idx(3,1,2,i)] -= ptc(1,2) * Xc;
                dydt[idx(3,2,2,i)] += ptc(1,2) * Xc;
            }
            if (k == 2) {
                dydt[idx(3,2,0,i)] -= (ptu(2,3) + ptu(2,4)) * Xu;
                dydt[idx(3,3,0,i)] += ptu(2,3) * Xu;
                dydt[idx(3,4,0,i)] += ptu(2,4) * Xu;

                dydt[idx(3,2,2,i)] -= (ptc(2,3) + ptc(2,4)) * Xc;
                dydt[idx(3,3,2,i)] += ptc(2,3) * Xc;
                dydt[idx(3,4,2,i)] += ptc(2,4) * Xc;
            }
            if (k == 3) {
                dydt[idx(3,3,0,i)] -= ptu(3,4) * Xu;
                dydt[idx(3,4,0,i)] += ptu(3,4) * Xu;

                dydt[idx(3,3,2,i)] -= ptc(3,4) * Xc;
                dydt[idx(3,4,2,i)] += ptc(3,4) * Xc;
            }
        }

    } // end age group loop

    return dydt;
}

// =============================================================================
// EULER INTEGRATOR (replace with RK4 for production runs)
// For calibration / exploratory use; switch to RK4 below if needed.
// =============================================================================
std::vector<double> euler_step(double t, double dt,
                               const std::vector<double>& y,
                               const Params& p) {
    auto dy = rhs(t, y, p);
    std::vector<double> y_new(640);
    for (int i = 0; i < 640; ++i)
        y_new[i] = y[i] + dt * dy[i];
    return y_new;
}

// RK4 integrator
std::vector<double> rk4_step(double t, double dt,
                              const std::vector<double>& y,
                              const Params& p) {
    auto k1 = rhs(t,        y,                                  p);
    std::vector<double> y2(640), y3(640), y4(640);
    for (int i = 0; i < 640; ++i) y2[i] = y[i] + 0.5*dt*k1[i];
    auto k2 = rhs(t+0.5*dt, y2,                                 p);
    for (int i = 0; i < 640; ++i) y3[i] = y[i] + 0.5*dt*k2[i];
    auto k3 = rhs(t+0.5*dt, y3,                                 p);
    for (int i = 0; i < 640; ++i) y4[i] = y[i] + dt*k3[i];
    auto k4 = rhs(t+dt,     y4,                                 p);

    std::vector<double> y_new(640);
    for (int i = 0; i < 640; ++i)
        y_new[i] = y[i] + (dt/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    return y_new;
}

// =============================================================================
// RCPP INTERFACE
// =============================================================================

// [[Rcpp::export]]
NumericMatrix run_sim(List params_r, List data_r) {

    // ─── unpack params ────────────────────────────────────────────────────
    Params p;

    p.q       = as<double>(params_r["q"]);
    p.kappa   = as<double>(params_r["kappa"]);
    p.iota1   = as<double>(params_r["iota1"]);
    p.iota2   = as<double>(params_r["iota2"]);
    p.rho     = as<double>(params_r["rho"]);

    // SVR12 rates (effective, pre-weighted in R before passing)
    p.alpha_NC     = as<double>(params_r["alpha_NC"]);
    p.alpha_DC_pos = as<double>(params_r["alpha_DC_pos"]);
    p.alpha_DC_neg = as<double>(params_r["alpha_DC_neg"]);
    p.alpha_HCC    = as<double>(params_r["alpha_HCC"]);

    // tau: treatment initiation rates (length-4 vector, one per stage)
    NumericVector tau_r = as<NumericVector>(params_r["tau"]);
    p.tau.assign(tau_r.begin(), tau_r.end());

    // baseline progression
    p.p_NC_CC   = as<double>(params_r["p_NC_CC"]);
    p.p_CC_DC   = as<double>(params_r["p_CC_DC"]);
    p.p_CC_HCC  = as<double>(params_r["p_CC_HCC"]);
    p.p_DC_HCC  = as<double>(params_r["p_DC_HCC"]);

    // GT3 relative risks
    p.r3_NC_CC   = as<double>(params_r["r3_NC_CC"]);
    p.r3_CC_DC   = as<double>(params_r["r3_CC_DC"]);
    p.r3_CC_HCC  = as<double>(params_r["r3_CC_HCC"]);
    p.r3_DC_HCC  = as<double>(params_r["r3_DC_HCC"]);

    // compute effective (genotype-weighted) progression rates
    p.ptc_NC_CC   = (p.rho * p.r3_NC_CC  + (1.0-p.rho)) * p.p_NC_CC;
    p.ptc_CC_DC   = (p.rho * p.r3_CC_DC  + (1.0-p.rho)) * p.p_CC_DC;
    p.ptc_CC_HCC  = (p.rho * p.r3_CC_HCC + (1.0-p.rho)) * p.p_CC_HCC;
    p.ptc_DC_HCC  = (p.rho * p.r3_DC_HCC + (1.0-p.rho)) * p.p_DC_HCC;

    // post-SVR modifiers
    p.phi_CC_DC   = as<double>(params_r["phi_CC_DC"]);
    p.phi_CC_HCC  = as<double>(params_r["phi_CC_HCC"]);
    p.phi_DC_HCC  = as<double>(params_r["phi_DC_HCC"]);

    // effective post-SVR progression
    p.ptu_CC_DC   = p.phi_CC_DC  * p.ptc_CC_DC;
    p.ptu_CC_HCC  = p.phi_CC_HCC * p.ptc_CC_HCC;
    p.ptu_DC_HCC  = p.phi_DC_HCC * p.ptc_DC_HCC;

    // mortality
    NumericVector mu_r = as<NumericVector>(params_r["mu"]);
    p.mu.assign(mu_r.begin(), mu_r.end());
    p.omega   = as<double>(params_r["omega"]);
    p.mu_DC   = as<double>(params_r["mu_DC"]);
    p.mu_HCC  = as<double>(params_r["mu_HCC"]);
    p.psi_DC  = as<double>(params_r["psi_DC"]);
    p.psi_HCC = as<double>(params_r["psi_HCC"]);

    // incarceration (lambda1 is pre-computed in R as lambda3 * c_true)
    NumericVector l1_r = as<NumericVector>(params_r["lambda1"]);
    NumericVector l2_r = as<NumericVector>(params_r["lambda2"]);
    NumericVector l3_r = as<NumericVector>(params_r["lambda3"]);
    p.lambda1.assign(l1_r.begin(), l1_r.end());
    p.lambda2.assign(l2_r.begin(), l2_r.end());
    p.lambda3.assign(l3_r.begin(), l3_r.end());
    p.pi_recid = as<double>(params_r["pi_recid"]);

    // contact rate and entry
    NumericMatrix C_r = as<NumericMatrix>(params_r["C_contact"]);
    p.C_contact = arma::mat(C_r.begin(), 10, 10, false); // no copy, column-major
    NumericVector beta_r = as<NumericVector>(params_r["beta"]);
    p.beta.assign(beta_r.begin(), beta_r.end());

    // ─── unpack simulation settings from data ────────────────────────────
    double t_start = as<double>(data_r["t_start"]);
    double t_end   = as<double>(data_r["t_end"]);
    double dt      = as<double>(data_r["dt"]);
    int    n_steps = static_cast<int>((t_end - t_start) / dt) + 1;

    NumericVector y0_r = as<NumericVector>(data_r["y0"]);
    std::vector<double> y(y0_r.begin(), y0_r.end());

    // ─── output matrix: rows = time steps, cols = 640 compartments + time
    NumericMatrix out(n_steps, 641);
    out(0, 0) = t_start;
    for (int c = 0; c < 640; ++c) out(0, c+1) = y[c];

    // ─── integrate ────────────────────────────────────────────────────────
    double t = t_start;
    for (int step = 1; step < n_steps; ++step) {
        y = rk4_step(t, dt, y, p);
        t += dt;
        // clamp negatives (numerical safety)
        for (int c = 0; c < 640; ++c) if (y[c] < 0.0) y[c] = 0.0;
        out(step, 0) = t;
        for (int c = 0; c < 640; ++c) out(step, c+1) = y[c];

        // Ageing process
        // For each age grp: y[i] loses y[i]/5 to y[i+1], last age receives inflow only.
        // if (std::fabs(t - std::round(t)) < 1e-9) {
            std::vector<double> y_new = y;

            for (int s = 0; s < 4; ++s) {
                for (int k = 1; k <= 4; ++k) {
                    for (int h = 0; h < 4; ++h) {
                        int base = idx(s, k, h, 0);

                        for (int i = 0; i < 9; ++i) {
                            double y_change = y[base + i] / 5.0 * dt;
                            y_new[base + i]     -= y_change;
                            y_new[base + i + 1] += y_change;
                        }
                    }
                }
            }

            y.swap(y_new);

            // Update output after ageing step
            for (int c = 0; c < 640; ++c) out(step, c+1) = y[c];
        // }
    }

    return out;
}
