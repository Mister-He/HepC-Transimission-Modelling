// =============================================================================
// sim.cpp — HCV Compartmental Model among PWID in Singapore (3-strata)
// =============================================================================
// Rcpp interface for a continuous-time ODE system.
// Compile via:  Rcpp::sourceCpp("sim.cpp")
//
// Compartment index convention
// ─────────────────────────────
//   Stratum  s ∈ {D=0, J=1, X=2}
//              D = never or formerly detained PWID (current PWID in community)
//              J = detained ex-PWID (currently incarcerated; no injecting)
//              X = formerly detained ex-PWID (no longer injecting; no re-arrest)
//   Stage    k ∈ {1,2,3,4}               (NC, CC, DC, HCC)
//   State    h ∈ {u=0, a=1, c=2, t=3, s=4}
//              u = susceptible (never infected)
//              a = acute infection
//              c = chronic infection
//              t = on treatment
//              s = seropositive, RNA-negative (spontaneously cleared or post-SVR)
//   Age grp  i ∈ {0,...,5}              (6 age groups)
//
// Flat index: idx(s,k,h,i) = s*4*5*6 + (k-1)*5*6 + h*6 + i
// Total compartments: 3*4*5*6 = 360
//
// STRUCTURAL JUSTIFICATION (revised calibration scope; see
// docs/calibration/DECISIONS.md and docs/calibration/model_audit.md):
// The source paper (Park et al., medRxiv 10.1101/2025.10.24.25338708) models
// THREE strata: D (never or formerly detained PWID), J (detained ex-PWID) and
// X (formerly detained ex-PWID). Upon release from J, individuals return to D
// with probability pi_recid (relapse to injecting) or move to X with
// probability 1-pi_recid (permanently quit). Only D (current community PWID)
// is at risk of acquiring or transmitting HCV through needle sharing.
// The previous four-stratum implementation kept a separate F (ever-
// incarcerated active PWID) pool with its own re-arrest loop (J -> F -> J).
// That recycling artificially concentrated infected recidivists and saturated
// young-age J seroprevalence near 100% regardless of transmission strength.
// The three-stratum structure matches the paper's model, removes that
// artifact, and lets D prevalence (and hence J prevalence) be driven by
// community transmission and entry of uninfected PWID.
//
// The s-state was added so the model can report anti-HCV seroprevalence
// (acute + chronic + treated + cleared), matching the serology-based prison
// screening target (Changi Prison universal screening, 2014-2016).
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <stdexcept>
using namespace Rcpp;

// ─── helper: flat compartment index ─────────────────────────────────────────
// s in {0,1,2}, k in {1,2,3,4}, h in {0,1,2,3,4}, i in {0,...,5}
inline int idx(int s, int k, int h, int i) {
    return s * 4 * 5 * 6 + (k - 1) * 5 * 6 + h * 6 + i;
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
    double alpha_NC;     // stages 1-2, RBV-free
    double alpha_DC_pos; // stage 3 with RBV+
    double alpha_DC_neg; // stage 3 RBV-free
    double alpha_HCC;    // stage 4

    // --- Treatment initiation rates (per scenario) -------------------------
    std::vector<double> tau; // length 4

    // --- Baseline progression rates (other GT, per year) -------------------
    double p_NC_CC;
    double p_CC_DC;
    double p_CC_HCC;
    double p_DC_HCC;

    // --- GT3 relative risks ------------------------------------------------
    double r3_NC_CC;
    double r3_CC_DC;
    double r3_CC_HCC;
    double r3_DC_HCC;

    // --- Effective progression rates (computed from above) -----------------
    double ptc_NC_CC;
    double ptc_CC_DC;
    double ptc_CC_HCC;
    double ptc_DC_HCC;

    // --- Post-SVR / cleared progression modifiers --------------------------
    double phi_CC_DC;
    double phi_CC_HCC;
    double phi_DC_HCC;

    // --- Effective post-SVR progression rates (computed) -------------------
    double ptu_CC_DC;
    double ptu_CC_HCC;
    double ptu_DC_HCC;

    // --- Background mortality (age-varying, length 6) ----------------------
    std::vector<double> mu;   // mu[i] for age group i (Singapore 2015 baseline)

    // --- Standardized mortality ratio for ever-PWIDs -----------------------
    double omega;  // scalar SMR applied as mu[i] * omega

    // --- Disease-specific excess mortality ---------------------------------
    double mu_DC;   // additional mortality in decompensated cirrhosis
    double mu_HCC;  // additional mortality in HCC

    // --- SVR mortality modifiers -------------------------------------------
    double psi_DC;  // relative mortality risk DC after SVR
    double psi_HCC; // relative mortality risk HCC after SVR

    // --- Seropositive background mortality multiplier (age-varying, len 6) --
    // Applied to cleared (s) individuals: mu_eff_s = mu[i]*omega*eta_s[i].
    std::vector<double> eta_s; // length 6, default all 1

    // --- Incarceration rates (age-varying, length 6 each) ------------------
    std::vector<double> lambda1; // arrest rate of community PWID (D -> J)
    std::vector<double> lambda2; // release rate (J -> D/X)
    std::vector<double> lambda3; // retained for R compatibility (unused)
    double pi_recid;             // prob released detainee relapses to injecting

    // --- Needle-sharing contact rate ---------------------------------------
    arma::mat C_contact; // 6x6

    // --- Population entry (age-varying, time-constant approximation) -------
    std::vector<double> beta; // entry into D_{u,1,i}

    // --- Time-varying transmission multiplier ------------------------------
    // m(t) = m_min + (m_max - m_min) / (1 + exp(-(t - t0)/tau))
    // Rising from m_min (historical, ~1970) to m_max (current, ~2015),
    // anchored so that the multiplier is essentially flat over the final 5+
    // years of the horizon (t0 = 25, tau = 3 with t_end = 45 => flat after
    // ~2005). Represents the historical growth of the PWID HCV epidemic in
    // Singapore and produces the cohort (age-hump) pattern observed in the
    // prison screening data. Defaults m_max = m_min = 1 recover constant
    // transmission.
    double m_min;
    double m_max;
    double t0;
    double m_slope;
};

// =============================================================================
// FORCE OF INFECTION
// gamma_i(t) = q * sum_j C_contact(i,j) * (infectious_D_j / active_D_j)
// Infectious = acute (a) + chronic (c) in D (community PWID).
// Active = all D compartments (u,a,c,t,s across liver stages).
// J and X do not share needles (prison / former PWID), matching the paper.
// =============================================================================
double forceOfInfection(int i, const std::vector<double>& y, const Params& p) {

    double lambda_i = 0.0;

    for (int j = 0; j < 6; ++j) {

        double infectious_j = 0.0;
        for (int k = 1; k <= 4; ++k) {
            infectious_j += y[idx(0,k,1,j)] + y[idx(0,k,2,j)];
        }

        double active_j = 0.0;
        for (int k = 1; k <= 4; ++k)
            for (int h = 0; h <= 4; ++h) {
                active_j += y[idx(0,k,h,j)];
            }

        if (active_j <= 0.0) continue;

        lambda_i += p.C_contact(i, j) * (infectious_j / active_j);
    }

    return p.q * lambda_i;
}

// Time-varying transmission multiplier
double transMult(double t, const Params& p) {
    double m = p.m_min + (p.m_max - p.m_min) /
        (1.0 + std::exp(-(t - p.t0) / p.m_slope));
    return m;
}

// =============================================================================
// ODE RIGHT-HAND SIDE (360 compartments)
// =============================================================================
std::vector<double> rhs(double t,
                        const std::vector<double>& y,
                        const Params& p) {

    std::vector<double> dydt(360, 0.0);

    // SVR12 rates by stage
    auto svr = [&](int k) -> double {
        if (k <= 2) return p.alpha_NC;
        if (k == 3) return p.alpha_DC_pos;
        return p.alpha_HCC;
    };

    // Effective mortality at stage k for state h
    // h=3 (treated): apply SVR mortality modifier for stages 3,4
    // h=4 (cleared/seropositive): background * eta_s[i]
    auto mu_eff = [&](int k, int h, int i) -> double {
        double base = p.mu[i] * p.omega;
        if (h == 4) return base * p.eta_s[i];
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

    // Post-SVR / cleared progression rates (applies to s compartments)
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
    for (int i = 0; i < 6; ++i) {

        double gam  = transMult(t, p) * forceOfInfection(i, y, p);
        double l1   = p.lambda1[i];
        double l2   = p.lambda2[i];
        double pi   = p.pi_recid;

        // =====================================================================
        //  STRATUM D — Never or formerly detained PWID (community, at risk)
        // =====================================================================
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);

            double Du = y[idx(0,k,0,i)];
            double Da = y[idx(0,k,1,i)];
            double Dc = y[idx(0,k,2,i)];
            double Dt = y[idx(0,k,3,i)];
            double Ds = y[idx(0,k,4,i)];
            double Ju = y[idx(1,k,0,i)];
            double Ja = y[idx(1,k,1,i)];
            double Jc = y[idx(1,k,2,i)];
            double Jt = y[idx(1,k,3,i)];
            double Js = y[idx(1,k,4,i)];

            double mu_u = mu_eff(k, 0, i);
            double mu_a = mu_eff(k, 1, i);
            double mu_c = mu_eff(k, 2, i);
            double mu_t = mu_eff(k, 3, i);
            double mu_s = mu_eff(k, 4, i);

            // u: infection, arrest, mortality; inflow beta at k=1; J releases
            dydt[idx(0,k,0,i)] +=
                  (k == 1 ? p.beta[i] : 0.0)
                + pi * l2 * Ju
                - (gam + l1 + mu_u) * Du;

            // a: infection of u and s; clearance to s; chronification to c
            dydt[idx(0,k,1,i)] +=
                  gam * (Du + Ds)
                + pi * l2 * Ja
                - (1.0/p.iota1 + l1 + mu_a) * Da;

            // c: from a (chronic) and t (non-SVR); progression, treatment, mortality
            double prog_c_out = p.tau[k-1];
            double prog_c_in  = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) { prog_c_out += ptc(2,3) + ptc(2,4); }
            if (k == 3) { prog_c_out += ptc(3,4); }
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(0,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(0,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(0,2,2,i)] + ptc(3,4) * y[idx(0,3,2,i)];

            dydt[idx(0,k,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Da
                + ((1.0-a) / p.iota2) * Dt
                + prog_c_in
                + pi * l2 * Jc
                - (prog_c_out + l1 + mu_c) * Dc;

            // t: treatment from c; SVR to s; non-SVR to c
            dydt[idx(0,k,3,i)] +=
                  p.tau[k-1] * Dc
                + pi * l2 * Jt
                - (1.0/p.iota2 + l1 + mu_t) * Dt;

            // s: spontaneous clearance + SVR; reinfection; progression (post-SVR)
            double prog_s_out = 0.0;
            double prog_s_in  = 0.0;
            if (k == 2) prog_s_out += ptu(2,3) + ptu(2,4);
            if (k == 3) prog_s_out += ptu(3,4);
            if (k == 2) prog_s_in += ptu(2,3) * y[idx(0,2,4,i)] + ptu(2,4) * y[idx(0,2,4,i)];
            if (k == 3) prog_s_in += ptu(3,4) * y[idx(0,3,4,i)];

            dydt[idx(0,k,4,i)] +=
                  (p.kappa / p.iota1) * Da
                + (a / p.iota2) * Dt
                + prog_s_in
                + pi * l2 * Js
                - (gam + l1 + mu_s + prog_s_out) * Ds;
        }

        // =====================================================================
        //  STRATUM J — Detained ex-PWID (no new infection; releases to D with
        //  prob pi and X with prob 1-pi)
        // =====================================================================
        for (int k = 1; k <= 4; ++k) {
            for (int h = 0; h <= 4; ++h) {
                double Jkh = y[idx(1,k,h,i)];
                double Dkh = y[idx(0,k,h,i)];
                double inflow = l1 * Dkh;
                double out = l2 + mu_eff(k,h,i);
                dydt[idx(1,k,h,i)] += inflow - out * Jkh;
            }
        }

        // J internal disease dynamics (mirrors D minus gamma/lambda1)
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);

            double Ja = y[idx(1,k,1,i)];
            double Jc = y[idx(1,k,2,i)];
            double Jt = y[idx(1,k,3,i)];

            // a: no new infection in J; recovery/chronification
            dydt[idx(1,k,1,i)] += - (1.0/p.iota1) * Ja;

            // c: from a and t(non-SVR); treatment; progression
            double prog_c_out = p.tau[k-1];
            double prog_c_in  = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) prog_c_out += ptc(2,3) + ptc(2,4);
            if (k == 3) prog_c_out += ptc(3,4);
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(1,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(1,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(1,2,2,i)] + ptc(3,4) * y[idx(1,3,2,i)];

            dydt[idx(1,k,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Ja
                + ((1.0-a) / p.iota2) * Jt
                + prog_c_in
                - prog_c_out * Jc;

            // t
            dydt[idx(1,k,3,i)] +=
                  p.tau[k-1] * Jc
                - (1.0/p.iota2) * Jt;

            // s: clearance + SVR; progression (post-SVR)
            double prog_s_out = 0.0;
            double prog_s_in  = 0.0;
            if (k == 2) prog_s_out += ptu(2,3) + ptu(2,4);
            if (k == 3) prog_s_out += ptu(3,4);
            if (k == 2) prog_s_in += ptu(2,3) * y[idx(1,2,4,i)] + ptu(2,4) * y[idx(1,2,4,i)];
            if (k == 3) prog_s_in += ptu(3,4) * y[idx(1,3,4,i)];

            dydt[idx(1,k,4,i)] +=
                  (p.kappa / p.iota1) * Ja
                + (a / p.iota2) * Jt
                + prog_s_in
                - prog_s_out * y[idx(1,k,4,i)];
        }

        // =====================================================================
        //  STRATUM X — Formerly detained ex-PWID (no infection risk, no arrest)
        // =====================================================================
        for (int k = 1; k <= 4; ++k) {
            for (int h = 0; h <= 4; ++h) {
                double Xkh = y[idx(2,k,h,i)];
                double Jkh = y[idx(1,k,h,i)];
                double inflow = (1.0 - pi) * l2 * Jkh;
                double out = mu_eff(k,h,i);
                dydt[idx(2,k,h,i)] += inflow - out * Xkh;
            }
        }

        // X internal disease dynamics
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);

            double Xa = y[idx(2,k,1,i)];
            double Xc = y[idx(2,k,2,i)];
            double Xt = y[idx(2,k,3,i)];

            // a: no new infection; recovery/chronification
            dydt[idx(2,k,1,i)] += - (1.0/p.iota1) * Xa;

            // c
            double prog_c_out = p.tau[k-1];
            double prog_c_in  = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) prog_c_out += ptc(2,3) + ptc(2,4);
            if (k == 3) prog_c_out += ptc(3,4);
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(2,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(2,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(2,2,2,i)] + ptc(3,4) * y[idx(2,3,2,i)];

            dydt[idx(2,k,2,i)] +=
                  ((1.0-p.kappa) / p.iota1) * Xa
                + ((1.0-a) / p.iota2) * Xt
                + prog_c_in
                - prog_c_out * Xc;

            // t
            dydt[idx(2,k,3,i)] +=
                  p.tau[k-1] * Xc
                - (1.0/p.iota2) * Xt;

            // s
            double prog_s_out = 0.0;
            double prog_s_in  = 0.0;
            if (k == 2) prog_s_out += ptu(2,3) + ptu(2,4);
            if (k == 3) prog_s_out += ptu(3,4);
            if (k == 2) prog_s_in += ptu(2,3) * y[idx(2,2,4,i)] + ptu(2,4) * y[idx(2,2,4,i)];
            if (k == 3) prog_s_in += ptu(3,4) * y[idx(2,3,4,i)];

            dydt[idx(2,k,4,i)] +=
                  (p.kappa / p.iota1) * Xa
                + (a / p.iota2) * Xt
                + prog_s_in
                - prog_s_out * y[idx(2,k,4,i)];
        }

    } // end age group loop

    return dydt;
}

// =============================================================================
// RK4 INTEGRATOR
// =============================================================================
std::vector<double> rk4_step(double t, double dt,
                             const std::vector<double>& y,
                             const Params& p) {
    auto k1 = rhs(t,        y,                                  p);
    std::vector<double> y2(360), y3(360), y4(360);
    for (int i = 0; i < 360; ++i) y2[i] = y[i] + 0.5*dt*k1[i];
    auto k2 = rhs(t+0.5*dt, y2,                                 p);
    for (int i = 0; i < 360; ++i) y3[i] = y[i] + 0.5*dt*k2[i];
    auto k3 = rhs(t+0.5*dt, y3,                                 p);
    for (int i = 0; i < 360; ++i) y4[i] = y[i] + dt*k3[i];
    auto k4 = rhs(t+dt,     y4,                                 p);

    std::vector<double> y_new(360);
    for (int i = 0; i < 360; ++i)
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

    p.alpha_NC     = as<double>(params_r["alpha_NC"]);
    p.alpha_DC_pos = as<double>(params_r["alpha_DC_pos"]);
    p.alpha_DC_neg = as<double>(params_r["alpha_DC_neg"]);
    p.alpha_HCC    = as<double>(params_r["alpha_HCC"]);

    NumericVector tau_r = as<NumericVector>(params_r["tau"]);
    p.tau.assign(tau_r.begin(), tau_r.end());

    p.p_NC_CC   = as<double>(params_r["p_NC_CC"]);
    p.p_CC_DC   = as<double>(params_r["p_CC_DC"]);
    p.p_CC_HCC  = as<double>(params_r["p_CC_HCC"]);
    p.p_DC_HCC  = as<double>(params_r["p_DC_HCC"]);

    p.r3_NC_CC   = as<double>(params_r["r3_NC_CC"]);
    p.r3_CC_DC   = as<double>(params_r["r3_CC_DC"]);
    p.r3_CC_HCC  = as<double>(params_r["r3_CC_HCC"]);
    p.r3_DC_HCC  = as<double>(params_r["r3_DC_HCC"]);

    p.ptc_NC_CC   = (p.rho * p.r3_NC_CC  + (1.0-p.rho)) * p.p_NC_CC;
    p.ptc_CC_DC   = (p.rho * p.r3_CC_DC  + (1.0-p.rho)) * p.p_CC_DC;
    p.ptc_CC_HCC  = (p.rho * p.r3_CC_HCC + (1.0-p.rho)) * p.p_CC_HCC;
    p.ptc_DC_HCC  = (p.rho * p.r3_DC_HCC + (1.0-p.rho)) * p.p_DC_HCC;

    p.phi_CC_DC   = as<double>(params_r["phi_CC_DC"]);
    p.phi_CC_HCC  = as<double>(params_r["phi_CC_HCC"]);
    p.phi_DC_HCC  = as<double>(params_r["phi_DC_HCC"]);

    p.ptu_CC_DC   = p.phi_CC_DC  * p.ptc_CC_DC;
    p.ptu_CC_HCC  = p.phi_CC_HCC * p.ptc_CC_HCC;
    p.ptu_DC_HCC  = p.phi_DC_HCC * p.ptc_DC_HCC;

    NumericVector mu_r = as<NumericVector>(params_r["mu"]);
    p.mu.assign(mu_r.begin(), mu_r.end());
    p.omega   = as<double>(params_r["omega"]);
    p.mu_DC   = as<double>(params_r["mu_DC"]);
    p.mu_HCC  = as<double>(params_r["mu_HCC"]);
    p.psi_DC  = as<double>(params_r["psi_DC"]);
    p.psi_HCC = as<double>(params_r["psi_HCC"]);

    if (params_r.containsElementNamed("eta_s")) {
        NumericVector eta_r = as<NumericVector>(params_r["eta_s"]);
        p.eta_s.assign(eta_r.begin(), eta_r.end());
    } else {
        p.eta_s.assign(6, 1.0);
    }

    p.m_min = params_r.containsElementNamed("m_min") ?
        as<double>(params_r["m_min"]) : 1.0;
    p.m_max = params_r.containsElementNamed("m_max") ?
        as<double>(params_r["m_max"]) : 1.0;
    p.t0    = params_r.containsElementNamed("m_t0") ?
        as<double>(params_r["m_t0"]) : 25.0;
    p.m_slope = params_r.containsElementNamed("m_tau") ?
        as<double>(params_r["m_tau"]) : 3.0;

    NumericVector l1_r = as<NumericVector>(params_r["lambda1"]);
    NumericVector l2_r = as<NumericVector>(params_r["lambda2"]);
    NumericVector l3_r = as<NumericVector>(params_r["lambda3"]);
    p.lambda1.assign(l1_r.begin(), l1_r.end());
    p.lambda2.assign(l2_r.begin(), l2_r.end());
    p.lambda3.assign(l3_r.begin(), l3_r.end());
    p.pi_recid = as<double>(params_r["pi_recid"]);

    NumericMatrix C_r = as<NumericMatrix>(params_r["C_contact"]);
    p.C_contact = arma::mat(C_r.begin(), 6, 6, false);
    NumericVector beta_r = as<NumericVector>(params_r["beta"]);
    p.beta.assign(beta_r.begin(), beta_r.end());

    // ─── unpack simulation settings from data ────────────────────────────
    double t_start = as<double>(data_r["t_start"]);
    double t_end   = as<double>(data_r["t_end"]);
    double dt      = as<double>(data_r["dt"]);
    int    n_steps = static_cast<int>((t_end - t_start) / dt) + 1;

    NumericVector y0_r = as<NumericVector>(data_r["y0"]);
    if (y0_r.size() != 360) stop("y0 must have length 360 (3 strata, 5 HCV states).");
    std::vector<double> y(y0_r.begin(), y0_r.end());

    // ─── output matrix: rows = time steps, cols = 360 compartments + time ──
    NumericMatrix out(n_steps, 361);
    out(0, 0) = t_start;
    for (int c = 0; c < 360; ++c) out(0, c+1) = y[c];

    // ─── integrate ────────────────────────────────────────────────────────
    double t = t_start;
    for (int step = 1; step < n_steps; ++step) {
        y = rk4_step(t, dt, y, p);
        t += dt;
        for (int c = 0; c < 360; ++c) if (y[c] < 0.0) y[c] = 0.0;
        out(step, 0) = t;
        for (int c = 0; c < 360; ++c) out(step, c+1) = y[c];

        // Ageing: y[i] loses y[i]/10 per year to y[i+1]; last group receives
        // inflow only (corrected 10-year adult band implementation, extended
        // to 3 strata and 5 states).
        std::vector<double> y_new = y;
        for (int s = 0; s < 3; ++s) {
            for (int k = 1; k <= 4; ++k) {
                for (int h = 0; h < 5; ++h) {
                    int base = idx(s, k, h, 0);
                    for (int i = 0; i < 5; ++i) {
                        double y_change = y[base + i] / 10.0 * dt;
                        y_new[base + i]     -= y_change;
                        y_new[base + i + 1] += y_change;
                    }
                }
            }
        }
        y.swap(y_new);

        for (int c = 0; c < 360; ++c) out(step, c+1) = y[c];
    }

    return out;
}
