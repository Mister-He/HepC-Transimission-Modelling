// =============================================================================
// sim.cpp — HCV Compartmental Model among PWID in Singapore (4-strata)
// =============================================================================
// Rcpp interface for a continuous-time ODE system.
// Compile via:  Rcpp::sourceCpp("sim.cpp")
//
// Compartment index convention
// ─────────────────────────────
//   Stratum  s ∈ {D=0, J=1, F=2, X=3}
//              D = never arrested, active PWID
//              J = currently arrested/detained
//              F = ever arrested, active PWID
//              X = former PWID (no injecting, no re-arrest)
//   Stage    k ∈ {1,2,3,4}               (NC, CC, DC, HCC)
//   State    h ∈ {u=0, a=1, c=2, t=3}
//              u = susceptible / post-SVR
//              a = acute infection
//              c = chronic infection
//              t = on treatment
//   Age grp  i ∈ {0,...,5}              (6 age groups)
//
// Flat index: idx(s,k,h,i) = s*4*4*6 + (k-1)*4*6 + h*6 + i
// Total compartments: 4*4*4*6 = 384
//
// Flows:
//   D --lambda1--> J (first arrest)
//   F --lambda3--> J (re-arrest)
//   J --lambda2--> F (pi_recid) or X (1-pi_recid)
//   beta enters D_{u,1,i}
//   D and F are active PWID at risk of infection; J and X are not.
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <stdexcept>
using namespace Rcpp;

inline int idx(int s, int k, int h, int i) {
    return s * 4 * 4 * 6 + (k - 1) * 4 * 6 + h * 6 + i;
}

struct Params {
    double q;
    double kappa;
    double iota1;
    double iota2;
    double rho;
    double alpha_NC;
    double alpha_DC_pos;
    double alpha_DC_neg;
    double alpha_HCC;
    std::vector<double> tau;
    std::vector<double> tau_stratum; // length 4: eligible strata (0/1)
    int tau_min_age;                 // 1-based minimum eligible age group
    double p_NC_CC, p_CC_DC, p_CC_HCC, p_DC_HCC;
    double r3_NC_CC, r3_CC_DC, r3_CC_HCC, r3_DC_HCC;
    double ptc_NC_CC, ptc_CC_DC, ptc_CC_HCC, ptc_DC_HCC;
    double phi_CC_DC, phi_CC_HCC, phi_DC_HCC;
    double ptu_CC_DC, ptu_CC_HCC, ptu_DC_HCC;
    std::vector<double> mu;
    double omega;
    double mu_DC, mu_HCC;
    double psi_DC, psi_HCC;
    std::vector<double> eta_s;
    std::vector<double> lambda1; // first arrest (D -> J)
    std::vector<double> lambda2; // release (J -> F/X)
    std::vector<double> lambda3; // re-arrest (F -> J)
    double pi_recid;
    arma::mat C_contact;
    std::vector<double> beta;
    double m_min, m_max, t0, m_slope;
};

double forceOfInfection(int i, const std::vector<double>& y, const Params& p) {
    double lambda_i = 0.0;
    for (int j = 0; j < 6; ++j) {
        double infectious_j = 0.0;
        double active_j = 0.0;
        for (int s = 0; s <= 3; ++s) {
            if (s == 1 || s == 3) continue; // J and X are not active
            for (int k = 1; k <= 4; ++k) {
                infectious_j += y[idx(s,k,1,j)] + y[idx(s,k,2,j)];
                for (int h = 0; h <= 3; ++h) active_j += y[idx(s,k,h,j)];
            }
        }
        if (active_j <= 0.0) continue;
        lambda_i += p.C_contact(i, j) * (infectious_j / active_j);
    }
    return p.q * lambda_i;
}

double transMult(double t, const Params& p) {
    return p.m_min + (p.m_max - p.m_min) /
        (1.0 + std::exp(-(t - p.t0) / p.m_slope));
}

std::vector<double> rhs(double t, const std::vector<double>& y, const Params& p) {
    std::vector<double> dydt(384, 0.0);

    auto svr = [&](int k) -> double {
        if (k <= 2) return p.alpha_NC;
        if (k == 3) return p.alpha_DC_pos;
        return p.alpha_HCC;
    };

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

    auto ptc = [&](int from_k, int to_k) -> double {
        if (from_k == 1 && to_k == 2) return p.ptc_NC_CC;
        if (from_k == 2 && to_k == 3) return p.ptc_CC_DC;
        if (from_k == 2 && to_k == 4) return p.ptc_CC_HCC;
        if (from_k == 3 && to_k == 4) return p.ptc_DC_HCC;
        return 0.0;
    };

    for (int i = 0; i < 6; ++i) {
        double gam = transMult(t, p) * forceOfInfection(i, y, p);
        double l1 = p.lambda1[i];
        double l2 = p.lambda2[i];
        double l3 = p.lambda3[i];
        double pi = p.pi_recid;
        auto treat = [&](int s) -> double {
            double elig = p.tau_stratum[s] * ((i + 1) >= p.tau_min_age ? 1.0 : 0.0);
            return elig;
        };

        // ---- D: never arrested, active ----
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);
            double Du = y[idx(0,k,0,i)], Da = y[idx(0,k,1,i)];
            double Dc = y[idx(0,k,2,i)], Dt = y[idx(0,k,3,i)];
            double mu_u = mu_eff(k,0,i), mu_a = mu_eff(k,1,i);
            double mu_c = mu_eff(k,2,i), mu_t = mu_eff(k,3,i);

            double treat_rate = p.tau[k-1] * treat(0);
            double prog_c_out = treat_rate;
            double prog_c_in = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) prog_c_out += ptc(2,3) + ptc(2,4);
            if (k == 3) prog_c_out += ptc(3,4);
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(0,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(0,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(0,2,2,i)] + ptc(3,4) * y[idx(0,3,2,i)];

            dydt[idx(0,k,0,i)] += (k == 1 ? p.beta[i] : 0.0)
                + (p.kappa / p.iota1) * Da + (a / p.iota2) * Dt
                - (gam + l1 + mu_u) * Du;
            dydt[idx(0,k,1,i)] += gam * Du
                - (1.0/p.iota1 + l1 + mu_a) * Da;
            dydt[idx(0,k,2,i)] += ((1.0-p.kappa)/p.iota1) * Da
                + ((1.0-a)/p.iota2) * Dt + prog_c_in
                - (prog_c_out + l1 + mu_c) * Dc;
            dydt[idx(0,k,3,i)] += treat_rate * Dc
                - (1.0/p.iota2 + l1 + mu_t) * Dt;
        }

        // ---- F: ever arrested, active ----
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);
            double Fu = y[idx(2,k,0,i)], Fa = y[idx(2,k,1,i)];
            double Fc = y[idx(2,k,2,i)], Ft = y[idx(2,k,3,i)];
            double Ju = y[idx(1,k,0,i)], Ja = y[idx(1,k,1,i)];
            double Jc = y[idx(1,k,2,i)], Jt = y[idx(1,k,3,i)];
            double mu_u = mu_eff(k,0,i), mu_a = mu_eff(k,1,i);
            double mu_c = mu_eff(k,2,i), mu_t = mu_eff(k,3,i);

            double treat_rate = p.tau[k-1] * treat(2);
            double prog_c_out = treat_rate;
            double prog_c_in = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) prog_c_out += ptc(2,3) + ptc(2,4);
            if (k == 3) prog_c_out += ptc(3,4);
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(2,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(2,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(2,2,2,i)] + ptc(3,4) * y[idx(2,3,2,i)];

            dydt[idx(2,k,0,i)] += pi * l2 * Ju
                + (p.kappa/p.iota1) * Fa + (a/p.iota2) * Ft
                - (gam + l3 + mu_u) * Fu;
            dydt[idx(2,k,1,i)] += pi * l2 * Ja + gam * Fu
                - (1.0/p.iota1 + l3 + mu_a) * Fa;
            dydt[idx(2,k,2,i)] += pi * l2 * Jc
                + ((1.0-p.kappa)/p.iota1) * Fa + ((1.0-a)/p.iota2) * Ft
                + prog_c_in - (prog_c_out + l3 + mu_c) * Fc;
            dydt[idx(2,k,3,i)] += pi * l2 * Jt + treat_rate * Fc
                - (1.0/p.iota2 + l3 + mu_t) * Ft;
        }

        // ---- J: currently arrested, no injecting ----
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);
            for (int h = 0; h <= 3; ++h) {
                double inflow = l1 * y[idx(0,k,h,i)] + l3 * y[idx(2,k,h,i)];
                double out = l2 + mu_eff(k,h,i);
                dydt[idx(1,k,h,i)] += inflow - out * y[idx(1,k,h,i)];
            }
            double Ja = y[idx(1,k,1,i)], Jc = y[idx(1,k,2,i)], Jt = y[idx(1,k,3,i)];
            double treat_rate = p.tau[k-1] * treat(1);
            double prog_c_out = treat_rate;
            double prog_c_in = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) prog_c_out += ptc(2,3) + ptc(2,4);
            if (k == 3) prog_c_out += ptc(3,4);
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(1,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(1,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(1,2,2,i)] + ptc(3,4) * y[idx(1,3,2,i)];
            dydt[idx(1,k,1,i)] += - (1.0/p.iota1) * Ja;
            dydt[idx(1,k,2,i)] += ((1.0-p.kappa)/p.iota1) * Ja
                + ((1.0-a)/p.iota2) * Jt + prog_c_in - prog_c_out * Jc;
            dydt[idx(1,k,3,i)] += treat_rate * Jc - (1.0/p.iota2) * Jt;
        }

        // ---- X: former, no injecting, no arrest ----
        for (int k = 1; k <= 4; ++k) {
            double a = svr(k);
            for (int h = 0; h <= 3; ++h) {
                double inflow = (1.0 - pi) * l2 * y[idx(1,k,h,i)];
                dydt[idx(3,k,h,i)] += inflow - mu_eff(k,h,i) * y[idx(3,k,h,i)];
            }
            double Xa = y[idx(3,k,1,i)], Xc = y[idx(3,k,2,i)], Xt = y[idx(3,k,3,i)];
            double treat_rate = p.tau[k-1] * treat(3);
            double prog_c_out = treat_rate;
            double prog_c_in = 0.0;
            if (k == 1) prog_c_out += ptc(1,2);
            if (k == 2) prog_c_out += ptc(2,3) + ptc(2,4);
            if (k == 3) prog_c_out += ptc(3,4);
            if (k == 2) prog_c_in += ptc(1,2) * y[idx(3,1,2,i)];
            if (k == 3) prog_c_in += ptc(2,3) * y[idx(3,2,2,i)];
            if (k == 4) prog_c_in += ptc(2,4) * y[idx(3,2,2,i)] + ptc(3,4) * y[idx(3,3,2,i)];
            dydt[idx(3,k,1,i)] += - (1.0/p.iota1) * Xa;
            dydt[idx(3,k,2,i)] += ((1.0-p.kappa)/p.iota1) * Xa
                + ((1.0-a)/p.iota2) * Xt + prog_c_in - prog_c_out * Xc;
            dydt[idx(3,k,3,i)] += treat_rate * Xc - (1.0/p.iota2) * Xt;
        }
    }
    return dydt;
}

std::vector<double> rk4_step(double t, double dt,
                             const std::vector<double>& y, const Params& p) {
    auto k1 = rhs(t, y, p);
    std::vector<double> y2(384), y3(384), y4(384);
    for (int i = 0; i < 384; ++i) y2[i] = y[i] + 0.5*dt*k1[i];
    auto k2 = rhs(t+0.5*dt, y2, p);
    for (int i = 0; i < 384; ++i) y3[i] = y[i] + 0.5*dt*k2[i];
    auto k3 = rhs(t+0.5*dt, y3, p);
    for (int i = 0; i < 384; ++i) y4[i] = y[i] + dt*k3[i];
    auto k4 = rhs(t+dt, y4, p);
    std::vector<double> y_new(384);
    for (int i = 0; i < 384; ++i)
        y_new[i] = y[i] + (dt/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    return y_new;
}

// [[Rcpp::export]]
NumericMatrix run_sim(List params_r, List data_r) {
    Params p;
    p.q = as<double>(params_r["q"]);
    p.kappa = as<double>(params_r["kappa"]);
    p.iota1 = as<double>(params_r["iota1"]);
    p.iota2 = as<double>(params_r["iota2"]);
    p.rho = as<double>(params_r["rho"]);
    p.alpha_NC = as<double>(params_r["alpha_NC"]);
    p.alpha_DC_pos = as<double>(params_r["alpha_DC_pos"]);
    p.alpha_DC_neg = as<double>(params_r["alpha_DC_neg"]);
    p.alpha_HCC = as<double>(params_r["alpha_HCC"]);
    NumericVector tau_r = as<NumericVector>(params_r["tau"]);
    p.tau.assign(tau_r.begin(), tau_r.end());
    if (params_r.containsElementNamed("tau_stratum")) {
        NumericVector ts = as<NumericVector>(params_r["tau_stratum"]);
        p.tau_stratum.assign(ts.begin(), ts.end());
    } else {
        p.tau_stratum.assign(4, 1.0);
    }
    p.tau_min_age = params_r.containsElementNamed("tau_min_age") ?
        as<int>(params_r["tau_min_age"]) : 1;
    p.p_NC_CC = as<double>(params_r["p_NC_CC"]);
    p.p_CC_DC = as<double>(params_r["p_CC_DC"]);
    p.p_CC_HCC = as<double>(params_r["p_CC_HCC"]);
    p.p_DC_HCC = as<double>(params_r["p_DC_HCC"]);
    p.r3_NC_CC = as<double>(params_r["r3_NC_CC"]);
    p.r3_CC_DC = as<double>(params_r["r3_CC_DC"]);
    p.r3_CC_HCC = as<double>(params_r["r3_CC_HCC"]);
    p.r3_DC_HCC = as<double>(params_r["r3_DC_HCC"]);
    p.ptc_NC_CC = (p.rho*p.r3_NC_CC + (1.0-p.rho))*p.p_NC_CC;
    p.ptc_CC_DC = (p.rho*p.r3_CC_DC + (1.0-p.rho))*p.p_CC_DC;
    p.ptc_CC_HCC = (p.rho*p.r3_CC_HCC + (1.0-p.rho))*p.p_CC_HCC;
    p.ptc_DC_HCC = (p.rho*p.r3_DC_HCC + (1.0-p.rho))*p.p_DC_HCC;
    p.phi_CC_DC = as<double>(params_r["phi_CC_DC"]);
    p.phi_CC_HCC = as<double>(params_r["phi_CC_HCC"]);
    p.phi_DC_HCC = as<double>(params_r["phi_DC_HCC"]);
    p.ptu_CC_DC = p.phi_CC_DC * p.ptc_CC_DC;
    p.ptu_CC_HCC = p.phi_CC_HCC * p.ptc_CC_HCC;
    p.ptu_DC_HCC = p.phi_DC_HCC * p.ptc_DC_HCC;
    NumericVector mu_r = as<NumericVector>(params_r["mu"]);
    p.mu.assign(mu_r.begin(), mu_r.end());
    p.omega = as<double>(params_r["omega"]);
    p.mu_DC = as<double>(params_r["mu_DC"]);
    p.mu_HCC = as<double>(params_r["mu_HCC"]);
    p.psi_DC = as<double>(params_r["psi_DC"]);
    p.psi_HCC = as<double>(params_r["psi_HCC"]);
    if (params_r.containsElementNamed("eta_s")) {
        NumericVector eta_r = as<NumericVector>(params_r["eta_s"]);
        p.eta_s.assign(eta_r.begin(), eta_r.end());
    } else {
        p.eta_s.assign(6, 1.0);
    }
    p.m_min = params_r.containsElementNamed("m_min") ? as<double>(params_r["m_min"]) : 1.0;
    p.m_max = params_r.containsElementNamed("m_max") ? as<double>(params_r["m_max"]) : 1.0;
    p.t0 = params_r.containsElementNamed("m_t0") ? as<double>(params_r["m_t0"]) : 25.0;
    p.m_slope = params_r.containsElementNamed("m_tau") ? as<double>(params_r["m_tau"]) : 3.0;
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

    double t_start = as<double>(data_r["t_start"]);
    double t_end = as<double>(data_r["t_end"]);
    double dt = as<double>(data_r["dt"]);
    int n_steps = static_cast<int>((t_end - t_start) / dt) + 1;
    NumericVector y0_r = as<NumericVector>(data_r["y0"]);
    if (y0_r.size() != 384) stop("y0 must have length 384 (4 strata, 4 HCV states).");
    std::vector<double> y(y0_r.begin(), y0_r.end());
    NumericMatrix out(n_steps, 385);
    out(0, 0) = t_start;
    for (int c = 0; c < 384; ++c) out(0, c+1) = y[c];

    double t = t_start;
    for (int step = 1; step < n_steps; ++step) {
        y = rk4_step(t, dt, y, p);
        t += dt;
        for (int c = 0; c < 384; ++c) if (y[c] < 0.0) y[c] = 0.0;
        out(step, 0) = t;
        for (int c = 0; c < 384; ++c) out(step, c+1) = y[c];
        std::vector<double> y_new = y;
        for (int s = 0; s < 4; ++s) {
            for (int k = 1; k <= 4; ++k) {
                for (int h = 0; h < 4; ++h) {
                    int base = idx(s, k, h, 0);
                    for (int i = 0; i < 5; ++i) {
                        double y_change = y[base + i] / 10.0 * dt;
                        y_new[base + i] -= y_change;
                        y_new[base + i + 1] += y_change;
                    }
                }
            }
        }
        y.swap(y_new);
        for (int c = 0; c < 384; ++c) out(step, c+1) = y[c];
    }
    return out;
}
