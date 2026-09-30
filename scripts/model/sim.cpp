// HCV transmission model: D/J/F/X x NC/CC/DC/HCC x u/a/c/t/r x six ages.
// r is antibody-reactive, susceptible to reinfection, and non-infectious.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

constexpr int N_STRATA = 4;
constexpr int N_STAGES = 4;
constexpr int N_STATES = 5;
constexpr int N_AGES = 6;
constexpr int N_COMP = N_STRATA * N_STAGES * N_STATES * N_AGES;

inline int idx(int s, int k, int h, int i) {
    return s * N_STAGES * N_STATES * N_AGES +
        (k - 1) * N_STATES * N_AGES + h * N_AGES + i;
}

struct Params {
    double q, kappa, iota1, iota2, rho;
    double alpha_NC, alpha_DC, alpha_DC_pos, alpha_DC_neg, alpha_HCC;
    std::vector<double> tau, tau_stratum, tau_by_stratum;
    int tau_min_age, tau_max_age;
    double p_NC_CC, p_CC_DC, p_CC_HCC, p_DC_HCC;
    double r3_NC_CC, r3_CC_DC, r3_CC_HCC, r3_DC_HCC;
    double ptc_NC_CC, ptc_CC_DC, ptc_CC_HCC, ptc_DC_HCC;
    double phi_CC_DC, phi_CC_HCC, phi_DC_HCC;
    double ptu_CC_DC, ptu_CC_HCC, ptu_DC_HCC;
    std::vector<double> mu;
    double omega, mu_DC, mu_HCC, psi_DC, psi_HCC;
    std::vector<double> lambda1, lambda2, lambda3;
    double pi_recid;
    arma::mat C_contact;
    std::vector<double> beta;
    std::vector<double> beta_base;
    std::vector<std::vector<double>> inflow_schedule;
    double inflow_schedule_start;
    bool has_inflow_schedule;
    double m_min, m_max, t0, m_slope;
};

double inflow_rate(int age, double time, const Params& p) {
    if (!p.has_inflow_schedule || time < p.inflow_schedule_start)
        return p.beta[age];
    // Model year 150 is the 2015 snapshot, so interval 150--151 is 2016.
    int row = static_cast<int>(std::floor(time - p.inflow_schedule_start + 1.0 + 1e-10));
    row = std::max(0, std::min(row, static_cast<int>(p.inflow_schedule.size()) - 1));
    double base = p.beta_base[age];
    double scale = base > 0.0 ? p.beta[age] / base : 0.0;
    return scale * p.inflow_schedule[row][age];
}

double force_of_infection(int age, const std::vector<double>& y,
                          const Params& p) {
    double lambda = 0.0;
    for (int j = 0; j < N_AGES; ++j) {
        double infectious = 0.0, active = 0.0;
        for (int s : {0, 2}) { // D and F inject; J and X do not.
            for (int k = 1; k <= N_STAGES; ++k) {
                infectious += y[idx(s, k, 1, j)] + y[idx(s, k, 2, j)];
                for (int h = 0; h < N_STATES; ++h)
                    active += y[idx(s, k, h, j)];
            }
        }
        if (active > 0.0) lambda += p.C_contact(age, j) * infectious / active;
    }
    return p.q * lambda;
}

double transmission_multiplier(double time, const Params& p) {
    return p.m_min + (p.m_max - p.m_min) /
        (1.0 + std::exp(-(time - p.t0) / p.m_slope));
}

std::vector<double> rhs(double time, const std::vector<double>& y,
                        const Params& p) {
    std::vector<double> dydt(N_COMP, 0.0);
    auto svr = [&](int k) {
        if (k <= 2) return p.alpha_NC;
        if (k == 3) return p.alpha_DC;
        return p.alpha_HCC;
    };
    auto mortality = [&](int k, int h, int age) {
        double base = p.mu[age] * p.omega;
        if (k == 3) return ((h == 3 || h == 4) ? p.psi_DC : 1.0) * (p.mu_DC + base);
        if (k == 4) return ((h == 3 || h == 4) ? p.psi_HCC : 1.0) * (p.mu_HCC + base);
        return base;
    };
    auto progression = [&](int from, int to) {
        if (from == 1 && to == 2) return p.ptc_NC_CC;
        if (from == 2 && to == 3) return p.ptc_CC_DC;
        if (from == 2 && to == 4) return p.ptc_CC_HCC;
        if (from == 3 && to == 4) return p.ptc_DC_HCC;
        return 0.0;
    };
    auto recovered_progression = [&](int from, int to) {
        if (from == 2 && to == 3) return p.ptu_CC_DC;
        if (from == 2 && to == 4) return p.ptu_CC_HCC;
        if (from == 3 && to == 4) return p.ptu_DC_HCC;
        return 0.0;
    };

    for (int age = 0; age < N_AGES; ++age) {
        const double gam = transmission_multiplier(time, p) *
            force_of_infection(age, y, p);
        const double l1 = p.lambda1[age], l2 = p.lambda2[age];
        const double l3 = p.lambda3[age], pi = p.pi_recid;
        auto treatment_rate = [&](int s, int k) {
            if (age + 1 < p.tau_min_age || age + 1 > p.tau_max_age) return 0.0;
            return p.tau_by_stratum[s * N_STAGES + k - 1];
        };

        // Active strata.  r can be reinfected; it remains antibody-reactive.
        for (int s : {0, 2}) {
            const double exit = (s == 0 ? l1 : l3);
            for (int k = 1; k <= N_STAGES; ++k) {
                const double alpha = svr(k);
                const double u = y[idx(s,k,0,age)], a = y[idx(s,k,1,age)];
                const double c = y[idx(s,k,2,age)], t = y[idx(s,k,3,age)];
                const double r = y[idx(s,k,4,age)];
                const double mu_u = mortality(k, 0, age), mu_a = mortality(k, 1, age);
                const double mu_c = mortality(k, 2, age), mu_t = mortality(k, 3, age);
                const double mu_r = mortality(k, 4, age);
                const double tau = treatment_rate(s, k);
                double c_out = tau, c_in = 0.0, r_out = 0.0, r_in = 0.0;
                if (k == 1) c_out += progression(1, 2);
                if (k == 2) c_out += progression(2, 3) + progression(2, 4);
                if (k == 3) c_out += progression(3, 4);
                if (k == 2) {
                    c_in += progression(1, 2) * y[idx(s,1,2,age)];
                    r_out = recovered_progression(2,3) + recovered_progression(2,4);
                }
                if (k == 3) {
                    c_in += progression(2, 3) * y[idx(s,2,2,age)];
                    r_in += recovered_progression(2,3) * y[idx(s,2,4,age)];
                    r_in += recovered_progression(3,4) * y[idx(s,3,4,age)];
                    r_out = recovered_progression(3,4);
                }
                if (k == 4) {
                    c_in += progression(2,4) * y[idx(s,2,2,age)] +
                        progression(3,4) * y[idx(s,3,2,age)];
                    r_in += recovered_progression(2,4) * y[idx(s,2,4,age)] +
                        recovered_progression(3,4) * y[idx(s,3,4,age)];
                }
                const double j_in = (s == 2 ? pi * l2 : 0.0);
                dydt[idx(s,k,0,age)] += (s == 0 && k == 1 ? inflow_rate(age, time, p) : 0.0) +
                    j_in * y[idx(1,k,0,age)] - (gam + exit + mu_u) * u;
                dydt[idx(s,k,1,age)] += j_in * y[idx(1,k,1,age)] +
                    gam * (u + r) - (1.0/p.iota1 + exit + mu_a) * a;
                dydt[idx(s,k,2,age)] += j_in * y[idx(1,k,2,age)] +
                    (1.0-p.kappa)/p.iota1*a + (1.0-alpha)/p.iota2*t + c_in -
                    (c_out + exit + mu_c) * c;
                dydt[idx(s,k,3,age)] += j_in * y[idx(1,k,3,age)] + tau * c -
                    (1.0/p.iota2 + exit + mu_t) * t;
                dydt[idx(s,k,4,age)] += j_in * y[idx(1,k,4,age)] +
                    p.kappa/p.iota1*a + alpha/p.iota2*t + r_in -
                    (gam + exit + mu_r + r_out) * r;
            }
        }

        // J: no infection while incarcerated, but all other state transitions continue.
        for (int k = 1; k <= N_STAGES; ++k) {
            const double alpha = svr(k);
            for (int h = 0; h < N_STATES; ++h) {
                double in = l1*y[idx(0,k,h,age)] + l3*y[idx(2,k,h,age)];
                dydt[idx(1,k,h,age)] += in -
                    (l2 + mortality(k,h,age))*y[idx(1,k,h,age)];
            }
            const double a = y[idx(1,k,1,age)], c = y[idx(1,k,2,age)];
            const double t = y[idx(1,k,3,age)], r = y[idx(1,k,4,age)];
            double c_out = treatment_rate(1,k), c_in = 0.0;
            double r_out = 0.0, r_in = 0.0;
            if (k == 1) c_out += progression(1,2);
            if (k == 2) {
                c_out += progression(2,3) + progression(2,4);
                c_in += progression(1,2)*y[idx(1,1,2,age)];
                r_out = recovered_progression(2,3)+recovered_progression(2,4);
            }
            if (k == 3) {
                c_out += progression(3,4);
                c_in += progression(2,3)*y[idx(1,2,2,age)];
                r_in += recovered_progression(2,3)*y[idx(1,2,4,age)];
                r_in += recovered_progression(3,4)*y[idx(1,3,4,age)];
                r_out = recovered_progression(3,4);
            }
            if (k == 4) {
                c_in += progression(2,4)*y[idx(1,2,2,age)] +
                    progression(3,4)*y[idx(1,3,2,age)];
                r_in += recovered_progression(2,4)*y[idx(1,2,4,age)] +
                    recovered_progression(3,4)*y[idx(1,3,4,age)];
            }
            dydt[idx(1,k,1,age)] += -a/p.iota1;
            dydt[idx(1,k,2,age)] += (1.0-p.kappa)/p.iota1*a +
                (1.0-alpha)/p.iota2*t + c_in - c_out*c;
            dydt[idx(1,k,3,age)] += treatment_rate(1,k)*c - t/p.iota2;
            dydt[idx(1,k,4,age)] += p.kappa/p.iota1*a + alpha/p.iota2*t +
                r_in - r_out*r;
        }

        // X: no infection and no re-arrest.
        for (int k = 1; k <= N_STAGES; ++k) {
            const double alpha = svr(k);
            for (int h = 0; h < N_STATES; ++h) {
                double in = (1.0-pi)*l2*y[idx(1,k,h,age)];
                dydt[idx(3,k,h,age)] += in -
                    mortality(k,h,age)*y[idx(3,k,h,age)];
            }
            const double a = y[idx(3,k,1,age)], c = y[idx(3,k,2,age)];
            const double t = y[idx(3,k,3,age)], r = y[idx(3,k,4,age)];
            double c_out = treatment_rate(3,k), c_in = 0.0;
            double r_out = 0.0, r_in = 0.0;
            if (k == 1) c_out += progression(1,2);
            if (k == 2) {
                c_out += progression(2,3) + progression(2,4);
                c_in += progression(1,2)*y[idx(3,1,2,age)];
                r_out = recovered_progression(2,3)+recovered_progression(2,4);
            }
            if (k == 3) {
                c_out += progression(3,4);
                c_in += progression(2,3)*y[idx(3,2,2,age)];
                r_in += recovered_progression(2,3)*y[idx(3,2,4,age)];
                r_in += recovered_progression(3,4)*y[idx(3,3,4,age)];
                r_out = recovered_progression(3,4);
            }
            if (k == 4) {
                c_in += progression(2,4)*y[idx(3,2,2,age)] +
                    progression(3,4)*y[idx(3,3,2,age)];
                r_in += recovered_progression(2,4)*y[idx(3,2,4,age)] +
                    recovered_progression(3,4)*y[idx(3,3,4,age)];
            }
            dydt[idx(3,k,1,age)] += -a/p.iota1;
            dydt[idx(3,k,2,age)] += (1.0-p.kappa)/p.iota1*a +
                (1.0-alpha)/p.iota2*t + c_in - c_out*c;
            dydt[idx(3,k,3,age)] += treatment_rate(3,k)*c - t/p.iota2;
            dydt[idx(3,k,4,age)] += p.kappa/p.iota1*a + alpha/p.iota2*t +
                r_in - r_out*r;
        }
    }
    return dydt;
}

std::vector<double> rk4_step(double time, double dt,
                             const std::vector<double>& y, const Params& p) {
    auto k1 = rhs(time, y, p);
    std::vector<double> y2(N_COMP), y3(N_COMP), y4(N_COMP);
    for (int i = 0; i < N_COMP; ++i) y2[i] = y[i] + 0.5*dt*k1[i];
    auto k2 = rhs(time+0.5*dt, y2, p);
    for (int i = 0; i < N_COMP; ++i) y3[i] = y[i] + 0.5*dt*k2[i];
    auto k3 = rhs(time+0.5*dt, y3, p);
    for (int i = 0; i < N_COMP; ++i) y4[i] = y[i] + dt*k3[i];
    auto k4 = rhs(time+dt, y4, p);
    std::vector<double> out(N_COMP);
    for (int i = 0; i < N_COMP; ++i)
        out[i] = y[i] + (dt/6.0)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    return out;
}

// [[Rcpp::export]]
NumericMatrix run_sim(List params_r, List data_r) {
    Params p;
    p.q = as<double>(params_r["q"]); p.kappa = as<double>(params_r["kappa"]);
    p.iota1 = as<double>(params_r["iota1"]); p.iota2 = as<double>(params_r["iota2"]);
    p.rho = as<double>(params_r["rho"]); p.alpha_NC = as<double>(params_r["alpha_NC"]);
    p.alpha_DC_pos = as<double>(params_r["alpha_DC_pos"]); p.alpha_DC_neg = as<double>(params_r["alpha_DC_neg"]);
    p.alpha_DC = params_r.containsElementNamed("alpha_DC") ? as<double>(params_r["alpha_DC"]) : p.alpha_DC_pos;
    p.alpha_HCC = as<double>(params_r["alpha_HCC"]);
    NumericVector tau = as<NumericVector>(params_r["tau"]); p.tau.assign(tau.begin(), tau.end());
    if (params_r.containsElementNamed("tau_stratum")) {
        NumericVector x = as<NumericVector>(params_r["tau_stratum"]); p.tau_stratum.assign(x.begin(), x.end());
    } else p.tau_stratum.assign(N_STRATA, 1.0);
    p.tau_min_age = params_r.containsElementNamed("tau_min_age") ? as<int>(params_r["tau_min_age"]) : 1;
    p.tau_max_age = params_r.containsElementNamed("tau_max_age") ? as<int>(params_r["tau_max_age"]) : N_AGES;
    if (params_r.containsElementNamed("tau_by_stratum")) {
        NumericVector x = as<NumericVector>(params_r["tau_by_stratum"]);
        if (x.size() != N_STRATA*N_STAGES) stop("tau_by_stratum must have length 16");
        p.tau_by_stratum.assign(x.begin(), x.end());
    } else {
        p.tau_by_stratum.resize(N_STRATA*N_STAGES);
        for (int s = 0; s < N_STRATA; ++s) for (int k = 0; k < N_STAGES; ++k)
            p.tau_by_stratum[s*N_STAGES+k] = p.tau[k] * p.tau_stratum[s];
    }
    p.p_NC_CC = as<double>(params_r["p_NC_CC"]); p.p_CC_DC = as<double>(params_r["p_CC_DC"]);
    p.p_CC_HCC = as<double>(params_r["p_CC_HCC"]); p.p_DC_HCC = as<double>(params_r["p_DC_HCC"]);
    p.r3_NC_CC = as<double>(params_r["r3_NC_CC"]); p.r3_CC_DC = as<double>(params_r["r3_CC_DC"]);
    p.r3_CC_HCC = as<double>(params_r["r3_CC_HCC"]); p.r3_DC_HCC = as<double>(params_r["r3_DC_HCC"]);
    p.ptc_NC_CC = (p.rho*p.r3_NC_CC + 1-p.rho)*p.p_NC_CC;
    p.ptc_CC_DC = (p.rho*p.r3_CC_DC + 1-p.rho)*p.p_CC_DC;
    p.ptc_CC_HCC = (p.rho*p.r3_CC_HCC + 1-p.rho)*p.p_CC_HCC;
    p.ptc_DC_HCC = (p.rho*p.r3_DC_HCC + 1-p.rho)*p.p_DC_HCC;
    p.phi_CC_DC = as<double>(params_r["phi_CC_DC"]); p.phi_CC_HCC = as<double>(params_r["phi_CC_HCC"]);
    p.phi_DC_HCC = as<double>(params_r["phi_DC_HCC"]);
    p.ptu_CC_DC = p.phi_CC_DC*p.ptc_CC_DC; p.ptu_CC_HCC = p.phi_CC_HCC*p.ptc_CC_HCC; p.ptu_DC_HCC = p.phi_DC_HCC*p.ptc_DC_HCC;
    NumericVector mu = as<NumericVector>(params_r["mu"]); p.mu.assign(mu.begin(), mu.end());
    p.omega = as<double>(params_r["omega"]); p.mu_DC = as<double>(params_r["mu_DC"]); p.mu_HCC = as<double>(params_r["mu_HCC"]);
    p.psi_DC = as<double>(params_r["psi_DC"]); p.psi_HCC = as<double>(params_r["psi_HCC"]);
    p.m_min = params_r.containsElementNamed("m_min") ? as<double>(params_r["m_min"]) : 1.0;
    p.m_max = params_r.containsElementNamed("m_max") ? as<double>(params_r["m_max"]) : 1.0;
    p.t0 = params_r.containsElementNamed("m_t0") ? as<double>(params_r["m_t0"]) : 25.0;
    p.m_slope = params_r.containsElementNamed("m_tau") ? as<double>(params_r["m_tau"]) : 3.0;
    NumericVector l1 = as<NumericVector>(params_r["lambda1"]), l2 = as<NumericVector>(params_r["lambda2"]), l3 = as<NumericVector>(params_r["lambda3"]);
    p.lambda1.assign(l1.begin(), l1.end()); p.lambda2.assign(l2.begin(), l2.end()); p.lambda3.assign(l3.begin(), l3.end());
    p.pi_recid = as<double>(params_r["pi_recid"]);
    NumericMatrix C = as<NumericMatrix>(params_r["C_contact"]); p.C_contact = arma::mat(C.begin(), N_AGES, N_AGES, false);
    NumericVector beta = as<NumericVector>(params_r["beta"]); p.beta.assign(beta.begin(), beta.end());
    if (params_r.containsElementNamed("beta_base")) {
        NumericVector beta_base = as<NumericVector>(params_r["beta_base"]);
        p.beta_base.assign(beta_base.begin(), beta_base.end());
    } else {
        p.beta_base = p.beta;
    }
    double t_start = as<double>(data_r["t_start"]), t_end = as<double>(data_r["t_end"]), dt = as<double>(data_r["dt"]);
    p.has_inflow_schedule = data_r.containsElementNamed("inflow_schedule");
    p.inflow_schedule_start = data_r.containsElementNamed("inflow_schedule_start") ?
        as<double>(data_r["inflow_schedule_start"]) : 0.0;
    if (p.has_inflow_schedule) {
        NumericMatrix schedule = as<NumericMatrix>(data_r["inflow_schedule"]);
        p.inflow_schedule.resize(schedule.nrow(), std::vector<double>(N_AGES));
        for (int r = 0; r < schedule.nrow(); ++r)
            for (int a = 0; a < N_AGES; ++a)
                p.inflow_schedule[r][a] = schedule(r, a);
    }
    int n_steps = static_cast<int>((t_end-t_start)/dt) + 1;
    NumericVector y0 = as<NumericVector>(data_r["y0"]);
    if (y0.size() != N_COMP) stop("y0 must have length 480 (five HCV states)");
    std::vector<double> y(y0.begin(), y0.end()); NumericMatrix out(n_steps, N_COMP+1);
    out(0,0) = t_start; for (int c = 0; c < N_COMP; ++c) out(0,c+1) = y[c];
    double time = t_start;
    for (int step = 1; step < n_steps; ++step) {
        y = rk4_step(time, dt, y, p); time += dt;
        for (double& value : y) if (value < 0.0) value = 0.0;
        std::vector<double> aged = y;
        for (int s = 0; s < N_STRATA; ++s) for (int k = 1; k <= N_STAGES; ++k)
            for (int h = 0; h < N_STATES; ++h) for (int age = 0; age < N_AGES-1; ++age) {
                double moved = y[idx(s,k,h,age)] / 10.0 * dt;
                aged[idx(s,k,h,age)] -= moved; aged[idx(s,k,h,age+1)] += moved;
            }
        y.swap(aged); out(step,0) = time;
        for (int c = 0; c < N_COMP; ++c) out(step,c+1) = y[c];
    }
    return out;
}
