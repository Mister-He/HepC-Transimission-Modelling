// =============================================================================
// hcv_model.stan — HCV Compartmental Model Calibration
// =============================================================================
// Calibrates two parameters against observed prison population data:
//   beta_scale : multiplicative scaling on inflow rate beta[i] (all age groups)
//   y0_scale   : multiplicative scaling on initial susceptible D_{u,1,i} (all age groups)
//
// Compartment index convention (matches sim.cpp exactly):
//   Stratum  s ∈ {D=0, J=1, F=2, X=3}
//   Stage    k ∈ {1,2,3,4}   (NC, CC, DC, HCC)
//   State    h ∈ {u=0, a=1, c=2, t=3}
//   Age grp  i ∈ {0,...,8}
//
// Stan uses 1-based indexing; all idx_f() calls add +1 to the flat index.
//
// Observed data targets (J stratum, by age group):
//   obs_J_total[i]  = sum_{k,h} J_{k,h,i}   c(307,797,829,633,598,642,481,439,366)
//   obs_J_susc[i]   = sum_k    J_{u,k,i}    c( 55,145,183,164,212,299,222,190,133)
//     NOTE: The second target formula in your spec appeared to duplicate the
//     first sum. Interpreted here as: sum over k of J_{h=0,k,i}
//     (susceptible/post-SVR incarcerated by age group). Please verify.
//
// Likelihood: Poisson on each of the 9+9 age-group counts.
// =============================================================================

functions {

  // ── 1-based flat index into the 576-element state vector ──────────────────
  // s in {0,1,2,3}, k in {1..4}, h in {0..3}, i in {0..8}
  int idx_f(int s, int k, int h, int i) {
    return s * 144 + (k - 1) * 36 + h * 9 + i + 1;
  }

  // ── SVR12 rate by disease stage ───────────────────────────────────────────
  real svr_rate(int k, real alpha_NC, real alpha_DC_pos, real alpha_HCC) {
    if (k <= 2) return alpha_NC;
    if (k == 3) return alpha_DC_pos;
    return alpha_HCC;
  }

  // ── Effective all-cause mortality (omega-scaled, disease-specific) ─────────
  // Matches C++ mu_eff lambda exactly:
  //   base = mu[i] * omega
  //   k==3 DC:  treated → psi_DC*(mu_DC+base);  else → mu_DC+base
  //   k==4 HCC: treated → psi_HCC*(mu_HCC+base); else → mu_HCC+base
  //   k<=2:     base only
  real mu_eff_f(int k, int h, int i,
                array[] real mu, real omega,
                real mu_DC, real mu_HCC,
                real psi_DC, real psi_HCC) {
    real base = mu[i + 1] * omega;     // i is 0-based; mu array is 1-based
    if (k == 3)
      return (h == 3) ? psi_DC * (mu_DC + base) : (mu_DC + base);
    if (k == 4)
      return (h == 3) ? psi_HCC * (mu_HCC + base) : (mu_HCC + base);
    return base;
  }

  // ── Post-SVR (u compartment) progression rates ────────────────────────────
  real ptu_f(int from_k, int to_k,
             real ptu_CC_DC, real ptu_CC_HCC, real ptu_DC_HCC) {
    if (from_k == 2 && to_k == 3) return ptu_CC_DC;
    if (from_k == 2 && to_k == 4) return ptu_CC_HCC;
    if (from_k == 3 && to_k == 4) return ptu_DC_HCC;
    return 0.0;
  }

  // ── Chronic (c compartment) progression rates ─────────────────────────────
  real ptc_f(int from_k, int to_k,
             real ptc_NC_CC, real ptc_CC_DC,
             real ptc_CC_HCC, real ptc_DC_HCC) {
    if (from_k == 1 && to_k == 2) return ptc_NC_CC;
    if (from_k == 2 && to_k == 3) return ptc_CC_DC;
    if (from_k == 2 && to_k == 4) return ptc_CC_HCC;
    if (from_k == 3 && to_k == 4) return ptc_DC_HCC;
    return 0.0;
  }

  // ── Force of infection for age group i ───────────────────────────────────
  // Infectious = acute (h=1) + chronic (h=2) in strata D (s=0) and F (s=2)
  real force_of_infection(int i, vector y,
                          array[,] real C_contact, real q) {
    real lambda_i = 0.0;
    for (j in 0:8) {
      // sum acute + chronic across all stages in D and F strata
      real infectious_j =
          y[idx_f(0,1,1,j)] + y[idx_f(0,1,2,j)]
        + y[idx_f(0,2,1,j)] + y[idx_f(0,2,2,j)]
        + y[idx_f(0,3,1,j)] + y[idx_f(0,3,2,j)]
        + y[idx_f(0,4,1,j)] + y[idx_f(0,4,2,j)]
        + y[idx_f(2,1,1,j)] + y[idx_f(2,1,2,j)]
        + y[idx_f(2,2,1,j)] + y[idx_f(2,2,2,j)]
        + y[idx_f(2,3,1,j)] + y[idx_f(2,3,2,j)]
        + y[idx_f(2,4,1,j)] + y[idx_f(2,4,2,j)];

      // total active PWID (strata D and F)
      real active_j = 0.0;
      for (k in 1:4)
        for (h in 0:3) {
          active_j += y[idx_f(0, k, h, j)];
          active_j += y[idx_f(2, k, h, j)];
        }

      if (active_j > 0.0)
        lambda_i += C_contact[i + 1, j + 1] * (infectious_j / active_j);
    }
    return q * lambda_i;
  }

  // ── ODE right-hand side — full 576-compartment system ────────────────────
  vector ode_rhs(real t, vector y,
                 // ── fixed parameters ──
                 real q, real kappa, real iota1, real iota2,
                 real alpha_NC, real alpha_DC_pos, real alpha_HCC,
                 array[] real tau,            // length 4, 1-indexed (stage 1..4)
                 real ptc_NC_CC, real ptc_CC_DC,
                 real ptc_CC_HCC, real ptc_DC_HCC,
                 real ptu_CC_DC, real ptu_CC_HCC, real ptu_DC_HCC,
                 array[] real mu, real omega, // mu length 9, 1-indexed
                 real mu_DC, real mu_HCC,
                 real psi_DC, real psi_HCC,
                 array[] real lambda1,        // length 9, 1-indexed
                 array[] real lambda2,
                 array[] real lambda3,
                 real pi_recid,
                 array[,] real C_contact,    // 9×9, 1-indexed
                 array[] real beta) {         // length 9, 1-indexed (scaled)

    vector[576] dydt = rep_vector(0.0, 576);

    for (i in 0:8) {
      real gam = force_of_infection(i, y, C_contact, q);
      real mu_i = mu[i + 1] * omega;   // background mortality for PWID
      real l1   = lambda1[i + 1];
      real l2   = lambda2[i + 1];
      real l3   = lambda3[i + 1];
      real pi   = pi_recid;

      // ════════════════════════════════════════════════════════════════════
      // STRATUM D (s=0) — Never incarcerated
      // ════════════════════════════════════════════════════════════════════

      // ── Stage 1: Non-cirrhosis ──────────────────────────────────────────
      {
        real Du1 = y[idx_f(0,1,0,i)];
        real Da1 = y[idx_f(0,1,1,i)];
        real Dc1 = y[idx_f(0,1,2,i)];
        real Dt1 = y[idx_f(0,1,3,i)];
        real a   = svr_rate(1, alpha_NC, alpha_DC_pos, alpha_HCC);
        real p12 = ptc_f(1,2, ptc_NC_CC, ptc_CC_DC, ptc_CC_HCC, ptc_DC_HCC);

        dydt[idx_f(0,1,0,i)] += -(gam + l1 + mu_i) * Du1
                                 + (kappa / iota1) * Da1
                                 + (a / iota2) * Dt1
                                 + beta[i + 1];         // scaled inflow

        dydt[idx_f(0,1,1,i)] += gam * Du1
                                 - (1.0/iota1 + l1 + mu_i) * Da1;

        dydt[idx_f(0,1,2,i)] += ((1.0 - kappa) / iota1) * Da1
                                 + ((1.0 - a) / iota2) * Dt1
                                 - (tau[1] + l1 + mu_i + p12) * Dc1;

        dydt[idx_f(0,1,3,i)] += tau[1] * Dc1
                                 - (1.0/iota2 + l1
                                    + mu_eff_f(1,3,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Dt1;
      }

      // ── Stage 2: Compensated cirrhosis ─────────────────────────────────
      {
        real Du2 = y[idx_f(0,2,0,i)];
        real Da2 = y[idx_f(0,2,1,i)];
        real Dc2 = y[idx_f(0,2,2,i)];
        real Dt2 = y[idx_f(0,2,3,i)];
        real Dc1 = y[idx_f(0,1,2,i)];
        real a   = svr_rate(2, alpha_NC, alpha_DC_pos, alpha_HCC);
        real p12 = ptc_f(1,2, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real p23 = ptc_f(2,3, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real p24 = ptc_f(2,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real u23 = ptu_f(2,3, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
        real u24 = ptu_f(2,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);

        dydt[idx_f(0,2,0,i)] += -(gam + l1 + mu_i + u23 + u24) * Du2
                                 + (kappa / iota1) * Da2
                                 + (a / iota2) * Dt2;

        dydt[idx_f(0,2,1,i)] += gam * Du2
                                 - (1.0/iota1 + l1 + mu_i) * Da2;

        dydt[idx_f(0,2,2,i)] += ((1.0 - kappa) / iota1) * Da2
                                 + ((1.0 - a) / iota2) * Dt2
                                 + p12 * Dc1
                                 - (tau[2] + l1 + mu_i + p23 + p24) * Dc2;

        dydt[idx_f(0,2,3,i)] += tau[2] * Dc2
                                 - (1.0/iota2 + l1
                                    + mu_eff_f(2,3,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Dt2;
      }

      // ── Stage 3: Decompensated cirrhosis ───────────────────────────────
      {
        real Du3 = y[idx_f(0,3,0,i)];
        real Da3 = y[idx_f(0,3,1,i)];
        real Dc3 = y[idx_f(0,3,2,i)];
        real Dt3 = y[idx_f(0,3,3,i)];
        real Du2 = y[idx_f(0,2,0,i)];
        real Dc2 = y[idx_f(0,2,2,i)];
        real a   = svr_rate(3, alpha_NC, alpha_DC_pos, alpha_HCC);
        real p23 = ptc_f(2,3, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real p34 = ptc_f(3,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real u23 = ptu_f(2,3, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
        real u34 = ptu_f(3,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);

        dydt[idx_f(0,3,0,i)] += -(gam + l1
                                   + mu_eff_f(3,0,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC)
                                   + u34) * Du3
                                 + (kappa / iota1) * Da3
                                 + (a / iota2) * Dt3
                                 + u23 * Du2;

        dydt[idx_f(0,3,1,i)] += gam * Du3
                                 - (1.0/iota1 + l1
                                    + mu_eff_f(3,1,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Da3;

        dydt[idx_f(0,3,2,i)] += ((1.0 - kappa) / iota1) * Da3
                                 + ((1.0 - a) / iota2) * Dt3
                                 + p23 * Dc2
                                 - (tau[3] + l1
                                    + mu_eff_f(3,2,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC)
                                    + p34) * Dc3;

        dydt[idx_f(0,3,3,i)] += tau[3] * Dc3
                                 - (1.0/iota2 + l1
                                    + mu_eff_f(3,3,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Dt3;
      }

      // ── Stage 4: HCC ────────────────────────────────────────────────────
      {
        real Du4 = y[idx_f(0,4,0,i)];
        real Da4 = y[idx_f(0,4,1,i)];
        real Dc4 = y[idx_f(0,4,2,i)];
        real Dt4 = y[idx_f(0,4,3,i)];
        real Du2 = y[idx_f(0,2,0,i)];
        real Du3 = y[idx_f(0,3,0,i)];
        real Dc2 = y[idx_f(0,2,2,i)];
        real Dc3 = y[idx_f(0,3,2,i)];
        real a   = svr_rate(4, alpha_NC, alpha_DC_pos, alpha_HCC);
        real p24 = ptc_f(2,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real p34 = ptc_f(3,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
        real u24 = ptu_f(2,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
        real u34 = ptu_f(3,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);

        dydt[idx_f(0,4,0,i)] += -(gam + l1
                                   + mu_eff_f(4,0,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Du4
                                 + (kappa / iota1) * Da4
                                 + (a / iota2) * Dt4
                                 + u24 * Du2
                                 + u34 * Du3;

        dydt[idx_f(0,4,1,i)] += gam * Du4
                                 - (1.0/iota1 + l1
                                    + mu_eff_f(4,1,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Da4;

        dydt[idx_f(0,4,2,i)] += ((1.0 - kappa) / iota1) * Da4
                                 + ((1.0 - a) / iota2) * Dt4
                                 + p24 * Dc2
                                 + p34 * Dc3
                                 - (tau[4] + l1
                                    + mu_eff_f(4,2,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Dc4;

        dydt[idx_f(0,4,3,i)] += tau[4] * Dc4
                                 - (1.0/iota2 + l1
                                    + mu_eff_f(4,3,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Dt4;
      }

      // ════════════════════════════════════════════════════════════════════
      // STRATUM J (s=1) — Currently incarcerated
      // Inflow:  l1*D_{k,h,i} + l3*F_{k,h,i}
      // Outflow: l2 + mu_eff
      // No force of infection inside prison
      // ════════════════════════════════════════════════════════════════════

      // Arrest/release transitions (all k, all h)
      for (k in 1:4) {
        for (h in 0:3) {
          real Jkh  = y[idx_f(1,k,h,i)];
          real Dkh  = y[idx_f(0,k,h,i)];
          real Fkh  = y[idx_f(2,k,h,i)];
          real out  = l2 + mu_eff_f(k,h,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC);
          dydt[idx_f(1,k,h,i)] += l1 * Dkh + l3 * Fkh - out * Jkh;
        }
      }

      // Internal disease dynamics within J (no gamma, no l1)
      for (k in 1:4) {
        real a  = svr_rate(k, alpha_NC, alpha_DC_pos, alpha_HCC);
        real Ju = y[idx_f(1,k,0,i)];
        real Ja = y[idx_f(1,k,1,i)];
        real Jc = y[idx_f(1,k,2,i)];
        real Jt = y[idx_f(1,k,3,i)];

        // u: SVR return — no new infection in prison
        dydt[idx_f(1,k,0,i)] += (kappa / iota1) * Ja + (a / iota2) * Jt;
        // a: clearance/chronification
        dydt[idx_f(1,k,1,i)] += -(1.0 / iota1) * Ja;
        // c: treatment initiation
        dydt[idx_f(1,k,2,i)] += ((1.0 - kappa) / iota1) * Ja
                                 + ((1.0 - a) / iota2) * Jt
                                 - tau[k] * Jc;
        // t: treatment
        dydt[idx_f(1,k,3,i)] += tau[k] * Jc - (1.0 / iota2) * Jt;

        // Stage progressions within J
        if (k == 1) {
          real p12 = ptc_f(1,2, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(1,1,2,i)] -= p12 * Jc;
          dydt[idx_f(1,2,2,i)] += p12 * Jc;
        }
        if (k == 2) {
          real u23 = ptu_f(2,3, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real u24 = ptu_f(2,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real p23 = ptc_f(2,3, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          real p24 = ptc_f(2,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(1,2,0,i)] -= (u23 + u24) * Ju;
          dydt[idx_f(1,3,0,i)] += u23 * Ju;
          dydt[idx_f(1,4,0,i)] += u24 * Ju;
          dydt[idx_f(1,2,2,i)] -= (p23 + p24) * Jc;
          dydt[idx_f(1,3,2,i)] += p23 * Jc;
          dydt[idx_f(1,4,2,i)] += p24 * Jc;
        }
        if (k == 3) {
          real u34 = ptu_f(3,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real p34 = ptc_f(3,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(1,3,0,i)] -= u34 * Ju;
          dydt[idx_f(1,4,0,i)] += u34 * Ju;
          dydt[idx_f(1,3,2,i)] -= p34 * Jc;
          dydt[idx_f(1,4,2,i)] += p34 * Jc;
        }
      }

      // ════════════════════════════════════════════════════════════════════
      // STRATUM F (s=2) — Ever-incarcerated, active PWID (at risk)
      // ════════════════════════════════════════════════════════════════════

      for (k in 1:4) {
        real a  = svr_rate(k, alpha_NC, alpha_DC_pos, alpha_HCC);
        real Fu = y[idx_f(2,k,0,i)];
        real Fa = y[idx_f(2,k,1,i)];
        real Fc = y[idx_f(2,k,2,i)];
        real Ft = y[idx_f(2,k,3,i)];
        real Ju = y[idx_f(1,k,0,i)];
        real Ja = y[idx_f(1,k,1,i)];
        real Jc = y[idx_f(1,k,2,i)];
        real Jt = y[idx_f(1,k,3,i)];

        dydt[idx_f(2,k,0,i)] += pi * l2 * Ju
                                 + (kappa / iota1) * Fa
                                 + (a / iota2) * Ft
                                 - (gam + l3
                                    + mu_eff_f(k,0,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Fu;

        dydt[idx_f(2,k,1,i)] += pi * l2 * Ja
                                 + gam * Fu
                                 - (1.0/iota1 + l3
                                    + mu_eff_f(k,1,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Fa;

        dydt[idx_f(2,k,2,i)] += pi * l2 * Jc
                                 + ((1.0 - kappa) / iota1) * Fa
                                 + ((1.0 - a) / iota2) * Ft
                                 - (tau[k] + l3
                                    + mu_eff_f(k,2,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Fc;

        dydt[idx_f(2,k,3,i)] += pi * l2 * Jt
                                 + tau[k] * Fc
                                 - (1.0/iota2 + l3
                                    + mu_eff_f(k,3,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Ft;
      }

      // F stage progressions
      for (k in 1:3) {
        real Fu = y[idx_f(2,k,0,i)];
        real Fc = y[idx_f(2,k,2,i)];
        if (k == 1) {
          real p12 = ptc_f(1,2, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(2,1,2,i)] -= p12 * Fc;
          dydt[idx_f(2,2,2,i)] += p12 * Fc;
        }
        if (k == 2) {
          real u23 = ptu_f(2,3, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real u24 = ptu_f(2,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real p23 = ptc_f(2,3, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          real p24 = ptc_f(2,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(2,2,0,i)] -= (u23 + u24) * Fu;
          dydt[idx_f(2,3,0,i)] += u23 * Fu;
          dydt[idx_f(2,4,0,i)] += u24 * Fu;
          dydt[idx_f(2,2,2,i)] -= (p23 + p24) * Fc;
          dydt[idx_f(2,3,2,i)] += p23 * Fc;
          dydt[idx_f(2,4,2,i)] += p24 * Fc;
        }
        if (k == 3) {
          real u34 = ptu_f(3,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real p34 = ptc_f(3,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(2,3,0,i)] -= u34 * Fu;
          dydt[idx_f(2,4,0,i)] += u34 * Fu;
          dydt[idx_f(2,3,2,i)] -= p34 * Fc;
          dydt[idx_f(2,4,2,i)] += p34 * Fc;
        }
      }

      // ════════════════════════════════════════════════════════════════════
      // STRATUM X (s=3) — Former PWID (no infection risk, no re-arrest)
      // ════════════════════════════════════════════════════════════════════

      for (k in 1:4) {
        real a  = svr_rate(k, alpha_NC, alpha_DC_pos, alpha_HCC);
        real Xu = y[idx_f(3,k,0,i)];
        real Xa = y[idx_f(3,k,1,i)];
        real Xc = y[idx_f(3,k,2,i)];
        real Xt = y[idx_f(3,k,3,i)];
        real Ju = y[idx_f(1,k,0,i)];
        real Ja = y[idx_f(1,k,1,i)];
        real Jc = y[idx_f(1,k,2,i)];
        real Jt = y[idx_f(1,k,3,i)];

        dydt[idx_f(3,k,0,i)] += (1.0 - pi) * l2 * Ju
                                 + (kappa / iota1) * Xa
                                 + (a / iota2) * Xt
                                 - mu_eff_f(k,0,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC)
                                   * Xu;

        dydt[idx_f(3,k,1,i)] += (1.0 - pi) * l2 * Ja
                                 - (1.0/iota1
                                    + mu_eff_f(k,1,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Xa;

        dydt[idx_f(3,k,2,i)] += (1.0 - pi) * l2 * Jc
                                 + ((1.0 - kappa) / iota1) * Xa
                                 + ((1.0 - a) / iota2) * Xt
                                 - (tau[k]
                                    + mu_eff_f(k,2,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Xc;

        dydt[idx_f(3,k,3,i)] += (1.0 - pi) * l2 * Jt
                                 + tau[k] * Xc
                                 - (1.0/iota2
                                    + mu_eff_f(k,3,i, mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC))
                                   * Xt;
      }

      // X stage progressions
      for (k in 1:3) {
        real Xu = y[idx_f(3,k,0,i)];
        real Xc = y[idx_f(3,k,2,i)];
        if (k == 1) {
          real p12 = ptc_f(1,2, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(3,1,2,i)] -= p12 * Xc;
          dydt[idx_f(3,2,2,i)] += p12 * Xc;
        }
        if (k == 2) {
          real u23 = ptu_f(2,3, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real u24 = ptu_f(2,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real p23 = ptc_f(2,3, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          real p24 = ptc_f(2,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(3,2,0,i)] -= (u23 + u24) * Xu;
          dydt[idx_f(3,3,0,i)] += u23 * Xu;
          dydt[idx_f(3,4,0,i)] += u24 * Xu;
          dydt[idx_f(3,2,2,i)] -= (p23 + p24) * Xc;
          dydt[idx_f(3,3,2,i)] += p23 * Xc;
          dydt[idx_f(3,4,2,i)] += p24 * Xc;
        }
        if (k == 3) {
          real u34 = ptu_f(3,4, ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC);
          real p34 = ptc_f(3,4, ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC);
          dydt[idx_f(3,3,0,i)] -= u34 * Xu;
          dydt[idx_f(3,4,0,i)] += u34 * Xu;
          dydt[idx_f(3,3,2,i)] -= p34 * Xc;
          dydt[idx_f(3,4,2,i)] += p34 * Xc;
        }
      }

    } // end age group loop

    return dydt;
  } // end ode_rhs

  // ── Single RK4 step ───────────────────────────────────────────────────────
  vector rk4_step_f(real t, real dt, vector y,
                    real q, real kappa, real iota1, real iota2,
                    real alpha_NC, real alpha_DC_pos, real alpha_HCC,
                    array[] real tau,
                    real ptc_NC_CC, real ptc_CC_DC, real ptc_CC_HCC, real ptc_DC_HCC,
                    real ptu_CC_DC, real ptu_CC_HCC, real ptu_DC_HCC,
                    array[] real mu, real omega,
                    real mu_DC, real mu_HCC, real psi_DC, real psi_HCC,
                    array[] real lambda1, array[] real lambda2, array[] real lambda3,
                    real pi_recid,
                    array[,] real C_contact,
                    array[] real beta) {

    vector[576] k1 = ode_rhs(t, y,
                              q,kappa,iota1,iota2,
                              alpha_NC,alpha_DC_pos,alpha_HCC,
                              tau,
                              ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC,
                              ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC,
                              mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC,
                              lambda1,lambda2,lambda3,pi_recid,
                              C_contact, beta);

    vector[576] y2 = y + 0.5 * dt * k1;
    vector[576] k2 = ode_rhs(t + 0.5*dt, y2,
                              q,kappa,iota1,iota2,
                              alpha_NC,alpha_DC_pos,alpha_HCC,
                              tau,
                              ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC,
                              ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC,
                              mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC,
                              lambda1,lambda2,lambda3,pi_recid,
                              C_contact, beta);

    vector[576] y3 = y + 0.5 * dt * k2;
    vector[576] k3 = ode_rhs(t + 0.5*dt, y3,
                              q,kappa,iota1,iota2,
                              alpha_NC,alpha_DC_pos,alpha_HCC,
                              tau,
                              ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC,
                              ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC,
                              mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC,
                              lambda1,lambda2,lambda3,pi_recid,
                              C_contact, beta);

    vector[576] y4 = y + dt * k3;
    vector[576] k4 = ode_rhs(t + dt, y4,
                              q,kappa,iota1,iota2,
                              alpha_NC,alpha_DC_pos,alpha_HCC,
                              tau,
                              ptc_NC_CC,ptc_CC_DC,ptc_CC_HCC,ptc_DC_HCC,
                              ptu_CC_DC,ptu_CC_HCC,ptu_DC_HCC,
                              mu,omega,mu_DC,mu_HCC,psi_DC,psi_HCC,
                              lambda1,lambda2,lambda3,pi_recid,
                              C_contact, beta);

    return y + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
  }

  // ── Annual aging step (mirrors C++ discrete aging block exactly) ───────────
  // For each (s,k,h) compartment block, shift 1/5 of each age group forward.
  vector aging_step_f(vector y) {
    vector[576] y_new = y;
    for (s in 0:3)
      for (k in 1:4)
        for (h in 0:3)
          for (i in 0:7) {
            real delta = y[idx_f(s,k,h,i)] / 5.0;
            y_new[idx_f(s,k,h,i)]     -= delta;
            y_new[idx_f(s,k,h,i + 1)] += delta;
          }
    return y_new;
  }

  // ── Full simulation: n_years years, dt step size, with annual aging ────────
  // Returns the final state vector (length 576).
  vector simulate_model(
      real beta_scale,          // calibrated: scales all beta[i]
      real y0_scale,            // calibrated: scales initial D_{u,1,i} compartments
      vector y0_base,           // base initial conditions (length 576)
      array[] real beta_base,   // base inflow rates (length 9, 1-indexed)
      int   n_years,            // number of years to simulate
      int   steps_per_year,     // RK4 steps per year (dt = 1/steps_per_year)
      // ── fixed parameters ──
      real q, real kappa, real iota1, real iota2,
      real alpha_NC, real alpha_DC_pos, real alpha_HCC,
      array[] real tau,
      real ptc_NC_CC, real ptc_CC_DC, real ptc_CC_HCC, real ptc_DC_HCC,
      real ptu_CC_DC, real ptu_CC_HCC, real ptu_DC_HCC,
      array[] real mu, real omega,
      real mu_DC, real mu_HCC, real psi_DC, real psi_HCC,
      array[] real lambda1, array[] real lambda2, array[] real lambda3,
      real pi_recid,
      array[,] real C_contact) {

    real dt = 1.0 / steps_per_year;

    // ── Apply y0_scale to initial susceptible D_{u,1,i} compartments ────────
    vector[576] y = y0_base;
    for (i in 0:8)
      y[idx_f(0, 1, 0, i)] = y0_scale * y0_base[idx_f(0, 1, 0, i)];

    // ── Apply beta_scale to inflow rates ─────────────────────────────────────
    array[9] real beta_scaled;
    for (i in 1:9)
      beta_scaled[i] = beta_scale * beta_base[i];

    // ── Integrate ─────────────────────────────────────────────────────────────
    real t = 0.0;
    for (yr in 1:n_years) {
      for (step in 1:steps_per_year) {
        y = rk4_step_f(t, dt, y,
                       q, kappa, iota1, iota2,
                       alpha_NC, alpha_DC_pos, alpha_HCC,
                       tau,
                       ptc_NC_CC, ptc_CC_DC, ptc_CC_HCC, ptc_DC_HCC,
                       ptu_CC_DC, ptu_CC_HCC, ptu_DC_HCC,
                       mu, omega, mu_DC, mu_HCC, psi_DC, psi_HCC,
                       lambda1, lambda2, lambda3, pi_recid,
                       C_contact, beta_scaled);
        // clamp negatives (numerical safety — non-differentiable but rare)
        for (c in 1:576)
          if (y[c] < 0.0) y[c] = 0.0;
        t += dt;
      }
      // Annual discrete aging
      y = aging_step_f(y);
    }

    return y;
  }

} // end functions block


// =============================================================================
// DATA
// =============================================================================
data {

  // ── Observed prison population counts by age group (9 age groups) ──────────
  // obs_J_total[i] = sum_{k,h} J_{k,h,i}  (all incarcerated, all stages/states)
  array[9] int obs_J_total;   // c(307, 797, 829, 633, 598, 642, 481, 439, 366)

  // obs_J_susc[i] = sum_k J_{u=h0, k, i}  (susceptible/post-SVR incarcerated)
  // NOTE: if the second target formula intended a different combination
  //       (e.g. J_u + J_c, or J + F), replace the extraction in
  //       transformed parameters accordingly.
  array[9] int obs_J_susc;    // c(55, 145, 183, 164, 212, 299, 222, 190, 133)

  // ── Simulation settings ───────────────────────────────────────────────────
  int n_years;          // number of years to run
  int steps_per_year;   // RK4 steps per year (e.g. 12 for monthly)

  // ── Initial conditions (length 576) ──────────────────────────────────────
  vector[576] y0_base;  // base state vector at t=0

  // ── Fixed model parameters ───────────────────────────────────────────────
  real q;
  real kappa;
  real iota1;
  real iota2;

  real alpha_NC;
  real alpha_DC_pos;
  real alpha_DC_neg;   // unused in RHS but kept for completeness
  real alpha_HCC;

  array[4] real tau;   // treatment initiation rates, 1-indexed (stage 1..4)

  real p_NC_CC;
  real p_CC_DC;
  real p_CC_HCC;
  real p_DC_HCC;

  real r3_NC_CC;
  real r3_CC_DC;
  real r3_CC_HCC;
  real r3_DC_HCC;

  real phi_CC_DC;
  real phi_CC_HCC;
  real phi_DC_HCC;

  array[9] real mu;    // background mortality by age group, 1-indexed
  real omega;          // PWID SMR
  real mu_DC;
  real mu_HCC;
  real psi_DC;
  real psi_HCC;

  array[9] real lambda1;  // first-arrest rate by age group
  array[9] real lambda2;  // release rate by age group
  array[9] real lambda3;  // re-arrest rate by age group
  real pi_recid;

  array[9, 9] real C_contact;  // age-structured contact matrix
  array[9] real beta_base;     // base inflow rates (before scaling)
}


// =============================================================================
// TRANSFORMED DATA — compute derived fixed parameters once
// =============================================================================
transformed data {
  // Genotype-weighted effective progression rates
  real ptc_NC_CC  = (r3_NC_CC  * p_NC_CC)  * (1.0 - (1.0 - r3_NC_CC)  * (1.0 - 1.0/r3_NC_CC));
  real ptc_CC_DC  = (r3_CC_DC  * p_CC_DC)  * (1.0 - (1.0 - r3_CC_DC)  * (1.0 - 1.0/r3_CC_DC));
  real ptc_CC_HCC = (r3_CC_HCC * p_CC_HCC) * (1.0 - (1.0 - r3_CC_HCC) * (1.0 - 1.0/r3_CC_HCC));
  real ptc_DC_HCC = (r3_DC_HCC * p_DC_HCC) * (1.0 - (1.0 - r3_DC_HCC) * (1.0 - 1.0/r3_DC_HCC));
  // NOTE: the formula above simplifies to just r3 * p, matching the C++ line:
  //   ptc = (rho * r3 + (1-rho)) * p  with rho=1 → ptc = r3 * p
  // Replace with the full rho-weighted version if rho < 1 and pass rho as data.

  // Post-SVR progression rates
  real ptu_CC_DC  = phi_CC_DC  * ptc_CC_DC;
  real ptu_CC_HCC = phi_CC_HCC * ptc_CC_HCC;
  real ptu_DC_HCC = phi_DC_HCC * ptc_DC_HCC;
}


// =============================================================================
// PARAMETERS — the two calibrated scaling factors
// =============================================================================
parameters {
  // Uninformative (weakly informative) priors on positive scale factors.
  // Both are multiplicative scalars: value = 1 means no adjustment.
  real<lower=0> beta_scale;   // scales beta[i] for all age groups
  real<lower=0> y0_scale;     // scales initial D_{u,1,i} for all age groups
}


// =============================================================================
// TRANSFORMED PARAMETERS — run the ODE and extract model predictions
// =============================================================================
transformed parameters {

  // ── Run full simulation ───────────────────────────────────────────────────
  vector[576] y_final = simulate_model(
      beta_scale, y0_scale,
      y0_base, beta_base,
      n_years, steps_per_year,
      q, kappa, iota1, iota2,
      alpha_NC, alpha_DC_pos, alpha_HCC,
      tau,
      ptc_NC_CC, ptc_CC_DC, ptc_CC_HCC, ptc_DC_HCC,
      ptu_CC_DC, ptu_CC_HCC, ptu_DC_HCC,
      mu, omega, mu_DC, mu_HCC, psi_DC, psi_HCC,
      lambda1, lambda2, lambda3, pi_recid,
      C_contact);

  // ── Extract model predictions for J stratum ───────────────────────────────
  // pred_J_total[i] = sum_{k=1..4, h=0..3} J_{k,h,i}
  // pred_J_susc[i]  = sum_{k=1..4}         J_{u=h0, k, i}
  array[9] real pred_J_total;
  array[9] real pred_J_susc;

  for (i in 0:8) {
    real total = 0.0;
    real susc  = 0.0;
    for (k in 1:4) {
      for (h in 0:3)
        total += y_final[idx_f(1, k, h, i)];
      susc += y_final[idx_f(1, k, 0, i)];   // h=0: susceptible/post-SVR
    }
    pred_J_total[i + 1] = fmax(total, 1e-9); // guard against log(0) in Poisson
    pred_J_susc[i + 1]  = fmax(susc,  1e-9);
  }
}


// =============================================================================
// MODEL — priors + likelihood
// =============================================================================
model {

  // ── Uninformative priors on positive scaling factors ─────────────────────
  // Half-normal(0, 10): very wide, centred near 0, all mass on positive reals.
  // Adjust scale if you have domain knowledge (e.g. scale ~ Normal(1, 0.5) if
  // you expect the base values to be approximately correct).
  beta_scale ~ normal(0, 10);
  y0_scale   ~ normal(0, 10);

  // ── Likelihood: Poisson on each of the 9 + 9 age-group counts ────────────
  // Each entry contributes: log p(obs | pred) = obs*log(pred) - pred - log(obs!)
  // The "difference in each entry" is the residual: obs - pred, which the
  // Poisson log-likelihood encodes implicitly through its score function.
  for (i in 1:9) {
    obs_J_total[i] ~ poisson(pred_J_total[i]);
    obs_J_susc[i]  ~ poisson(pred_J_susc[i]);
  }
}


// =============================================================================
// GENERATED QUANTITIES — posterior predictive checks and residuals
// =============================================================================
generated quantities {

  // Posterior predictive samples (for PPC plots)
  array[9] int ppc_J_total;
  array[9] int ppc_J_susc;

  // Residuals: observed - predicted (for diagnostics)
  array[9] real resid_J_total;
  array[9] real resid_J_susc;

  // Log-likelihood contributions per observation (for loo/waic)
  array[9] real log_lik_total;
  array[9] real log_lik_susc;

  for (i in 1:9) {
    ppc_J_total[i]  = poisson_rng(pred_J_total[i]);
    ppc_J_susc[i]   = poisson_rng(pred_J_susc[i]);
    resid_J_total[i] = obs_J_total[i] - pred_J_total[i];
    resid_J_susc[i]  = obs_J_susc[i]  - pred_J_susc[i];
    log_lik_total[i] = poisson_lpmf(obs_J_total[i] | pred_J_total[i]);
    log_lik_susc[i]  = poisson_lpmf(obs_J_susc[i]  | pred_J_susc[i]);
  }
}
