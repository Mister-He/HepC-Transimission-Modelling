# =============================================================================
# facts.py — every number the manuscript and S1 Text use, read from saved outputs (no
# simulation). Each fact keeps the output file it came from, so the build can bind every
# sentence to its evidence. `python3 scripts/paper/facts.py [ROOT]` prints them.
# =============================================================================
import csv, json, math, os, re, subprocess, sys

F = float


class Facts:
    def __init__(self, root):
        self.root = root
        r = lambda p: list(csv.DictReader(open(os.path.join(root, p))))
        j = lambda p: json.load(open(os.path.join(root, p)))
        S = "output/scenarios_2026/"
        C = "output/calibration/"
        self.files = {  # output key -> repository path (evidence sources)
            "summary": S + "scenario_2030_summary.csv", "traj": S + "scenario_trajectories.csv",
            "notx": S + "no_treatment_trajectory.csv", "annual": S + "scenario_annual_services_costs.csv",
            "costev": S + "cost_evidence.csv", "config": "data/scenario_policy_config.json",
            "draws_valid": S + "draw_validity.txt", "targets": "data/calibration_targets.json",
            "cnb": "data/singstat_new_drug_abusers.csv", "prison": "output/prison_2015_counts.csv",
            "setup": "scripts/model/setup.R", "sim": "scripts/model/sim.cpp",
            "likelihood": "scripts/calibration/likelihood.R", "equilibrium": "scripts/calibration/equilibrium.R",
            "scen_code": "scripts/scenarios_2026.R", "rates": "output/rearrest_rate_fit.rds",
            "release": "output/release_rate_fit.rds", "nm": C + "r_20261027_agepop/solutions.csv",
            "nm_cfg": C + "r_20261027_agepop/run_config.json", "nm_pred": C + "r_20261027_agepop/predictions.csv",
            "mcmc_cfg": C + "mcmc_reference/mcmc_run_config.json", "mcmc_diag": C + "mcmc_reference/mcmc_diagnostics.csv",
            "mcmc_chain": C + "mcmc_reference/mcmc_chain_summary.csv", "mcmc_par": C + "mcmc_reference/mcmc_parameter_summary.csv",
            "strata": C + "mcmc_reference/model_prediction_2015_strata.csv", "overlap": "figures/calibration_interval_overlap.csv",
            "mcmc_sbc": C + "npe_vs_mcmc/mcmc_sbc/mcmc_sbc_summary.csv", "npe_sim": C + "npe_lhs/simulation_summary.json",
            "npe_sbc": C + "npe_lhs/sbc_summary.csv", "npe_acc": C + "npe_vs_mcmc/accuracy_summary.csv",
            "npe_pol": C + "npe_vs_mcmc/policy_comparison.csv", "npe_support": C + "npe_vs_mcmc/npe_support.txt",
            "npe_cfg": C + "npe_lhs/npe_config.json", "stratum_tot": S + "scenario_stratum_totals.csv",
            "sps": "data/sps.csv and output/cleaned_data.rds (record-level, not distributed)",
            "stan_release": "scripts/model/release_rate.stan", "stan_rearrest": "scripts/model/rearrest_rate.stan",
            "preprocess": "scripts/model/preprocess.R", "figures": "scripts/figures.py",
        }
        self.summ = {x["scenario_id"]: x for x in r(self.files["summary"])}
        self.traj = r(self.files["traj"]); self.notx = r(self.files["notx"]); self.ann = r(self.files["annual"])
        self.costev = r(self.files["costev"]); self.cfg = j(self.files["config"]); self.tg = j(self.files["targets"])
        self.cnb = [x for x in r(self.files["cnb"]) if x["year"] == "2015"]
        self.prison = r(self.files["prison"]); self.nm = r(self.files["nm"]); self.nm_cfg = j(self.files["nm_cfg"])
        self.nm_pred = r(self.files["nm_pred"]); self.mcfg = j(self.files["mcmc_cfg"])
        self.mdiag = [x for x in r(self.files["mcmc_diag"]) if x["subset"] == "all_retained"]
        self.mchain = r(self.files["mcmc_chain"]); self.mpar = {x["parameter"]: x for x in r(self.files["mcmc_par"])}
        self.strata = {x["measure"]: x for x in r(self.files["strata"])}; self.ovl = r(self.files["overlap"])
        self.msbc = {x["parameter"]: x for x in r(self.files["mcmc_sbc"])}; self.npe_sim = j(self.files["npe_sim"])
        self.npe_sbc = r(self.files["npe_sbc"]); self.npe_acc = {x["run"]: x for x in r(self.files["npe_acc"])}
        self.npe_pol = r(self.files["npe_pol"]); self.npe_cfg = j(self.files["npe_cfg"])
        self.npe_support = dict(l.split() for l in open(os.path.join(root, self.files["npe_support"])))
        self.setup_txt = open(os.path.join(root, self.files["setup"])).read()
        self.valid = dict(l.split() for l in open(os.path.join(root, self.files["draws_valid"])))
        self._rates()
        self._derive()

    # ---- helpers ----
    def setup(self, name):
        v = re.search(rf"^\s*{name}\s*=\s*(c\([^)]*\)|[0-9.eE+-]+)", self.setup_txt, re.M).group(1)
        x = [F(a) for a in re.findall(r"[0-9.]+(?:[eE][+-]?[0-9]+)?", v)]
        return x if v.startswith("c(") else x[0]

    def tr(self, sid, yr, rows=None):
        return next(x for x in (rows or self.traj) if x["scenario_id"] == sid and int(F(x["year"])) == yr)

    def an(self, sid, yr, st):
        return next(x for x in self.ann if x["scenario_id"] == sid and int(F(x["year"])) == yr and x["stratum"] == st)

    def _rates(self):
        # Posterior mean and 95% interval of the Stan rate fits (aggregate summaries only).
        code = ('suppressMessages(library(posterior)); g<-function(p) as_draws_df(readRDS(p)[[2]]); '
                'r<-g("output/release_rate_fit.rds"); a<-g("output/rearrest_rate_fit.rds"); '
                'q<-function(v) sprintf("%.6f,%.6f,%.6f", mean(v), quantile(v,.025), quantile(v,.975)); '
                'for(i in 1:6) cat("lambda2",i,q(exp(r[[sprintf("log_lambda[%d]",i)]])),"\\n"); '
                'for(i in 1:6) cat("lambda3",i,q(exp(a[[sprintf("log_lambda[%d]",i)]])),"\\n"); '
                'cat("pi",1,q(plogis(a$pi_raw)),"\\n"); '
                'd<-readRDS("output/cleaned_data.rds"); cat("sps", length(unique(d$IDENTIFER)), nrow(d), format(min(d$start), "%Y"), format(max(d$start), "%Y"), "\\n"); '
                'for (f in c("output/release_rate_fit.rds","output/rearrest_rate_fit.rds")) '
                'cat("rhat",f,max(readRDS(f)[[1]]$rhat),"\\n")')
        out = subprocess.run(["Rscript", "-e", code], cwd=self.root, capture_output=True, text=True, check=True).stdout
        self.rate_ci = {}; self.rate_rhat = []
        for line in out.splitlines():
            p = line.split()
            if p[0] == "rhat": self.rate_rhat.append(F(p[2])); continue
            if p[0] == "sps":  # aggregate counts only (people, stays, first and last admission year)
                self.sps_people, self.sps_stays, self.sps_first, self.sps_last = int(p[1]), int(p[2]), p[3], p[4]; continue
            self.rate_ci[(p[0], int(p[1]))] = [F(v) for v in p[2].split(",")]

    def _derive(self):
        s = self.summ; self.ids = list(s)
        self.finite = [k for k in s if s[k]["scenario_type"] == "finite"]
        self.inc_pass = [k for k in self.finite if s[k]["who_incidence_pass"] == "TRUE"]
        self.mort_pass = [k for k in self.finite if s[k]["who_mortality_pass"] == "TRUE"]
        self.ndraw = int(F(s["J_only"]["n_posterior_draws"]))
        cF, cX = self.an("J_only", 2026, "F"), self.an("J_only", 2026, "X")
        self.f_eligible = F(cF["screening_eligible"]); self.x_eligible = F(cX["screening_eligible"])
        self.f_chronic = F(cF["chronic_stock_start"]); self.x_chronic = F(cX["chronic_stock_start"])
        self.x_chronic_share = self.x_chronic / (self.x_chronic + self.f_chronic)
        j26, j30 = self.tr("J_only", 2026), self.tr("J_only", 2030)
        self.fx_death_share_2026 = F(j26["fx_mortality"]) / F(j26["mortality"])
        self.fx_death_share_2030 = F(j30["fx_mortality"]) / F(j30["mortality"])
        self.d_inc_share_sinf = 1 - F(s["SInf"]["fx_incidence"]) / F(s["SInf"]["incidence"])
        self.d_inc_share_j = 1 - F(s["J_only"]["fx_incidence"]) / F(s["J_only"]["incidence"])
        b = self.tr("J_only", 2015)
        self.base_inc, self.base_mort = F(b["incidence"]), F(b["mortality"])
        self.base_fx_death_share = F(b["fx_mortality"]) / F(b["mortality"])
        quota = lambda sid, st, c: sum(F(x[c]) for x in self.ann if x["scenario_id"] == sid and x["stratum"] == st)
        self.quota_use = {k: {st: quota(k, st, "actual_screened") / quota(k, st, "screening_quota")
                              for st in ("F", "X")} for k in self.finite}
        rr = [x for x in self.ann if x["scenario_id"] == "S5000_F25"]
        self.drug_share_5000 = sum(F(x["discounted_cost_daa_sgd"]) for x in rr) / sum(F(x["discounted_cost_total_sgd"]) for x in rr)
        self.cnb_counts = [int(x["new_drug_abusers"]) for x in self.cnb]
        self.j_obs_2015 = sum(int(x["n_in_prison"]) for x in self.prison)
        self.screened_2015 = sum(int(x["screened"]) for x in self.prison)
        self.positive_2015 = sum(int(x["positive"]) for x in self.prison)
        # Seroprevalence targets against positive/screened in the current preprocessing output.
        self.sero_recomputed = [int(x["positive"]) / int(x["screened"]) for x in self.prison]
        self.sero_max_diff = max(abs(a - b) for a, b in zip(self.sero_recomputed, self.tg["seroprevalence"]))
        nm = sorted(self.nm, key=lambda x: F(x["objective"]))
        self.nm_best = nm[0]
        self.nm_same = [x for x in nm if F(x["objective"]) - F(nm[0]["objective"]) < 1]
        self.nm_other = [x for x in nm if x not in self.nm_same]
        self.c_post = [math.exp(F(self.mpar["theta13"][k])) for k in ("median", "q025", "q975")]
        self.entry_post = {i: [math.exp(F(self.mpar[f"theta{i}"][k])) for k in ("median", "q025", "q975")] for i in range(7, 13)}
        self.contact_post = {i: [math.exp(F(self.mpar[f"theta{i}"][k])) for k in ("median", "q025", "q975")] for i in range(1, 7)}
        self.prior_tsd = [trunc_sd(m, sd, lo, hi) for m, sd, lo, hi in
                          zip(self.mcfg["prior_mean"], self.mcfg["prior_sd"], self.mcfg["theta_lo"], self.mcfg["theta_hi"])]
        self.sd_ratio = [F(self.mpar[f"theta{i}"]["sd"]) / self.prior_tsd[i - 1] for i in range(1, 14)]
        self.prior_dominated = [i for i in range(1, 14) if self.sd_ratio[i - 1] > 0.7]
        self.rhat_max = max(F(x["rhat"]) for x in self.mdiag)
        self.ess_bulk_min = min(F(x["ess_bulk"]) for x in self.mdiag)
        self.ess_tail_min = min(F(x["ess_tail"]) for x in self.mdiag)
        self.mcse_max = max(F(x["mcse_mean"]) / F(x["sd"]) for x in self.mdiag)
        self.acc = [F(x["acceptance"]) for x in self.mchain]
        self.n_ovl = sum(x["MCMC_overlap"] == "TRUE" for x in self.ovl)
        self.npe_bad_shape = [x["parameter"] for x in self.npe_sbc if x["rank_shape"] != "uniform"]
        self.npe_pass_eq = int(self.npe_support["npe_draws_equilibrium_pass"]); self.npe_checked = int(self.npe_support["npe_draws_checked"])
        self.inc_order_ok = all(
            [F(s[f"S{cap}_F{x}"]["incidence_reduction"]) for x in (25, 50, 75)] ==
            sorted(F(s[f"S{cap}_F{x}"]["incidence_reduction"]) for x in (25, 50, 75)) for cap in self.cfg["annual_screening_totals"])
        self.mort_order_ok = all(
            [F(s[f"S{cap}_F{x}"]["mortality_reduction"]) for x in (25, 50, 75)] ==
            sorted((F(s[f"S{cap}_F{x}"]["mortality_reduction"]) for x in (25, 50, 75)), reverse=True)
            for cap in self.cfg["annual_screening_totals"])
        self.tau_j = -math.log1p(-self.cfg["j_annual_start_probability"])
        self.tau_fx = -math.log1p(-self.cfg["fx_annual_start_probability"])


def trunc_sd(m, sd, lo, hi):  # SD of Normal(m, sd) truncated to [lo, hi]
    phi = lambda z: math.exp(-z * z / 2) / math.sqrt(2 * math.pi)
    Phi = lambda z: 0.5 * (1 + math.erf(z / math.sqrt(2)))
    a, b = (lo - m) / sd, (hi - m) / sd; Z = Phi(b) - Phi(a)
    return sd * math.sqrt(1 + (a * phi(a) - b * phi(b)) / Z - ((phi(a) - phi(b)) / Z) ** 2)


if __name__ == "__main__":
    f = Facts(os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "."))
    for k in ["ndraw", "inc_pass", "mort_pass", "f_eligible", "x_eligible", "f_chronic", "x_chronic", "x_chronic_share",
              "fx_death_share_2026", "fx_death_share_2030", "d_inc_share_sinf", "d_inc_share_j", "base_inc", "base_mort",
              "base_fx_death_share", "drug_share_5000", "j_obs_2015", "screened_2015", "positive_2015", "sero_max_diff",
              "c_post", "prior_dominated", "sd_ratio", "rhat_max", "ess_bulk_min", "ess_tail_min", "mcse_max", "acc", "n_ovl",
              "npe_bad_shape", "npe_pass_eq", "npe_checked", "inc_order_ok", "mort_order_ok", "tau_j", "tau_fx", "rate_ci",
              "rate_rhat", "entry_post", "contact_post"]:
        print(k, getattr(f, k))
    print("nm", [(x["start_id"], round(F(x["objective"]), 2)) for x in f.nm])
    print("quota_use", {k: {a: round(b, 4) for a, b in v.items()} for k, v in f.quota_use.items()})
    for k, x in f.summ.items():
        print(k, x["incidence_reduction_text"], x["mortality_reduction_text"], x["joint_pass_probability"],
              x["total_cost_text"], x["incremental_cost_text"], x["fx_starts_text"], x["j_starts_text"])
    for y in (2019, 2025, 2026, 2030):
        t = f.tr("J_only", y); n = f.tr("no_treatment", y, f.notx)
        print("J", y, t["incidence_ratio_2015"], t["mortality_ratio_2015"], "notx", n["incidence_ratio_2015"], n["incidence_ratio_2015_lo"], n["incidence_ratio_2015_hi"], n["mortality_ratio_2015"])
