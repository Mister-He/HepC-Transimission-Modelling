#!/usr/bin/env python3
# =============================================================================
# figures.py — redraws every project figure from saved outputs with the
# scientific-visualization skill (style "default": 8 pt sans, Okabe-Ito colours;
# export_figure: PNG at 300 dpi + SVG with editable text). Widths follow the
# skill's PLOS snapshot (full 190.5 mm, text column 132 mm).
# Inputs: output/scenarios_2026/*.csv (scenarios_2026.R), output/calibration/
# {mcmc_reference,npe_lhs,npe_vs_mcmc}/ (run_mcmc_reference.R, run_npe.R,
# plot_calibration.R, predict_2015_posterior.R, compare_npe_mcmc.R),
# data/calibration_targets.json. Nothing is re-simulated.
# Usage: python scripts/figures.py [ROOT] [--skill-dir DIR] [--only NAME ...]
# =============================================================================
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--skill-dir", default=str(Path.home() / ".claude/skills/scientific-visualization"))
ap.add_argument("--only", nargs="*", help="figure names (file stems) to draw; default all")
args = ap.parse_args()
ROOT = Path(args.root).resolve()
SKILL = Path(args.skill_dir).expanduser()
if not (SKILL / "scripts/figure_export.py").exists():
    sys.exit(f"scientific-visualization skill not found at {SKILL} (pass --skill-dir)")
sys.path[:0] = [str(SKILL / "scripts"), str(SKILL / "assets")]

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch  # noqa: E402
from matplotlib.ticker import FuncFormatter, PercentFormatter  # noqa: E402

from figure_export import export_figure  # noqa: E402
from style_presets import style_context  # noqa: E402

SC = ROOT / "output/scenarios_2026"
CAL = ROOT / "output/calibration"
MM = 1 / 25.4
FULL, COL = 190.5, 132.0  # PLOS widths (mm), skill snapshot accessed 2026-10-01
STRATUM = {"D": "#595959", "J": "#0072B2", "F": "#D55E00", "X": "#009E73"}
ITEM = {"SOF/VEL": "#0072B2", "HCV antibody": "#D55E00", "HCV RNA": "#009E73", "Consultation": "#CC79A7"}
ITEM_COLS = {"SOF/VEL": "discounted_cost_daa_sgd", "HCV antibody": "discounted_cost_antibody_sgd",
             "HCV RNA": "discounted_cost_rna_sgd", "Consultation": "discounted_cost_consultation_sgd"}
STATUS = {True: "#0072B2", False: "#D55E00"}  # pass / fail; always paired with text or a threshold line
SHARE_LS = {0.25: "-", 0.5: (0, (4, 2)), 0.75: (0, (1, 1.5))}
CAP_MARK = ["o", "s", "^", "D", "v"]
pct = PercentFormatter(1.0, decimals=0)
comma = FuncFormatter(lambda v, _: f"{v:,.0f}")
TH = [f"theta{i}" for i in range(1, 8)]  # six log contact scales + log c
written = []


def want(name):
    return not args.only or name in args.only


def save(fig, path, sources):
    export_figure(fig, path, formats=["png", "svg"], dpi=300, overwrite=True,
                  provenance={"script": "scripts/figures.py", "sources": sources})
    plt.close(fig)
    written.append(str(Path(path).relative_to(ROOT)))


def panel_label(ax, text):
    ax.set_title(text, loc="left", fontweight="bold")


# ---------------------------------------------------------------- shared data
defs = pd.read_csv(SC / "scenario_definitions.csv")
summ = pd.read_csv(SC / "scenario_2030_summary.csv")
caps = sorted(defs.loc[defs.scenario_type == "finite", "screening_total"].astype(float).astype(int).unique())
CAP_COL = dict(zip(caps, plt.get_cmap("plasma")(np.linspace(0.0, 0.7, len(caps)))))
order = list(defs.scenario_id)
label = dict(zip(defs.scenario_id, defs.scenario_label))
label.update(J_only="J only [reference]", SInf="No cap [reference]")
finite_ids = list(defs.loc[defs.scenario_type == "finite", "scenario_id"])
meta = defs.set_index("scenario_id")


def cap_legend(ax_or_fig, **kw):
    h = [Line2D([], [], color=CAP_COL[c], marker=CAP_MARK[i], label=f"{c:,}") for i, c in enumerate(caps)]
    return ax_or_fig.legend(handles=h, title="Screens per year", ncols=len(caps), **kw)


# ---------------------------------------------------------------- Fig 1 schematic
if want("paper_fig1_schematic"):
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(FULL * MM, 95 * MM), layout="constrained")
        ax.set(xlim=(-0.4, 9.1), ylim=(-0.2, 4.8)); ax.axis("off")
        boxes = {"D": (1, 2, "Active PWID,\nnever incarcerated"), "J": (4, 2, "Currently\nincarcerated"),
                 "F": (7.6, 3.4, "Released,\nstill injecting"), "X": (7.6, 0.6, "Released,\nceased injecting")}
        for k, (x, y, txt) in boxes.items():
            ax.add_patch(FancyBboxPatch((x - 1.05, y - 0.55), 2.1, 1.1, boxstyle="round,pad=0.02",
                                        fc=matplotlib.colors.to_rgba(STRATUM[k], 0.15), ec=STRATUM[k], lw=1.5))
            ax.text(x, y + 0.22, k, ha="center", va="center", fontsize=11, fontweight="bold", color=STRATUM[k])
            ax.text(x, y - 0.15, txt, ha="center", va="center", fontsize=8)

        def arrow(a, b, rad=0.0):
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=10, lw=0.9,
                                         color="0.2", connectionstyle=f"arc3,rad={rad}"))
        arrow((2.1, 2), (2.9, 2)); arrow((5.1, 2.35), (6.5, 3.15)); arrow((5.1, 1.65), (6.5, 0.85))
        arrow((6.5, 3.75), (4.3, 2.6), rad=0.35)
        ax.text(2.5, 2.2, "First arrest", ha="center", va="bottom")
        ax.text(5.7, 3.0, "Release,\nreturn to injecting", ha="right", va="center")
        ax.text(5.7, 1.0, "Release,\nstop injecting", ha="right", va="center")
        ax.text(4.9, 4.05, "Re-arrest", ha="center")
        note = dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", lw=0.6)
        ax.text(4, 0.3, "Treatment in J from 2019", ha="center", color=STRATUM["J"], fontweight="bold", bbox=note)
        ax.text(7.6, 4.5, "Screening quotas for F and X from 2026", ha="center", color=STRATUM["F"],
                fontweight="bold", bbox=note)
        ax.text(1, 3.2, "Transmission in D and F only", ha="center", color="0.15", bbox=note)
        save(fig, SC / "paper_fig1_schematic", ["model structure in scripts/model/sim.cpp"])

# ---------------------------------------------------------------- Fig 2 trajectories
if want("paper_fig2_trajectories"):
    tr = pd.read_csv(SC / "scenario_trajectories.csv")
    nt = pd.read_csv(SC / "no_treatment_trajectory.csv")
    pick = {"no_treatment": ("No treatment", "0.35", (0, (1, 1.5))),
            "J_only": ("Prison treatment only (J)", STRATUM["J"], "-"),
            "S8000_F25": ("J + 8,000 screens/yr, F 25%", STRATUM["F"], (0, (4, 2))),
            "S8000_F75": ("J + 8,000 screens/yr, F 75%", STRATUM["X"], (0, (6, 1.5, 1, 1.5)))}
    d = pd.concat([nt, tr])
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(1, 2, figsize=(FULL * MM, 85 * MM), layout="constrained", sharey=True)
        for ax, (m, title, thr) in zip(axes, [("incidence_ratio_2015", "A  HCV incidence", 0.2),
                                              ("mortality_ratio_2015", "B  HCV-related mortality", 0.35)]):
            for sid, (lab, col, ls) in pick.items():
                s = d[d.scenario_id == sid]
                if sid.startswith("S8000"):  # screening scenarios share the J-only history
                    s = s[s.year >= 2025]
                ax.fill_between(s.year, s[m + "_lo"], s[m + "_hi"], color=col, alpha=0.12, lw=0)
                ax.plot(s.year, s[m], color=col, ls=ls, lw=1.3, label=lab)
            ax.axhline(thr, color="black", ls="--", lw=0.7)
            ax.text(2015.2, thr, f"2030 target ({thr:.0%} of 2015)", va="bottom", fontsize=7)
            for x in (2019, 2026):
                ax.axvline(x, color="0.6", ls=":", lw=0.7)
            ax.set(xticks=[2015, 2019, 2026, 2030], ylim=(0, 1.4), xlabel="Calendar year")
            ax.yaxis.set_major_formatter(pct); panel_label(ax, title)
        axes[0].set_ylabel("Annual events relative to 2015")
        fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncols=2)
        save(fig, SC / "paper_fig2_trajectories", ["scenario_trajectories.csv", "no_treatment_trajectory.csv"])

# ---------------------------------------------------------------- Fig 3 allocation
fin = summ[summ.scenario_type == "finite"].copy()
fin["cap"] = fin.screening_total.astype(float).astype(int)
if want("paper_fig3_allocation"):
    ref = summ.set_index("scenario_id")
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(1, 2, figsize=(FULL * MM, 90 * MM), layout="constrained")
        for ax, (m, title, thr) in zip(axes, [("incidence_reduction", "A  Incidence reduction, 2030", 0.8),
                                              ("mortality_reduction", "B  Mortality reduction, 2030", 0.65)]):
            ax.axhline(thr, color="black", ls="--", lw=0.7)
            ax.text(23, thr, f"Target {thr:.0%}", va="bottom", fontsize=7)
            ax.axhline(ref.loc["J_only", m], color="0.4", ls=":", lw=0.8)
            ax.axhline(ref.loc["SInf", m], color="0.4", ls="-", lw=0.8)
            for i, c in enumerate(caps):
                s = fin[fin.cap == c].sort_values("f_share")
                x = 100 * s.f_share + (i - 2) * 2.0  # dodge
                ax.plot(x, s[m], color=CAP_COL[c], lw=0.8)
                ax.errorbar(x, s[m], yerr=[s[m] - s[m + "_lo"], s[m + "_hi"] - s[m]], fmt=CAP_MARK[i],
                            color=CAP_COL[c], ms=4, elinewidth=0.8, capsize=0)
            ax.set(xticks=[25, 50, 75], xlabel="Share of screening quota allocated to F",
                   ylabel="Reduction versus 2015")
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
            ax.yaxis.set_major_formatter(pct); panel_label(ax, title)
        refs = [Line2D([], [], color="0.4", ls=":", label="J only"), Line2D([], [], color="0.4", label="No screening cap")]
        cap_legend(fig, loc="outside lower left")
        fig.legend(handles=refs, title="Reference", loc="outside lower right", ncols=2)
        save(fig, SC / "paper_fig3_allocation", ["scenario_2030_summary.csv (medians, 95% UI)"])

# ---------------------------------------------------------------- Fig S5 F-pool quota use
ann = pd.read_csv(SC / "scenario_annual_services_costs.csv")
if want("paper_fig4_f_pool"):
    a = ann[ann.stratum.isin(["F", "X"]) & ann.scenario_id.isin(finite_ids)]
    tot = a.groupby(["scenario_id", "stratum"])[["screening_quota", "actual_screened"]].sum().reset_index()
    tot = tot.join(meta[["screening_total", "f_share"]], on="scenario_id")
    tot["use"] = tot.actual_screened / tot.screening_quota
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(1, 2, figsize=(FULL * MM, 80 * MM), layout="constrained", sharey=True)
        for ax, (k, title) in zip(axes, [("F", "F (released, still injecting)"), ("X", "X (released, ceased injecting)")]):
            for i, c in enumerate(caps):
                s = tot[(tot.stratum == k) & (tot.screening_total.astype(float) == c)].sort_values("f_share")
                ax.plot(100 * s.f_share, s.use, color=CAP_COL[c], marker=CAP_MARK[i], ms=4)
            ax.set(xticks=[25, 50, 75], ylim=(0, 1.02), xlabel="Share of screening quota allocated to F")
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
            ax.yaxis.set_major_formatter(pct); panel_label(ax, title)
        axes[0].set_ylabel("Share of quota used, 2026–2030")
        cap_legend(fig, loc="outside lower center")
        save(fig, SC / "paper_fig4_f_pool", ["scenario_annual_services_costs.csv (sum of annual medians)"])

# ---------------------------------------------------------------- calibration
if want("calibration_agepop"):
    tg = json.loads((ROOT / "data/calibration_targets.json").read_text())
    ages = tg["age_groups"]; y = np.array(tg["seroprevalence"]); phi = tg["beta_precision"]
    # Seroprevalence is the only target; its interval is the 95% range of the Beta likelihood.
    target = (y, stats.beta.ppf(0.025, y * phi, (1 - y) * phi), stats.beta.ppf(0.975, y * phi, (1 - y) * phi))
    srcs = {"Target (Beta likelihood, 95% range)": (None, "black", "^"),
            "Reference MCMC (production)": (CAL / "mcmc_reference/model_prediction_2015_age_draws.csv", "#0072B2", "s"),
            "Single-round NPE (50,000 LHS)": (CAL / "npe_lhs/model_prediction_2015_age_draws.csv", "#D55E00", "o")}
    chk = pd.read_csv(ROOT / "figures/calibration_interval_overlap.csv")
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(COL * MM, 95 * MM), layout="constrained")
        for j, (name, (f, col, mk)) in enumerate(srcs.items()):
            if f is None:
                est, lo, hi = target
            else:
                dr = pd.read_csv(f); dr = dr[dr.equilibrium_pass.astype(bool)]
                v = dr[[f"p{i}" for i in range(1, 7)]].to_numpy()
                est, lo, hi = np.median(v, 0), np.quantile(v, .025, 0), np.quantile(v, .975, 0)
                who = "MCMC" if "MCMC" in name else "NPE"
                assert np.allclose(lo, chk[f"{who}_lower"], rtol=1e-6) and np.allclose(hi, chk[f"{who}_upper"], rtol=1e-6), \
                    "interval drift"
            xs = np.arange(6) + (j - 1) * 0.22
            ax.errorbar(xs, est, yerr=[est - lo, hi - est], fmt=mk, color=col, ms=4.5, elinewidth=1,
                        capsize=2, label=name)
        ax.set(xticks=range(6), xticklabels=ages, xlabel="Age group (years)", ylim=(0, None),
               ylabel="HCV seroprevalence in prison, 2015")
        ax.yaxis.set_major_formatter(pct)
        fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncols=1)
        save(fig, ROOT / "figures/calibration_agepop",
             ["data/calibration_targets.json", "*/model_prediction_2015_age_draws.csv (median, 2.5-97.5%)"])

# ---------------------------------------------------------------- Fig S2 targets
traj = pd.read_csv(SC / "scenario_trajectories.csv")
ylab = [label[s] for s in order][::-1]  # first scenario at the top
pos = {s: i for i, s in enumerate(order[::-1])}
S = summ.set_index("scenario_id").loc[order]

if want("scenario_target_trajectories"):
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig = plt.figure(figsize=(FULL * MM, 215 * MM), layout="constrained")
        top, bot = fig.subfigures(2, 1, height_ratios=[1, 1.35])
        ta = top.subplots(1, 2, sharey=True); ba = bot.subplots(1, 2, sharey=True)
        for ax, (m, name, thr) in zip(ta, [("incidence_ratio_2015", "Incidence", .2), ("mortality_ratio_2015", "Mortality", .35)]):
            ax.axhline(1, color="0.35", ls=":", lw=0.6); ax.axhline(thr, color="black", ls="--", lw=0.7)
            ax.axvline(2019, color="0.65", ls=":", lw=0.6); ax.axvline(2026, color="0.65", ls="--", lw=0.6)
            h = traj[(traj.scenario_id == "J_only") & (traj.year <= 2025)]
            ax.plot(h.year, h[m], color="0.15", lw=1.3)
            for sid in order:
                s = traj[(traj.scenario_id == sid) & (traj.year >= 2025)]
                if sid in finite_ids:
                    col, ls = CAP_COL[int(meta.loc[sid, "screening_total"])], SHARE_LS[meta.loc[sid, "f_share"]]
                else:
                    col, ls = ("0.35", "-") if sid == "J_only" else ("black", "-")
                ax.plot(s.year, s[m], color=col, ls=ls, lw=1.0)
            ax.set(xticks=[2015, 2019, 2026, 2030], xlim=(2015, 2030), ylim=(0, 1.2), xlabel="Calendar year")
            ax.yaxis.set_major_formatter(pct); ax.set_title(name)
        ta[0].set_ylabel("Annual events relative to 2015")
        top.suptitle("A  WHO impact indicators, 2015–2030 (medians). Dotted 2019: J treatment starts; "
                     "dashed 2026: F/X screening starts", x=0.01, ha="left", fontweight="bold", fontsize=8)
        h1 = [Line2D([], [], color=CAP_COL[c], label=f"{c:,}/yr") for c in caps] + \
             [Line2D([], [], color="0.35", label="J only"), Line2D([], [], color="black", label="No cap")]
        h2 = [Line2D([], [], color="0.3", ls=SHARE_LS[f], label=f"F {f:.0%}") for f in (.25, .5, .75)]
        top.legend(handles=h1, title="Annual F/X screening", loc="outside lower left", ncols=4, fontsize=7)
        top.legend(handles=h2, title="F allocation", loc="outside lower right", ncols=3, fontsize=7)
        for ax, (k, name, thr) in zip(ba, [("incidence", "Incidence", .8), ("mortality", "Mortality", .65)]):
            m = S[f"{k}_reduction"]; ok = S[f"who_{k}_pass"].fillna(False).astype(bool)  # missing = not met
            y = np.array([pos[s] for s in S.index])
            ax.barh(y, m, color=[STATUS[v] for v in ok], height=0.7, edgecolor="0.2", lw=0.3)
            ax.errorbar(m, y, xerr=[m - S[f"{k}_reduction_lo"], S[f"{k}_reduction_hi"] - m], fmt="none",
                        ecolor="0.15", elinewidth=0.6, capsize=1.5)
            ax.axvline(thr, color="black", ls="--", lw=0.8)
            ax.set(xlim=(0, 1), xlabel="Reduction versus 2015", title=f"{name} (target {thr:.0%})")
            ax.xaxis.set_major_formatter(pct)
        ba[0].set(yticks=range(len(order)), yticklabels=ylab)
        bot.suptitle("B  2030 attainment: bars = median, whiskers = 95% UI, dashed = target",
                     x=0.01, ha="left", fontweight="bold", fontsize=8)
        bot.legend(handles=[Patch(fc=STATUS[True], ec="0.2", label="Median meets target"),
                            Patch(fc=STATUS[False], ec="0.2", label="Median misses target")],
                   loc="outside lower center", ncols=2)
        save(fig, SC / "scenario_target_trajectories", ["scenario_trajectories.csv", "scenario_2030_summary.csv"])

if want("scenario_target_heatmap"):
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, ax = plt.subplots(figsize=(COL * MM, 150 * MM), layout="constrained")
        cols = [("Incidence", "incidence_reduction", "who_incidence_pass"),
                ("Mortality", "mortality_reduction", "who_mortality_pass"), ("Joint", None, "who_joint_pass")]
        for j, (name, v, p) in enumerate(cols):
            for sid, r in S.iterrows():
                ok = bool(r[p]) if pd.notna(r[p]) else False; y = pos[sid]  # missing = not met
                ax.add_patch(plt.Rectangle((j, y - .5), 1, 1, fc=matplotlib.colors.to_rgba(STATUS[ok], .3),
                                           ec="white", lw=1.5))
                txt = ("PASS" if ok else "FAIL") if v is None else f"{100 * r[v]:.1f}%  {'pass' if ok else 'fail'}"
                ax.text(j + .5, y, txt, ha="center", va="center", fontsize=7,
                        fontweight="bold" if v is None else "normal")
        ax.set(xlim=(0, 3), ylim=(-.5, len(order) - .5), xticks=[.5, 1.5, 2.5], xticklabels=[c[0] for c in cols],
               yticks=range(len(order)), yticklabels=ylab)
        ax.xaxis.tick_top(); ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xlabel("2030 reduction versus 2015 (median)\nTargets: incidence ≥80%, mortality ≥65%; joint needs both")
        save(fig, SC / "scenario_target_heatmap", ["scenario_2030_summary.csv"])

# ---------------------------------------------------------------- Fig S7 cost bars
if want("scenario_cumulative_costs"):
    st = pd.read_csv(SC / "scenario_stratum_totals.csv")
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(1, 2, figsize=(FULL * MM, 130 * MM), layout="constrained", sharey=True)
        specs = [("discounted_cost_total_sgd", "total_discounted_cost_2026_2030_sgd", 1e6,
                  "Discounted cost (2026 S$ millions)", "A  2026–2030 discounted costs (3%/yr to 2026)", "{:.1f}"),
                 ("treatment_starts", "total_treatment_starts", 1, "Treatment starts (courses, undiscounted)",
                  "B  2026–2030 treatment starts", "{:,.0f}")]
        y = np.array([pos[s] for s in order])
        for ax, (col, tcol, div, xl, title, f) in zip(axes, specs):
            left = np.zeros(len(order))
            for k in ["F", "J", "X"]:
                v = st[st.stratum == k].set_index("scenario_id").reindex(order)[col].fillna(0).to_numpy() / div
                ax.barh(y, v, left=left, color=STRATUM[k], height=0.7, label=k); left += v
            t = S[tcol] / div; lo = S[tcol + "_lo"] / div; hi = S[tcol + "_hi"] / div
            ax.errorbar(t, y, xerr=[t - lo, hi - t], fmt="none", ecolor="0.1", elinewidth=0.6, capsize=1.5)
            for yy, h, tt in zip(y, hi, t):
                ax.text(h, yy, " " + f.format(tt), va="center", fontsize=6.5)
            ax.set_xlim(0, hi.max() * 1.2); ax.set_xlabel(xl); panel_label(ax, title)
            ax.xaxis.set_major_formatter(comma)
        axes[0].set(yticks=range(len(order)), yticklabels=ylab)
        fig.legend(*axes[0].get_legend_handles_labels(), title="Stratum (segment = stratum median); "
                   "label = median total; whisker = 95% UI", loc="outside lower center", ncols=3)
        save(fig, SC / "scenario_cumulative_costs", ["scenario_stratum_totals.csv", "scenario_2030_summary.csv"])

# ---------------------------------------------------------------- F/X annual grids
shares = [.25, .5, .75]
fx = ann[ann.stratum.isin(["F", "X"]) & ann.scenario_id.isin(finite_ids)].join(
    meta[["screening_total", "f_share"]], on="scenario_id")


def grid(name, ylabel, draw, sharey="row", legend=None, sources=("scenario_annual_services_costs.csv",)):
    if not want(name):
        return
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(3, 5, figsize=(FULL * MM, 130 * MM), layout="constrained",
                                 sharex=True, sharey=sharey)
        for r, f in enumerate(shares):
            for c, cap in enumerate(caps):
                ax = axes[r, c]
                ax.set_xticks([2026, 2028, 2030]); ax.yaxis.set_major_formatter(comma)
                draw(ax, fx[(fx.f_share == f) & (fx.screening_total.astype(float) == cap)])
                if r == 0:
                    ax.set_title(f"{cap:,} screens/yr")
                if c == 4:
                    ax.yaxis.set_label_position("right")
                    ax.set_ylabel(f"F {f:.0%} / X {1 - f:.0%}", rotation=270, va="bottom")
        fig.supxlabel("Calendar year", fontsize=8); fig.supylabel(ylabel, fontsize=8)
        fig.legend(handles=legend, loc="outside upper center", ncols=len(legend))
        save(fig, SC / name, list(sources))


SFX = [Line2D([], [], color=STRATUM[k], label=f"{k}") for k in "FX"]


def draw_starts(ax, d):
    for k in "FX":
        s = d[d.stratum == k].sort_values("year")
        ax.fill_between(s.year, s.treatment_starts_lo, s.treatment_starts_hi, color=STRATUM[k], alpha=.18, lw=0)
        ax.plot(s.year, s.treatment_starts, color=STRATUM[k], marker="o", ms=2.5, mfc="white", lw=1)
        ax.plot(s.year, s.treatment_capacity, color=STRATUM[k], ls="--", lw=0.9)


same = np.allclose(fx.treatment_capacity, fx.treatment_starts, rtol=1e-9)
grid("scenario_treatment_starts", "Treatment starts or quota (courses)", draw_starts,
     legend=SFX + [Line2D([], [], color="0.3", marker="o", mfc="white", label="Actual starts (median, 95% UI band)"),
                   Line2D([], [], color="0.3", ls="--", label="Screen-derived quota" +
                          (" (equals starts; hidden)" if same else ""))])
grid("scenario_cumulative_treatment", "Cumulative treatment starts (courses)",
     lambda ax, d: [ax.plot(s.year, s.treatment_starts.cumsum(), color=STRATUM[k], ls="-" if k == "F" else "--", lw=1)
                    for k in "FX" for s in [d[d.stratum == k].sort_values("year")]],
     legend=[Line2D([], [], color=STRATUM["F"], label="F"), Line2D([], [], color=STRATUM["X"], ls="--", label="X")])


def draw_util(ax, d):
    for k in "FX":
        s = d[(d.stratum == k) & (d.treatment_capacity > 0)].sort_values("year")
        ax.plot(s.year, s.treatment_starts / s.treatment_capacity, color=STRATUM[k], ls="-" if k == "F" else "--", lw=1)
    ax.set_ylim(0, 1.01); ax.yaxis.set_major_formatter(pct)


grid("scenario_capacity_utilization", "Screen-derived treatment quota used", draw_util, sharey=True,
     legend=[Line2D([], [], color=STRATUM["F"], label="F"), Line2D([], [], color=STRATUM["X"], ls="--", label="X")])
SERV = {"actual_screened": "Screened", "rna_tests": "RNA tests", "rna_positive": "RNA positive",
        "chronic_identified": "Chronic identified", "treatment_starts": "Treatment starts"}
SERV_COL = dict(zip(SERV, plt.get_cmap("plasma")(np.linspace(0.0, 0.7, 5))))
SERV_MK = dict(zip(SERV, ["o", "s", "^", "D", "v"]))
grid("scenario_care_cascade", "Expected F+X events (sum of stratum medians)",
     lambda ax, d: [ax.plot(g.index, g[k], color=SERV_COL[k], marker=SERV_MK[k], ms=2.5, lw=1)
                    for g in [d.groupby("year")[list(SERV)].sum()] for k in SERV],
     legend=[Line2D([], [], color=SERV_COL[k], marker=SERV_MK[k], label=v) for k, v in SERV.items()])

# ---------------------------------------------------------------- annual costs
if want("scenario_annual_costs"):
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(4, 5, figsize=(FULL * MM, 200 * MM), layout="constrained", sharex=True)
        for ax, sid in zip(axes.flat, order):
            g = ann[ann.scenario_id == sid].groupby("year")[list(ITEM_COLS.values())].sum() / 1e6
            bottom = np.zeros(len(g))
            for it, c in ITEM_COLS.items():
                ax.bar(g.index, g[c], bottom=bottom, color=ITEM[it], width=0.75, label=it); bottom += g[c].to_numpy()
            ax.set_title(label[sid], fontsize=7); ax.set_xticks([2026, 2028, 2030])
        for ax in axes.flat[len(order):]:
            ax.axis("off")
        fig.supxlabel("Calendar year", fontsize=8)
        fig.supylabel("Discounted cost (2026 S$ millions; 3%/yr; sum of stratum medians)", fontsize=8)
        fig.legend(*axes.flat[0].get_legend_handles_labels(), title="Cost item", loc="outside upper center", ncols=4)
        save(fig, SC / "scenario_annual_costs", ["scenario_annual_services_costs.csv"])

# ---------------------------------------------------------------- Fig S6 pies
if want("scenario_cost_composition_pies"):
    it = pd.read_csv(SC / "scenario_item_totals.csv").set_index("scenario_id")
    pass_ids = [s for s in order if s != "J_only" and summ.set_index("scenario_id").loc[s, "target_status"] == "PASS"]
    ids = ["J_only"] + pass_ids
    with style_context("default", palette_name="okabe_ito_on_white"):
        # Grid of at most 3 pies per row with one shared legend; each pie lists its shares below
        # the title, so no per-pie legend can overlap its neighbour.
        ncol = min(3, len(ids)); nrow = -(-len(ids) // ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=((COL if ncol == 1 else FULL) * MM, (25 + 55 * nrow) * MM),
                                 layout="constrained", squeeze=False)
        for ax in axes.flat[len(ids):]:
            ax.set_visible(False)
        for ax, sid in zip(axes.flat, ids):
            v = it.loc[sid, list(ITEM_COLS.values())].to_numpy(float)
            fr = v / v.sum()
            ax.pie(fr, colors=list(ITEM.values()), startangle=90, counterclock=False,
                   wedgeprops=dict(ec="white", lw=0.6))
            name = "Status quo (J only)" if sid == "J_only" else label[sid]
            sh = [f"{k} {100 * x:.1f}%" for k, x in zip(ITEM, fr)]
            shares = ", ".join(sh[:2]) + "\n" + ", ".join(sh[2:])
            ax.set_title(f"{name}: S${S.loc[sid, 'total_discounted_cost_2026_2030_sgd'] / 1e6:.1f}m\n{shares}",
                         fontsize=6.5)
        fig.legend([Patch(fc=c) for c in ITEM.values()], list(ITEM), loc="outside lower center",
                   ncols=len(ITEM), frameon=False)
        fig.suptitle("2026–2030 discounted cost composition\n(2026 S$, 3%/yr; shares of item medians)",
                     fontsize=8, fontweight="bold")
        save(fig, SC / "scenario_cost_composition_pies", ["scenario_item_totals.csv", "scenario_2030_summary.csv"])

# ---------------------------------------------------------------- calibration diagnostics
if want("mcmc_traceplots"):
    dr = pd.read_csv(CAL / "mcmc_reference/mcmc_draws.csv")
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(3, 3, figsize=(FULL * MM, 150 * MM), layout="constrained", sharex=True)
        for ax in axes.flat[len(TH):]:
            ax.set_visible(False)
        chain_col = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
        for ax, t in zip(axes.flat, TH):
            for i, (ch, g) in enumerate(dr.groupby("chain")):
                ax.plot(g.iteration, g[t], color=chain_col[i], lw=0.25, alpha=0.8, label=f"Chain {ch}")
            ax.set_title(f"θ{t[5:]}", fontsize=8)
        thin = json.loads((CAL / "mcmc_reference/mcmc_run_config.json").read_text())["thin_storage"]
        fig.supxlabel(f"Retained iteration (every {thin}th stored)", fontsize=8)
        fig.supylabel("θ (log scale)", fontsize=8)
        fig.legend(*axes.flat[0].get_legend_handles_labels(), loc="outside upper center", ncols=4)
        save(fig, CAL / "mcmc_reference/mcmc_traceplots", ["mcmc_reference/mcmc_draws.csv"])

if want("model_prediction_2015_strata"):
    p = pd.read_csv(CAL / "mcmc_reference/model_prediction_2015_strata.csv").set_index("measure")
    groups = ["D", "J", "F", "X", "current_PWID", "ever_PWID"]
    names = ["D", "J", "F", "X", "Current\nPWID\n(D+F)", "Ever-\nPWID\n(D+J+F+X)"]
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(1, 2, figsize=(FULL * MM, 80 * MM), layout="constrained")
        for ax, (pre, title, col, mk) in zip(axes, [("population_", "A  Total population", "#0072B2", "o"),
                                                    ("active_hcv_", "B  Active HCV cases (a+c+t)", "#D55E00", "^")]):
            r = p.loc[[pre + g for g in groups]]
            ax.errorbar(range(6), r["median"], yerr=[r["median"] - r.lower_95, r.upper_95 - r["median"]],
                        fmt=mk, color=col, ms=4.5, capsize=2)
            ax.set(xticks=range(6), xticklabels=names, ylim=(0, None), ylabel="People, 2015")
            ax.yaxis.set_major_formatter(comma); panel_label(ax, title)
        fig.supxlabel(f"Posterior median and 95% credible interval ({int(p.n_draws_used.iloc[0]):,} reference MCMC draws)",
                      fontsize=8)
        save(fig, CAL / "mcmc_reference/model_prediction_2015_strata", ["model_prediction_2015_strata.csv"])

if want("sbc_rank_histograms"):
    ranks = np.load(CAL / "npe_lhs/sbc_ranks.npy"); n, nb = len(ranks), 20
    lo, hi = stats.binom.ppf([0.005, 0.995], n, 1 / nb)
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(2, 4, figsize=(FULL * MM, 85 * MM), layout="constrained", sharex=True, sharey=True)
        for ax in axes.flat[ranks.shape[1]:]:
            ax.set_visible(False)
        for i, ax in zip(range(ranks.shape[1]), axes.flat):
            ax.axhspan(lo, hi, color="#56B4E9", alpha=0.3, lw=0)
            ax.hist(ranks[:, i], bins=nb, range=(0, 1), color="0.45", edgecolor="white", lw=0.5)
            ax.set_title(f"θ{i + 1}", fontsize=8)
        fig.supxlabel(f"Normalised rank of the true parameter ({n} SBC datasets, {nb} bins); "
                      "band = 99% range under uniformity", fontsize=8)
        fig.supylabel("Datasets per bin", fontsize=8)
        save(fig, CAL / "npe_lhs/sbc_rank_histograms", ["npe_lhs/sbc_ranks.npy"])

CMP = CAL / "npe_vs_mcmc/figures"
if want("marginals") or want("pairs"):
    ref = pd.read_csv(CAL / "mcmc_reference/mcmc_draws.csv")[TH]
    npe = pd.read_csv(CAL / "npe_lhs/posterior_samples_npe.csv")[TH]
    rng = np.random.default_rng(1)
    LEG = [Patch(fc="0.8", label="Reference MCMC"), Line2D([], [], color="#D55E00", label="NPE, 50,000 LHS")]

if want("marginals"):
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(2, 4, figsize=(FULL * MM, 90 * MM), layout="constrained")
        for ax in axes.flat[len(TH):]:
            ax.set_visible(False)
        for ax, t in zip(axes.flat, TH):
            a, b = ref[t].to_numpy(), npe[t].to_numpy()
            g = np.linspace(min(a.min(), b.min()), max(a.max(), b.max()), 300)
            ax.fill_between(g, stats.gaussian_kde(a, "silverman")(g), color="0.8", lw=0)
            ax.plot(g, stats.gaussian_kde(b, "silverman")(g), color="#D55E00", lw=1)
            ax.set_title(f"θ{t[5:]}", fontsize=8); ax.set_yticks([])
        fig.supxlabel("θ (log scale); Gaussian KDE, Silverman bandwidth", fontsize=8)
        fig.supylabel("Density", fontsize=8)
        fig.legend(handles=LEG, loc="outside upper center", ncols=2)
        save(fig, CMP / "marginals", ["mcmc_reference/mcmc_draws.csv", "npe_lhs/posterior_samples_npe.csv"])

if want("accuracy"):
    pc = pd.read_csv(CAL / "npe_vs_mcmc/parameter_comparison.csv").set_index("parameter").loc[TH[::-1]]
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(1, 2, figsize=(FULL * MM, 80 * MM), layout="constrained", sharey=True)
        for ax, (c, band, xl) in zip(axes, [("smd", (-.1, .1), "Standardised mean difference vs reference"),
                                            ("width95_ratio", (.9, 1.1), "95% interval width ratio vs reference")]):
            ax.axvspan(*band, color="0.85", lw=0, label="Tolerance")
            ax.plot(pc[c], range(len(TH)), "o", color="#D55E00", ms=4, label="NPE, 50,000 LHS")
            ax.set_xlabel(xl)
        axes[0].set(yticks=range(len(TH)), yticklabels=[f"θ{t[5:]}" for t in pc.index])
        fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncols=2)
        save(fig, CMP / "accuracy", ["npe_vs_mcmc/parameter_comparison.csv"])

if want("pairs"):
    with style_context("default", palette_name="okabe_ito_on_white"):
        fig, axes = plt.subplots(2, 2, figsize=(COL * MM, 125 * MM), layout="constrained")
        for ax, (u, v) in zip(axes.flat, [("theta1", "theta7"), ("theta6", "theta7"), ("theta2", "theta3"),
                                          ("theta4", "theta5")]):
            sets = [(ref, "0.2", "-"), (npe, "#D55E00", "--")]
            xy = [s[[u, v]].to_numpy()[rng.choice(len(s), 5000, replace=False)].T for s, _, _ in sets]
            allp = np.hstack(xy)
            gx, gy = np.meshgrid(*[np.linspace(allp[i].min(), allp[i].max(), 120) for i in (0, 1)])
            for (s, col, ls), pts in zip(sets, xy):
                zz = stats.gaussian_kde(pts)(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
                ax.contour(gx, gy, zz, levels=6, colors=col, linestyles=ls, linewidths=0.8)
            ax.set(xlabel=f"θ{u[5:]}", ylabel=f"θ{v[5:]}")
        fig.legend(handles=[Line2D([], [], color="0.2", label="Reference MCMC"),
                            Line2D([], [], color="#D55E00", ls="--", label="NPE, 50,000 LHS")],
                   loc="outside lower center", ncols=2)
        save(fig, CMP / "pairs", ["mcmc_reference/mcmc_draws.csv", "npe_lhs/posterior_samples_npe.csv",
                                  "5,000 draws each; 2D Gaussian KDE, 6 levels"])

print("\n".join(f"wrote {w}.png/.svg" for w in written))
