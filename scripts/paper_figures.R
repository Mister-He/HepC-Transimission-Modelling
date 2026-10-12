#!/usr/bin/env Rscript
# =============================================================================
# paper_figures.R — manuscript figures (model schematic, trajectories, allocation and
# F-quota use) built from output/scenarios_2026/*.csv; run scenarios_2026.R first.
# =============================================================================

suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
ROOT <- normalizePath(if (length(args)) args[1] else getwd(), mustWork = TRUE)
OUT <- file.path(ROOT, "output/scenarios_2026")
theme_paper <- function() theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
cols <- c(D = "#7F7F7F", J = "#0072B2", F = "#D55E00", X = "#009E73")

# Figure 1: model structure -------------------------------------------------
box <- data.frame(id = c("D", "J", "F", "X"), x = c(1, 4, 7.6, 7.6), y = c(2, 2, 3.4, 0.6),
  label = c("D\nActive PWID,\nnever incarcerated", "J\nCurrently\nincarcerated",
            "F\nReleased,\nstill injecting", "X\nReleased,\nceased injecting"))
arr <- data.frame(x = c(2.1, 5.1, 5.1), xend = c(2.9, 6.5, 6.5), y = c(2, 2.4, 1.6), yend = c(2, 3.2, 0.8))
p1 <- ggplot() +
  geom_tile(data = box, aes(x, y, fill = id), width = 2.2, height = 1.2, alpha = .9) +
  geom_text(data = box, aes(x, y, label = label), colour = "white", size = 3.3, fontface = "bold") +
  geom_segment(data = arr, aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(.12, "in"), type = "closed"), linewidth = .6, colour = "grey30") +
  geom_curve(data = data.frame(x = 6.5, xend = 4.3, y = 3.75, yend = 2.6), aes(x = x, xend = xend, y = y, yend = yend),
    curvature = .35, arrow = arrow(length = unit(.12, "in"), type = "closed"), linewidth = .6, colour = "grey30") +
  annotate("text", x = 2.5, y = 2.25, label = "First arrest", size = 3) +
  annotate("text", x = 5.75, y = 3.1, label = "Release,\nreturn to injecting", size = 2.9, hjust = 1, lineheight = .9) +
  annotate("text", x = 5.75, y = 0.95, label = "Release,\nstop injecting", size = 2.9, hjust = 1, lineheight = .9) +
  annotate("text", x = 4.9, y = 4.05, label = "Re-arrest", size = 3) +
  annotate("label", x = 4, y = 0.3, size = 3.1, fill = "white",
    label = "Treatment in J from 2019", colour = cols[["J"]], fontface = "bold") +
  annotate("label", x = 7.6, y = 4.45, size = 3.1, fill = "white",
    label = "Screening quotas for F and X from 2026", colour = cols[["F"]], fontface = "bold") +
  annotate("label", x = 1, y = 3.2, size = 3.1, fill = "white",
    label = "Transmission in D and F only", colour = "grey20") +
  scale_fill_manual(values = cols) + coord_cartesian(xlim = c(-0.3, 9), ylim = c(-0.1, 4.7)) +
  theme_void() + theme(legend.position = "none")
ggsave(file.path(OUT, "paper_fig1_schematic.png"), p1, width = 8, height = 4.2, dpi = 300, bg = "white")

# Figure 2: incidence and mortality trajectories -----------------------------
tr <- read.csv(file.path(OUT, "scenario_trajectories.csv"))
nt <- read.csv(file.path(OUT, "no_treatment_trajectory.csv"))
pick <- c(no_treatment = "No treatment", J_only = "Prison treatment only (J)",
          S8000_F25 = "J + 8,000 screens/yr, F 25%", S8000_F75 = "J + 8,000 screens/yr, F 75%")
d <- rbind(nt[, names(nt) %in% names(tr)], tr[tr$scenario_id %in% names(pick), names(tr) %in% names(nt)])
d <- d[d$scenario_id %in% names(pick), ]
# Screening scenarios share the J-only history; draw them from 2025 only.
d <- d[!(grepl("^S8000", d$scenario_id) & d$year < 2025), ]
d$label <- factor(pick[d$scenario_id], levels = pick)
long <- rbind(
  data.frame(d[c("label", "year")], panel = "A  HCV incidence", m = d$incidence_ratio_2015,
             lo = d$incidence_ratio_2015_lo, hi = d$incidence_ratio_2015_hi),
  data.frame(d[c("label", "year")], panel = "B  HCV-related mortality", m = d$mortality_ratio_2015,
             lo = d$mortality_ratio_2015_lo, hi = d$mortality_ratio_2015_hi))
thr <- data.frame(panel = c("A  HCV incidence", "B  HCV-related mortality"), y = c(.2, .35))
line_cols <- c("grey35", "#0072B2", "#D55E00", "#009E73")
names(line_cols) <- pick
p2 <- ggplot(long, aes(year, m, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .12, colour = NA) + geom_line(linewidth = .8) +
  geom_hline(data = thr, aes(yintercept = y), linetype = 2, colour = "black", linewidth = .4) +
  geom_vline(xintercept = c(2019, 2026), linetype = 3, colour = "grey60") +
  facet_wrap(~panel) + scale_colour_manual(values = line_cols) + scale_fill_manual(values = line_cols) +
  scale_x_continuous(breaks = c(2015, 2019, 2026, 2030)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.4)) +
  labs(x = "Calendar year", y = "Annual events relative to 2015", colour = NULL, fill = NULL) +
  guides(colour = guide_legend(nrow = 2)) + theme_paper()
ggsave(file.path(OUT, "paper_fig2_trajectories.png"), p2, width = 9, height = 4.6, dpi = 300, bg = "white")

# Figure 3: 2030 reductions by capacity and F share ---------------------------
s <- read.csv(file.path(OUT, "scenario_2030_summary.csv"))
f <- s[s$scenario_type == "finite", ]
f$cap <- factor(format(f$screening_total, big.mark = ","), levels = format(sort(unique(f$screening_total)), big.mark = ","))
f$fshare <- 100 * f$f_share
mk <- function(prefix, label) data.frame(f[c("cap", "fshare")], panel = label,
  m = f[[prefix]], lo = f[[paste0(prefix, "_lo")]], hi = f[[paste0(prefix, "_hi")]])
g <- rbind(mk("incidence_reduction", "A  Incidence reduction, 2030"), mk("mortality_reduction", "B  Mortality reduction, 2030"))
thr3 <- data.frame(panel = unique(g$panel), y = c(.8, .65))
ref <- s[s$scenario_id %in% c("J_only", "SInf"), ]
refd <- rbind(data.frame(panel = unique(g$panel)[1], ref[c("scenario_id")], y = ref$incidence_reduction),
              data.frame(panel = unique(g$panel)[2], ref[c("scenario_id")], y = ref$mortality_reduction))
refd$lab <- ifelse(refd$scenario_id == "J_only", "J only", "No screening cap")
p3 <- ggplot(g, aes(fshare, m, colour = cap, group = cap)) +
  geom_hline(data = thr3, aes(yintercept = y), linetype = 2, linewidth = .4) +
  geom_hline(data = refd, aes(yintercept = y, linetype = lab), colour = "grey40", linewidth = .4) +
  geom_line(position = position_dodge(width = 4), linewidth = .6) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 4), size = .35) +
  facet_wrap(~panel) + scale_colour_viridis_d(option = "D", end = .85) +
  scale_linetype_manual(values = c(`J only` = 3, `No screening cap` = 1)) +
  scale_x_continuous(breaks = c(25, 50, 75), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Share of screening quota allocated to F", y = "Reduction versus 2015",
       colour = "Screens per year", linetype = "Reference") + theme_paper()
ggsave(file.path(OUT, "paper_fig3_allocation.png"), p3, width = 9, height = 4.6, dpi = 300, bg = "white")

# Figure 4: why a larger F share can waste quota ------------------------------
a <- read.csv(file.path(OUT, "scenario_annual_services_costs.csv"))
a <- a[a$stratum %in% c("F", "X") & !a$scenario_id %in% c("J_only", "SInf"), ]
tot <- aggregate(cbind(screening_quota, actual_screened) ~ scenario_id + stratum, a, sum)
meta <- unique(s[s$scenario_type == "finite", c("scenario_id", "screening_total", "f_share")])
tot <- merge(tot, meta)
tot$use <- tot$actual_screened / tot$screening_quota
tot$cap <- factor(format(tot$screening_total, big.mark = ","), levels = levels(f$cap))
tot$fshare <- 100 * tot$f_share
p4 <- ggplot(tot, aes(fshare, use, colour = cap, group = cap)) + geom_line() + geom_point() +
  facet_wrap(~stratum, labeller = labeller(stratum = c(F = "F (released, still injecting)", X = "X (released, ceased injecting)"))) +
  scale_colour_viridis_d(option = "D", end = .85) +
  scale_x_continuous(breaks = c(25, 50, 75), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.02)) +
  labs(x = "Share of screening quota allocated to F", y = "Share of quota used, 2026-2030",
       colour = "Screens per year") + theme_paper()
ggsave(file.path(OUT, "paper_fig4_f_pool.png"), p4, width = 8, height = 4, dpi = 300, bg = "white")
cat("Wrote paper figures to", OUT, "\n")
