# ==============================================================================
# 02_replot.R — Re-plot all figures at publication pointsize
#
# NOTE FOR READERS: you do not need to run this script to reproduce the
# analysis or figures.  Running 01_analysis.R is sufficient.
#
# What this script does: it loads the saved analysis objects from
# output/analysis_checkpoint.RData (written by 01_analysis.R) and redraws
# every figure with pointsize = 18, which produces larger axis labels and
# annotations suitable for journal submission without having to re-run the
# full analysis (which includes several simulation-intensive steps).
# It is provided for transparency and in case you wish to adjust figure
# text size.  The checkpoint file is not included in the repository; run
# 01_analysis.R once to generate it.
#
# Run order: 01_analysis.R → 02_replot.R (optional) → 03_composites.R (optional)
# ==============================================================================

# ==============================================================================
# CONFIGURATION — paths derive automatically from the project root (.Rproj)
# ==============================================================================
library(here)
OUTPUT_DIR <- here("output")
# ==============================================================================

library(spatstat)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(sparr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

cat("Loading analysis objects...\n")
load(file.path(OUTPUT_DIR, "analysis_checkpoint.RData"))
setwd(OUTPUT_DIR)
cat("Objects loaded. Working dir:", getwd(), "\n")

# Colour palettes
pal_br  <- c(TYPE_A="gold2", TYPE_B="purple", TYPE_C="black", uncovered="deepskyblue3")
pal_mat <- c(BONE="blue", LITHIC="deeppink")
PS <- 18   # base pointsize for all pdf() calls — 50% larger than default 12


# ==============================================================================
# SECTION 3 FIGURES — Q1: Cementation Gradient — Plan View (Combined)
# ==============================================================================

cat("\n=== Re-plotting Section 3 figures (Q1 cementation, plan view) ===\n")

# Fig1a: Spatial distribution
pdf("Fig1a_BrecciaDist_plan.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 1.5, 2), cex.axis = 1.2, cex.lab = 1.2)
plot(ori_br, cols = pal_br[levels(marks(ori_br))],
     pch = 19, cex = 1.2, main = "", axes = TRUE)
dev.off()
cat("Fig1a written.\n")

# Fig1b: Density curves along x and y axes
db_br_clean <- subset(db_combined, !is.na(breccia_type))
p_denx <- ggplot(db_br_clean, aes(x = x, color = breccia_type, fill = breccia_type)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = pal_br) + scale_fill_manual(values = pal_br) +
  labs(title = "Cementation density — X axis", x = "X (m)", y = "Density",
       color = "Type", fill = "Type") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(size = 12))
p_deny <- ggplot(db_br_clean, aes(x = y, color = breccia_type, fill = breccia_type)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = pal_br) + scale_fill_manual(values = pal_br) +
  labs(title = "Cementation density — Y axis", x = "Y (m)", y = "Density",
       color = "Type", fill = "Type") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(size = 12))
pdf("Fig1b_BrecciaDensityX.pdf", height = 8, width = 10, pointsize = PS); print(p_denx); dev.off()
pdf("Fig1b_BrecciaDensityY.pdf", height = 8, width = 10, pointsize = PS); print(p_deny); dev.off()
cat("Fig1b written.\n")

# Fig1c: Overall kernel intensity
pdf("Fig1c_BrecciaIntensity_all.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 1.5, 2), cex.axis = 1.2, cex.lab = 1.2)
plot(den_br_all, axes = TRUE, main = "Kernel Smoothed Intensity -- all cementation")
dev.off()
cat("Fig1c written.\n")

# Fig1d: Kernel intensity per cementation type (4-panel)
pdf("Fig1d_BrecciaIntensity_types.pdf", height = 12, width = 12, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(dA_im, main = "Intensity TYPE_A",    axes = TRUE); contour(dA_im, add = TRUE)
plot(dB_im, main = "Intensity TYPE_B",    axes = TRUE); contour(dB_im, add = TRUE)
plot(dC_im, main = "Intensity TYPE_C",    axes = TRUE); contour(dC_im, add = TRUE)
plot(dU_im, main = "Intensity uncovered", axes = TRUE); contour(dU_im, add = TRUE)
dev.off()
cat("Fig1d written.\n")

# Fig1e: Hotspot maps — scanLRTS (4-panel)
pdf("Fig1e_BrecciaHotspots.pdf", height = 12, width = 12, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(LR_A, main = "Hot spots TYPE_A",    pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
plot(LR_B, main = "Hot spots TYPE_B",    pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
plot(LR_C, main = "Hot spots TYPE_C",    pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
plot(LR_U, main = "Hot spots uncovered", pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
dev.off()
cat("Fig1e written.\n")

# Fig1f: Adaptive relative risk
pdf("Fig1f_BrecciaRelRisk.pdf", height = 12, width = 12, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(rr_BA, main = "B vs A")
plot(rr_CA, main = "C vs A")
plot(rr_CB, main = "C vs B")
plot(rr_UA, main = "uncovered vs A")
dev.off()
cat("Fig1f written.\n")

# Fig1g: Dominant cementation type map
pdf("Fig1g_BrecciaDominantType.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 2, 4), cex.axis = 1.2, cex.lab = 1.2)
plot(dom_br,
     col  = pal_br[c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")],
     main = "Dominant cementation type", axes = TRUE)
dev.off()
cat("Fig1g written.\n")

# Supplementary: Kinhom per cementation type (re-use saved envelopes from Eb objects)
# Note: Kinhom envelopes are recomputed here because the loop objects are not individually
# named in the workspace — regenerated with the same seed for reproducibility.
set.seed(42)
for (tp in c("A", "B", "C", "U")) {
  obj   <- get(paste0("ori_br", tp))
  lab   <- if (tp == "U") "uncovered" else paste0("TYPE_", tp)
  fname <- paste0("FigS_BrecciaKinhom_", tp, ".pdf")
  den_tp <- density(obj, bw.ppl(obj), positive = TRUE)
  pdf(fname, height = 9, width = 9, pointsize = PS)
  par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
  Eb <- envelope(obj, modKinhom,
                 simulate = expression(rpoispp(den_tp)), nsim = 39, global = TRUE)
  plot(Eb, lwd = 2, main = paste("Kinhom —", lab))
  dev.off()
  cat(fname, "written.\n")
}

# Fig1h: nnequal — cementation (rlabel envelope)
set.seed(42)
pdf("Fig1h_BrecciaNNequal.pdf", height = 9, width = 9, pointsize = PS)
par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
nneq_br <- envelope(ori_br, nnequal,
                    simulate = expression(rlabel(ori_br)), nsim = 39)
plot(nneq_br, lwd = 2, main = "Nearest neighbour equality — cementation types")
dev.off()
cat("Fig1h written.\n")


# ==============================================================================
# SECTION 4 FIGURES — Q1: Cementation Gradient — Section View xz (Combined)
# ==============================================================================

cat("\n=== Re-plotting Section 4 figures (Q1 cementation, section view xz) ===\n")

db_br_clean <- subset(db_combined, !is.na(breccia_type))

# FigS: Spatial distribution xz
pdf("FigS_BrecciaDist_section.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 1.5, 2), cex.axis = 1.2, cex.lab = 1.2)
plot(ori_br_xz, cols = pal_br[levels(marks(ori_br_xz))],
     pch = 19, cex = 1.2, main = "", axes = TRUE)
dev.off()
cat("FigS_BrecciaDist_section written.\n")

# FigS: Density along z axis (elevation)
p_denz <- ggplot(db_br_clean, aes(x = z, color = breccia_type, fill = breccia_type)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = pal_br) + scale_fill_manual(values = pal_br) +
  labs(title = "Cementation density — Z axis (elevation)", x = "Z / elevation (m)",
       y = "Density", color = "Type", fill = "Type") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(size = 12))
pdf("FigS_BrecciaDensityZ.pdf", height = 8, width = 10, pointsize = PS); print(p_denz); dev.off()
cat("FigS_BrecciaDensityZ written.\n")

# FigS: Kernel intensity per type — section view
pdf("FigS_BrecciaIntensity_section.pdf", height = 12, width = 12, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(dA_xz, main = "Intensity TYPE_A (xz)",    axes = TRUE); contour(dA_xz, add = TRUE)
plot(dB_xz, main = "Intensity TYPE_B (xz)",    axes = TRUE); contour(dB_xz, add = TRUE)
plot(dC_xz, main = "Intensity TYPE_C (xz)",    axes = TRUE); contour(dC_xz, add = TRUE)
plot(dU_xz, main = "Intensity uncovered (xz)", axes = TRUE); contour(dU_xz, add = TRUE)
dev.off()
cat("FigS_BrecciaIntensity_section written.\n")

# FigS: Dominant type — section view
pdf("FigS_BrecciaDominantType_section.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 2, 4), cex.axis = 1.2, cex.lab = 1.2)
plot(dom_br_xz,
     col  = pal_br[c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")],
     main = "Dominant cementation type (xz section)", axes = TRUE)
dev.off()
cat("FigS_BrecciaDominantType_section written.\n")


# ==============================================================================
# SECTION 5 FIGURES — Q2: Bones vs Lithics — Plan View (Combined)
# ==============================================================================

cat("\n=== Re-plotting Section 5 figures (Q2 bones vs lithics, plan view) ===\n")

db_mat_clean <- subset(db_combined, !is.na(code))

# Fig2a: Spatial distribution
pdf("Fig2a_BoneLithicDist.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 1.5, 2), cex.axis = 1.2, cex.lab = 1.2)
plot(ori_mat, cols = pal_mat[levels(marks(ori_mat))],
     pch = 19, cex = 1.2, main = "", axes = TRUE)
dev.off()
cat("Fig2a written.\n")

# Fig2b: Density curves along x and y axes
p_matx <- ggplot(db_mat_clean, aes(x = x, color = code, fill = code)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = pal_mat) + scale_fill_manual(values = pal_mat) +
  labs(title = "Bone/Lithic density — X axis", x = "X (m)", y = "Density",
       color = "Material", fill = "Material") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(size = 12))
p_maty <- ggplot(db_mat_clean, aes(x = y, color = code, fill = code)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  scale_color_manual(values = pal_mat) + scale_fill_manual(values = pal_mat) +
  labs(title = "Bone/Lithic density — Y axis", x = "Y (m)", y = "Density",
       color = "Material", fill = "Material") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", plot.title = element_text(size = 12))
pdf("Fig2b_BoneLithicDensityX.pdf", height = 8, width = 10, pointsize = PS); print(p_matx); dev.off()
pdf("Fig2b_BoneLithicDensityY.pdf", height = 8, width = 10, pointsize = PS); print(p_maty); dev.off()
cat("Fig2b written.\n")

# Fig2c: Overall kernel intensity
pdf("Fig2c_BoneLithicIntensity_all.pdf", height = 10, width = 10, pointsize = PS)
par(mar = c(4, 4, 2, 4), cex.axis = 1.2, cex.lab = 1.2)
plot(den_mat_all, axes = TRUE, main = "Kernel smoothed intensity — all materials")
dev.off()
cat("Fig2c written.\n")

# Fig2d: Kernel intensity per material type
pdf("Fig2d_BoneLithicIntensity_types.pdf", height = 12, width = 8, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(c(1, 2), nrow = 2, ncol = 1))
plot(density(ori_bone,   sigma = bw.ppl(ori_bone),   eps = 0.01),
     main = "Intensity BONE",   axes = TRUE)
plot(density(ori_lithic, sigma = bw.ppl(ori_lithic), eps = 0.01),
     main = "Intensity LITHIC", axes = TRUE)
dev.off()
cat("Fig2d written.\n")

# Fig2e: Hotspot maps
pdf("Fig2e_BoneLithicHotspots.pdf", height = 12, width = 8, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(c(1, 2), nrow = 2))
plot(LR_bone,   main = "Hot spots BONE",   pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
plot(LR_lithic, main = "Hot spots LITHIC", pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
dev.off()
cat("Fig2e written.\n")

# Fig2f: Adaptive relative risk
pdf("Fig2f_BoneLithicRelRisk.pdf", height = 12, width = 8, pointsize = PS)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
layout(matrix(c(1, 2), nrow = 2))
plot(rr_BL, main = "BONE vs LITHIC")
plot(rr_LB, main = "LITHIC vs BONE")
dev.off()
cat("Fig2f written.\n")

# Supplementary: Kinhom per material type (recomputed with same seed)
set.seed(42)
for (tp in c("bone", "lithic")) {
  obj    <- get(paste0("ori_", tp))
  lab    <- toupper(tp)
  fname  <- paste0("FigS_BoneLithicKinhom_", lab, ".pdf")
  den_tp <- density(obj, bw.ppl(obj), positive = TRUE)
  pdf(fname, height = 10, width = 10, pointsize = PS)
  par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
  Eb <- envelope(obj, modKinhom,
                 simulate = expression(rpoispp(den_tp)), nsim = 39, global = TRUE)
  plot(Eb, lwd = 4, main = paste("Kinhom --", lab))
  dev.off()
  cat(fname, "written.\n")
}

# Fig2g: Cross-K random labelling
set.seed(42)
lev_mat     <- levels(marks(ori_mat))
myKcross_RL <- function(Y, i, j, ...) {
  lam <- density.ppp(Y, sigma = bw.ppl, adjust = 2)
  Kcross.inhom(Y, i = i, j = j, lambdaI = lam, lambdaJ = lam, ...)
}
sim_RL_mat  <- function(Y, ...) {
  Z <- rlabel(Y); marks(Z) <- factor(marks(Z), levels = lev_mat); Z
}
pdf("Fig2g_CrossK_RL.pdf", height = 9, width = 9, pointsize = PS)
par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
Ecross_RL <- envelope(ori_mat, fun = myKcross_RL, nsim = 39,
                      i = "BONE", j = "LITHIC", simulate = sim_RL_mat)
plot(Ecross_RL, lwd = 2,
     main = "Kcross BONE–LITHIC (random labelling null: mixing)")
dev.off()
cat("Fig2g written.\n")

# FigS: Cross-K independence
set.seed(42)
W_mat  <- Window(ori_mat)
nB     <- npoints(ori_bone)
nL     <- npoints(ori_lithic)
lamB      <- density.ppp(ori_bone,   sigma = bw.ppl(ori_bone),   adjust = 2)
lamB$v    <- pmax(lamB$v, 0)
lamL      <- density.ppp(ori_lithic, sigma = bw.ppl(ori_lithic), adjust = 2)
lamL$v    <- pmax(lamL$v, 0)
myKcross_indep <- function(Y, i, j, ...) {
  lamI <- if (i == "BONE") lamB else lamL
  lamJ <- if (j == "BONE") lamB else lamL
  Kcross.inhom(Y, i = i, j = j, lambdaI = lamI, lambdaJ = lamJ, ...)
}
sim_indep_mat <- function(Y, ...) {
  YB <- rpoint(nB, f = lamB, win = W_mat)
  YL <- rpoint(nL, f = lamL, win = W_mat)
  Z  <- superimpose(BONE = YB, LITHIC = YL, W = W_mat)
  marks(Z) <- factor(marks(Z), levels = lev_mat); Z
}
pdf("FigS_CrossK_indep.pdf", height = 9, width = 9, pointsize = PS)
par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
Ecross_indep <- envelope(ori_mat, fun = myKcross_indep, nsim = 39,
                         i = "BONE", j = "LITHIC", simulate = sim_indep_mat)
plot(Ecross_indep, lwd = 2, main = "Kcross BONE–LITHIC (independent processes null)")
dev.off()
cat("FigS_CrossK_indep written.\n")

# Fig2h: nnequal — material type (rlabel envelope)
set.seed(42)
pdf("Fig2h_BoneLithicNNequal.pdf", height = 9, width = 9, pointsize = PS)
par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
nneq_mat <- envelope(ori_mat, nnequal,
                     simulate = expression(rlabel(ori_mat)), nsim = 39)
plot(nneq_mat, lwd = 2, main = "Nearest neighbour equality — BONE vs LITHIC")
dev.off()
cat("Fig2h written.\n")


# ==============================================================================
# SECTION 7 FIGURES — Per-Level Analyses Q1 + Q2 (Supplementary)
# ==============================================================================

cat("\n=== Re-plotting Section 7 figures (per-level supplementary) ===\n")

level_list  <- list("2" = db2, "3" = db3, "4" = db4)
level_short <- list("2" = "Lv2", "3" = "Lv3", "4" = "Lv4")

for (lv in names(level_list)) {
  dl     <- level_list[[lv]]
  lshort <- level_short[[lv]]
  cat("\n--- Level", lv, "(n =", nrow(dl), ") ---\n")

  win_xz_lv <- make_win_xz_lv(dl)
  pp <- make_level_ppp(dl, win_xy)

  # ==== Q1: Breccia per level ====
  if (npoints(pp$brA) >= 5 && npoints(pp$brB) >= 5) {

    pdf(paste0("FigS_", lshort, "_BrecciaDist.pdf"), height = 10, width = 10, pointsize = PS)
    par(mar = c(4, 4, 2, 2), cex.axis = 1.2, cex.lab = 1.2)
    plot(pp$br, cols = pal_br[levels(marks(pp$br))], pch = 19, cex = 1.2,
         main = paste("Cementation —", lv), axes = TRUE)
    dev.off()

    dl_br_c <- subset(dl, !is.na(breccia_type))
    p_bxy <- (
      ggplot(dl_br_c, aes(x = x, color = breccia_type, fill = breccia_type)) +
        geom_density(alpha = 0.15, linewidth = 0.8) +
        scale_color_manual(values = pal_br) + scale_fill_manual(values = pal_br) +
        labs(title = paste("Cementation density X —", lv), x = "X (m)", y = "Density",
             color = "Type", fill = "Type") +
        theme_bw(base_size = 13) +
        theme(legend.position = "bottom", plot.title = element_text(size = 11))
    ) / (
      ggplot(dl_br_c, aes(x = y, color = breccia_type, fill = breccia_type)) +
        geom_density(alpha = 0.15, linewidth = 0.8) +
        scale_color_manual(values = pal_br) + scale_fill_manual(values = pal_br) +
        labs(title = paste("Cementation density Y —", lv), x = "Y (m)", y = "Density",
             color = "Type", fill = "Type") +
        theme_bw(base_size = 13) +
        theme(legend.position = "bottom", plot.title = element_text(size = 11))
    )
    pdf(paste0("FigS_", lshort, "_BrecciaDensityXY.pdf"), height = 10, width = 10, pointsize = PS)
    print(p_bxy); dev.off()

    pdf(paste0("FigS_", lshort, "_BrecciaIntensity.pdf"), height = 12, width = 12, pointsize = PS)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
    layout(matrix(1:4, nrow = 2, byrow = TRUE))
    for (tp in c("brA","brB","brC","brU")) {
      obj_tp <- pp[[tp]]
      lab_tp <- switch(tp, brA="TYPE_A", brB="TYPE_B", brC="TYPE_C", brU="uncovered")
      if (npoints(obj_tp) >= 3) {
        plot(density(obj_tp, sigma = bw.ppl(obj_tp), eps = 0.01),
             main = paste("Intensity", lab_tp, lv), axes = TRUE)
      } else { plot.new(); title(paste("Too few points:", lab_tp)) }
    }; dev.off()

    pdf(paste0("FigS_", lshort, "_BrecciaHotspots.pdf"), height = 12, width = 12, pointsize = PS)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
    layout(matrix(1:4, nrow = 2, byrow = TRUE))
    for (tp in c("brA","brB","brC","brU")) {
      obj_tp <- pp[[tp]]
      lab_tp <- switch(tp, brA="TYPE_A", brB="TYPE_B", brC="TYPE_C", brU="uncovered")
      if (npoints(obj_tp) >= 5) {
        LR_tp <- scanLRTS(obj_tp, r = bw.ppl(obj_tp))
        plot(LR_tp, main = paste("Hot spots", lab_tp, lv), pch = 19, cex = 0.2, axes = TRUE)
        plot(win_xy, add = TRUE, border = "white")
      } else { plot.new(); title(paste("Too few:", lab_tp)) }
    }; dev.off()

    if (npoints(pp$brA) >= 10 && npoints(pp$brB) >= 10) {
      h0_lv <- OS(pp$br, nstar = "geometric")
      pdf(paste0("FigS_", lshort, "_BrecciaRelRisk.pdf"), height = 14, width = 7, pointsize = PS)
      par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
      layout(matrix(1:3, nrow = 3))
      for (pair in list(c("brB","brA","B vs A"), c("brC","brA","C vs A"), c("brC","brB","C vs B"))) {
        tryCatch({
          rr <- risk(pp[[pair[1]]], pp[[pair[2]]], h0 = h0_lv, adapt = TRUE, tolerate = TRUE,
                     hp = OS(pp$br) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
          plot(rr, main = paste(pair[3], "--", lv))
        }, error = function(e) { plot.new(); title(paste("Error:", pair[3])) })
      }; dev.off()
    }

    pdf(paste0("FigS_", lshort, "_BrecciaNNequal.pdf"), height = 9, width = 9, pointsize = PS)
    par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
    tryCatch({
      nneq_br_lv <- envelope(pp$br, nnequal,
                             simulate = expression(rlabel(pp$br)), nsim = 39)
      plot(nneq_br_lv, lwd = 2, main = paste("NNequal cementation —", lv))
    }, error = function(e) { plot.new(); title(paste("NNequal error:", lv)) })
    dev.off()

    cat("  Q1 figures done for level", lv, "\n")
  } else cat("  Skipping Q1 level", lv, "-- insufficient points.\n")

  # ==== Q2: Bones vs Lithics per level ====
  if (npoints(pp$bone) >= 5 && npoints(pp$lithic) >= 5) {

    pdf(paste0("FigS_", lshort, "_BoneLithicDist.pdf"), height = 10, width = 10, pointsize = PS)
    par(mar = c(4, 4, 2, 2), cex.axis = 1.2, cex.lab = 1.2)
    plot(pp$mat, cols = pal_mat[levels(marks(pp$mat))], pch = 19, cex = 1.2,
         main = paste("Bone/Lithic —", lv), axes = TRUE); dev.off()

    dl_mc <- subset(dl, !is.na(code))
    p_mxy <- (
      ggplot(dl_mc, aes(x = x, color = code, fill = code)) +
        geom_density(alpha = 0.15, linewidth = 0.8) +
        scale_color_manual(values = pal_mat) + scale_fill_manual(values = pal_mat) +
        labs(title = paste("Bone/Lithic density X —", lv), x = "X (m)", y = "Density",
             color = "Material", fill = "Material") +
        theme_bw(base_size = 13) +
        theme(legend.position = "bottom", plot.title = element_text(size = 11))
    ) / (
      ggplot(dl_mc, aes(x = y, color = code, fill = code)) +
        geom_density(alpha = 0.15, linewidth = 0.8) +
        scale_color_manual(values = pal_mat) + scale_fill_manual(values = pal_mat) +
        labs(title = paste("Bone/Lithic density Y —", lv), x = "Y (m)", y = "Density",
             color = "Material", fill = "Material") +
        theme_bw(base_size = 13) +
        theme(legend.position = "bottom", plot.title = element_text(size = 11))
    )
    pdf(paste0("FigS_", lshort, "_BoneLithicDensityXY.pdf"), height = 10, width = 10, pointsize = PS)
    print(p_mxy); dev.off()

    pdf(paste0("FigS_", lshort, "_BoneLithicIntensity.pdf"), height = 12, width = 8, pointsize = PS)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
    layout(matrix(c(1, 2), nrow = 2))
    plot(density(pp$bone,   sigma = bw.ppl(pp$bone),   eps = 0.01),
         main = paste("Intensity BONE --",   lv), axes = TRUE)
    plot(density(pp$lithic, sigma = bw.ppl(pp$lithic), eps = 0.01),
         main = paste("Intensity LITHIC --", lv), axes = TRUE)
    dev.off()

    pdf(paste0("FigS_", lshort, "_BoneLithicHotspots.pdf"), height = 12, width = 8, pointsize = PS)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
    layout(matrix(c(1, 2), nrow = 2))
    LR_b <- scanLRTS(pp$bone,   r = bw.ppl(pp$bone))
    LR_l <- scanLRTS(pp$lithic, r = bw.ppl(pp$lithic))
    plot(LR_b, main = paste("Hot spots BONE --",   lv), pch = 19, cex = 0.2, axes = TRUE)
    plot(win_xy, add = TRUE, border = "white")
    plot(LR_l, main = paste("Hot spots LITHIC --", lv), pch = 19, cex = 0.2, axes = TRUE)
    plot(win_xy, add = TRUE, border = "white")
    dev.off()

    if (npoints(pp$bone) >= 10 && npoints(pp$lithic) >= 10) {
      h0_ml <- OS(pp$mat, nstar = "geometric")
      pdf(paste0("FigS_", lshort, "_BoneLithicRelRisk.pdf"), height = 12, width = 8, pointsize = PS)
      par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1, cex.axis = 1.2, cex.lab = 1.2)
      layout(matrix(c(1, 2), nrow = 2))
      tryCatch({
        rr <- risk(pp$bone, pp$lithic, h0=h0_ml, adapt=TRUE, tolerate=TRUE,
                   hp=OS(pp$mat)/2, pilot.symmetry="pooled", davies.baddeley=0.05)
        plot(rr, main = paste("BONE vs LITHIC --", lv))
      }, error = function(e) { plot.new(); title("Error B-L") })
      tryCatch({
        rr <- risk(pp$lithic, pp$bone, h0=h0_ml, adapt=TRUE, tolerate=TRUE,
                   hp=OS(pp$mat)/2, pilot.symmetry="pooled", davies.baddeley=0.05)
        plot(rr, main = paste("LITHIC vs BONE --", lv))
      }, error = function(e) { plot.new(); title("Error L-B") })
      dev.off()
    }

    # Cross-K random labelling — only for levels with n_lithic >= 30
    # (Level 3 has only 19 lithics and is therefore skipped; no CrossK PDF will
    #  be produced for it, and the composite figure will use a 3×2 layout
    #  without a CrossK panel — see 03_composites.R)
    if (npoints(pp$lithic) >= 30) {
      lev_ml <- levels(marks(pp$mat))
      myKcross_RL_lv <- function(Y, i, j, ...) {
        lam <- density.ppp(Y, sigma = bw.ppl, adjust = 2)
        Kcross.inhom(Y, i=i, j=j, lambdaI=lam, lambdaJ=lam, ...)
      }
      sim_RL_lv <- function(Y, ...) {
        Z <- rlabel(Y); marks(Z) <- factor(marks(Z), levels = lev_ml); Z
      }
      pdf(paste0("FigS_", lshort, "_CrossK_RL.pdf"), height = 9, width = 9, pointsize = PS)
      par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
      tryCatch({
        Ecrl <- envelope(pp$mat, fun=myKcross_RL_lv, nsim=39,
                         i="BONE", j="LITHIC", simulate=sim_RL_lv)
        plot(Ecrl, lwd=2, main = paste("Kcross BONE–LITHIC (RL) —", lv))
      }, error = function(e) { plot.new(); title(paste("CrossK error:", lv)) })
      dev.off()
    } else {
      cat("  NOTE: CrossK_RL skipped for level", lv,
          "— too few lithics (n =", npoints(pp$lithic), "< 30).",
          "No CrossK PDF will be generated for this level.\n")
    }

    pdf(paste0("FigS_", lshort, "_BoneLithicNNequal.pdf"), height = 9, width = 9, pointsize = PS)
    par(mar = c(4, 4, 2.5, 2), cex.axis = 1.2, cex.lab = 1.2)
    tryCatch({
      nneq_lv <- envelope(pp$mat, nnequal,
                          simulate = expression(rlabel(pp$mat)), nsim = 39)
      plot(nneq_lv, lwd = 2, main = paste("NNequal BONE/LITHIC —", lv))
    }, error = function(e) { plot.new(); title(paste("NNequal error:", lv)) })
    dev.off()

    cat("  Q2 figures done for level", lv, "\n")
  } else cat("  Skipping Q2 level", lv, "-- insufficient points.\n")
}

cat("\n*** 02_replot.R complete — all figures replotted with PS =", PS, "***\n")
