# ==============================================================================
#   Spatial Analysis — Gruta da Companheira (GdC), Chamber 3
#   Script: 01_analysis.R
#
#   Alzate-Casallas et al. (in prep.)
#   "Spatial point pattern analysis of archaeological finds and cementation
#    types from Chamber 3, Gruta da Companheira (Portugal)"
#
# ------------------------------------------------------------------------------
#   HOW TO USE THIS SCRIPT
#
#   1. Open GDC_Spatial_Analysis.Rproj in RStudio (or set your working
#      directory to the project root).  No path editing is required.
#
#   2. Run the script in full (source("code/01_analysis.R") or open it in
#      RStudio and use Session > Source File).  All figures and statistics
#      tables are saved to output/.
#
#   3. Required packages (install once if needed):
#        install.packages(c("here",
#                            "spatstat", "spatstat.geom", "spatstat.explore",
#                            "spatstat.random", "sparr", "ggplot2", "patchwork",
#                            "dplyr", "tidyr", "flextable", "dixon", "readxl"))
#      For exact package versions: install.packages("renv"); renv::restore()
#
#   Expected run time: 10–20 minutes (simulation-intensive steps use nsim = 999
#   for the Hopkins–Skellam and scan tests, nsim = 39 for K-function envelopes).
#   See README.md for full instructions and output file descriptions.
# ------------------------------------------------------------------------------
#
#   Research questions:
#     Q1 — Cementation gradient (TYPE_A / TYPE_B / TYPE_C / uncovered)
#     Q2 — Material type distribution (BONE / LITHIC)
#     Q3 — Cross-tabulation cementation × material type
#
#   Sections:
#     0 — Packages and setup
#     1 — Data loading and cleaning
#     2 — Spatial windows, ppp objects, helper functions
#     3 — Q1: Cementation gradient, plan view (combined)
#     4 — Q1: Cementation gradient, section view xz (combined)
#     5 — Q2: Bones vs Lithics, plan view (combined)
#     6 — Q3: Cross-tabulation cementation × material type
#     7 — Per-level analyses Q1 + Q2 (supplementary)
#     8 — Summary statistics tables
#
#   Outputs (saved to output/):
#     Fig1a–h_*.pdf  — Q1 cementation (plan view, combined levels)
#     Fig2a–h_*.pdf  — Q2 bone/lithic (plan view, combined levels)
#     Fig3_*.pdf     — Q3 contingency bar chart
#     FigS_*.pdf     — Supplementary (section view, per-level analyses)
#     Table1_CombinedStats.csv
#     Table2_Contingency.csv
#     TableS_PerLevel_Stats.csv
# ==============================================================================


# ==============================================================================
# SECTION 0 — Packages and Setup
# ==============================================================================

library(here)
library(spatstat)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(sparr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(flextable)
library(dixon)
library(readxl)

# ==============================================================================
# CONFIGURATION — paths derive automatically from the project root (.Rproj)
# ==============================================================================
DATA_DIR   <- here("data")
OUTPUT_DIR <- here("output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
# ==============================================================================

options(max.print = .Machine$integer.max)
setwd(DATA_DIR)


# ==============================================================================
# SECTION 1 — Data Loading and Cleaning
# ==============================================================================

gdc_raw <- read.table(
  file             = "GDC_SPATIAL_FEB26_archmaterials.txt",
  sep              = "\t",
  header           = TRUE,
  stringsAsFactors = FALSE
)

cat("Raw dimensions:", nrow(gdc_raw), "rows,", ncol(gdc_raw), "cols\n")
cat("Raw levels:\n");        print(table(gdc_raw$level, useNA = "always"))
cat("Raw breccia types:\n"); print(table(gdc_raw$breccia_type, useNA = "always"))

# Fix 1: strip trailing/leading spaces from breccia_type
gdc_raw$breccia_type <- trimws(gdc_raw$breccia_type)

# Fix 3: rename NA / empty breccia entries → "uncovered"
gdc_raw$breccia_type[is.na(gdc_raw$breccia_type) |
                       gdc_raw$breccia_type == "NA" |
                       gdc_raw$breccia_type == ""] <- "uncovered"

cat("\nLevels after cleaning:\n");            print(table(gdc_raw$level))
cat("Breccia types after cleaning:\n");      print(table(gdc_raw$breccia_type, useNA = "always"))

# Ordered factor levels — consistent colours across all figures
br_levels  <- c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")
mat_levels <- c("BONE", "LITHIC")

gdc_raw$breccia_type <- factor(gdc_raw$breccia_type, levels = br_levels)
gdc_raw$code         <- factor(gdc_raw$code,         levels = mat_levels)

# Exclude Level 1 (late intrusive sediment — not part of primary deposit)
gdc <- subset(gdc_raw, level != "1")
cat("\nRows after excluding Level 1:", nrow(gdc), "\n")

# Create db_combined: levels 104_e + 2 + 3 + 4
# NOTE: the field level "104E" is stored as "104_e" in the data file
# Level 3 = former Level 3 + former Level 3/1B (merged by field team)
# Level 4 = former Level 1B (renamed by field team)
valid_levels <- c("104_e", "2", "3", "4")
db_combined  <- subset(gdc, level %in% valid_levels)

cat("db_combined n =", nrow(db_combined), "\n")
cat("Level breakdown:\n");   print(table(db_combined$level))
cat("Breccia types:\n");     print(table(db_combined$breccia_type, useNA = "always"))
cat("Material codes:\n");    print(table(db_combined$code,         useNA = "always"))

# Per-level subsets (exclude 104_e from per-level analyses)
db2 <- subset(db_combined, level == "2")
db3 <- subset(db_combined, level == "3")
db4 <- subset(db_combined, level == "4")

cat("\nPer-level n: Lv2 =", nrow(db2),
    "| Lv3 =", nrow(db3),
    "| Lv4 =", nrow(db4), "\n")


# ==============================================================================
# SECTION 2 — Spatial Windows, ppp Objects, Helper Functions
# ==============================================================================

# ==============================================================================
# Observation windows — defined from excavation boundary reference data
# ==============================================================================
#
# Source: GDC_CLEAN_DATABASE_FEB_2026.xlsx
#
# Boundary reference points used:
#   MM5, MM6, MM8, MM9, MM10, MM11, MM12 — micromorphology block front faces
#     at deposit level (z < 97 m). MM2 (z=98.2) and MM3 (z=97.6) excluded:
#     taken from the upper loose sediment during Level 1 surface cleaning, not
#     representative of the main deposit boundary.
#   Cave_topo1–8 — total-station points on the northern cave wall
#     (y ≈ 101.4–101.5 m), forming the NE boundary of the excavation.
#   north_profile, northwest_profile, east_profile, south_profile,
#   Limpeza perfil N, Limpeza_perfil — profile-cleaning points defining
#     the excavated section faces.
#
# Plan view window (win_xy):
#   Convex hull of (all finds + deposit-level boundary references).
#   Area: 7.99 m² vs bounding rectangle 10.89 m² (27% reduction).
#   All 1012 analysis points verified inside the window.
#
# Section view window (win_xz):
#   Convex hull of the combined find cloud in x–z projection.
#   Area: 4.04 m² vs bounding rectangle 7.00 m² (42% reduction).
#   Using the finds-only hull is appropriate because the section view
#   captures the vertical profile of the deposit; the boundary references
#   (MM, profiles) define the plan footprint, not a distinct xz constraint.
#
# Per-level section view windows (make_win_xz_lv):
#   Level-specific convex hull in xz — correctly bounds each level to its
#   own observed depth range. Intensity is reported per m², so different
#   window sizes between levels do not prevent comparison.

db_boundary <- read_excel("GDC_CLEAN_DATABASE_FEB_2026.xlsx")

mm_dep_codes    <- c("MM5","MM6","MM8","MM9","MM10","MM11","MM12")
cave_topo_codes <- paste0("Cave_topo", 1:8)
profile_levels  <- c("north_profile","northwest_profile","east_profile",
                     "south_profile","Limpeza perfil N","Limpeza_perfil")

bnd_pts <- subset(db_boundary,
                  code  %in% c(mm_dep_codes, cave_topo_codes) |
                  level %in% profile_levels)

cat("Boundary reference points loaded:", nrow(bnd_pts), "\n")

# Plan view: convex hull of finds + boundary references
win_xy <- convexhull.xy(
  x = c(db_combined$x, bnd_pts$x),
  y = c(db_combined$y, bnd_pts$y)
)

# Section view: convex hull of finds in xz
win_xz <- convexhull.xy(
  x = db_combined$x,
  y = db_combined$z
)

# Per-level section view window builder
make_win_xz_lv <- function(db) convexhull.xy(db$x, db$z)

cat(sprintf("win_xy area: %.3f m2  (rectangle was %.3f m2, -27%%)\n",
            area.owin(win_xy),
            diff(range(db_combined$x)) * diff(range(db_combined$y))))
cat(sprintf("win_xz area: %.3f m2  (rectangle was %.3f m2, -42%%)\n",
            area.owin(win_xz),
            diff(range(db_combined$x)) * diff(range(db_combined$z))))

# Colour palettes — consistent across all figures in the paper
pal_br  <- c(TYPE_A = "gold2", TYPE_B = "purple", TYPE_C = "black",
             uncovered = "deepskyblue3")
pal_mat <- c(BONE = "blue", LITHIC = "deeppink")

# --- Q1: cementation (plan view) ---
ori_br  <- unique(ppp(db_combined$x, db_combined$y, win_xy,
                      marks = db_combined$breccia_type))

db_A    <- subset(db_combined, breccia_type == "TYPE_A")
db_B    <- subset(db_combined, breccia_type == "TYPE_B")
db_C    <- subset(db_combined, breccia_type == "TYPE_C")
db_U    <- subset(db_combined, breccia_type == "uncovered")

ori_brA <- unique(ppp(db_A$x, db_A$y, win_xy))
ori_brB <- unique(ppp(db_B$x, db_B$y, win_xy))
ori_brC <- unique(ppp(db_C$x, db_C$y, win_xy))
ori_brU <- unique(ppp(db_U$x, db_U$y, win_xy))

cat("Q1 ppp (plan): A =", npoints(ori_brA), "| B =", npoints(ori_brB),
    "| C =", npoints(ori_brC), "| uncovered =", npoints(ori_brU),
    "| combined =", npoints(ori_br), "\n")

# --- Q1: cementation (section view xz) ---
ori_br_xz  <- unique(ppp(db_combined$x, db_combined$z, win_xz,
                         marks = db_combined$breccia_type))
ori_brA_xz <- unique(ppp(db_A$x, db_A$z, win_xz))
ori_brB_xz <- unique(ppp(db_B$x, db_B$z, win_xz))
ori_brC_xz <- unique(ppp(db_C$x, db_C$z, win_xz))
ori_brU_xz <- unique(ppp(db_U$x, db_U$z, win_xz))

# --- Q2: material type (plan view) ---
ori_mat    <- unique(ppp(db_combined$x, db_combined$y, win_xy,
                         marks = db_combined$code))
ori_bone   <- unique(ppp(subset(db_combined, code == "BONE")$x,
                         subset(db_combined, code == "BONE")$y, win_xy))
ori_lithic <- unique(ppp(subset(db_combined, code == "LITHIC")$x,
                         subset(db_combined, code == "LITHIC")$y, win_xy))

cat("Q2 ppp: BONE =", npoints(ori_bone),
    "| LITHIC =", npoints(ori_lithic),
    "| combined =", npoints(ori_mat), "\n")

# --- nnequal function ---
# Nearest-neighbour equality function
# Developed by Adrian Baddeley (personal communication)
# Used here with permission
nnequal <- function(X, ..., kmax = 20, ratio = TRUE, cumulative = TRUE) {
  stopifnot(is.ppp(X))
  stopifnot(is.multitype(X))
  N        <- nnwhich(X, k = 1:kmax)
  D        <- nndist(X,  k = 1:kmax)
  B        <- bdist.points(X)
  observed <- (D <= B)
  marx     <- marks(X)
  mI       <- matrix(marx[row(N)], ncol = kmax)
  mJ       <- matrix(marx[N],      ncol = kmax)
  counted  <- (mI == mJ)
  if (cumulative) {
    numer <- rowSums(apply(observed & counted, 1, cumsum))
    denom <- rowSums(apply(observed,           1, cumsum))
  } else {
    numer <- colSums(observed & counted)
    denom <- colSums(observed)
  }
  estimate <- ifelse(denom > 0, numer / denom, 0)
  m        <- as.integer(table(marx))
  n        <- npoints(X)
  pequal   <- sum(m * (m - 1) / (n * (n - 1)))
  df       <- data.frame(k = 1:kmax, theo = pequal, bord = estimate)
  desc     <- c("neighbour order k",
                "theoretical %s",
                "border-corrected estimate of %s")
  labl     <- c("k", "%s[theo](k)", "hat(%s)[bord](k)")
  dendf    <- data.frame(k = 1:kmax, theo = denom, bord = denom)
  yexp     <- ylab <- quote(E(k))
  fname    <- "E"
  ratfv(df, NULL, dendf,
        argu  = "k", ylab = ylab, valu = "bord", fmla = . ~ k,
        alim  = c(1, kmax), labl = labl, desc = desc,
        fname = fname, yexp = yexp, ratio = ratio,
        unitname = c("neighbour step", "neighbour steps"))
}

# --- modKinhom helper ---
modKinhom <- function(X, ...) {
  den <- density(X, bw.ppl)
  Kinhom(X, den, ...)
}

# --- make_level_ppp: build all ppp objects for one stratigraphic level ---
make_level_ppp <- function(db, win) {
  db$breccia_type <- factor(db$breccia_type, levels = br_levels)
  db$code         <- factor(db$code,         levels = mat_levels)
  list(
    mat    = unique(ppp(db$x, db$y, win, marks = db$code)),
    bone   = unique(ppp(subset(db, code == "BONE")$x,
                        subset(db, code == "BONE")$y, win)),
    lithic = unique(ppp(subset(db, code == "LITHIC")$x,
                        subset(db, code == "LITHIC")$y, win)),
    br     = unique(ppp(db$x, db$y, win, marks = db$breccia_type)),
    brA    = unique(ppp(subset(db, breccia_type == "TYPE_A")$x,
                        subset(db, breccia_type == "TYPE_A")$y, win)),
    brB    = unique(ppp(subset(db, breccia_type == "TYPE_B")$x,
                        subset(db, breccia_type == "TYPE_B")$y, win)),
    brC    = unique(ppp(subset(db, breccia_type == "TYPE_C")$x,
                        subset(db, breccia_type == "TYPE_C")$y, win)),
    brU    = unique(ppp(subset(db, breccia_type == "uncovered")$x,
                        subset(db, breccia_type == "uncovered")$y, win))
  )
}

# --- dominant_type_map: argmax over a list of density images ---
# Uses pixel matrix approach (sapply over im objects fails in spatstat 3.x)
dominant_type_map <- function(im_list, labels) {
  # im_list: named list of spatstat im objects (same grid)
  mats   <- lapply(im_list, `[[`, "v")       # extract pixel matrices
  arr    <- array(unlist(mats), dim = c(dim(mats[[1]]), length(mats)))
  idx    <- apply(arr, c(1, 2), function(x) {
    x[is.na(x)] <- -Inf
    if (all(!is.finite(x))) return(NA_integer_)
    which.max(x)
  })
  fac_v  <- factor(as.vector(idx), levels = seq_along(labels), labels = labels)
  ref    <- im_list[[1]]
  im(matrix(as.integer(fac_v), nrow = nrow(ref$v), ncol = ncol(ref$v)),
     xcol = ref$xcol, yrow = ref$yrow)
}


# ==============================================================================
# SECTION 3 — Q1: Cementation Gradient — Plan View (Combined)
# ==============================================================================

# Redirect all output to output/ folder
setwd(OUTPUT_DIR)
cat("Output directory:", getwd(), "\n")

cat("\n=== Section 3: Q1 Cementation gradient (plan view, combined) ===\n")

main_stats  <- data.frame()   # collector for Table 2
ori_br_u    <- unmark(ori_br) # unmarked combined breccia pattern
db_br_clean <- subset(db_combined, !is.na(breccia_type))

# Fig 2a: Spatial distribution
pdf("Fig1a_BrecciaDist_plan.pdf", height = 10, width = 10)
par(mar = c(4, 4, 1.5, 2))
plot(ori_br, cols = pal_br[levels(marks(ori_br))],
     pch = 19, cex = 1.2, main = "", axes = TRUE)
dev.off()
cat("Fig1a written.\n")

# Fig 2b: Density curves along x and y axes
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
pdf("Fig1b_BrecciaDensityX.pdf", height = 8, width = 10); print(p_denx); dev.off()
pdf("Fig1b_BrecciaDensityY.pdf", height = 8, width = 10); print(p_deny); dev.off()
cat("Fig1b written.\n")

# Fig 2c: Overall kernel intensity (bandwidth by likelihood cross-validation)
den_br_all <- density(ori_br_u, bw.ppl(ori_br_u), eps = 0.01, positive = TRUE)
pdf("Fig1c_BrecciaIntensity_all.pdf", height = 10, width = 10)
plot(den_br_all, axes = TRUE, main = "Kernel Smoothed Intensity -- all cementation")
dev.off()
cat("Fig1c written.\n")

# Fig 2d: Kernel intensity per cementation type (4-panel)
dA_im <- density(ori_brA, sigma = bw.ppl(ori_brA), eps = 0.01)
dB_im <- density(ori_brB, sigma = bw.ppl(ori_brB), eps = 0.01)
dC_im <- density(ori_brC, sigma = bw.ppl(ori_brC), eps = 0.01)
dU_im <- density(ori_brU, sigma = bw.ppl(ori_brU), eps = 0.01)
pdf("Fig1d_BrecciaIntensity_types.pdf", height = 12, width = 12)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(dA_im, main = "Intensity TYPE_A",    axes = TRUE); contour(dA_im, add = TRUE)
plot(dB_im, main = "Intensity TYPE_B",    axes = TRUE); contour(dB_im, add = TRUE)
plot(dC_im, main = "Intensity TYPE_C",    axes = TRUE); contour(dC_im, add = TRUE)
plot(dU_im, main = "Intensity uncovered", axes = TRUE); contour(dU_im, add = TRUE)
dev.off()
cat("Fig1d written.\n")

# Hopkins-Skellam test (-> Table 2)
# H_A << 1 indicates clustering relative to CSR
set.seed(42)
cat("--- Hopkins-Skellam tests Q1 (nsim=999) ---\n")
hop_brA <- hopskel.test(ori_brA, alternative = "clustered", method = "MonteCarlo", nsim = 999)
hop_brB <- hopskel.test(ori_brB, alternative = "clustered", method = "MonteCarlo", nsim = 999)
hop_brC <- hopskel.test(ori_brC, alternative = "clustered", method = "MonteCarlo", nsim = 999)
hop_brU <- hopskel.test(ori_brU, alternative = "clustered", method = "MonteCarlo", nsim = 999)
hop_br  <- hopskel.test(ori_br,  alternative = "clustered", method = "MonteCarlo", nsim = 999)
cat("TYPE_A:    A =", hop_brA$statistic, "p =", hop_brA$p.value, "\n")
cat("TYPE_B:    A =", hop_brB$statistic, "p =", hop_brB$p.value, "\n")
cat("TYPE_C:    A =", hop_brC$statistic, "p =", hop_brC$p.value, "\n")
cat("uncovered: A =", hop_brU$statistic, "p =", hop_brU$p.value, "\n")
cat("combined:  A =", hop_br$statistic,  "p =", hop_br$p.value,  "\n")

# Average intensity (pieces/m2)
int_A      <- intensity(ori_brA)
int_B      <- intensity(ori_brB)
int_C      <- intensity(ori_brC)
int_U      <- intensity(ori_brU)
int_br_all <- intensity(ori_br_u)
cat("Intensity (pieces/m2): A =", int_A, "| B =", int_B,
    "| C =", int_C, "| uncovered =", int_U, "| all =", int_br_all, "\n")

# Fig 2e: Hotspot maps — scanLRTS (4-panel)
pdf("Fig1e_BrecciaHotspots.pdf", height = 12, width = 12)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
LR_A <- scanLRTS(ori_brA, r = bw.ppl(ori_brA))
plot(LR_A, main = "Hot spots TYPE_A",    pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
LR_B <- scanLRTS(ori_brB, r = bw.ppl(ori_brB))
plot(LR_B, main = "Hot spots TYPE_B",    pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
LR_C <- scanLRTS(ori_brC, r = bw.ppl(ori_brC))
plot(LR_C, main = "Hot spots TYPE_C",    pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
LR_U <- scanLRTS(ori_brU, r = bw.ppl(ori_brU))
plot(LR_U, main = "Hot spots uncovered", pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
dev.off()
cat("Fig1e written.\n")

# Scan test (-> Table 2)
set.seed(42)
sc_brA <- scan.test(ori_brA, r = bw.ppl(ori_brA), nsim = 999, verbose = FALSE)
sc_brB <- scan.test(ori_brB, r = bw.ppl(ori_brB), nsim = 999, verbose = FALSE)
sc_brC <- scan.test(ori_brC, r = bw.ppl(ori_brC), nsim = 999, verbose = FALSE)
sc_brU <- scan.test(ori_brU, r = bw.ppl(ori_brU), nsim = 999, verbose = FALSE)
cat("Scan LR: A =", sc_brA$statistic, "p =", sc_brA$p.value,
    "| B =", sc_brB$statistic, "p =", sc_brB$p.value,
    "| C =", sc_brC$statistic, "p =", sc_brC$p.value,
    "| U =", sc_brU$statistic, "p =", sc_brU$p.value, "\n")

# Fig 2f: Adaptive relative risk
h0_br <- OS(ori_br, nstar = "geometric")
pdf("Fig1f_BrecciaRelRisk.pdf", height = 12, width = 12)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
rr_BA <- risk(ori_brB, ori_brA, h0 = h0_br, adapt = TRUE, tolerate = TRUE,
              hp = OS(ori_br) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
plot(rr_BA, main = "B vs A")
rr_CA <- risk(ori_brC, ori_brA, h0 = h0_br, adapt = TRUE, tolerate = TRUE,
              hp = OS(ori_br) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
plot(rr_CA, main = "C vs A")
rr_CB <- risk(ori_brC, ori_brB, h0 = h0_br, adapt = TRUE, tolerate = TRUE,
              hp = OS(ori_br) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
plot(rr_CB, main = "C vs B")
rr_UA <- risk(ori_brU, ori_brA, h0 = h0_br, adapt = TRUE, tolerate = TRUE,
              hp = OS(ori_br) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
plot(rr_UA, main = "uncovered vs A")
dev.off()
cat("Fig1f written.\n")

# Segregation test (-> Table 2)
set.seed(42)
ss_br <- segregation.test(ori_br, nsim = 99, verbose = FALSE, sigma = bw.ppl(ori_br_u))
cat("Segregation T =", ss_br$statistic, "p =", ss_br$p.value, "\n")

# Fig 2g: Dominant cementation type map
# Uses pixel-matrix argmax; listof + sapply is unreliable with im objects
pdf("Fig1g_BrecciaDominantType.pdf", height = 10, width = 10)
par(mar = c(4, 4, 2, 4))
dom_br <- dominant_type_map(
  list(dA_im, dB_im, dC_im, dU_im),
  c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")
)
plot(dom_br,
     col  = pal_br[c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")],
     main = "Dominant cementation type", axes = TRUE)
dev.off()
cat("Fig1g written.\n")

# Supplementary: Kinhom per cementation type (nsim=39)
set.seed(42)
for (tp in c("A", "B", "C", "U")) {
  obj   <- get(paste0("ori_br", tp))
  lab   <- if (tp == "U") "uncovered" else paste0("TYPE_", tp)
  fname <- paste0("FigS_BrecciaKinhom_", tp, ".pdf")
  den_tp <- density(obj, bw.ppl(obj), positive = TRUE)
  pdf(fname, height = 9, width = 9)
  par(mar = c(4, 4, 2.5, 2))
  Eb <- envelope(obj, modKinhom,
                 simulate = expression(rpoispp(den_tp)), nsim = 39, global = TRUE)
  plot(Eb, lwd = 2, main = paste("Kinhom —", lab))
  dev.off()
  cat(fname, "written.\n")
}

# Fig 2h: nnequal — cementation (rlabel envelope, nsim=39)
set.seed(42)
pdf("Fig1h_BrecciaNNequal.pdf", height = 9, width = 9)
par(mar = c(4, 4, 2.5, 2))
nneq_br <- envelope(ori_br, nnequal,
                    simulate = expression(rlabel(ori_br)), nsim = 39)
plot(nneq_br, lwd = 2, main = "Nearest neighbour equality — cementation types")
dev.off()
cat("Fig1h written.\n")

# Collect Q1 combined stats
main_stats <- rbind(main_stats, data.frame(
  Question     = "Q1",
  Category     = c("TYPE_A","TYPE_B","TYPE_C","uncovered","combined"),
  n            = c(npoints(ori_brA),npoints(ori_brB),npoints(ori_brC),npoints(ori_brU),npoints(ori_br)),
  intensity_m2 = c(int_A, int_B, int_C, int_U, int_br_all),
  HopSkel_A    = c(hop_brA$statistic,hop_brB$statistic,hop_brC$statistic,hop_brU$statistic,hop_br$statistic),
  HopSkel_p    = c(hop_brA$p.value,  hop_brB$p.value,  hop_brC$p.value,  hop_brU$p.value, hop_br$p.value),
  ScanLR       = c(sc_brA$statistic, sc_brB$statistic, sc_brC$statistic, sc_brU$statistic, NA),
  Scan_p       = c(sc_brA$p.value,   sc_brB$p.value,   sc_brC$p.value,   sc_brU$p.value,  NA),
  Segreg_T     = c(rep(NA, 4), ss_br$statistic),
  Segreg_p     = c(rep(NA, 4), ss_br$p.value),
  stringsAsFactors = FALSE
))
cat("Section 3 complete.\n")


# ==============================================================================
# SECTION 4 — Q1: Cementation Gradient — Section View xz (Combined)
# ==============================================================================

cat("\n=== Section 4: Q1 Cementation gradient (section view xz) ===\n")

# FigS: Spatial distribution xz
pdf("FigS_BrecciaDist_section.pdf", height = 10, width = 10)
par(mar = c(4, 4, 1.5, 2))
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
pdf("FigS_BrecciaDensityZ.pdf", height = 8, width = 10); print(p_denz); dev.off()
cat("FigS_BrecciaDensityZ written.\n")

# FigS: Kernel intensity per type — section view
dA_xz <- density(ori_brA_xz, sigma = bw.ppl(ori_brA_xz), eps = 0.01)
dB_xz <- density(ori_brB_xz, sigma = bw.ppl(ori_brB_xz), eps = 0.01)
dC_xz <- density(ori_brC_xz, sigma = bw.ppl(ori_brC_xz), eps = 0.01)
dU_xz <- density(ori_brU_xz, sigma = bw.ppl(ori_brU_xz), eps = 0.01)
pdf("FigS_BrecciaIntensity_section.pdf", height = 12, width = 12)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
plot(dA_xz, main = "Intensity TYPE_A (xz)",    axes = TRUE); contour(dA_xz, add = TRUE)
plot(dB_xz, main = "Intensity TYPE_B (xz)",    axes = TRUE); contour(dB_xz, add = TRUE)
plot(dC_xz, main = "Intensity TYPE_C (xz)",    axes = TRUE); contour(dC_xz, add = TRUE)
plot(dU_xz, main = "Intensity uncovered (xz)", axes = TRUE); contour(dU_xz, add = TRUE)
dev.off()
cat("FigS_BrecciaIntensity_section written.\n")

# FigS: Dominant type — section view
pdf("FigS_BrecciaDominantType_section.pdf", height = 10, width = 10)
par(mar = c(4, 4, 2, 4))
dom_br_xz <- dominant_type_map(
  list(dA_xz, dB_xz, dC_xz, dU_xz),
  c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")
)
plot(dom_br_xz,
     col  = pal_br[c("TYPE_A", "TYPE_B", "TYPE_C", "uncovered")],
     main = "Dominant cementation type (xz section)", axes = TRUE)
dev.off()
cat("FigS_BrecciaDominantType_section written.\n")
cat("Section 4 complete.\n")


# ==============================================================================
# SECTION 5 — Q2: Bones vs Lithics — Plan View (Combined)
# ==============================================================================

cat("\n=== Section 5: Q2 Bones vs Lithics (plan view, combined) ===\n")

ori_mat_u    <- unmark(ori_mat)
db_mat_clean <- subset(db_combined, !is.na(code))

# Fig 1a: Spatial distribution
pdf("Fig2a_BoneLithicDist.pdf", height = 10, width = 10)
par(mar = c(4, 4, 1.5, 2))
plot(ori_mat, cols = pal_mat[levels(marks(ori_mat))],
     pch = 19, cex = 1.2, main = "", axes = TRUE)
dev.off()
cat("Fig2a written.\n")

# Fig 1b: Density curves along x and y axes
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
pdf("Fig2b_BoneLithicDensityX.pdf", height = 8, width = 10); print(p_matx); dev.off()
pdf("Fig2b_BoneLithicDensityY.pdf", height = 8, width = 10); print(p_maty); dev.off()
cat("Fig2b written.\n")

# Fig 1c: Overall kernel intensity
den_mat_all <- density(ori_mat_u, bw.ppl(ori_mat_u), eps = 0.01, positive = TRUE)
pdf("Fig2c_BoneLithicIntensity_all.pdf", height = 10, width = 10)
par(mar = c(4, 4, 2, 4))
plot(den_mat_all, axes = TRUE, main = "Kernel smoothed intensity — all materials")
dev.off()
cat("Fig2c written.\n")

# Fig 1d: Kernel intensity per material type
pdf("Fig2d_BoneLithicIntensity_types.pdf", height = 12, width = 8)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(c(1, 2), nrow = 2, ncol = 1))
plot(density(ori_bone,   sigma = bw.ppl(ori_bone),   eps = 0.01),
     main = "Intensity BONE",   axes = TRUE)
plot(density(ori_lithic, sigma = bw.ppl(ori_lithic), eps = 0.01),
     main = "Intensity LITHIC", axes = TRUE)
dev.off()
cat("Fig2d written.\n")

# Hopkins-Skellam (-> Table 2)
set.seed(42)
cat("--- Hopkins-Skellam Q2 (nsim=999) ---\n")
hop_bone   <- hopskel.test(ori_bone,   alternative = "clustered", method = "MonteCarlo", nsim = 999)
hop_lithic <- hopskel.test(ori_lithic, alternative = "clustered", method = "MonteCarlo", nsim = 999)
hop_mat    <- hopskel.test(ori_mat,    alternative = "clustered", method = "MonteCarlo", nsim = 999)
cat("BONE:     A =", hop_bone$statistic,   "p =", hop_bone$p.value, "\n")
cat("LITHIC:   A =", hop_lithic$statistic, "p =", hop_lithic$p.value, "\n")
cat("combined: A =", hop_mat$statistic,    "p =", hop_mat$p.value, "\n")

int_bone    <- intensity(ori_bone)
int_lithic  <- intensity(ori_lithic)
int_mat_all <- intensity(ori_mat_u)
cat("Intensity: BONE =", int_bone, "| LITHIC =", int_lithic,
    "| all =", int_mat_all, "\n")

# Fig 1e: Hotspot maps
pdf("Fig2e_BoneLithicHotspots.pdf", height = 12, width = 8)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(c(1, 2), nrow = 2))
LR_bone   <- scanLRTS(ori_bone,   r = bw.ppl(ori_bone))
LR_lithic <- scanLRTS(ori_lithic, r = bw.ppl(ori_lithic))
plot(LR_bone,   main = "Hot spots BONE",   pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
plot(LR_lithic, main = "Hot spots LITHIC", pch = 19, cex = 0.2, axes = TRUE)
plot(win_xy, add = TRUE, border = "white")
dev.off()
cat("Fig2e written.\n")

# Scan tests (-> Table 2)
set.seed(42)
sc_bone   <- scan.test(ori_bone,   r = bw.ppl(ori_bone),   nsim = 999, verbose = FALSE)
sc_lithic <- scan.test(ori_lithic, r = bw.ppl(ori_lithic), nsim = 999, verbose = FALSE)
cat("Scan LR: BONE =", sc_bone$statistic, "p =", sc_bone$p.value,
    "| LITHIC =", sc_lithic$statistic, "p =", sc_lithic$p.value, "\n")

# Fig 1f: Adaptive relative risk
h0_mat <- OS(ori_mat, nstar = "geometric")
pdf("Fig2f_BoneLithicRelRisk.pdf", height = 12, width = 8)
par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
layout(matrix(c(1, 2), nrow = 2))
rr_BL <- risk(ori_bone, ori_lithic, h0 = h0_mat, adapt = TRUE, tolerate = TRUE,
              hp = OS(ori_mat) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
plot(rr_BL, main = "BONE vs LITHIC")
rr_LB <- risk(ori_lithic, ori_bone, h0 = h0_mat, adapt = TRUE, tolerate = TRUE,
              hp = OS(ori_mat) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
plot(rr_LB, main = "LITHIC vs BONE")
dev.off()
cat("Fig2f written.\n")

# Segregation test (-> Table 2)
set.seed(42)
ss_mat <- segregation.test(ori_mat, nsim = 99, verbose = FALSE, sigma = bw.ppl(ori_mat_u))
cat("Segregation T =", ss_mat$statistic, "p =", ss_mat$p.value, "\n")

# Supplementary: Kinhom per material type (nsim=39)
set.seed(42)
for (tp in c("bone", "lithic")) {
  obj    <- get(paste0("ori_", tp))
  lab    <- toupper(tp)
  fname  <- paste0("FigS_BoneLithicKinhom_", lab, ".pdf")
  den_tp <- density(obj, bw.ppl(obj), positive = TRUE)
  pdf(fname, height = 10, width = 10)
  Eb <- envelope(obj, modKinhom,
                 simulate = expression(rpoispp(den_tp)), nsim = 39, global = TRUE)
  plot(Eb, lwd = 4, main = paste("Kinhom --", lab))
  dev.off()
  cat(fname, "written.\n")
}

# Fig 1g: Cross-K random labelling
# Null hypothesis: labels are randomly distributed across observed locations (mixing null).
# Tests whether BONE and LITHIC are interchangeable with respect to spatial pattern.
lev_mat        <- levels(marks(ori_mat))
myKcross_RL    <- function(Y, i, j, ...) {
  lam <- density.ppp(Y, sigma = bw.ppl, adjust = 2)
  Kcross.inhom(Y, i = i, j = j, lambdaI = lam, lambdaJ = lam, ...)
}
sim_RL_mat     <- function(Y, ...) {
  Z <- rlabel(Y); marks(Z) <- factor(marks(Z), levels = lev_mat); Z
}
set.seed(42)
pdf("Fig2g_CrossK_RL.pdf", height = 9, width = 9)
par(mar = c(4, 4, 2.5, 2))
Ecross_RL <- envelope(ori_mat, fun = myKcross_RL, nsim = 39,
                      i = "BONE", j = "LITHIC", simulate = sim_RL_mat)
plot(Ecross_RL, lwd = 2,
     main = "Kcross BONE–LITHIC (random labelling null: mixing)")
dev.off()
cat("Fig2g written.\n")

# FigS: Cross-K independence
# Null hypothesis: bones and lithics are spatially independent Poisson processes,
# each with its own intensity surface (estimated separately from the data).
W_mat  <- Window(ori_mat)
nB     <- npoints(ori_bone)
nL     <- npoints(ori_lithic)
# Clamp density images to >= 0 (rpoint requires non-negative probabilities)
# Note: eval.im(pmax(density.ppp(...), 0)) fails in spatstat >= 3.x when the
# density call is nested — compute image first, then clamp pixel matrix directly.
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
set.seed(42)
pdf("FigS_CrossK_indep.pdf", height = 9, width = 9)
par(mar = c(4, 4, 2.5, 2))
Ecross_indep <- envelope(ori_mat, fun = myKcross_indep, nsim = 39,
                         i = "BONE", j = "LITHIC", simulate = sim_indep_mat)
plot(Ecross_indep, lwd = 2, main = "Kcross BONE–LITHIC (independent processes null)")
dev.off()
cat("FigS_CrossK_indep written.\n")

# Fig 1h: nnequal — material type (rlabel envelope, nsim=39)
set.seed(42)
pdf("Fig2h_BoneLithicNNequal.pdf", height = 9, width = 9)
par(mar = c(4, 4, 2.5, 2))
nneq_mat <- envelope(ori_mat, nnequal,
                     simulate = expression(rlabel(ori_mat)), nsim = 39)
plot(nneq_mat, lwd = 2, main = "Nearest neighbour equality — BONE vs LITHIC")
dev.off()
cat("Fig2h written.\n")

# Collect Q2 combined stats
main_stats <- rbind(main_stats, data.frame(
  Question     = "Q2",
  Category     = c("BONE","LITHIC","combined"),
  n            = c(npoints(ori_bone),npoints(ori_lithic),npoints(ori_mat)),
  intensity_m2 = c(int_bone, int_lithic, int_mat_all),
  HopSkel_A    = c(hop_bone$statistic,  hop_lithic$statistic, hop_mat$statistic),
  HopSkel_p    = c(hop_bone$p.value,    hop_lithic$p.value,   hop_mat$p.value),
  ScanLR       = c(sc_bone$statistic,   sc_lithic$statistic,  NA),
  Scan_p       = c(sc_bone$p.value,     sc_lithic$p.value,    NA),
  Segreg_T     = c(NA, NA, ss_mat$statistic),
  Segreg_p     = c(NA, NA, ss_mat$p.value),
  stringsAsFactors = FALSE
))
cat("Section 5 complete.\n")

# Checkpoint: save all analysis objects so the plot-only script can reload them
save.image(file.path(OUTPUT_DIR, "analysis_checkpoint.RData"))
cat("Checkpoint saved: analysis_checkpoint.RData\n")


# ==============================================================================
# SECTION 6 — Q3: Cross-tabulation Cementation x Material Type
# ==============================================================================

cat("\n=== Section 6: Q3 cross-tabulation ===\n")

ct <- table(db_combined$code, db_combined$breccia_type)
cat("\nContingency table (code x breccia_type):\n"); print(ct)
cat("\nRow proportions (within material type):\n")
print(round(prop.table(ct, margin = 1), 3))
cat("\nColumn proportions (within cementation type):\n")
print(round(prop.table(ct, margin = 2), 3))

# Chi-square or Fisher's exact test
exp_counts <- chisq.test(ct)$expected
if (any(exp_counts < 5)) {
  cat("\nSome expected counts < 5; using Fisher's exact test.\n")
  ct_test <- fisher.test(ct, simulate.p.value = TRUE, B = 9999)
  print(ct_test)
} else {
  ct_test <- chisq.test(ct)
  cat("\nChi-square test:\n"); print(ct_test)
}

# Write Table 3
ct_df       <- as.data.frame.matrix(ct)
ct_df$code  <- rownames(ct_df)
write.csv(ct_df, "Table2_Contingency.csv", row.names = FALSE)
cat("Table2_Contingency.csv written.\n")

# Fig 3: Contingency bar chart (proportions stacked)
db_ct <- db_combined %>%
  filter(!is.na(code), !is.na(breccia_type)) %>%
  group_by(code, breccia_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(code) %>%
  mutate(prop = n / sum(n))

p_ct <- ggplot(db_ct, aes(x = code, y = prop, fill = breccia_type)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = pal_br, name = "Cementation type") +
  labs(x = "Material type", y = "Proportion",
       title = "Cementation type by material type") +
  theme_bw(base_size = 14)
pdf("Fig3_ContingencyBarplot.pdf", height = 8, width = 8); print(p_ct); dev.off()
cat("Fig3_ContingencyBarplot written.\n")

# Per-level cross-tabulations (supplementary)
cat("\n--- Per-level contingency tables ---\n")
for (lv in c("2", "3", "4")) {
  dl <- subset(db_combined, level == lv)
  cat("\nLevel", lv, ":\n"); print(table(dl$code, dl$breccia_type))
}
cat("Section 6 complete.\n")


# ==============================================================================
# SECTION 7 — Per-Level Analyses Q1 + Q2 (Supplementary)
# ==============================================================================

cat("\n=== Section 7: Per-level analyses (supplementary) ===\n")

supp_stats  <- data.frame()
level_list  <- list("2" = db2, "3" = db3, "4" = db4)
level_short <- list("2" = "Lv2", "3" = "Lv3", "4" = "Lv4")

for (lv in names(level_list)) {
  dl     <- level_list[[lv]]
  lshort <- level_short[[lv]]
  cat("\n--- Level", lv, "(n =", nrow(dl), ") ---\n")
  # Per-level plan view window: shared win_xy (same excavation footprint all levels)
  # Per-level section view window: level-specific convex hull in xz
  win_xz_lv <- make_win_xz_lv(dl)
  pp <- make_level_ppp(dl, win_xy)

  # ==== Q1: Breccia per level ====
  cat("  Q1: A =", npoints(pp$brA), "B =", npoints(pp$brB),
      "C =", npoints(pp$brC), "U =", npoints(pp$brU), "\n")

  if (npoints(pp$brA) >= 5 && npoints(pp$brB) >= 5) {

    pdf(paste0("FigS_", lshort, "_BrecciaDist.pdf"), height = 10, width = 10)
    par(mar = c(4, 4, 2, 2))
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
    pdf(paste0("FigS_", lshort, "_BrecciaDensityXY.pdf"), height = 10, width = 10)
    print(p_bxy); dev.off()

    pdf(paste0("FigS_", lshort, "_BrecciaIntensity.pdf"), height = 12, width = 12)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
    layout(matrix(1:4, nrow = 2, byrow = TRUE))
    for (tp in c("brA","brB","brC","brU")) {
      obj_tp <- pp[[tp]]
      lab_tp <- switch(tp, brA="TYPE_A", brB="TYPE_B", brC="TYPE_C", brU="uncovered")
      if (npoints(obj_tp) >= 3) {
        plot(density(obj_tp, sigma = bw.ppl(obj_tp), eps = 0.01),
             main = paste("Intensity", lab_tp, lv), axes = TRUE)
      } else { plot.new(); title(paste("Too few points:", lab_tp)) }
    }; dev.off()

    pdf(paste0("FigS_", lshort, "_BrecciaHotspots.pdf"), height = 12, width = 12)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
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

    set.seed(42)
    hop_A_lv <- tryCatch(hopskel.test(pp$brA, alternative="clustered", method="MonteCarlo", nsim=999), error=function(e) NULL)
    hop_B_lv <- tryCatch(hopskel.test(pp$brB, alternative="clustered", method="MonteCarlo", nsim=999), error=function(e) NULL)
    hop_C_lv <- tryCatch(hopskel.test(pp$brC, alternative="clustered", method="MonteCarlo", nsim=999), error=function(e) NULL)
    sc_A_lv  <- tryCatch(scan.test(pp$brA, r=bw.ppl(pp$brA), nsim=999, verbose=FALSE), error=function(e) NULL)
    sc_B_lv  <- tryCatch(scan.test(pp$brB, r=bw.ppl(pp$brB), nsim=999, verbose=FALSE), error=function(e) NULL)
    sc_C_lv  <- tryCatch(scan.test(pp$brC, r=bw.ppl(pp$brC), nsim=999, verbose=FALSE), error=function(e) NULL)

    if (npoints(pp$brA) >= 10 && npoints(pp$brB) >= 10) {
      h0_lv <- OS(pp$br, nstar = "geometric")
      pdf(paste0("FigS_", lshort, "_BrecciaRelRisk.pdf"), height = 14, width = 7)
      par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
      layout(matrix(1:3, nrow = 3))
      for (pair in list(c("brB","brA","B vs A"), c("brC","brA","C vs A"), c("brC","brB","C vs B"))) {
        tryCatch({
          rr <- risk(pp[[pair[1]]], pp[[pair[2]]], h0 = h0_lv, adapt = TRUE, tolerate = TRUE,
                     hp = OS(pp$br) / 2, pilot.symmetry = "pooled", davies.baddeley = 0.05)
          plot(rr, main = paste(pair[3], "--", lv))
        }, error = function(e) { plot.new(); title(paste("Error:", pair[3])) })
      }; dev.off()
    }

    ss_br_lv <- tryCatch(segregation.test(pp$br, nsim=99, verbose=FALSE), error=function(e) NULL)

    pdf(paste0("FigS_", lshort, "_BrecciaNNequal.pdf"), height = 9, width = 9)
    par(mar = c(4, 4, 2.5, 2))
    tryCatch({
      nneq_br_lv <- envelope(pp$br, nnequal,
                             simulate = expression(rlabel(pp$br)), nsim = 39)
      plot(nneq_br_lv, lwd = 2, main = paste("NNequal cementation —", lv))
    }, error = function(e) { plot.new(); title(paste("NNequal error:", lv)) })
    dev.off()

    for (tp in c("A", "B", "C")) {
      hop_tp <- get(paste0("hop_", tp, "_lv"))
      sc_tp  <- get(paste0("sc_",  tp, "_lv"))
      supp_stats <- rbind(supp_stats, data.frame(
        Level    = lv, Question = "Q1", Category = paste0("TYPE_", tp),
        n        = npoints(pp[[paste0("br", tp)]]),
        HopSkel_A = if (!is.null(hop_tp)) hop_tp$statistic else NA,
        HopSkel_p = if (!is.null(hop_tp)) hop_tp$p.value   else NA,
        ScanLR    = if (!is.null(sc_tp))  sc_tp$statistic  else NA,
        Scan_p    = if (!is.null(sc_tp))  sc_tp$p.value    else NA,
        Segreg_T  = if (tp == "C" && !is.null(ss_br_lv)) ss_br_lv$statistic else NA,
        Segreg_p  = if (tp == "C" && !is.null(ss_br_lv)) ss_br_lv$p.value   else NA,
        stringsAsFactors = FALSE
      ))
    }
    cat("  Q1 done for level", lv, "\n")
  } else cat("  Skipping Q1 level", lv, "-- insufficient points.\n")

  # ==== Q2: Bones vs Lithics per level ====
  cat("  Q2: BONE =", npoints(pp$bone), "LITHIC =", npoints(pp$lithic), "\n")

  if (npoints(pp$bone) >= 5 && npoints(pp$lithic) >= 5) {

    pdf(paste0("FigS_", lshort, "_BoneLithicDist.pdf"), height = 10, width = 10)
    par(mar = c(4, 4, 2, 2))
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
    pdf(paste0("FigS_", lshort, "_BoneLithicDensityXY.pdf"), height = 10, width = 10)
    print(p_mxy); dev.off()

    pdf(paste0("FigS_", lshort, "_BoneLithicIntensity.pdf"), height = 12, width = 8)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
    layout(matrix(c(1, 2), nrow = 2))
    plot(density(pp$bone,   sigma = bw.ppl(pp$bone),   eps = 0.01),
         main = paste("Intensity BONE --",   lv), axes = TRUE)
    plot(density(pp$lithic, sigma = bw.ppl(pp$lithic), eps = 0.01),
         main = paste("Intensity LITHIC --", lv), axes = TRUE)
    dev.off()

    pdf(paste0("FigS_", lshort, "_BoneLithicHotspots.pdf"), height = 12, width = 8)
    par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
    layout(matrix(c(1, 2), nrow = 2))
    LR_b <- scanLRTS(pp$bone,   r = bw.ppl(pp$bone))
    LR_l <- scanLRTS(pp$lithic, r = bw.ppl(pp$lithic))
    plot(LR_b, main = paste("Hot spots BONE --",   lv), pch = 19, cex = 0.2, axes = TRUE)
    plot(win_xy, add = TRUE, border = "white")
    plot(LR_l, main = paste("Hot spots LITHIC --", lv), pch = 19, cex = 0.2, axes = TRUE)
    plot(win_xy, add = TRUE, border = "white")
    dev.off()

    set.seed(42)
    hop_bone_lv   <- tryCatch(hopskel.test(pp$bone,   alternative="clustered", method="MonteCarlo", nsim=999), error=function(e) NULL)
    hop_lithic_lv <- tryCatch(hopskel.test(pp$lithic, alternative="clustered", method="MonteCarlo", nsim=999), error=function(e) NULL)
    sc_bone_lv    <- tryCatch(scan.test(pp$bone,   r=bw.ppl(pp$bone),   nsim=999, verbose=FALSE), error=function(e) NULL)
    sc_lithic_lv  <- tryCatch(scan.test(pp$lithic, r=bw.ppl(pp$lithic), nsim=999, verbose=FALSE), error=function(e) NULL)

    if (npoints(pp$bone) >= 10 && npoints(pp$lithic) >= 10) {
      h0_ml <- OS(pp$mat, nstar = "geometric")
      pdf(paste0("FigS_", lshort, "_BoneLithicRelRisk.pdf"), height = 12, width = 8)
      par(mar = c(3, 3, 2, 3.5), oma = c(0.3, 0.3, 0.3, 0.3), cex.main = 1)
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

    ss_mat_lv <- tryCatch(segregation.test(pp$mat, nsim=99, verbose=FALSE), error=function(e) NULL)

    # Cross-K random labelling — only for levels with n_lithic >= 30
    if (npoints(pp$lithic) >= 30) {
      lev_ml <- levels(marks(pp$mat))
      myKcross_RL_lv <- function(Y, i, j, ...) {
        lam <- density.ppp(Y, sigma = bw.ppl, adjust = 2)
        Kcross.inhom(Y, i=i, j=j, lambdaI=lam, lambdaJ=lam, ...)
      }
      sim_RL_lv <- function(Y, ...) {
        Z <- rlabel(Y); marks(Z) <- factor(marks(Z), levels = lev_ml); Z
      }
      pdf(paste0("FigS_", lshort, "_CrossK_RL.pdf"), height = 9, width = 9)
      par(mar = c(4, 4, 2.5, 2))
      tryCatch({
        Ecrl <- envelope(pp$mat, fun=myKcross_RL_lv, nsim=39,
                         i="BONE", j="LITHIC", simulate=sim_RL_lv)
        plot(Ecrl, lwd=2, main = paste("Kcross BONE–LITHIC (RL) —", lv))
      }, error = function(e) { plot.new(); title(paste("CrossK error:", lv)) })
      dev.off()
    }

    pdf(paste0("FigS_", lshort, "_BoneLithicNNequal.pdf"), height = 9, width = 9)
    par(mar = c(4, 4, 2.5, 2))
    tryCatch({
      nneq_lv <- envelope(pp$mat, nnequal,
                          simulate = expression(rlabel(pp$mat)), nsim = 39)
      plot(nneq_lv, lwd = 2, main = paste("NNequal BONE/LITHIC —", lv))
    }, error = function(e) { plot.new(); title(paste("NNequal error:", lv)) })
    dev.off()

    for (tp in c("bone", "lithic")) {
      hop_tp <- if (tp == "bone") hop_bone_lv else hop_lithic_lv
      sc_tp  <- if (tp == "bone") sc_bone_lv  else sc_lithic_lv
      supp_stats <- rbind(supp_stats, data.frame(
        Level    = lv, Question = "Q2", Category = toupper(tp),
        n        = npoints(pp[[tp]]),
        HopSkel_A = if (!is.null(hop_tp)) hop_tp$statistic else NA,
        HopSkel_p = if (!is.null(hop_tp)) hop_tp$p.value   else NA,
        ScanLR    = if (!is.null(sc_tp))  sc_tp$statistic  else NA,
        Scan_p    = if (!is.null(sc_tp))  sc_tp$p.value    else NA,
        Segreg_T  = if (tp == "lithic" && !is.null(ss_mat_lv)) ss_mat_lv$statistic else NA,
        Segreg_p  = if (tp == "lithic" && !is.null(ss_mat_lv)) ss_mat_lv$p.value   else NA,
        stringsAsFactors = FALSE
      ))
    }
    cat("  Q2 done for level", lv, "\n")
  } else cat("  Skipping Q2 level", lv, "-- insufficient points.\n")
}

write.csv(supp_stats, "TableS_PerLevel_Stats.csv", row.names = FALSE)
cat("\nTableS_PerLevel_Stats.csv written (", nrow(supp_stats), "rows )\n")
cat("Section 7 complete.\n")


# ==============================================================================
# SECTION 8 — Main Text Statistics Tables
# ==============================================================================

cat("\n=== Section 8: Compile main text tables ===\n")
write.csv(main_stats, "Table1_CombinedStats.csv", row.names = FALSE)
cat("Table1_CombinedStats.csv written.\n")
print(main_stats)
cat("\n*** GDC_SPATIAL_ANALYSIS_MAIN.R complete ***\n")

# Save session information for reproducibility records
sink(file.path(OUTPUT_DIR, "session_info.txt"))
cat("Session info for GDC Spatial Analysis\n")
cat("Generated:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()
cat("session_info.txt written to output/\n")
