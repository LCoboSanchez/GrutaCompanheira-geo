# ==============================================================================
# 03_composites.R — Build all composite figures from individual PDFs
# Reads PDFs from output/, assembles multi-panel PNG composites, and
# saves them to output/figures/composite/.
# Run after 02_replot.R (or 01_analysis.R) has completed.
# ==============================================================================

# ==============================================================================
# CONFIGURATION — paths derive automatically from the project root (.Rproj)
# ==============================================================================
library(here)
OUTPUT_DIR    <- here("output")
COMPOSITE_DIR <- here("output", "figures", "composite")
dir.create(COMPOSITE_DIR, showWarnings = FALSE, recursive = TRUE)
# ==============================================================================

library(magick)
library(cowplot)
library(ggplot2)

# ── Helper: read one PDF page, trim white border, add generous padding ───────
panel <- function(fname) {
  path <- file.path(OUTPUT_DIR, fname)
  if (!file.exists(path)) {
    cat("MISSING:", fname, "\n")
    return(ggdraw() + ggplot2::theme(plot.background =
             ggplot2::element_rect(fill = "grey90", colour = NA)))
  }
  img <- magick::image_read_pdf(path, density = 200)[1]
  img <- magick::image_trim(img, fuzz = 8)
  img <- magick::image_border(img, "white", "28x28")
  ggdraw() + draw_image(img) +
    ggplot2::theme(plot.background =
                     ggplot2::element_rect(fill = "white", colour = NA))
}

# ── Helper: assemble and save composite ─────────────────────────────────────
make_grid <- function(panels, nc, nr, labels, w, h, fname, lsize = 40) {
  fig <- plot_grid(plotlist = panels, ncol = nc, nrow = nr,
                   labels = labels, label_size = lsize,
                   label_fontface = "bold", label_colour = "black",
                   hjust = -0.2, vjust = 1.3)
  fig <- ggdraw(fig) +
    ggplot2::theme(plot.background =
                     ggplot2::element_rect(fill = "white", colour = NA))
  out <- file.path(COMPOSITE_DIR, fname)
  ggsave(out, fig, width = w, height = h, units = "in", dpi = 150, bg = "white")
  cat("Saved:", fname, "\n")
}


# ==============================================================================
# MAIN FIGURE 1 — Cementation gradient, plan view (3×3, A–I)
# ==============================================================================
cat("\n=== Building Fig1_cementation ===\n")
panels1 <- lapply(c(
  "Fig1a_BrecciaDist_plan.pdf",
  "Fig1b_BrecciaDensityX.pdf",
  "Fig1b_BrecciaDensityY.pdf",
  "Fig1c_BrecciaIntensity_all.pdf",
  "Fig1d_BrecciaIntensity_types.pdf",
  "Fig1e_BrecciaHotspots.pdf",
  "Fig1f_BrecciaRelRisk.pdf",
  "Fig1g_BrecciaDominantType.pdf",
  "Fig1h_BrecciaNNequal.pdf"
), panel)
make_grid(panels1, 3, 3, LETTERS[1:9], 36, 36, "Fig1_cementation.png")


# ==============================================================================
# MAIN FIGURE 2 — Bones vs Lithics, plan view (3×3, A–I)
# ==============================================================================
cat("\n=== Building Fig2_materials ===\n")
panels2 <- lapply(c(
  "Fig2a_BoneLithicDist.pdf",
  "Fig2b_BoneLithicDensityX.pdf",
  "Fig2b_BoneLithicDensityY.pdf",
  "Fig2c_BoneLithicIntensity_all.pdf",
  "Fig2d_BoneLithicIntensity_types.pdf",
  "Fig2e_BoneLithicHotspots.pdf",
  "Fig2f_BoneLithicRelRisk.pdf",
  "Fig2g_CrossK_RL.pdf",
  "Fig2h_BoneLithicNNequal.pdf"
), panel)
make_grid(panels2, 3, 3, LETTERS[1:9], 36, 36, "Fig2_materials.png")


# ==============================================================================
# SUPPLEMENTARY — Cementation sections (2×2)
# ==============================================================================
cat("\n=== Building FigS_cementation_sections ===\n")
panelsS_sec <- lapply(c(
  "FigS_BrecciaDist_section.pdf",
  "FigS_BrecciaDensityZ.pdf",
  "FigS_BrecciaIntensity_section.pdf",
  "FigS_BrecciaDominantType_section.pdf"
), panel)
make_grid(panelsS_sec, 2, 2, LETTERS[1:4], 24, 24, "FigS_cementation_sections.png")


# ==============================================================================
# SUPPLEMENTARY — Density lines, cementation (3×1: X, Y, Z)
# ==============================================================================
cat("\n=== Building FigS_cementation_density ===\n")
panelsS_dc <- lapply(c(
  "Fig1b_BrecciaDensityX.pdf",
  "Fig1b_BrecciaDensityY.pdf",
  "FigS_BrecciaDensityZ.pdf"
), panel)
make_grid(panelsS_dc, 3, 1, LETTERS[1:3], 30, 11, "FigS_cementation_density.png")


# ==============================================================================
# SUPPLEMENTARY — Density lines, Bone/Lithic (2×1: X, Y)
# ==============================================================================
cat("\n=== Building FigS_materials_density ===\n")
panelsS_dm <- lapply(c(
  "Fig2b_BoneLithicDensityX.pdf",
  "Fig2b_BoneLithicDensityY.pdf"
), panel)
make_grid(panelsS_dm, 2, 1, LETTERS[1:2], 20, 11, "FigS_materials_density.png")


# ==============================================================================
# SUPPLEMENTARY — Cementation Kinhom (2×2: A, B, C, U)
# ==============================================================================
cat("\n=== Building FigS_cementation_Kinhom ===\n")
panelsS_bk <- lapply(c(
  "FigS_BrecciaKinhom_A.pdf",
  "FigS_BrecciaKinhom_B.pdf",
  "FigS_BrecciaKinhom_C.pdf",
  "FigS_BrecciaKinhom_U.pdf"
), panel)
make_grid(panelsS_bk, 2, 2, LETTERS[1:4], 24, 24, "FigS_cementation_Kinhom.png")


# ==============================================================================
# SUPPLEMENTARY — Bone/Lithic Kinhom (2×1: BONE, LITHIC)
# ==============================================================================
cat("\n=== Building FigS_materials_Kinhom ===\n")
panelsS_mk <- lapply(c(
  "FigS_BoneLithicKinhom_BONE.pdf",
  "FigS_BoneLithicKinhom_LITHIC.pdf"
), panel)
make_grid(panelsS_mk, 2, 1, LETTERS[1:2], 24, 13, "FigS_materials_Kinhom.png")


# ==============================================================================
# SUPPLEMENTARY — Cross-K independence (standalone)
# ==============================================================================
cat("\n=== Building FigS_CrossK_independence ===\n")
panelsS_cx <- list(panel("FigS_CrossK_indep.pdf"))
make_grid(panelsS_cx, 1, 1, "", 14, 14, "FigS_CrossK_independence.png", lsize = 0)


# ==============================================================================
# SUPPLEMENTARY — Per-level analyses (4 levels × 2 questions)
# ==============================================================================

levels_info <- list(
  list(short = "Lv2", label = "Level 2", out = "Level2"),
  list(short = "Lv3", label = "Level 3", out = "Level3"),
  list(short = "Lv4", label = "Level 4", out = "Level4")
)

for (li in levels_info) {
  lshort <- li$short
  lout   <- li$out

  # ── Cementation per level (3×2: Dist, DensityXY, Intensity, Hotspots, RelRisk, NNequal)
  cat("\n=== Building FigS_", lout, "_cementation ===\n", sep = "")
  panels_br <- lapply(paste0("FigS_", lshort, "_",
                             c("BrecciaDist", "BrecciaDensityXY", "BrecciaIntensity",
                               "BrecciaHotspots", "BrecciaRelRisk", "BrecciaNNequal"),
                             ".pdf"), panel)
  make_grid(panels_br, 3, 2, LETTERS[1:6], 36, 24,
            paste0("FigS_", lout, "_cementation.png"))

  # ── Materials per level
  # Layout depends on whether CrossK_RL was produced (requires n_lithic >= 30).
  # If CrossK exists: 4×2 grid (7 panels + blank), panels A–G.
  # If CrossK missing: 3×2 grid (6 panels), panels A–F — no blank area.
  cat("=== Building FigS_", lout, "_materials ===\n", sep = "")
  crossk_file <- paste0("FigS_", lshort, "_CrossK_RL.pdf")
  has_crossk  <- file.exists(file.path(OUTPUT_DIR, crossk_file))
  if (!has_crossk) {
    cat("  NOTE: CrossK_RL not found for", lshort,
        "— using 3x2 layout (too few lithics for that level).\n")
  }
  if (has_crossk) {
    panels_bl <- lapply(paste0("FigS_", lshort, "_",
                               c("BoneLithicDist", "BoneLithicDensityXY",
                                 "BoneLithicIntensity", "BoneLithicHotspots",
                                 "BoneLithicRelRisk", "CrossK_RL",
                                 "BoneLithicNNequal"),
                               ".pdf"), panel)
    # blank 8th cell to complete 4×2 grid
    panels_bl[[8]] <- ggdraw() +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    make_grid(panels_bl, 4, 2, c(LETTERS[1:7], ""), 40, 22,
              paste0("FigS_", lout, "_materials.png"))
  } else {
    # 6-panel 3×2 grid — skip CrossK panel entirely
    panels_bl <- lapply(paste0("FigS_", lshort, "_",
                               c("BoneLithicDist", "BoneLithicDensityXY",
                                 "BoneLithicIntensity", "BoneLithicHotspots",
                                 "BoneLithicRelRisk", "BoneLithicNNequal"),
                               ".pdf"), panel)
    make_grid(panels_bl, 3, 2, LETTERS[1:6], 36, 24,
              paste0("FigS_", lout, "_materials.png"))
  }
}

cat("\n=== ALL COMPOSITE FIGURES COMPLETE ===\n")
cat("Output folder:", COMPOSITE_DIR, "\n")
