# >>> deposit header
# Purpose: Fit the six primary models of AMPKα C-terminal length and mean predicted disorder against AMPKγ
#          absence, AMPKβ absence and the retained AMPKβ module score, with Holm adjustment across the six, and
#          their declared one-dimension sensitivities, by phylogenetic least squares.
# Methods: Phylogenetic regressions
# Inputs:  `ctdx_lib.R` and the ten inputs it checks; `03_reproduction/GAMMA_ALPHA_REPRODUCTION_GATE.tsv` and
#          `03_reproduction/gamma_alpha_reproduction.tsv` (written by the reproduction step that precedes these fits)
# Outputs: `04_cross_ampk/CTD_cross_AMPK_primary.tsv`, `gamma_alpha_extension_primary.tsv`,
#          `CTD_cross_AMPK_sensitivity.tsv`, `sample_flow_cross_ampk.tsv` and `encoding_checks_cross_ampk.tsv`
# Tools:   R with ape and phylolm 2.6.5 (Pagel's lambda, Brownian motion as a sensitivity); versions: see envs/
# <<< deposit header
# Ran as: ctdx_cross_ampk.R
## deposit repo root: begin
.REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(f[1]) else "."
  normalizePath(file.path(d, ".."))
})
## deposit repo root: end
## The six cross-subunit models of the AMPKα C-terminal region (CTD), specified before the analysis: tail length
## and mean predicted disorder regressed on AMPKγ absence, AMPKβ absence and the retained AMPKβ module score
## (CORE5), with Holm adjustment across exactly these six; their declared sensitivities, each changing one
## dimension; and the declared sensitivities of the two reproduced models of tail length and disorder on the
## AMPKγ contact score (REPRODUCED_PRIMARY_EXTENSION, outside every Holm family).
source(file.path(.REPO_ROOT, "5_phylogenetic_analyses", "lib", "tail_model_lib.R"))
S <- "ctdx_cross_ampk.R"
g <- read.delim(file.path(CTDX, "03_reproduction/GAMMA_ALPHA_REPRODUCTION_GATE.tsv"), comment.char = "#")
if (g$STATUS != "PASS") stop("REFUSED: gamma-alpha reproduction did not pass")
M <- build_frames()
F6 <- "SIX_PRESPECIFIED_CTD_CROSS_AMPK_MODELS"
L <- "CTD_LENGTH_MEDIAN"; D <- "CTD_MEAN_DISORDER_MEDIAN"
prim <- list(
  spec("XA_GLOSS_LEN", "CROSS_AMPK_GAMMA_WHOLE_LOSS", "PRIMARY", L, "median_primary", "AMPKg_absent", "callable_whole_subunit_loss", holm = F6),
  spec("XA_GLOSS_DIS", "CROSS_AMPK_GAMMA_WHOLE_LOSS", "PRIMARY", D, "median_primary", "AMPKg_absent", "callable_whole_subunit_loss", holm = F6),
  spec("XA_BLOSS_LEN", "CROSS_AMPK_BETA_WHOLE_LOSS", "PRIMARY", L, "median_primary", "AMPKb_absent", "callable_whole_subunit_loss", holm = F6),
  spec("XA_BLOSS_DIS", "CROSS_AMPK_BETA_WHOLE_LOSS", "PRIMARY", D, "median_primary", "AMPKb_absent", "callable_whole_subunit_loss", holm = F6),
  spec("XA_CBM_LEN", "CROSS_AMPK_BETA_CBM_RETAINED", "PRIMARY", L, "median_primary", "CBM_CORE5_MEDIAN", "core5_median_primary", holm = F6),
  spec("XA_CBM_DIS", "CROSS_AMPK_BETA_CBM_RETAINED", "PRIMARY", D, "median_primary", "CBM_CORE5_MEDIAN", "core5_median_primary", holm = F6))
ext <- list(
  spec("GA_CONTACT_LEN", "GAMMA_ALPHA_RETAINED_CONTACT", "REPRODUCED_PRIMARY_EXTENSION", L, "median_primary", "SCORE_SITE1_MEDIAN", "site1_median_primary"),
  spec("GA_CONTACT_DIS", "GAMMA_ALPHA_RETAINED_CONTACT", "REPRODUCED_PRIMARY_EXTENSION", D, "median_primary", "SCORE_SITE1_MEDIAN", "site1_median_primary"))
sens <- list()
for (p in c(prim, ext)) sens <- c(sens, alpha_sensitivities(p))
for (p in prim[5:6]) sens <- c(sens, list(
  predictor_swap(p, "CBM_CORE5_MEDIAN_JOINT", "core5_joint_alignment", "CBM_JOINT", "BETA_PREDICTOR_JOINT"),
  predictor_swap(p, "CBM_STRUCTURAL7_MEDIAN", "structural7_siteset", "CBM_STRUCT7", "BETA_PREDICTOR_STRUCTURAL7"),
  predictor_swap(p, "CBM_CORE5_MEDIAN_HIGHPROV", "core5_high_provenance", "CBM_HIGHPROV", "BETA_PREDICTOR_HIGH_PROVENANCE"),
  predictor_swap(p, "CBM_CORE5_MEDIAN_NO_MODEL_SENSITIVE", "core5_no_model_sensitive", "CBM_NOMODELSENS", "BETA_PREDICTOR_NO_MODEL_SENSITIVE"),
  predictor_swap(p, "CBM_CONTACT_MINUS_BACKGROUND_MEDIAN", "contact_minus_background", "CBM_CMB", "BETA_PREDICTOR_CONTACT_MINUS_BACKGROUND"),
  predictor_swap(p, "CBM_CORE5_MEDIAN", "core5_median_primary", "CBM_NONCONTACT_COV", "BETA_NONCONTACT_IDENTITY_COVARIATE",
                 cov = "CBM_NONCONTACT_IDENTITY_MEDIAN")))
for (p in ext) sens <- c(sens, list(
  predictor_swap(p, "SCORE_SITE1_MEDIAN_HIGHPROV", "site1_high_provenance", "G_HIGHPROV", "GAMMA_PREDICTOR_HIGH_PROVENANCE"),
  predictor_swap(p, "SCORE_SITE1_MEDIAN_NO_MODEL_SENSITIVE", "site1_no_model_sensitive", "G_NOMODELSENS", "GAMMA_PREDICTOR_NO_MODEL_SENSITIVE"),
  predictor_swap(p, "SCORE_SITE1_MEDIAN_JOINT", "site1_joint_alignment", "G_JOINT", "GAMMA_PREDICTOR_JOINT"),
  predictor_swap(p, "SCORE_SITE3_MEDIAN", "site3_siteset", "G_SITE3", "GAMMA_PREDICTOR_SITE3"),
  predictor_swap(p, "SCORE_SITE4_MEDIAN", "site4_siteset", "G_SITE4", "GAMMA_PREDICTOR_SITE4"),
  predictor_swap(p, "CONTACT_MINUS_BACKGROUND_MEDIAN", "contact_minus_background", "G_CMB", "GAMMA_PREDICTOR_CONTACT_MINUS_BACKGROUND"),
  predictor_swap(p, "SCORE_SITE1_MEDIAN", "site1_median_primary", "G_NONCONTACT_COV", "GAMMA_NONCONTACT_IDENTITY_COVARIATE",
                 cov = "CBS_NONCONTACT_IDENTITY_MEDIAN")))
res <- run_catalogue(M, c(prim, ext, sens))
res <- apply_holm(res, F6, 6)
pid <- vapply(prim, `[[`, "", "MODEL_ID"); eid <- vapply(ext, `[[`, "", "MODEL_ID"); sid <- vapply(sens, `[[`, "", "MODEL_ID")

## the new data frame must give the same estimates for the two AMPKγ contact-score models as the reproduction
## step (same species, same model code)
rz <- read.delim(file.path(CTDX, "03_reproduction/gamma_alpha_reproduction.tsv"), comment.char = "#", colClasses = "character")
for (k in 1:2) {
  want <- rz$REPRODUCED[rz$ROW == c(1, 5)[k] & rz$FIELD == "COEF"]
  got <- sprintf("%.4f", res[[eid[k]]]$raw$coef)
  if (!identical(want, got)) stop("REFUSED: the cross-feature frame does not reproduce ", eid[k], ": ", got, " vs ", want)
}

N1 <- c("Cross-species phylogenetic associations only; no causal or ordering reading. phylolm lambda is primary.",
        "Whole-subunit loss is the frozen binary state (PRESENT 0, callable ABSENT 1, otherwise NA); retained-subunit",
        "sequence scores exist only for retaining species and are never set to zero for a lost subunit.",
        "HOLM_P uses a fixed family size equal to the declared family; an inestimable model keeps HOLM_P NA.")
write_tsv("04_cross_ampk/CTD_cross_AMPK_primary.tsv", bind_rows(res, pid), N1, S)
write_tsv("04_cross_ampk/gamma_alpha_extension_primary.tsv", bind_rows(res, eid), c(N1,
  "REPRODUCED_PRIMARY_EXTENSION: the frozen RBX1 gamma-alpha primary rows refit here only to carry declared sensitivities and LOCO."), S)
write_tsv("04_cross_ampk/CTD_cross_AMPK_sensitivity.tsv", sens_table(res, sid), c(N1,
  "One dimension changed per row. DELTA_BETA across a changed response or predictor scale is not a magnitude comparison."), S)
write_tsv("04_cross_ampk/sample_flow_cross_ampk.tsv", bind_flow(res, c(pid, eid, sid)), "N_FINAL equals the engine's fitted N.", S)
write_tsv("04_cross_ampk/encoding_checks_cross_ampk.tsv", enc_table(), "Counts feed the z5 gates.", S)
cat("z2 complete\n")
