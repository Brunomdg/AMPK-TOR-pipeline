# >>> deposit header
# Purpose: Fit the 18 phylogenetic logistic regressions of AMPK subunit absence (alpha, beta, gamma) on TORC1 or
#          TORC2 core disruption: six prespecified models with a Holm adjustment and two sensitivity sets, each
#          fitted by a penalised-likelihood multistart search.
# Methods: Phylogenetic regressions
# Inputs:  `TORC_core_states.tsv` (under `16_control_repair/refit/rbx2`; a core is intact when every component is
#          present, disrupted when at least one is absent where absence can be assessed); `rbx1_species_traits.tsv`
#          (proteome sizes); `species_chronogram_209.treefile`
# Outputs: `core_AMPK_phyloglm.tsv` and one .rds per fit, in `02_derived/core18_penalized_20260918`
# Tools:   R with ape, stats (glm, p.adjust) and phylolm 2.6.5 (phyloglm, logistic_MPLE, btol = 30); sources
#          `core_multistart_helper.R`; other versions: see envs/
# <<< deposit header
## deposit path defaults: begin
.TA210_SCRATCH <- local({v <- Sys.getenv("TA210_SCRATCH"); if (nzchar(v)) v else "/gpfs2/scratch/bdigenov"})
## deposit path defaults: end
## deposit repo root: begin
.REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(f[1]) else "."
  normalizePath(file.path(d, ".."))
})
## deposit repo root: end
## Phylogenetic models on the 209-tip chronogram, which was fixed before these fits.
## Binary state models use phylolm::phyloglm (logistic_MPLE, btol = 30), exactly as the regressions of the
## publication tables for candidate replacements and for lifestyle associations do.
## The continuous-trait function below (fit_cont, not called in this file) uses phylolm::phylolm with
## model = "lambda" (primary) and model = "BM" (sensitivity).
## The Holm adjustment is applied only across the six prespecified core-state models (family A below) and
## is a secondary column. The AMPK gamma contact and alpha C-terminal tables are read but not used here.
suppressMessages({library(ape); library(phylolm)})
SEED <- 488761235
set.seed(SEED)
## The output folder of these fits. It was the first line of core_multistart_helper.R until the two
## multistart helpers were merged; the path is unchanged.
D <- file.path(.TA210_SCRATCH, "ms210_phase2/TOR-AMPK_orthologs_10target_R1",
               "20_port_68taxa_analyses/02_derived", "core18_penalized_20260918")
source(file.path(.REPO_ROOT, "5_phylogenetic_analyses", "lib", "multistart_search.R"))
X <- Sys.getenv("RBX1_ROOT")
CONCORD <- paste0(.TA210_SCRATCH, "/concord/ms210")
TREE <- file.path(CONCORD, "P04b_chronogram", "species_chronogram_209.treefile")
COMPONENTS <- c("TOR", "RAPTOR", "RICTOR")
SUBUNITS <- c("AMPKa", "AMPKb", "AMPKg")
Z <- 1.959964
rd <- function(p) read.delim(p, comment.char = "#", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
tree <- read.tree(TREE)
sp <- rd(file.path(X, "05_component_loss_context", "rbx1_species_traits.tsv"))
gam <- rd(file.path(X, "03_ampkg_contacts", "AMPKg_contact_species.tsv"))
alp <- rd(file.path(X, "04_ampka_ctd", "AMPKa_CTD_species.tsv"))
num <- function(v) suppressWarnings(as.numeric(ifelse(v %in% c("NA", "-", ""), NA, v)))

banner <- function(what) c(
  sprintf("## RBX1 %s. phylolm %s on the %d-tip frozen chronogram %s. Seed %d.",
          what, as.character(packageVersion("phylolm")), Ntip(tree), basename(TREE), SEED),
  "## Binary responses: phyloglm logistic_MPLE btol=30, as Release B. Continuous: phylolm lambda (primary) and BM (sensitivity).",
  "## Every row carries its own N, coefficient, SE, 95% Wald interval and nominal p. NONPHYLO_* is the same model",
  "## without the phylogenetic correction, reported beside it so that what accounting for relatedness did is visible.",
  "## A p value is nominal unless the row names a declared testing family in HOLM_FAMILY.")

fit_binary <- function(dat, resp, pred, cov) {
  d <- dat[, c("SPECIES", resp, pred, cov), drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  d <- d[d$SPECIES %in% tree$tip.label, , drop = FALSE]
  if (nrow(d) < 10 || length(unique(d[[resp]])) < 2 || length(unique(d[[pred]])) < 2)
    return(list(status = "NOT_ESTIMABLE", n = nrow(d)))
  t2 <- keep.tip(tree, d$SPECIES); d <- d[match(t2$tip.label, d$SPECIES), ]; rownames(d) <- d$SPECIES
  f <- as.formula(paste(resp, "~", paste(c(pred, cov), collapse = " + ")))
  sep <- sum(d[[pred]] == 1 & d[[resp]] == 1) == 0 || sum(d[[pred]] == 1 & d[[resp]] == 0) == 0 || sum(d[[pred]] == 0 & d[[resp]] == 1) == 0 || sum(d[[pred]] == 0 & d[[resp]] == 0) == 0   ## quasi-separation: any empty cell of the 2 x 2 table (criterion chosen on 17 September 2026)
  pg <- fit_multiple(f,d,t2)
  gl <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
  if (is.null(pg)) return(list(status = "MODEL_FAILED_TO_FIT", n = nrow(d), sep = sep))
  s <- summary(pg)$coefficients
  list(status = "OK", n = nrow(d), coef = s[pred, 1], se = s[pred, 2], p = s[pred, 4],
       alpha = pg$alpha, sep = sep,
       np_coef = if (!is.null(gl)) coef(summary(gl))[pred, 1] else NA,
       np_p = if (!is.null(gl)) coef(summary(gl))[pred, 4] else NA)
}

fit_cont <- function(dat, resp, pred, cov, phylomodel) {
  d <- dat[, c("SPECIES", resp, pred, cov), drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  d <- d[d$SPECIES %in% tree$tip.label, , drop = FALSE]
  if (nrow(d) < 10 || length(unique(d[[pred]])) < 2) return(list(status = "NOT_ESTIMABLE", n = nrow(d)))
  t2 <- keep.tip(tree, d$SPECIES); d <- d[match(t2$tip.label, d$SPECIES), ]; rownames(d) <- d$SPECIES
  f <- as.formula(paste(resp, "~", paste(c(pred, cov), collapse = " + ")))
  pl <- tryCatch(phylolm(f, data = d, phy = t2, model = phylomodel), error = function(e) NULL)
  lm0 <- tryCatch(lm(f, data = d), error = function(e) NULL)
  if (is.null(pl)) return(list(status = "MODEL_FAILED_TO_FIT", n = nrow(d)))
  s <- summary(pl)$coefficients
  list(status = "OK", n = nrow(d), coef = s[pred, 1], se = s[pred, 2], p = s[pred, 4],
       alpha = if (!is.null(pl$optpar)) pl$optpar else NA,
       np_coef = if (!is.null(lm0)) coef(summary(lm0))[pred, 1] else NA,
       np_p = if (!is.null(lm0)) coef(summary(lm0))[pred, 4] else NA)
}

row_of <- function(fam, resp, variant, pred, cov, pm, r, family = "-") {
  c(fam, resp, variant, pred, if (nzchar(cov)) cov else "none", pm, r$n,
    if (!is.null(r$coef)) sprintf("%.4f", r$coef) else "NA",
    if (!is.null(r$se)) sprintf("%.4f", r$se) else "NA",
    if (!is.null(r$coef)) sprintf("%.4f", r$coef - Z * r$se) else "NA",
    if (!is.null(r$coef)) sprintf("%.4f", r$coef + Z * r$se) else "NA",
    if (!is.null(r$p)) sprintf("%.4g", r$p) else "NA",
    if (!is.null(r$alpha) && !is.na(r$alpha)) sprintf("%.4g", r$alpha) else "NA",
    if (!is.null(r$np_coef) && !is.na(r$np_coef)) sprintf("%.4f", r$np_coef) else "NA",
    if (!is.null(r$np_p) && !is.na(r$np_p)) sprintf("%.4g", r$np_p) else "NA",
    if (!is.null(r$sep) && isTRUE(r$sep)) "QUASI_SEPARATED" else "no", r$status, family)
}
HDR <- c("MODEL_FAMILY", "RESPONSE", "RESPONSE_VARIANT", "PREDICTOR", "COVARIATE", "PHYLO_MODEL", "N",
         "COEF", "SE", "CI95_LO", "CI95_HI", "NOMINAL_P", "PHYLO_PARAM", "NONPHYLO_COEF", "NONPHYLO_P",
         "QUASI_SEPARATION", "STATUS", "HOLM_FAMILY")
wr <- function(path, rows, what, extra = character(0)) {
  path <- file.path(D,basename(path))
  if(file.exists(path))stop("Refuse to overwrite staged table")
  con <- file(path, "w")
  writeLines(c(banner(what), extra), con)
  writeLines(paste(HDR, collapse = "\t"), con)
  for (r in rows) writeLines(paste(r, collapse = "\t"), con)
  close(con)
  cat(sprintf("wrote %s rows=%d\n", basename(path), length(rows)))
}


X1 <- X
X2 <- paste0(.TA210_SCRATCH, "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx2")
got <- "49cb4591031090d4e0c1e78967ff4418039a79f6d3fcf3d0aab6ccb6715c6faa"
Z <- 1.959964
num <- function(v) suppressWarnings(as.numeric(ifelse(v %in% c("NA", "-", ""), NA, v)))
HDR2 <- c("MODEL_FAMILY", "RESPONSE", "RESPONSE_VARIANT", "PREDICTOR", "COVARIATE", "PHYLO_MODEL", "N",
          "COEF", "SE", "CI95_LO", "CI95_HI", "NOMINAL_P", "PHYLO_PARAM", "NONPHYLO_COEF", "NONPHYLO_P",
          "QUASI_SEPARATION", "STATUS", "HOLM_FAMILY", "HOLM_P")
wr2 <- function(path, rows, what, extra = character(0)) {
  path <- file.path(D,basename(path))
  stopifnot(!file.exists(path))
  con <- file(path, "w")
  writeLines(c(sprintf("## RBX2 %s. Penalized multistart audit derived from rbx2_models.R (sha256 %s): phyloglm",
                       what, substr(got, 1, 16)),
               "## logistic_MPLE btol=30 for binary responses and phylolm lambda (BM as sensitivity) for continuous,",
               sprintf("## on the frozen %d-tip chronogram, seed %d, phylolm %s.", Ntip(tree), SEED,
                       as.character(packageVersion("phylolm"))),
               "## Every row carries its own N, coefficient, SE, 95% Wald interval and nominal p. HOLM_P is filled",
               "## only for rows inside a declared primary family, named in HOLM_FAMILY; every other p is nominal.",
               extra), con)
  writeLines(paste(HDR2, collapse = "\t"), con)
  for (r in rows) writeLines(paste(r, collapse = "\t"), con)
  close(con)
  cat(sprintf("wrote %s rows=%d\n", basename(path), length(rows)))
}
mk <- function(fam, resp, variant, pred, cov, pm, r, family = "-", holm = "-") {
  c(fam, resp, variant, pred, if (nzchar(cov)) cov else "none", pm, r$n,
    if (!is.null(r$coef)) sprintf("%.4f", r$coef) else "NA",
    if (!is.null(r$se)) sprintf("%.4f", r$se) else "NA",
    if (!is.null(r$coef)) sprintf("%.4f", r$coef - Z * r$se) else "NA",
    if (!is.null(r$coef)) sprintf("%.4f", r$coef + Z * r$se) else "NA",
    if (!is.null(r$p)) sprintf("%.4g", r$p) else "NA",
    if (!is.null(r$alpha) && !is.na(r$alpha)) sprintf("%.4g", r$alpha) else "NA",
    if (!is.null(r$np_coef) && !is.na(r$np_coef)) sprintf("%.4f", r$np_coef) else "NA",
    if (!is.null(r$np_p) && !is.na(r$np_p)) sprintf("%.4g", r$np_p) else "NA",
    if (!is.null(r$sep) && isTRUE(r$sep)) "QUASI_SEPARATED" else "no", r$status, family, holm)
}

MODE <- (function(a) if (length(a)) a[1] else "core")(commandArgs(trailingOnly = TRUE))
cat(sprintf("mode=%s\n", MODE))

## ---- data ----
core <- rd(file.path(X2, "03_torc_core/TORC_core_states.tsv"))
tr <- rd(file.path(X1, "05_component_loss_context/rbx1_species_traits.tsv"))
gam <- rd(file.path(X1, "03_ampkg_contacts/AMPKg_contact_species.tsv"))
alp <- rd(file.path(X1, "04_ampka_ctd/AMPKa_CTD_species.tsv"))
core$TORC1_disrupted <- ifelse(core$TORC1_CORE_STATE == "DISRUPTED", 1,
                        ifelse(core$TORC1_CORE_STATE == "INTACT", 0, NA))
core$TORC2_disrupted <- ifelse(core$TORC2_CORE_STATE == "DISRUPTED", 1,
                        ifelse(core$TORC2_CORE_STATE == "INTACT", 0, NA))
for (s in SUBUNITS) core[[paste0(s, "_absent")]] <- num(core[[paste0(s, "_absent")]])
core$log10_n_proteins <- num(tr$log10_n_proteins[match(core$SPECIES, tr$SPECIES)])
PREDS <- c("TORC1_disrupted", "TORC2_disrupted")

if (MODE == "core") {
## ---- family A: the six prespecified models of AMPK subunit absence on core disruption ----
rows <- list(); praw <- c()
for (p in PREDS) for (s in SUBUNITS) {
  r <- fit_binary(core, paste0(s, "_absent"), p, character(0))
  rows[[length(rows) + 1]] <- mk("CORE_STATE", paste0(s, "_absent"), "all_callable_species", p, "",
                                 "phyloglm_logistic_MPLE", r, "six_prespecified_core_state_models")
  praw <- c(praw, if (!is.null(r$p)) r$p else NA)
}
.holm_n <- 6L  ## Holm denominator fixed at the declared family size of six, so a model without a p value does not shrink it
if (length(praw) != .holm_n) stop("REFUSED: family six_prespecified_core_state_models has ", length(praw),
    " members, declared ", .holm_n)
holm <- p.adjust(praw, method = "holm", n = .holm_n)
for (i in seq_along(rows)) rows[[i]][length(HDR2)] <- sprintf("%.4g", holm[i])
for (p in PREDS) for (s in SUBUNITS) {
  r <- fit_binary(core, paste0(s, "_absent"), p, "log10_n_proteins")
  rows[[length(rows) + 1]] <- mk("CORE_STATE", paste0(s, "_absent"), "all_callable_species", p,
                                 "log10_n_proteins", "phyloglm_logistic_MPLE", r)
  hp <- core[core[[paste0(s, "_R66_ONLY")]] != "yes", ]
  r2 <- fit_binary(hp, paste0(s, "_absent"), p, character(0))
  rows[[length(rows) + 1]] <- mk("CORE_STATE", paste0(s, "_absent"), "high_provenance_only", p, "",
                                 "phyloglm_logistic_MPLE", r2)
}
wr2(file.path(X2, "03_torc_core/core_AMPK_phyloglm.tsv"), rows, "canonical core integrity versus AMPK subunit loss",
    "## The six no-covariate models are the prespecified family; every other row is a declared sensitivity.")


}
