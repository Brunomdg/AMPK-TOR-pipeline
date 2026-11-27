#!/usr/bin/env Rscript
# >>> deposit header
# Purpose: The model engine of the sequence-phenotype regressions: the frozen chronogram, the seed, the binary
#          and continuous fitting functions (`fit_binary`, `fit_cont`), the output row format (`row_of`, `HDR`,
#          `wr`) and the banner. Its preamble is what the AMPKa tail models load; its own four analyses follow
#          the preamble and are superseded (see the note below).
# Methods: Phylogenetic regressions
# Inputs:  the 209-tip dated species tree `species_chronogram_209.treefile`; `rbx1_species_traits.tsv`,
#          `AMPKg_contact_species.tsv` and `AMPKa_CTD_species.tsv` under the root given in `RBX1_ROOT`
# Outputs: when run whole: `component_loss_AMPK_phyloglm.tsv`, `AMPKg_contact_models.tsv`,
#          `AMPKa_CTD_models.tsv` and `gamma_alpha_covariation.tsv`. When loaded by `ctdx_lib.R`, nothing: only
#          the preamble above the `## ---------------- Analysis 1` marker is evaluated
# Tools:   R with ape and phylolm 2.6.5 (phyloglm logistic_MPLE btol = 30; phylolm lambda primary, BM
#          sensitivity); versions: see envs/
# <<< deposit header
# Ran as: rbx1_models.R
# Two things a reader needs to know about this file, both recorded in `docs/deposit_gaps.md`.
#
# 1. This is not the revision that produced the published tail models. `ctdx_lib.R` pins the engine at sha256
#    526dc52ae74a87bc46a374ec55c4b83325d6ea8a9706b85fb0be49055f3abe71; this file is ef3350ad8eb85b33bbead953
#    add1f94c18f94f197aa1b8a27539b7efc384bec7, a later revision. The pinned revision was not preserved. The
#    preamble that the tail models load is unchanged in substance: the seed, the chronogram, `fit_cont` and
#    `row_of` all still satisfy the checks `ctdx_lib.R` makes of them.
# 2. Analysis 1 below is superseded. The nine prespecified state models were refitted by the audited penalised
#    multistart pipeline, and the guard on its first line stops this file from overwriting them. Analyses 2 to 4
#    are the engine's own AMPKg contact, AMPKa tail and covariation fits.
# deposit path defaults: begin
.TA210_SCRATCH <- local({v <- Sys.getenv("TA210_SCRATCH"); if (nzchar(v)) v else "/gpfs2/scratch/bdigenov"})
# deposit path defaults: end
## GOVERNED BY: RBX1_EXECUTION_CURRENT.md job x4. Phylogenetic models on the frozen 209-tip chronogram.
## Binary state models use phylolm::phyloglm logistic_MPLE btol=30, exactly as Release B Q0 and Q6 do.
## Continuous sequence models use phylolm::phylolm with model="lambda" (primary) and model="BM" (sensitivity).
## Holm adjustment is applied only across the nine prespecified state models and is a secondary column.
suppressMessages({library(ape); library(phylolm)})
SEED <- 488761235
set.seed(SEED)
X <- Sys.getenv("RBX1_ROOT")
CONCORD <- file.path(.TA210_SCRATCH, "concord/ms210")
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
  sep <- sum(d[[pred]] == 1 & d[[resp]] == 1) == 0 || sum(d[[pred]] == 1 & d[[resp]] == 0) == 0 || sum(d[[pred]] == 0 & d[[resp]] == 1) == 0 || sum(d[[pred]] == 0 & d[[resp]] == 0) == 0   ## any empty cell: Bruno 2026-09-17, QUASI_SEPARATION_LITERATURE_BASIS.md
  pg <- tryCatch(phyloglm(f, data = d, phy = t2, method = "logistic_MPLE", btol = 30), error = function(e) NULL)
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
  con <- file(path, "w")
  writeLines(c(banner(what), extra), con)
  writeLines(paste(HDR, collapse = "\t"), con)
  for (r in rows) writeLines(paste(r, collapse = "\t"), con)
  close(con)
  cat(sprintf("wrote %s rows=%d\n", basename(path), length(rows)))
}

## ---------------- Analysis 1: the nine prespecified state models ----------------
if(file.exists(paste0(.TA210_SCRATCH, "/ms210_phase2/TOR-AMPK_orthologs_10target_R1/20_port_68taxa_analyses/02_derived/component27_penalized_20260918/INSTALLATION.json")))stop("REFUSED: component fits are governed by the audited penalized multistart pipeline; use component27_penalized.R and validation/installation workflow before replacing them.")

for (tg in c(COMPONENTS, SUBUNITS)) sp[[paste0(tg, "_absent")]] <- num(sp[[paste0(tg, "_absent")]])
sp$log10_n_proteins <- num(sp$log10_n_proteins)
rows <- list(); praw <- c(); pidx <- c()
for (comp in COMPONENTS) for (sub in SUBUNITS) {
  r <- fit_binary(sp, paste0(sub, "_absent"), paste0(comp, "_absent"), character(0))
  rows[[length(rows) + 1]] <- row_of("STATE", paste0(sub, "_absent"), "all_callable_species",
                                     paste0(comp, "_absent"), "", "phyloglm_logistic_MPLE", r,
                                     "nine_prespecified_state_models")
  praw <- c(praw, if (!is.null(r$p)) r$p else NA); pidx <- c(pidx, length(rows))
}
.holm_n <- 9L  ## declared size of nine_prespecified_state_models; see HOLM_DENOMINATOR_FIX.tsv
if (length(praw) != .holm_n) stop("REFUSED: family nine_prespecified_state_models has ", length(praw),
    " members, declared ", .holm_n)
holm <- p.adjust(praw, method = "holm", n = .holm_n)
for (i in seq_along(pidx)) rows[[pidx[i]]] <- c(rows[[pidx[i]]], sprintf("%.4g", holm[i]))
HDR <- c(HDR, "HOLM_P_NINE_STATE_MODELS")
for (i in seq_along(rows)) if (length(rows[[i]]) < length(HDR)) rows[[i]] <- c(rows[[i]], "-")
## declared sensitivities: proteome size, and dropping species whose only ortholog is R66-resolved
for (comp in COMPONENTS) for (sub in SUBUNITS) {
  r <- fit_binary(sp, paste0(sub, "_absent"), paste0(comp, "_absent"), "log10_n_proteins")
  rows[[length(rows) + 1]] <- c(row_of("STATE", paste0(sub, "_absent"), "all_callable_species",
                                       paste0(comp, "_absent"), "log10_n_proteins", "phyloglm_logistic_MPLE", r), "-")
  hp <- sp[!(sp[[paste0(sub, "_R66_ONLY")]] == "yes" | sp[[paste0(comp, "_R66_ONLY")]] == "yes"), ]
  r2 <- fit_binary(hp, paste0(sub, "_absent"), paste0(comp, "_absent"), character(0))
  rows[[length(rows) + 1]] <- c(row_of("STATE", paste0(sub, "_absent"), "high_provenance_only",
                                       paste0(comp, "_absent"), "", "phyloglm_logistic_MPLE", r2), "-")
}
wr(file.path(X, "05_component_loss_context", "component_loss_AMPK_phyloglm.tsv"), rows,
   "component loss versus AMPK subunit loss",
   "## The nine no-covariate models are the prespecified family carrying HOLM_P_NINE_STATE_MODELS; every other row is a declared sensitivity and its p is nominal.")

## ---------------- Analysis 2: AMPKg AMP-contact scores ----------------
g <- merge(gam, sp, by = "SPECIES")
for (cn in grep("^SCORE_|^CBS_|^CONTACT_MINUS", names(g), value = TRUE)) g[[cn]] <- num(g[[cn]])
rows <- list()
variants <- list(c("SCORE_SITE1_MEDIAN", "site1_median_primary"), c("SCORE_SITE1_MIN", "site1_min_duplicate_sensitivity"),
                 c("SCORE_SITE1_MAX", "site1_max_duplicate_sensitivity"), c("SCORE_SITE1_MEDIAN_HIGHPROV", "site1_high_provenance"),
                 c("SCORE_SITE1_MEDIAN_NO_MODEL_SENSITIVE", "site1_no_model_sensitive"),
                 c("SCORE_SITE3_MEDIAN", "site3_median_sensitivity_siteset"), c("SCORE_SITE4_MEDIAN", "site4_median_sensitivity_siteset"),
                 c("SCORE_SITE1_MEDIAN_JOINT", "site1_median_joint_alignment_sensitivity"))
for (v in variants) for (comp in COMPONENTS) for (pm in c("lambda", "BM")) {
  r <- fit_cont(g, v[1], paste0(comp, "_absent"), character(0), pm)
  rows[[length(rows) + 1]] <- c(row_of("GAMMA_CONTACT", v[1], v[2], paste0(comp, "_absent"), "", paste0("phylolm_", pm), r), "-")
}
for (comp in COMPONENTS) for (cov in c("CBS_NONCONTACT_IDENTITY_MEDIAN", "log10_n_proteins")) {
  r <- fit_cont(g, "SCORE_SITE1_MEDIAN", paste0(comp, "_absent"), cov, "lambda")
  rows[[length(rows) + 1]] <- c(row_of("GAMMA_CONTACT", "SCORE_SITE1_MEDIAN", "site1_median_primary",
                                       paste0(comp, "_absent"), cov, "phylolm_lambda", r), "-")
}
wr(file.path(X, "03_ampkg_contacts", "AMPKg_contact_models.tsv"), rows,
   "AMP-contact score versus component loss among species retaining AMPKg",
   "## site1 is the primary frozen site set; site3 and site4 are predeclared sensitivity site sets and are never pooled with it.")

## ---------------- Analysis 3: AMPKa CTD length and disorder ----------------
a <- merge(alp, sp, by = "SPECIES")
for (cn in grep("^CTD_", names(a), value = TRUE)) a[[cn]] <- num(a[[cn]])
a$log1p_CTD_LENGTH_MEDIAN <- log1p(a$CTD_LENGTH_MEDIAN)
rows <- list()
resp <- list(c("CTD_LENGTH_MEDIAN", "median_primary"), c("CTD_MEAN_DISORDER_MEDIAN", "median_primary"),
             c("CTD_LENGTH_MIN", "min_duplicate_sensitivity"), c("CTD_LENGTH_MAX", "max_duplicate_sensitivity"),
             c("CTD_MEAN_DISORDER_MIN", "min_duplicate_sensitivity"), c("CTD_MEAN_DISORDER_MAX", "max_duplicate_sensitivity"),
             c("CTD_LENGTH_MEDIAN_HIGHPROV", "high_provenance"), c("CTD_MEAN_DISORDER_MEDIAN_HIGHPROV", "high_provenance"),
             c("CTD_LENGTH_MEDIAN_NO_MODEL_SENSITIVE", "no_model_sensitive"),
             c("CTD_MEAN_DISORDER_MEDIAN_NO_MODEL_SENSITIVE", "no_model_sensitive"),
             c("log1p_CTD_LENGTH_MEDIAN", "log1p_sensitivity"), c("CTD_LENGTH_DOMAIN_END_MEDIAN", "domain_end_definition_sensitivity"))
for (v in resp) for (comp in COMPONENTS) for (pm in c("lambda", "BM")) {
  if (pm == "BM" && v[2] != "median_primary") next
  r <- fit_cont(a, v[1], paste0(comp, "_absent"), character(0), pm)
  rows[[length(rows) + 1]] <- c(row_of("ALPHA_CTD", v[1], v[2], paste0(comp, "_absent"), "", paste0("phylolm_", pm), r), "-")
}
for (v in c("CTD_LENGTH_MEDIAN", "CTD_MEAN_DISORDER_MEDIAN")) for (comp in COMPONENTS) {
  r <- fit_cont(a, v, paste0(comp, "_absent"), "log10_n_proteins", "lambda")
  rows[[length(rows) + 1]] <- c(row_of("ALPHA_CTD", v, "median_primary", paste0(comp, "_absent"),
                                       "log10_n_proteins", "phylolm_lambda", r), "-")
}
wr(file.path(X, "04_ampka_ctd", "AMPKa_CTD_models.tsv"), rows,
   "alpha CTD length and disorder versus component loss among species retaining AMPKa",
   "## Only genome-derived records with a callable residue-290 anchor enter the species medians; raw length is primary and log1p is a declared sensitivity.")

## ---------------- Analysis 4: gamma-alpha covariation ----------------
ga <- merge(merge(gam, alp, by = "SPECIES"), sp, by = "SPECIES")
for (cn in grep("^SCORE_|^CTD_|^CBS_", names(ga), value = TRUE)) ga[[cn]] <- num(ga[[cn]])
rows <- list()
for (v in c("CTD_LENGTH_MEDIAN", "CTD_MEAN_DISORDER_MEDIAN")) {
  r <- fit_cont(ga, v, "SCORE_SITE1_MEDIAN", character(0), "lambda")
  rows[[length(rows) + 1]] <- c(row_of("GAMMA_ALPHA", v, "median_primary", "SCORE_SITE1_MEDIAN", "", "phylolm_lambda", r), "-")
  for (comp in COMPONENTS) {
    r2 <- fit_cont(ga, v, "SCORE_SITE1_MEDIAN", paste0(comp, "_absent"), "lambda")
    rows[[length(rows) + 1]] <- c(row_of("GAMMA_ALPHA", v, "median_primary", "SCORE_SITE1_MEDIAN",
                                         paste0(comp, "_absent"), "phylolm_lambda", r2), "-")
  }
}
wr(file.path(X, "07_integrated_models", "gamma_alpha_covariation.tsv"), rows,
   "alpha CTD versus gamma AMP-contact score among species retaining both",
   "## Covariation across species, not a coupling test and not a temporal statement. Interaction terms are not fitted: the state combinations do not carry independent phylogenetic replication.")
cat("x4 models complete\n")
