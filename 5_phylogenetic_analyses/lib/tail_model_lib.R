# >>> deposit header
# Purpose: The shared library of the AMPKa tail models: builds the model frames, encodes every state once,
#          applies the Holm adjustment within each declared family, runs the model catalogue and writes the
#          result, sensitivity and sample-flow tables. Loads the preamble of the model engine.
# Methods: Phylogenetic analyses
# Inputs:  `lib/model_engine.R` (its preamble only, up to the Analysis 1 marker); the dated chronogram; the
#          ten-target species-state table; the AMPKa tail, AMPKg contact, AMPKb module, upstream-kinase and
#          architecture species tables
# Outputs: none of its own; the two calling scripts write through `write_tsv()`
# Tools:   R with ape and phylolm 2.6.5; versions: see environment/
# <<< deposit header
# Ran as: ctdx_lib.R
# It carried no header in the first deposit because other scripts verified its checksum. Those checks are
# gone, so it is documented like every other file.
# Ran as: ctdx_lib.R
## deposit repo root: begin
.REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(f[1]) else "."
  normalizePath(file.path(d, "..", ".."))
})
## deposit repo root: end
## TEN-TARGET REFIT of ctdx_lib.R (Bruno 2026-09-16). Inputs are the sister's ten-target T2 and the rebuilt
## RBX1, RBX2 and RBX3 layers; the engine and the chronogram are the frozen ones and their pins are unchanged.
## GOVERNED BY: CTD_CROSSFEATURE_EXECUTION_CURRENT.md, shared by z1 to z4. Loads the frozen RBX1 engine preamble
## (everything before "## ---------------- Analysis 1") after a SHA256 check, exactly as rbx2c_loco.R does, so
## fit_cont, row_of, tree, sp, gam and alp are the frozen objects. None of the old analyses is executed and
## nothing is written outside the new release root. Frames encode states once; no score is recomputed.
suppressMessages({library(ape); library(phylolm)})
CTDX <- normalizePath(Sys.getenv("CTDX_ROOT"), mustWork = TRUE)
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("REFUSED: not in Slurm")
if (normalizePath(Sys.getenv("RBX1_ROOT")) != "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1")
  stop("REFUSED: RBX1_ROOT is not the rebuilt RBX1 refit root")
sha256 <- function(p) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE))
SHA <- list(
  ENGINE = "526dc52ae74a87bc46a374ec55c4b83325d6ea8a9706b85fb0be49055f3abe71",
  TREE = "300dcecaa8ef78d7e0583560a9428769472c1204854ba0cbe809c494bd906060",
  T2 = "e4e4946abac3d6cb41ed3c3477e8bd6d9d70f03820ef1e3930f1b1e9023e65d4",
  CTD = "06bf0a3e2070d92bb49c588b102a5c17149cc02a87e23d05f262ece787257241",
  GAM = "3bfc67076c79e04d59099534d7ffea03d5a0816004e6fa618ed758bff576a40d",
  TRAITS = "4ed9781d6ca3d4b71b946a89b4b34cccb390ce7b57590cc60d59a206842b5bb3",
  GA = "e555088b90254c71f225905c93878056c6e713b686e0377b2867f175ecc1d3ca",
  CBM = "5a2e606cc2687b6f5437c942bf9d37a11dc188e918fb03f9930da83bd336c192",
  UP2 = "70cc7fda7b81555728d9d774886fcc337ba4c948fd5cdfc0160223db96a2fe38",
  ARCH3 = "7e2f7bfc2b78fc37979602a02bd6f7ddcd02ea1ac053863e8667e9e8866092f6")
B <- "/gpfs2/scratch/bdigenov/ms210_phase2"
PATHS <- list(
  ENGINE = file.path(.REPO_ROOT, "5_phylogenetic_analyses", "lib", "model_engine.R"),
  TREE = "/gpfs2/scratch/bdigenov/concord/ms210/P04b_chronogram/species_chronogram_209.treefile",
  T2 = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/07_findings/publication_tables/T2_presence_long.tsv",
  CTD = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1/04_ampka_ctd/AMPKa_CTD_species.tsv",
  GAM = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1/03_ampkg_contacts/AMPKg_contact_species.tsv",
  TRAITS = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1/05_component_loss_context/rbx1_species_traits.tsv",
  GA = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx1/07_integrated_models/gamma_alpha_covariation.tsv",
  CBM = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx2/05_ampkb_cbm/AMPKb_CBM_species.tsv",
  UP2 = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx2/04_ampka_activation/upstream_kinase_states.tsv",
  ARCH3 = "/gpfs2/scratch/bdigenov/ms210_phase2/TOR-AMPK_orthologs_10target_R1/16_control_repair/refit/rbx3/05_upstream_architecture/CAMKK_LKB1_architecture.tsv")
for (k in names(SHA)) if (sha256(PATHS[[k]]) != SHA[[k]]) stop("REFUSED: sha256 mismatch for ", k, " ", PATHS[[k]])

## ---- the frozen engine preamble ----
src <- readLines(PATHS$ENGINE)
cut <- grep("^## ---------------- Analysis 1", src)[1]
if (is.na(cut)) stop("REFUSED: engine preamble marker not found")
eval(parse(text = paste(src[seq_len(cut - 1)], collapse = "\n")), envir = globalenv())
for (fn in c("fit_cont", "row_of")) if (!exists(fn)) stop("REFUSED: the engine lacks ", fn)
if (!grepl('phylolm\\(f, data = d, phy = t2, model = phylomodel\\)', paste(deparse(fit_cont), collapse = " ")))
  stop("REFUSED: fit_cont is not the frozen phylolm implementation")
if (SEED != 488761235) stop("REFUSED: engine seed changed")
if (normalizePath(TREE) != normalizePath(PATHS$TREE) || Ntip(tree) != 209) stop("REFUSED: engine tree is not the frozen chronogram")
if (as.character(packageVersion("phylolm")) != "2.6.5") stop("REFUSED: phylolm is not 2.6.5")

## ---- write guard ----
JOBTAG <- Sys.getenv("SLURM_JOB_ID")
guard <- function(p, script = "R") {
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  rp <- file.path(normalizePath(dirname(p), mustWork = TRUE), basename(p))
  if (file.exists(p)) rp <- normalizePath(p)
  inside <- startsWith(rp, paste0(CTDX, "/"))
  lg <- file.path(CTDX, "09_logs", "write_guard.tsv")
  if (!file.exists(lg)) cat("TIMESTAMP\tSLURM_JOB_ID\tSCRIPT\tPATH\tREALPATH\tINSIDE_ROOT\n", file = lg)
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), JOBTAG, script, p, rp,
              if (inside) "yes" else "no"), file = lg, append = TRUE)
  if (!inside) stop("REFUSED: write to ", p, " resolves outside ", CTDX)
  rp
}
write_tsv <- function(rel, df, notes = character(0), script = "R") {
  p <- guard(file.path(CTDX, rel), script)
  con <- file(p, "w")
  writeLines(c(sprintf("## CTD-CROSSFEATURE-R1-20260915 | %s | %s job %s | phylolm %s | seed %d | tree %s",
                       rel, script, JOBTAG, as.character(packageVersion("phylolm")), SEED, basename(PATHS$TREE)),
               paste("##", notes)), con)
  close(con)
  suppressWarnings(write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE, na = "NA"))
  cat(sprintf("wrote %s rows=%d\n", rel, nrow(df)))
}
chr0 <- function(v) { v <- as.character(v); v[is.na(v)] <- "NA"; v }
fmt <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "NA" else sprintf("%.6g", x)
uniq_or_stop <- function(df, what) if (anyDuplicated(df$SPECIES)) stop("REFUSED: duplicate SPECIES in ", what)

## ---- frames (encoding happens here and only here) ----
ENC <- list()   # gate counters
build_frames <- function() {
  uniq_or_stop(alp, "CTD"); uniq_or_stop(gam, "gamma"); uniq_or_stop(sp, "traits")
  cbm <- rd(PATHS$CBM); up <- rd(PATHS$UP2); ar <- rd(PATHS$ARCH3)
  uniq_or_stop(cbm, "CBM"); uniq_or_stop(up, "upstream"); uniq_or_stop(ar, "architecture")

  a <- alp
  for (cn in grep("^CTD_", names(a), value = TRUE)) a[[cn]] <- num(a[[cn]])
  a$log1p_CTD_LENGTH_MEDIAN <- log1p(a$CTD_LENGTH_MEDIAN)
  s <- sp[, c("SPECIES", "CLADE", "ON_CHRONOGRAM", "AMPKa_STATE", "AMPKg_STATE", "AMPKg_ABSENCE_CALLABLE", "AMPKg_absent",
              "AMPKb_STATE", "AMPKb_ABSENCE_CALLABLE", "AMPKb_absent", "log10_n_proteins")]
  s$AMPKg_absent <- num(s$AMPKg_absent); s$AMPKb_absent <- num(s$AMPKb_absent); s$log10_n_proteins <- num(s$log10_n_proteins)
  ## the frozen coding is used as-is; this only counts cells that would contradict it
  bad_bin <- function(state, callable, code) sum(
    (!is.na(code) & code == 0 & state != "PRESENT") |
    (!is.na(code) & code == 1 & !(state == "ABSENT" & callable == "yes")) |
    (is.na(code) & (state == "PRESENT" | (state == "ABSENT" & callable == "yes"))))
  ENC$N_TRAITS_BINARY_CONTRADICTIONS <<- bad_bin(s$AMPKg_STATE, s$AMPKg_ABSENCE_CALLABLE, s$AMPKg_absent) +
                                         bad_bin(s$AMPKb_STATE, s$AMPKb_ABSENCE_CALLABLE, s$AMPKb_absent)
  if (!all(a$SPECIES %in% s$SPECIES)) stop("REFUSED: CTD species missing from traits")
  if (any(s$AMPKa_STATE[match(a$SPECIES, s$SPECIES)] != "PRESENT")) stop("REFUSED: CTD row for a species without AMPKa PRESENT")
  M <- merge(a, s, by = "SPECIES", all.x = TRUE, sort = TRUE)

  gcols <- c("SCORE_SITE1_MEDIAN", "SCORE_SITE1_MEDIAN_HIGHPROV", "SCORE_SITE1_MEDIAN_NO_MODEL_SENSITIVE",
             "SCORE_SITE1_MEDIAN_JOINT", "SCORE_SITE3_MEDIAN", "SCORE_SITE4_MEDIAN", "CBS_NONCONTACT_IDENTITY_MEDIAN",
             "CONTACT_MINUS_BACKGROUND_MEDIAN")
  g <- gam[, c("SPECIES", gcols)]; for (cn in gcols) g[[cn]] <- num(g[[cn]])
  M <- merge(M, g, by = "SPECIES", all.x = TRUE, sort = TRUE)
  bcols <- c("CBM_CORE5_MEDIAN", "CBM_CORE5_MEDIAN_JOINT", "CBM_STRUCTURAL7_MEDIAN", "CBM_CORE5_MEDIAN_HIGHPROV",
             "CBM_CORE5_MEDIAN_NO_MODEL_SENSITIVE", "CBM_NONCONTACT_IDENTITY_MEDIAN", "CBM_CONTACT_MINUS_BACKGROUND_MEDIAN")
  cb <- cbm[, c("SPECIES", bcols)]; for (cn in bcols) cb[[cn]] <- num(cb[[cn]])
  M <- merge(M, cb, by = "SPECIES", all.x = TRUE, sort = TRUE)
  ## a whole-subunit loss must never carry a sequence score, and no score is ever set to zero
  gscore <- rowSums(!is.na(M[, gcols])) > 0; bscore <- rowSums(!is.na(M[, bcols])) > 0
  ENC$N_GAMMA_LOSS_ENCODED_AS_CONTACT_SCORE_ZERO <<- sum(gscore & !is.na(M$AMPKg_absent) & M$AMPKg_absent == 1)
  ENC$N_BETA_LOSS_ENCODED_AS_CBM_ZERO <<- sum(bscore & !is.na(M$AMPKb_absent) & M$AMPKb_absent == 1)
  ENC$N_GAMMA_SCORE_WITHOUT_RETAINED_GAMMA <<- sum(gscore & (is.na(M$AMPKg_absent) | M$AMPKg_absent != 0))
  ENC$N_BETA_SCORE_WITHOUT_RETAINED_BETA <<- sum(bscore & (is.na(M$AMPKb_absent) | M$AMPKb_absent != 0))
  if (ENC$N_GAMMA_LOSS_ENCODED_AS_CONTACT_SCORE_ZERO + ENC$N_BETA_LOSS_ENCODED_AS_CBM_ZERO > 0)
    stop("REFUSED: a whole-subunit loss carries a sequence score")

  ## RBX2 upstream, the governing source
  u <- up[, c("SPECIES", "CAMKK_STATE", "CAMKK_CALLABLE", "LKB1_STATE", "LKB1_CALLABLE", "UPSTREAM_KINASE_STATE",
              "NONE_CAPTURED_BINARY")]
  for (cn in names(u)[-1]) u[[cn]] <- chr0(u[[cn]])
  u$NONE_CAPTURED_ASGIVEN <- ifelse(u$NONE_CAPTURED_BINARY == "NONE_CAPTURED", 1,
                                    ifelse(u$NONE_CAPTURED_BINARY == "ANY_CAPTURED", 0, NA))
  u$NONE_CAPTURED <- ifelse(u$UPSTREAM_KINASE_STATE == "UNKNOWN", NA, u$NONE_CAPTURED_ASGIVEN)
  binabs <- function(state, callable) ifelse(state == "PRESENT", 0, ifelse(state == "ABSENT" & callable == "yes", 1, NA))
  u$CAMKK_absent <- binabs(u$CAMKK_STATE, u$CAMKK_CALLABLE)
  u$LKB1_absent <- binabs(u$LKB1_STATE, u$LKB1_CALLABLE)
  M <- merge(M, u, by = "SPECIES", all.x = TRUE, sort = TRUE)
  if (any(is.na(M$UPSTREAM_KINASE_STATE))) stop("REFUSED: CTD species missing from RBX2 upstream states")

  ## RBX3 sensitivities, AMPKa-present rows only
  r3 <- ar[chr0(ar$AMPKa_STATE) == "PRESENT", ]
  for (cn in c("LKB1_USABLE_SENSITIVITY_A", "LKB1_USABLE_SENSITIVITY_B", "UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED",
               "LKB1_STATE_FROZEN")) r3[[cn]] <- chr0(r3[[cn]])
  usA <- chr0(ar$LKB1_USABLE_SENSITIVITY_A); usB <- chr0(ar$LKB1_USABLE_SENSITIVITY_B)
  ## Re-pinned twice. The frozen design was 204/198. The ten-target rebuild first gave 209/209, against the
  ## architecture table of 2026-09-16 22:51. It now gives 196/189, against the table r4 rebuilt on
  ## 2026-09-17 from the masked-primary LKB1 tree that L57 ruled for. LKB1 is present in 44 of the 210
  ## species; of those, 10 are declared positives, 13 are strongly supported, 7 compatible and 14
  ## unresolved, so tier A holds 30 of the 44 and tier B 23, and every species outside a tier is one where
  ## LKB1 is present and no retained copy reaches it. The 30 reaching tier A are the same 30 of 42 that L57
  ## recorded for the masked tree.
  if (sum(usA == "yes") != 196 || sum(usB == "yes") != 189) stop("REFUSED: LKB1 tier masks differ from the rebuilt design")
  if (any(usB == "yes" & usA != "yes")) stop("REFUSED: LKB1 tier B is not a subset of tier A")
  r3$NONE_CAPTURED_SEQ_SUPPORTED <- ifelse(r3$UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED == "NEITHER_CAPTURED", 1,
    ifelse(r3$UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED %in% c("BOTH_CAPTURED", "CAMKK_CATEGORY_ONLY", "LKB1_ONLY"), 0, NA))
  M <- merge(M, r3[, c("SPECIES", "LKB1_STATE_FROZEN", "LKB1_USABLE_SENSITIVITY_A", "LKB1_USABLE_SENSITIVITY_B",
                       "UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED", "NONE_CAPTURED_SEQ_SUPPORTED")],
             by = "SPECIES", all.x = TRUE, sort = TRUE)
  ## the literal RBX3 tier line
  M$LKB1_absent_A <- ifelse(chr0(M$LKB1_USABLE_SENSITIVITY_A) == "yes", M$LKB1_absent, NA)
  M$LKB1_absent_B <- ifelse(chr0(M$LKB1_USABLE_SENSITIVITY_B) == "yes", M$LKB1_absent, NA)
  if (any(!is.na(M$LKB1_absent) & M$LKB1_absent == 1 & (is.na(M$LKB1_absent_A) | is.na(M$LKB1_absent_B))))
    stop("REFUSED: a callable LKB1 absence was removed by a sequence tier")

  ## gate counters for the encodings actually used
  unk <- M$UPSTREAM_KINASE_STATE == "UNKNOWN"
  ENC$N_UNKNOWN_WITH_NONE_CAPTURED_CODED <<- sum(unk & !is.na(M$NONE_CAPTURED))
  ENC$N_UNKNOWN_ADMITTED_BY_ASGIVEN_SENSITIVITY <<- sum(unk & !is.na(M$NONE_CAPTURED_ASGIVEN))
  ENC$N_UNKNOWN_SEQ_SUPPORTED_CODED <<- sum(M$UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED %in% c("UNKNOWN", "UNKNOWN_LKB1_SUPPORT") &
                                            !is.na(M$NONE_CAPTURED_SEQ_SUPPORTED))
  unresolved <- chr0(M$LKB1_STATE_FROZEN) == "PRESENT"
  ENC$N_LKB1_UNRESOLVED_POSITIVES_RECODED_AS_ABSENT <<-
    sum(unresolved & chr0(M$LKB1_USABLE_SENSITIVITY_A) != "yes" & !is.na(M$LKB1_absent_A) & M$LKB1_absent_A == 1) +
    sum(unresolved & chr0(M$LKB1_USABLE_SENSITIVITY_B) != "yes" & !is.na(M$LKB1_absent_B) & M$LKB1_absent_B == 1) +
    sum(M$UPSTREAM_ARCHITECTURE_SEQ_SUPPORTED %in% "UNKNOWN_LKB1_SUPPORT" & !is.na(M$NONE_CAPTURED_SEQ_SUPPORTED) &
        M$NONE_CAPTURED_SEQ_SUPPORTED == 1)
  nc_bin <- function(code, state, callable) sum(!is.na(code) & !(state == "PRESENT" | (state == "ABSENT" & callable == "yes")))
  ENC$N_NONCALLABLE_BINARY_STATES_USED_AS_0_OR_1 <<- ENC$N_TRAITS_BINARY_CONTRADICTIONS +
    nc_bin(M$CAMKK_absent, M$CAMKK_STATE, M$CAMKK_CALLABLE) + nc_bin(M$LKB1_absent, M$LKB1_STATE, M$LKB1_CALLABLE)
  M
}

## ---- clade map from R1 T2, never relabelled ----
clade_map <- function() {
  t2 <- rd(PATHS$T2); t2 <- unique(t2[, c("SPECIES", "CLADE")])
  if (anyDuplicated(t2$SPECIES)) stop("REFUSED: a species carries two clades in T2")
  setNames(t2$CLADE, t2$SPECIES)
}

## ---- one model row ----
MODEL_COLS <- c("MODEL_ID", "MODEL_FAMILY", "ANALYSIS_CLASS", "RESPONSE", "RESPONSE_VARIANT", "PREDICTOR",
                "PREDICTOR_VARIANT", "COVARIATE", "PHYLO_MODEL", "N", "COEF", "SE", "CI95_LO", "CI95_HI", "NOMINAL_P",
                "HOLM_FAMILY", "HOLM_P", "PHYLO_PARAM", "NONPHYLO_COEF", "NONPHYLO_P", "STATUS",
                "SOURCE_RESPONSE_SHA256", "SOURCE_PREDICTOR_SHA256", "TREE_SHA256", "ENGINE_SHA256")
pred_source <- function(pred) {
  if (grepl("^SCORE_|^CBS_|^CONTACT_MINUS", pred)) SHA$GAM
  else if (grepl("^CBM_", pred)) SHA$CBM
  else if (pred %in% c("AMPKg_absent", "AMPKb_absent")) SHA$TRAITS
  else if (pred %in% c("NONE_CAPTURED_SEQ_SUPPORTED", "LKB1_absent_A", "LKB1_absent_B")) paste0(SHA$UP2, "+", SHA$ARCH3)
  else SHA$UP2
}
fit_row <- function(M, m) {
  cov <- if (m$COVARIATE == "none") character(0) else m$COVARIATE
  pm <- sub("^phylolm_", "", m$PHYLO_MODEL)
  r <- fit_cont(M, m$RESPONSE, m$PREDICTOR, cov, pm)
  ontree <- M$SPECIES %in% tree$tip.label
  cols <- c(m$RESPONSE, m$PREDICTOR, cov)
  cc <- ontree & complete.cases(M[, cols, drop = FALSE])
  flow <- data.frame(MODEL_ID = m$MODEL_ID, N_START = nrow(M), N_ON_TREE = sum(ontree),
    N_RESPONSE_CALLABLE = sum(ontree & !is.na(M[[m$RESPONSE]])), N_PREDICTOR_CALLABLE = sum(ontree & !is.na(M[[m$PREDICTOR]])),
    N_COMPLETE_CASES = sum(cc), N_FINAL = r$n,
    EXCLUSION_REASONS = sprintf("not_on_tree=%d;response_NA=%d;predictor_NA=%d;covariate_NA=%d", sum(!ontree),
      sum(ontree & is.na(M[[m$RESPONSE]])), sum(ontree & !is.na(M[[m$RESPONSE]]) & is.na(M[[m$PREDICTOR]])),
      if (length(cov)) sum(ontree & !is.na(M[[m$RESPONSE]]) & !is.na(M[[m$PREDICTOR]]) & is.na(M[[cov]])) else 0L),
    stringsAsFactors = FALSE)
  if (r$n != sum(cc)) stop("REFUSED: sample flow does not reconcile for ", m$MODEL_ID)
  ok <- identical(r$status, "OK")
  row <- data.frame(MODEL_ID = m$MODEL_ID, MODEL_FAMILY = m$MODEL_FAMILY, ANALYSIS_CLASS = m$ANALYSIS_CLASS,
    RESPONSE = m$RESPONSE, RESPONSE_VARIANT = m$RESPONSE_VARIANT, PREDICTOR = m$PREDICTOR,
    PREDICTOR_VARIANT = m$PREDICTOR_VARIANT, COVARIATE = m$COVARIATE, PHYLO_MODEL = m$PHYLO_MODEL, N = r$n,
    COEF = if (ok) fmt(r$coef) else "NA", SE = if (ok) fmt(r$se) else "NA",
    CI95_LO = if (ok) fmt(r$coef - Z * r$se) else "NA", CI95_HI = if (ok) fmt(r$coef + Z * r$se) else "NA",
    NOMINAL_P = if (ok) fmt(r$p) else "NA", HOLM_FAMILY = m$HOLM_FAMILY, HOLM_P = "NA",
    PHYLO_PARAM = if (ok) fmt(r$alpha) else "NA", NONPHYLO_COEF = if (ok) fmt(r$np_coef) else "NA",
    NONPHYLO_P = if (ok) fmt(r$np_p) else "NA", STATUS = r$status,
    SOURCE_RESPONSE_SHA256 = SHA$CTD, SOURCE_PREDICTOR_SHA256 = pred_source(m$PREDICTOR),
    TREE_SHA256 = SHA$TREE, ENGINE_SHA256 = SHA$ENGINE, stringsAsFactors = FALSE)
  list(row = row, flow = flow, raw = r)
}
spec <- function(id, family, class, resp, rvar, pred, pvar, cov = "none", pm = "phylolm_lambda", holm = "-",
                 parent = "-", stype = "PRIMARY")
  list(MODEL_ID = id, MODEL_FAMILY = family, ANALYSIS_CLASS = class, RESPONSE = resp, RESPONSE_VARIANT = rvar,
       PREDICTOR = pred, PREDICTOR_VARIANT = pvar, COVARIATE = cov, PHYLO_MODEL = pm, HOLM_FAMILY = holm,
       PARENT = parent, SENSITIVITY_TYPE = stype)

## declared one-dimension sensitivities of a primary model
alpha_sensitivities <- function(p) {
  len <- p$RESPONSE == "CTD_LENGTH_MEDIAN"
  base <- if (len) "CTD_LENGTH" else "CTD_MEAN_DISORDER"
  mk <- function(suffix, stype, resp = p$RESPONSE, rvar = p$RESPONSE_VARIANT, cov = "none", pm = "phylolm_lambda")
    spec(paste0(p$MODEL_ID, "__", suffix), p$MODEL_FAMILY, "DECLARED_SENSITIVITY", resp, rvar, p$PREDICTOR,
         p$PREDICTOR_VARIANT, cov, pm, "-", p$MODEL_ID, stype)
  out <- list(mk("BM", "PHYLO_MODEL_BM", pm = "phylolm_BM"),
              mk("PROTEOME", "PROTEOME_SIZE_COVARIATE", cov = "log10_n_proteins"),
              mk("HIGHPROV", "ALPHA_HIGH_PROVENANCE", paste0(base, "_MEDIAN_HIGHPROV"), "high_provenance"),
              mk("NOMODELSENS", "ALPHA_NO_MODEL_SENSITIVE", paste0(base, "_MEDIAN_NO_MODEL_SENSITIVE"), "no_model_sensitive"),
              mk("DUPMIN", "DUPLICATE_MIN", paste0(base, "_MIN"), "min_duplicate_sensitivity"),
              mk("DUPMAX", "DUPLICATE_MAX", paste0(base, "_MAX"), "max_duplicate_sensitivity"))
  if (len) out <- c(out, list(
              mk("LOG1P", "LENGTH_DEFINITION_LOG1P", "log1p_CTD_LENGTH_MEDIAN", "log1p_sensitivity"),
              mk("DOMAINEND", "LENGTH_DEFINITION_DOMAIN_END", "CTD_LENGTH_DOMAIN_END_MEDIAN", "domain_end_definition_sensitivity")))
  out
}
predictor_swap <- function(p, pred, pvar, suffix, stype, cov = "none")
  spec(paste0(p$MODEL_ID, "__", suffix), p$MODEL_FAMILY, "DECLARED_SENSITIVITY", p$RESPONSE, p$RESPONSE_VARIANT,
       pred, pvar, cov, "phylolm_lambda", "-", p$MODEL_ID, stype)

## ---- run a catalogue, apply fixed-size Holm, attach sensitivity comparisons ----
run_catalogue <- function(M, specs) {
  res <- lapply(specs, function(m) { x <- fit_row(M, m); x$spec <- m; x })
  names(res) <- vapply(specs, `[[`, "", "MODEL_ID")
  res
}
apply_holm <- function(res, family, n) {
  ids <- names(res)[vapply(res, function(x) x$spec$HOLM_FAMILY == family, TRUE)]
  if (length(ids) != n) stop("REFUSED: family ", family, " has ", length(ids), " models, declared ", n)
  p <- vapply(ids, function(i) if (identical(res[[i]]$raw$status, "OK")) res[[i]]$raw$p else NA_real_, 0)
  h <- rep(NA_real_, n); okp <- !is.na(p)
  if (any(okp)) h[okp] <- p.adjust(p[okp], method = "holm", n = n)
  for (k in seq_along(ids)) res[[ids[k]]]$row$HOLM_P <- fmt(h[k])
  res
}
bind_rows <- function(res, ids) do.call(rbind, lapply(res[ids], `[[`, "row"))
bind_flow <- function(res, ids) do.call(rbind, lapply(res[ids], `[[`, "flow"))
sens_table <- function(res, ids) {
  do.call(rbind, lapply(ids, function(i) {
    x <- res[[i]]; pa <- res[[x$spec$PARENT]]
    ok <- identical(x$raw$status, "OK") && identical(pa$raw$status, "OK")
    lo <- function(r) r$coef - Z * r$se; hi <- function(r) r$coef + Z * r$se
    cbind(x$row, PARENT_MODEL_ID = x$spec$PARENT,
      DIRECTION_AGREES_PRIMARY = if (ok) (if (sign(x$raw$coef) == sign(pa$raw$coef)) "yes" else "no") else "NA",
      DELTA_BETA = if (ok) fmt(x$raw$coef - pa$raw$coef) else "NA",
      CI_OVERLAP_PRIMARY = if (ok) (if (lo(x$raw) <= hi(pa$raw) && lo(pa$raw) <= hi(x$raw)) "yes" else "no") else "NA",
      SENSITIVITY_TYPE = x$spec$SENSITIVITY_TYPE, stringsAsFactors = FALSE)
  }))
}
enc_table <- function() data.frame(CHECK = names(ENC), VALUE = unlist(ENC), stringsAsFactors = FALSE)
