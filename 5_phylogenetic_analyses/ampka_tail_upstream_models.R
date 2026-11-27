# >>> deposit header
# Purpose: Fit the two primary models of AMPKα C-terminal length and disorder on the absence of both captured
#          upstream kinase categories and the four secondary models on CAMKK-category and LKB1 absence, each family
#          with its own Holm adjustment, with their declared sensitivities, and tabulate the C-terminal measures
#          by upstream category.
# Methods: Phylogenetic regressions
# Inputs:  `ctdx_lib.R` and the inputs it checks; `03_reproduction/GAMMA_ALPHA_REPRODUCTION_GATE.tsv`
# Outputs: `05_upstream/CTD_upstream_primary.tsv`, `CTD_upstream_component_secondary.tsv`,
#          `CTD_upstream_sensitivity.tsv`, `sample_flow_upstream.tsv` and `encoding_checks_upstream.tsv`;
#          `07_tables/CTD_BY_UPSTREAM_ARCHITECTURE_DESCRIPTIVE.tsv`
# Tools:   R with ape and phylolm 2.6.5 (Pagel's lambda, Brownian motion as a sensitivity); versions: see envs/
# <<< deposit header
## deposit repo root: begin
.REPO_ROOT <- local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(f[1]) else "."
  normalizePath(file.path(d, ".."))
})
## deposit repo root: end
## Two upstream-repertoire models specified before the analysis (AMPKα C-terminal length and disorder regressed
## on the absence of both captured upstream kinase categories, Holm across the two), four secondary models guided
## by that hypothesis (length and disorder on CAMKK-category absence and on LKB1 absence, Holm across the four),
## their declared sensitivities, and a descriptive table of the C-terminal region (CTD) by the four known upstream
## categories. Species of unknown upstream state are meant never to enter a model of the primary or secondary
## families; the count below checks this.
source(file.path(.REPO_ROOT, "5_phylogenetic_analyses", "lib", "tail_model_lib.R"))
S <- "ctdx_upstream.R"
g <- read.delim(file.path(CTDX, "03_reproduction/GAMMA_ALPHA_REPRODUCTION_GATE.tsv"), comment.char = "#")
if (g$STATUS != "PASS") stop("REFUSED: gamma-alpha reproduction did not pass")
M <- build_frames()
F2 <- "TWO_PRESPECIFIED_CTD_UPSTREAM_REPERTOIRE_MODELS"; F4 <- "FOUR_SECONDARY_CTD_UPSTREAM_COMPONENT_MODELS"
L <- "CTD_LENGTH_MEDIAN"; D <- "CTD_MEAN_DISORDER_MEDIAN"
prim <- list(
  spec("UP_NONE_LEN", "UPSTREAM_REPERTOIRE", "PRIMARY", L, "median_primary", "NONE_CAPTURED", "rbx2_neither_captured_known_only", holm = F2),
  spec("UP_NONE_DIS", "UPSTREAM_REPERTOIRE", "PRIMARY", D, "median_primary", "NONE_CAPTURED", "rbx2_neither_captured_known_only", holm = F2))
secd <- list(
  spec("UP_CAMKK_LEN", "UPSTREAM_COMPONENT", "HYPOTHESIS_GUIDED_SECONDARY", L, "median_primary", "CAMKK_absent", "rbx2_callable_camkk_category_absence", holm = F4),
  spec("UP_CAMKK_DIS", "UPSTREAM_COMPONENT", "HYPOTHESIS_GUIDED_SECONDARY", D, "median_primary", "CAMKK_absent", "rbx2_callable_camkk_category_absence", holm = F4),
  spec("UP_LKB1_LEN", "UPSTREAM_COMPONENT", "HYPOTHESIS_GUIDED_SECONDARY", L, "median_primary", "LKB1_absent", "rbx2_callable_lkb1_absence", holm = F4),
  spec("UP_LKB1_DIS", "UPSTREAM_COMPONENT", "HYPOTHESIS_GUIDED_SECONDARY", D, "median_primary", "LKB1_absent", "rbx2_callable_lkb1_absence", holm = F4))
sens <- list()
for (p in prim) sens <- c(sens, alpha_sensitivities(p), list(
  predictor_swap(p, "NONE_CAPTURED_SEQ_SUPPORTED", "rbx3_sequence_supported_architecture", "SEQSUP", "UPSTREAM_SEQ_SUPPORTED_ARCHITECTURE"),
  predictor_swap(p, "NONE_CAPTURED_ASGIVEN", "rbx2_frozen_binary_as_given_admits_4_unknown", "ASGIVEN", "UPSTREAM_FROZEN_BINARY_AS_GIVEN")))
for (p in secd[3:4]) sens <- c(sens, list(
  predictor_swap(p, "LKB1_absent_A", "rbx3_lkb1_tier_A", "LKB1_TIER_A", "LKB1_SEQUENCE_TIER_A"),
  predictor_swap(p, "LKB1_absent_B", "rbx3_lkb1_tier_B", "LKB1_TIER_B", "LKB1_SEQUENCE_TIER_B")))
res <- run_catalogue(M, c(prim, secd, sens))
res <- apply_holm(res, F2, 2); res <- apply_holm(res, F4, 4)
pid <- vapply(prim, `[[`, "", "MODEL_ID"); cid <- vapply(secd, `[[`, "", "MODEL_ID"); sid <- vapply(sens, `[[`, "", "MODEL_ID")
## count the species of unknown upstream state that enter a primary or secondary fit; the count, expected to be
## zero, is written to the table of encoding checks (the script does not stop on it)
used_unknown <- 0
for (i in c(pid, cid)) {
  m <- res[[i]]$spec
  cc <- M$SPECIES %in% tree$tip.label & complete.cases(M[, c(m$RESPONSE, m$PREDICTOR)])
  used_unknown <- used_unknown + sum(cc & M$UPSTREAM_KINASE_STATE == "UNKNOWN")
}
ENC$N_UNKNOWN_UPSTREAM_STATES_USED_IN_MODELS <- used_unknown

N1 <- c("Cross-species phylogenetic associations only. NEITHER_CAPTURED means neither recognizable CAMKK-category nor",
        "LKB1 machinery was captured; it is not a statement that no AMPK kinase exists.",
        "RBX2 upstream states govern primary and secondary models; RBX3 enters only named sensitivities.",
        "HOLM_P uses a fixed family size equal to the declared family; an inestimable model keeps HOLM_P NA.")
write_tsv("05_upstream/CTD_upstream_primary.tsv", bind_rows(res, pid), N1, S)
write_tsv("05_upstream/CTD_upstream_component_secondary.tsv", bind_rows(res, cid), c(N1,
  "HYPOTHESIS_GUIDED_SECONDARY: a decomposition of the primary repertoire contrast, Holm within its own four."), S)
write_tsv("05_upstream/CTD_upstream_sensitivity.tsv", sens_table(res, sid), c(N1,
  "ASGIVEN admits the four UNKNOWN species that the frozen binary codes ANY_CAPTURED; it is a declared sensitivity only.",
  "LKB1 tiers set a frozen LKB1-positive species that fails the tier to NA, never to ABSENT."), S)
write_tsv("05_upstream/sample_flow_upstream.tsv", bind_flow(res, c(pid, cid, sid)), "N_FINAL equals the engine's fitted N.", S)
write_tsv("05_upstream/encoding_checks_upstream.tsv", enc_table(), "Counts feed the z5 gates.", S)

## descriptive table only, no inference
q <- function(x, pr) if (length(x)) sprintf("%.6g", quantile(x, pr, names = FALSE)) else "NA"
cats <- c("BOTH", "CAMKK_CATEGORY_ONLY", "LKB1_ONLY", "NEITHER_CAPTURED")
desc <- do.call(rbind, lapply(cats, function(k) {
  s <- M[M$UPSTREAM_KINASE_STATE == k & M$SPECIES %in% tree$tip.label, ]
  l <- s$CTD_LENGTH_MEDIAN[!is.na(s$CTD_LENGTH_MEDIAN)]; d <- s$CTD_MEAN_DISORDER_MEDIAN[!is.na(s$CTD_MEAN_DISORDER_MEDIAN)]
  data.frame(UPSTREAM_KINASE_STATE = k, N_SPECIES = nrow(s), N_LENGTH = length(l),
             CTD_LENGTH_MEDIAN = q(l, .5), CTD_LENGTH_Q1 = q(l, .25), CTD_LENGTH_Q3 = q(l, .75),
             CTD_LENGTH_IQR = if (length(l)) sprintf("%.6g", IQR(l)) else "NA", N_DISORDER = length(d),
             CTD_DISORDER_MEDIAN = q(d, .5), CTD_DISORDER_Q1 = q(d, .25), CTD_DISORDER_Q3 = q(d, .75),
             CTD_DISORDER_IQR = if (length(d)) sprintf("%.6g", IQR(d)) else "NA", stringsAsFactors = FALSE)
}))
write_tsv("07_tables/CTD_BY_UPSTREAM_ARCHITECTURE_DESCRIPTIVE.tsv", desc, c(
  "Descriptive only: species-level CTD medians by frozen RBX2 upstream category, on-tree species, UNKNOWN excluded.",
  "No p value is attached to these four categories; the inferential contrast is the NONE_CAPTURED model."), S)
cat("z3 complete\n")
