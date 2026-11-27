# >>> deposit header
# Purpose: Assigns each candidate tip of a CAMKK or LKB1 kinase-domain tree, and each reference tip by leave-one-out,
#          to a class from the smallest split holding it with another reference; tests that rerooting changes no call.
# Methods: Orthology calls
# Inputs:  arguments: an IQ-TREE `.treefile`, the kinase-domain table `<target>_kinase_domain_manifest.tsv`, the
#          output path and the target (CAMKK or LKB1)
# Outputs: the tree-call table (for example `CAMKK_tree_calls_masked.tsv`) and `<target>_reroot_invariance.tsv`
# Tools:   R with ape; versions: see envs/
# <<< deposit header
## Classifies each candidate and reference tip of a CAMKK or LKB1 kinase-domain tree by a fixed rule, set before
## any tree was examined. IQ-TREE produces an unrooted maximum-likelihood tree, so the rule is defined on
## bipartitions (splits), never on a displayed root: for each tip, take the smallest side of any split that holds
## the tip and at least one other reference sequence. The smallest such split is taken whatever its support; its
## ultrafast bootstrap support then grades the call (95 or more supported, at least 70 provisional, below 70
## unresolved), and a class is named only if the split holds references of one class and no negative sequence.
## Bipartitions do not change on rerooting, so neither does the call, and the rerooting test below checks this
## rather than assuming it. Nothing is classified by taxonomy or by an annotation name.
suppressMessages(library(ape))
a <- commandArgs(trailingOnly = TRUE)
TREEFILE <- a[1]; MANIFEST <- a[2]; OUT <- a[3]; LANE <- a[4]
SUPPORTED <- 95; PROVISIONAL <- 70

rd <- function(p) read.delim(p, comment.char = "#", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
man <- rd(MANIFEST)
## A Newick label cannot contain a colon, and the tree writer replaces one with an underscore. Three CAMKK
## candidates have a colon in their protein identifier, so a name built by replacing only "." would not match the
## tree, and the intersect() calls below would drop it silently. Every character a Newick label cannot hold is
## therefore replaced in the same way here, so the two sides agree; a name without such characters is unchanged.
newick_safe <- function(s) gsub("[^A-Za-z0-9_-]", "_", s)
man$TIP <- newick_safe(paste0(man$SPECIES, "__", gsub("\\.", "_", man$PROTEIN_ID)))
phy0 <- read.tree(TREEFILE)
if (is.null(phy0)) stop("REFUSED: could not read ", TREEFILE)

neg <- intersect(man$TIP[man$ROLE == "NEGATIVE"], phy0$tip.label)
ref <- intersect(man$TIP[grepl("REFERENCE", man$ROLE)], phy0$tip.label)
cand <- intersect(man$TIP[man$ROLE == "CANDIDATE"], phy0$tip.label)
## intersect() silently drops a name it cannot find, so every tip of the kinase-domain table must be in the tree;
## if one is not, the script stops here rather than writing a shorter table.
absent <- setdiff(man$TIP, phy0$tip.label)
if (length(absent)) {
  stop(sprintf("REFUSED: %d of %d manifest tips are not in the tree (%s%s); every copy must reach the tree",
               length(absent), nrow(man), paste(head(absent, 3), collapse = ", "),
               if (length(absent) > 3) ", ..." else ""))
}
if (length(cand) != sum(man$ROLE == "CANDIDATE")) {
  stop(sprintf("REFUSED: %d candidates in the manifest but %d matched in the tree",
               sum(man$ROLE == "CANDIDATE"), length(cand)))
}
cls <- setNames(man$CLASS_L1, man$TIP); cls2 <- setNames(man$CLASS_L2, man$TIP)
rid <- setNames(man$REFERENCE_ID, man$TIP)

## ---- every non-trivial bipartition of the unrooted tree, with the support value of its edge ----
bipartitions <- function(phy) {
  phy <- reorder(phy, "postorder")
  ntip <- length(phy$tip.label)
  supp <- suppressWarnings(as.numeric(phy$node.label))
  desc <- vector("list", ntip + phy$Nnode)
  for (i in seq_len(ntip)) desc[[i]] <- phy$tip.label[i]
  for (e in seq_len(nrow(phy$edge))) {
    p <- phy$edge[e, 1]; ch <- phy$edge[e, 2]
    desc[[p]] <- c(desc[[p]], desc[[ch]])
  }
  all_tips <- phy$tip.label
  root <- ntip + 1L
  out <- list(); seen <- character(0)
  for (nd in (ntip + 1L):(ntip + phy$Nnode)) {
    if (nd == root) next
    A <- desc[[nd]]
    if (length(A) < 2 || length(A) > ntip - 2) next
    B <- setdiff(all_tips, A)
    small <- if (length(A) <= length(B)) A else B
    key <- paste(sort(small), collapse = "|")
    if (key %in% seen) next
    seen <- c(seen, key)
    out[[length(out) + 1]] <- list(A = A, B = B, key = key, support = supp[nd - ntip])
  }
  out
}

## Support belongs to an edge, but ape stores it as a node label and does not move node labels when a tree is
## rerooted. On this dataset, rerooting left the set of bipartitions unchanged (179 of 179 shared) but moved the
## support value on up to 8 of them. The support table is therefore built once, from the tree as IQ-TREE wrote
## it, keyed by the bipartition itself, and never re-read from a rerooted copy; this makes the assignment
## independent of the root.
support_map <- function(bips) {
  ## Keyed by the bipartition. The table is held in two parallel vectors rather than in an R environment,
  ## because environment names are R symbols, limited to 10,000 bytes, and a bipartition key exceeds that in a
  ## large tree (196 tips is within the limit, 683 is not). Keys are unique, because the builder above skips a
  ## bipartition it has already seen, so the first match is the only match.
  list(keys = vapply(bips, function(b) b$key, ""),
       supp = lapply(bips, function(b) b$support))
}
lookup <- function(m, key) { i <- match(key, m$keys); if (is.na(i)) NA else m$supp[[i]] }

classify <- function(bips, smap) {
  res <- list()
  for (tp in c(cand, ref)) {
    best <- NULL
    for (b in bips) {
      side <- if (tp %in% b$A) b$A else b$B
      r_in <- setdiff(intersect(side, ref), tp)
      if (length(r_in) == 0) next
      if (is.null(best) || length(side) < length(best$side)) {
        best <- list(side = side, refs = r_in, negs = intersect(side, neg),
                     support = lookup(smap, b$key))
      }
    }
    if (is.null(best)) {
      res[[tp]] <- list(node = NA, support = NA, refs = character(0), negs = character(0),
                        l1 = "SUBFAMILY_UNRESOLVED", l2 = "-", conf = "UNRESOLVED",
                        basis = "no supported split places this tip with any reference", n = NA)
      next
    }
    sup <- best$support
    conf <- if (is.na(sup)) "UNRESOLVED" else if (sup >= SUPPORTED) "SUPPORTED" else
            if (sup >= PROVISIONAL) "PROVISIONAL_PLACEMENT" else "UNRESOLVED"
    l1s <- unique(cls[best$refs]); l2s <- unique(cls2[best$refs]); l2s <- l2s[l2s != "-"]
    basis <- "smallest supported split containing this tip and a reference"
    l1 <- if (length(best$negs) > 0) {
      basis <- "split contains a negative-family reference"; "SUBFAMILY_UNRESOLVED"
    } else if (length(l1s) != 1) {
      basis <- "split contains references of more than one constituent family"; "SUBFAMILY_UNRESOLVED"
    } else if (conf == "UNRESOLVED") {
      basis <- sprintf("support %s below the provisional floor of %d",
                       ifelse(is.na(sup), "NA", sprintf("%.1f", sup)), PROVISIONAL)
      "SUBFAMILY_UNRESOLVED"
    } else l1s[1]
    l2 <- if (l1 == "SUBFAMILY_UNRESOLVED" || LANE != "CAMKK") "-" else
          if (length(l2s) == 1) l2s[1] else
          if (l1 == "SNF1_SNRK1_ACTIVATOR_LIKE") "SNF1_ACTIVATOR_SUBTYPE_UNRESOLVED" else "-"
    res[[tp]] <- list(node = length(best$side), support = sup, refs = best$refs, negs = best$negs,
                      l1 = l1, l2 = l2, conf = conf, basis = basis, n = length(best$side))
  }
  res
}

base_bips <- bipartitions(phy0)
SMAP <- support_map(base_bips)          # the support table, read once from the tree as written
base <- classify(base_bips, SMAP)

## ---- rerooting test: the call must not depend on where the tree is rooted for display ----
set.seed(488761235)
probe <- unique(c(head(phy0$tip.label, 1), sample(phy0$tip.label, min(5, length(phy0$tip.label)))))
changes <- 0; tested <- 0; setdiffs <- 0; ape_relabel <- 0
basekeys <- sort(vapply(base_bips, function(b) b$key, ""))
for (tp in probe) {
  rr <- tryCatch(root(unroot(phy0), outgroup = tp, resolve.root = TRUE), error = function(e) NULL)
  if (is.null(rr)) next
  tested <- tested + 1
  rb <- bipartitions(rr)
  if (!identical(sort(vapply(rb, function(b) b$key, "")), basekeys)) setdiffs <- setdiffs + 1
  for (b in rb) {
    s0 <- lookup(SMAP, b$key)
    if (!is.na(s0) && !is.na(b$support) && !identical(s0, b$support)) ape_relabel <- ape_relabel + 1
  }
  alt <- classify(rb, SMAP)             # the same support table, on the rerooted topology
  for (k in names(base)) {
    if (is.null(alt[[k]]) || !identical(alt[[k]]$l1, base[[k]]$l1) ||
        !identical(alt[[k]]$l2, base[[k]]$l2) || !identical(alt[[k]]$conf, base[[k]]$conf)) {
      changes <- changes + 1
    }
  }
}
cat(sprintf("REROOT_INVARIANCE: %d rerootings, N_BIPARTITION_SET_DIFFERENCES = %d, N_REROOT_CLASSIFICATION_CHANGES = %d\n",
            tested, setdiffs, changes))
cat(sprintf("APE_NODE_LABEL_REASSIGNMENTS_OBSERVED = %d (an ape artefact, not a property of the rule)\n",
            ape_relabel))

res <- do.call(rbind, lapply(names(base), function(tp) {
  b <- base[[tp]]
  data.frame(TIP = tp,
             ROLE = ifelse(tp %in% ref, "REFERENCE_LEAVE_ONE_OUT", "CANDIDATE"),
             VALIDATION_TYPE = ifelse(tp %in% ref, "REFERENCE_LOO_SELF_CONSISTENCY", "CANDIDATE_CLASSIFICATION"),
             SPLIT_SIZE = ifelse(is.na(b$n), "NA", as.character(b$n)),
             TREE_SUPPORT = ifelse(is.na(b$support), "NA", sprintf("%.1f", b$support)),
             N_REFERENCES_IN_SPLIT = length(b$refs), N_NEGATIVES_IN_SPLIT = length(b$negs),
             REFERENCE_IDS = if (length(b$refs)) paste(sort(unique(ifelse(rid[b$refs] == "-", b$refs,
                             rid[b$refs]))), collapse = ";") else "-",
             NEGATIVE_TIPS = if (length(b$negs)) paste(head(sort(b$negs), 8), collapse = ";") else "-",
             CLASS_L1 = b$l1, CLASS_L2 = b$l2, CALL_CONFIDENCE = b$conf, CALL_BASIS = b$basis,
             stringsAsFactors = FALSE)
}))
hdr <- c(sprintf("## RBX3 %s assignment from supported bipartitions of the unrooted tree; no displayed root", LANE),
         sprintf("## enters the rule. REROOT_INVARIANCE tested at %d rerootings, N_REROOT_CLASSIFICATION_CHANGES=%d.",
                 tested, changes),
         "## Support is ultrafast bootstrap from IQ-TREE -B 1000. SUPPORTED >=95 is the conventional",
         "## confidence cut; 70 to 94.999 is PROVISIONAL_PLACEMENT and is never a validated family call.",
         "## Declared references are assigned leave-one-out and their VALIDATION_TYPE reads",
         "## REFERENCE_LOO_SELF_CONSISTENCY, which is self-consistency and not independent validation.",
         sprintf("## tree=%s tips=%d candidates=%d references=%d negatives=%d bipartitions=%d",
                 basename(TREEFILE), length(phy0$tip.label), length(cand), length(ref), length(neg),
                 length(base_bips)))
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
writeLines(hdr, OUT)
write.table(res, OUT, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
writeLines(c("QUANTITY\tVALUE\tNOTE",
             sprintf("N_REROOTINGS_TESTED\t%d\t-", tested),
             sprintf("N_BIPARTITION_SET_DIFFERENCES\t%d\tthe topology content must be identical under rerooting", setdiffs),
             sprintf("N_REROOT_CLASSIFICATION_CHANGES\t%d\tthe gate; must be zero", changes),
             sprintf("APE_NODE_LABEL_REASSIGNMENTS_OBSERVED\t%d\tape carries edge support as node labels and does not re-map them on rerooting; support is therefore keyed by bipartition and read once from the tree as IQ-TREE wrote it", ape_relabel),
             sprintf("RESULT\t%s\t-", ifelse(changes == 0 && setdiffs == 0, "PASS", "FAIL"))),
           file.path(dirname(OUT), sprintf("%s_reroot_invariance.tsv", LANE)))
cat(sprintf("%s: %d tips classified (%d candidates, %d leave-one-out references)\n",
            LANE, nrow(res), sum(res$ROLE == "CANDIDATE"), sum(res$ROLE != "CANDIDATE")))
print(table(res$CLASS_L1, res$CALL_CONFIDENCE))
if (changes != 0 || setdiffs != 0) stop("REFUSED: classification is not invariant to rerooting")
