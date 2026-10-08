#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(digest)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript fit_activity_expression_loeuf.R INPUT_TSV OUTPUT_DIR")
input <- normalizePath(args[[1]], mustWork = TRUE)
OUT <- normalizePath(args[[2]], mustWork = FALSE)
for (d in c("03_activity_primary", "07_logs")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
d <- fread(input, na.strings = c("NA", ""), check.names = FALSE)
if (nrow(d) != 1267L) stop("Expected frozen exact-match Activity expression complete cases N=1,267; observed ", nrow(d))
if (anyDuplicated(d$gene) || !identical(as.integer(d$model_row_order), seq_len(nrow(d)))) stop("Model input row order/unique gene QA failed")
d[, `:=`(hit = as.integer(hit), exposure = as.integer(exposure), bait_status = as.integer(bait_status),
          guide_count_z = as.numeric(guide_count_z), loeuf_z = as.numeric(loeuf_z), baseline_expression_z = as.numeric(baseline_expression_z))]
d[, selection_stratum := factor(selection_stratum)]
if (anyNA(d[, .(hit, exposure, bait_status, guide_count_z, loeuf_z, baseline_expression_z)])) stop("Model complete-case input contains missing covariates")
if (!all(d$hit %in% c(0L, 1L)) || !all(d$exposure %in% c(0L, 1L))) stop("Binary variable encoding invalid")

formulae <- c(
  reduced = "hit ~ exposure + bait_status + selection_stratum + guide_count_z + loeuf_z",
  full = "hit ~ exposure + bait_status + selection_stratum + guide_count_z + loeuf_z + baseline_expression_z"
)
gene_hash <- function(x) digest(paste(sort(unique(as.character(x))), collapse = "\n"), algo = "sha256", serialize = FALSE)
log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
firth_eval <- function(beta, X, y, want_score = TRUE) {
  eta <- as.vector(X %*% beta); pr <- plogis(eta); w <- pmax(pr * (1 - pr), 1e-12)
  info <- crossprod(X, X * w)
  deti <- try(determinant(info, logarithm = TRUE), silent = TRUE)
  if (inherits(deti, "try-error") || deti$sign <= 0 || !is.finite(as.numeric(deti$modulus))) return(list(ll = -Inf, score = rep(NA_real_, ncol(X))))
  ll <- sum(y * eta - log1pexp(eta)) + 0.5 * as.numeric(deti$modulus)
  if (!want_score) return(list(ll = ll, score = NULL))
  inv <- try(solve(info), silent = TRUE)
  if (inherits(inv, "try-error")) return(list(ll = -Inf, score = rep(NA_real_, ncol(X))))
  h <- w * rowSums((X %*% inv) * X)
  list(ll = ll, score = as.vector(crossprod(X, y - pr + h * (0.5 - pr))))
}
fit_firth <- function(X, y) {
  X <- as.matrix(X); storage.mode(X) <- "double"; y <- as.numeric(y)
  qrx <- qr(X)
  if (qrx$rank < ncol(X)) return(list(status = "RANK_DEFICIENT", rank = qrx$rank))
  opt <- try(optim(rep(0, ncol(X)), fn = function(b) -firth_eval(b, X, y)$ll,
    gr = function(b) -firth_eval(b, X, y)$score, method = "BFGS", control = list(maxit = 3000, reltol = 1e-10)), silent = TRUE)
  if (inherits(opt, "try-error") || !is.finite(opt$value) || opt$convergence != 0L) return(list(status = "OPTIMIZATION_FAILED", rank = qrx$rank))
  b <- opt$par; eta <- as.vector(X %*% b); w <- pmax(plogis(eta) * (1 - plogis(eta)), 1e-12)
  V <- try(solve(crossprod(X, X * w)), silent = TRUE)
  if (inherits(V, "try-error")) return(list(status = "COVARIANCE_FAILED", rank = qrx$rank))
  list(status = "COMPLETE", beta = b, V = V, max_ll = -opt$value, rank = qrx$rank, X = X, y = y)
}
profile_exposure <- function(fit, term = "exposure") {
  X <- fit$X; y <- fit$y; target <- match(term, colnames(X))
  if (is.na(target)) return(list(status = "EXPOSURE_TERM_MISSING"))
  profile_ll <- function(value) {
    start <- fit$beta[-target]
    objective <- function(u) { b <- numeric(ncol(X)); b[target] <- value; b[-target] <- u; -firth_eval(b, X, y, FALSE)$ll }
    gradient <- function(u) { b <- numeric(ncol(X)); b[target] <- value; b[-target] <- u; -firth_eval(b, X, y)$score[-target] }
    z <- optim(start, objective, gr = gradient, method = "BFGS", control = list(maxit = 2000, reltol = 1e-9))
    if (!is.finite(z$value)) -Inf else -z$value
  }
  ll0 <- profile_ll(0); p <- pchisq(max(0, 2 * (fit$max_ll - ll0)), 1, lower.tail = FALSE); cutoff <- qchisq(.95, 1)
  bound <- function(direction) {
    center <- fit$beta[target]; step <- 1; edge <- center + direction * step; val <- 2 * (fit$max_ll - profile_ll(edge))
    while (is.finite(val) && val < cutoff && step < 40) { step <- step * 1.6; edge <- center + direction * step; val <- 2 * (fit$max_ll - profile_ll(edge)) }
    if (!is.finite(val) || val < cutoff) return(NA_real_)
    tryCatch(uniroot(function(z) 2 * (fit$max_ll - profile_ll(z)) - cutoff, min(center, edge), max(center, edge), tol = 1e-6)$root, error = function(e) NA_real_)
  }
  list(status = "COMPLETE", beta = fit$beta[target], ci_low = bound(-1), ci_high = bound(1), p = p)
}

mle_fits <- list(); warnings_by <- list(); bad <- FALSE
for (id in names(formulae)) {
  warns <- character()
  fit <- withCallingHandlers(try(glm(as.formula(formulae[[id]]), data = as.data.frame(d), family = binomial()), silent = TRUE),
    warning = function(w) {warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")})
  warnings_by[[id]] <- warns
  ok <- !inherits(fit, "try-error") && !length(warns) && isTRUE(fit$converged)
  if (ok) {
    cc <- coef(fit); vv <- try(vcov(fit), silent = TRUE)
    ok <- !inherits(vv, "try-error") && "exposure" %in% names(cc) && is.finite(cc["exposure"]) &&
      is.finite(vv["exposure", "exposure"]) && sqrt(vv["exposure", "exposure"]) > 0 && abs(cc["exposure"]) <= 15
  }
  if (!ok) bad <- TRUE
  mle_fits[[id]] <- fit
}
estimator <- if (bad) "Firth; custom Jeffreys-penalized likelihood; uniform fallback for paired same-case reduced/full models" else "Maximum-likelihood logistic regression; Wald CI/P"
coeff_rows <- list(); model_rows <- list(); diagnostic_rows <- list(); firth_fits <- list()
if (bad) {
  for (id in names(formulae)) {
    X <- model.matrix(as.formula(formulae[[id]]), data = as.data.frame(d))
    firth_fits[[id]] <- fit_firth(X, d$hit)
    if (firth_fits[[id]]$status != "COMPLETE") stop("Uniform Firth fallback failed for ", id, ": ", firth_fits[[id]]$status)
  }
}
for (id in names(formulae)) {
  if (!bad) {
    fit <- mle_fits[[id]]; ss <- coef(summary(fit)); V <- vcov(fit)
    for (term in rownames(ss)) {
      b <- ss[term, "Estimate"]; se <- ss[term, "Std. Error"]; p <- ss[term, "Pr(>|z|)"]
      coeff_rows[[length(coeff_rows)+1L]] <- data.table(model_id=id, term=term, estimate_log_odds=b, standard_error=se,
        CI_low_log_odds=b-1.96*se, CI_high_log_odds=b+1.96*se, OR=exp(b), OR_CI_low=exp(b-1.96*se), OR_CI_high=exp(b+1.96*se), P_raw=p,
        inference=if (term=="exposure") "Wald z" else "Wald z")
    }
    b <- coef(fit)["exposure"]; se <- sqrt(V["exposure", "exposure"]); p <- 2*pnorm(-abs(b/se)); lo <- b-1.96*se; hi <- b+1.96*se
    status <- "COMPLETE_MLE"
  } else {
    ff <- firth_fits[[id]]; bvec <- ff$beta; V <- ff$V; terms <- names(bvec); profile <- profile_exposure(ff)
    for (term in terms) {
      b0 <- bvec[term]; se0 <- sqrt(V[term, term]); lo0 <- b0-1.96*se0; hi0 <- b0+1.96*se0; p0 <- 2*pnorm(-abs(b0/se0)); inf <- "Wald approximation from penalized-information covariance"
      if (term == "exposure") { b0 <- profile$beta; lo0 <- profile$ci_low; hi0 <- profile$ci_high; p0 <- profile$p; inf <- "profile likelihood CI and likelihood-ratio P" }
      coeff_rows[[length(coeff_rows)+1L]] <- data.table(model_id=id, term=term, estimate_log_odds=b0, standard_error=se0,
        CI_low_log_odds=lo0, CI_high_log_odds=hi0, OR=exp(b0), OR_CI_low=exp(lo0), OR_CI_high=exp(hi0), P_raw=p0, inference=inf)
    }
    b <- profile$beta; lo <- profile$ci_low; hi <- profile$ci_high; p <- profile$p; status <- "COMPLETE_FIRTH_UNIFORM_FALLBACK"
  }
  coef_table <- rbindlist(coeff_rows)[model_id == id]
  get_term <- function(term, col) { term_key <- term; z <- coef_table[coef_table$term == term_key]; if (!nrow(z)) NA_real_ else as.numeric(z[[col]][1]) }
  model_rows[[length(model_rows)+1L]] <- data.table(model_id=id, status=status, N=nrow(d), events=sum(d$hit), prey_N=sum(d$exposure==1L),
    non_PPI_nonbait_N=sum(d$exposure==0L & d$bait_status==0L), bait_N=sum(d$exposure==0L & d$bait_status==1L), exposure_zero_N=sum(d$exposure==0L),
    beta_exposure=b, OR=exp(b), CI_low_OR=exp(lo), CI_high_OR=exp(hi), P_raw=p,
    beta_loeuf=get_term("loeuf_z", "estimate_log_odds"), loeuf_CI_low_OR=get_term("loeuf_z", "OR_CI_low"), loeuf_CI_high_OR=get_term("loeuf_z", "OR_CI_high"), loeuf_P_raw=get_term("loeuf_z", "P_raw"),
    beta_baseline_expression=get_term("baseline_expression_z", "estimate_log_odds"), expression_CI_low_OR=get_term("baseline_expression_z", "OR_CI_low"), expression_CI_high_OR=get_term("baseline_expression_z", "OR_CI_high"), expression_P_raw=get_term("baseline_expression_z", "P_raw"),
    estimator=estimator, formula=formulae[[id]], population_hash=gene_hash(d$gene))
  diagnostic_rows[[length(diagnostic_rows)+1L]] <- data.table(model_id=id, estimator=estimator, mle_converged=if (bad) "NOT_USED_UNIFORM_FIRTH" else isTRUE(mle_fits[[id]]$converged),
    mle_warning_count=length(warnings_by[[id]]), mle_warnings=paste(warnings_by[[id]], collapse="; "),
    exposure_abs_beta_gt_15=abs(b)>15, separation_warning_detected=any(grepl("0 or 1|separation|did not converge", warnings_by[[id]], ignore.case=TRUE)),
    firth_status=if (bad) firth_fits[[id]]$status else "NOT_REQUIRED", model_matrix_rank=if (bad) firth_fits[[id]]$rank else NA_integer_, rows=NROW(d), population_hash=gene_hash(d$gene))
}
models <- rbindlist(model_rows, fill=TRUE); coefs <- rbindlist(coeff_rows, fill=TRUE); diagnostics <- rbindlist(diagnostic_rows, fill=TRUE)
if (length(unique(models$population_hash)) != 1L || length(unique(models$N)) != 1L) stop("Reduced/full model population hash or N differs")
if (!identical(as.integer(d$model_row_order), seq_len(nrow(d)))) stop("The shared model input order changed")
reduced <- models[model_id=="reduced"]; full <- models[model_id=="full"]
attenuation <- if (is.finite(reduced$beta_exposure) && reduced$beta_exposure != 0) 100*(reduced$beta_exposure-full$beta_exposure)/reduced$beta_exposure else NA_real_
or_change <- 100*(full$OR-reduced$OR)/reduced$OR
primary <- data.table(
  comparison_id="ACTIVITY_BASELINE_EXPRESSION_PLUS_LOEUF", same_complete_case_population=TRUE, N=full$N, events=full$events,
  prey_N=full$prey_N, non_PPI_nonbait_N=full$non_PPI_nonbait_N, bait_N=full$bait_N, exposure_zero_N=full$exposure_zero_N,
  reduced_formula=reduced$formula, reduced_prey_OR=reduced$OR, reduced_prey_CI_low_OR=reduced$CI_low_OR, reduced_prey_CI_high_OR=reduced$CI_high_OR, reduced_prey_P_raw=reduced$P_raw,
  full_formula=full$formula, full_prey_OR=full$OR, full_prey_CI_low_OR=full$CI_low_OR, full_prey_CI_high_OR=full$CI_high_OR, full_prey_P_raw=full$P_raw,
  logOR_attenuation_percent=attenuation, OR_change_percent=or_change,
  full_model_baseline_expression_beta=full$beta_baseline_expression, full_model_baseline_expression_OR_CI_low=full$expression_CI_low_OR, full_model_baseline_expression_OR_CI_high=full$expression_CI_high_OR, full_model_baseline_expression_P_raw=full$expression_P_raw,
  full_model_loeuf_beta=full$beta_loeuf, full_model_loeuf_OR_CI_low=full$loeuf_CI_low_OR, full_model_loeuf_OR_CI_high=full$loeuf_CI_high_OR, full_model_loeuf_P_raw=full$loeuf_P_raw,
  estimator=estimator, population_hash=full$population_hash)
fwrite(primary, file.path(OUT,"03_activity_primary/ACTIVITY_PRIMARY_BASELINE_EXPR_LOEUF_RESULTS.tsv"), sep="\t", na="NA", quote=FALSE)
fwrite(coefs, file.path(OUT,"03_activity_primary/ACTIVITY_PRIMARY_MODEL_COEFFICIENTS.tsv"), sep="\t", na="NA", quote=FALSE)
fwrite(models, file.path(OUT,"03_activity_primary/ACTIVITY_PAIRED_MODEL_ESTIMATES.tsv"), sep="\t", na="NA", quote=FALSE)
fwrite(diagnostics, file.path(OUT,"03_activity_primary/ACTIVITY_MODEL_DIAGNOSTICS.tsv"), sep="\t", na="NA", quote=FALSE)
counts <- data.table(population="same_expression_complete_case_population", N=nrow(d), events=sum(d$hit), prey_N=sum(d$exposure==1L), non_PPI_nonbait_N=sum(d$exposure==0L & d$bait_status==0L), bait_N=sum(d$exposure==0L & d$bait_status==1L), exposure_zero_N=sum(d$exposure==0L), population_hash=gene_hash(d$gene))
fwrite(counts, file.path(OUT,"03_activity_primary/ACTIVITY_COMPLETE_CASE_COUNTS.tsv"), sep="\t", na="NA", quote=FALSE)
fwrite(data.table(model_row_order=d$model_row_order, gene=d$gene), file.path(OUT,"03_activity_primary/ACTIVITY_R_MODEL_POPULATION_ORDER.tsv"), sep="\t", na="NA", quote=FALSE)
session <- capture.output(sessionInfo())
writeLines(c(paste("R:", R.version.string), paste("data.table:", as.character(packageVersion("data.table"))), paste("digest:", as.character(packageVersion("digest"))), "", session), file.path(OUT,"07_logs/R_SESSION_INFO.txt"))
writeLines(c("# Activity model QA", "", paste0("Same complete-case population was used for both models (N=",nrow(d),", events=",sum(d$hit),"). Population hash: ",gene_hash(d$gene),"."),
  paste0("Reduced formula: `",formulae[["reduced"]],"`."), paste0("Full formula: `",formulae[["full"]],"`."),
  paste0("Estimator: ",estimator,"."), paste0("Reduced prey OR ",signif(reduced$OR,5)," (95% CI ",signif(reduced$CI_low_OR,5),"–",signif(reduced$CI_high_OR,5),"), P=" ,signif(reduced$P_raw,5),"."),
  paste0("Full prey OR ",signif(full$OR,5)," (95% CI ",signif(full$CI_low_OR,5),"–",signif(full$CI_high_OR,5),"), P=" ,signif(full$P_raw,5),"."),
  paste0("Descriptive log(OR) attenuation = ",signif(attenuation,5),"%; OR change = ",signif(or_change,5),"%."),
  paste0("Baseline expression beta ",signif(full$beta_baseline_expression,5),"; LOEUF beta ",signif(full$beta_loeuf,5)," in the full model."),
  "Covariate and model uncertainty is reported in `ACTIVITY_PRIMARY_MODEL_COEFFICIENTS.tsv`; inferential scale for the prey term is OR with OR-scale 95% CI.",
  "The baseline was estimated from a closely matched same-study Day-21 WTC11 iPSC-neuron CROP-seq NTC reference; this does not establish exact CaMPARI2 batch matching.",
  "No P-dependent source, population, transform, covariate or threshold choices were made."), file.path(OUT,"03_activity_primary/ACTIVITY_MODEL_QA.md"))
cat(sprintf("Activity models complete: N=%d events=%d reduced OR=%.6f P=%.6g full OR=%.6f CI=[%.6f, %.6f] P=%.6g estimator=%s\n", nrow(d),sum(d$hit),reduced$OR,reduced$P_raw,full$OR,full$CI_low_OR,full$CI_high_OR,full$P_raw,estimator))
