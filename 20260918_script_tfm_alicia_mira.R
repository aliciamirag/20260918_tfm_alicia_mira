# TFM — CONSTRUCTING SUBJECTIVITIES AROUND REPRODUCTIVE LABOUR IN SPAIN
# (ESGE 2024 – ISSP Family and Changing Gender Roles)

# Alicia Mira Guirao
# 100453622


# PURPOSE OF THIS SCRIPT: This single file reproduces the whole empirical analysis of the thesis, from
# the raw CIS microdata to every table and figure reported in the text and in
# the appendices. I have merged what were originally seven working scripts into
# one so that the analysis can be executed end to end in a single run and
# audited without having to reconstruct the order in which the pieces were
# written.

# The script is organised in nine blocks, which must be run in order because
# each one depends on the objects created by the previous ones:
#
#   BLOCK 1  Data preparation: loading, recoding and construction of the
#            analysis file and of the two long-format files used in H1.
#   BLOCK 2  Measurement: ordinal CFA of the gender-attitude battery,
#            measurement invariance by sex and construction of the scales.
#   BLOCK 3  Multiple imputation of household income and composition.
#   BLOCK 4  Descriptive analysis, including the exploratory descriptives run
#            as a preliminary step for each hypothesis.
#   BLOCK 5  H1a: measurability of care relative to domestic work.
#   BLOCK 6  H1b: self-reports versus proxy-informant reports.
#   BLOCK 7  H2: family-to-work interference, material position and sex.
#   BLOCK 8  H3a and H3b: naturalising discourse and perceived proportionality.
#   BLOCK 9  Table manifest and session information.


# HYPOTHESES:

#   H1a  Care is harder to quantify in units of time than domestic work.
#   H1b  The hours of reproductive labour attributed to women are lower when
#        reported by a male proxy informant than when women report their own.
#   H2   At equal working hours women perceive greater interference of
#        reproductive labour in paid work, and this worsens with lower
#        equivalised income.
#   H3a  The naturalising discourse operates, in both sexes, as an incremental
#        factor in perceiving one's own contribution as proportionate.
#   H3b  Being male increases that perception independently of that discourse.
#

# INPUT AND OUTPUT

# Input   3506_num.csv — CIS Study 3506 microdata, semicolon-separated, comma
#         decimal mark, UTF-8. The file is expected in the working directory.
# Output  output/tables/<section>/<name>.html   every table, written with
#                                               stargazer for direct pasting
#                                               into the manuscript
#         output/figures/<section>/<name>.png   the two figures of the thesis
#         output/tables/_manifest.html          index of every table produced
#
# I export all tables as HTML rather than as CSV because the final deliverable
# is a Word document: stargazer output can be opened in a browser and pasted
# into the manuscript with its structure intact.
#
# Approximate running time: fifteen to twenty minutes, dominated by the
# multiple imputation (Block 3), the RIF cluster bootstrap (Block 6) and the
# Oaxaca-Blinder bootstrap (Block 7).


# =============================================================================
# BLOCK 0. SETUP
# =============================================================================
# I load every package used anywhere in the pipeline at the outset, so that a
# missing dependency fails immediately rather than twenty minutes into the run.

library(tidyverse)
library(survey)
library(psych)
library(lavaan)
library(mice)
library(survival)   # conditional logit for H1a
library(sandwich)   # heteroskedasticity-robust and clustered covariances
library(lmtest)     # coefficient tests with alternative covariances
library(quantreg)   # median regression, H1b sensitivity
library(MASS)       # multivariate normal draws for the marginal effects
library(stargazer)

# MASS::select masks dplyr::select, so I restore the tidyverse verbs explicitly
# rather than relying on the attachment order.
select <- dplyr::select

# -----------------------------------------------------------------------------
# 0.1 Global parameters
# -----------------------------------------------------------------------------

RAW_PATH   <- "3506_num.csv"
OUT_ROOT   <- "output"
REF_YEAR   <- 2025L        # fieldwork was conducted between March and May 2025

SEED_PREP   <- 20250912L   # preparation, measurement and imputation
SEED_MODELS <- 20250913L   # resampling-based inference in the model blocks

N_IMPUTATIONS   <- 20L     # the fraction of missing information is near a third
N_ITERATIONS    <- 10L     # mice converges well before the tenth iteration here
N_BOOT_RIF      <- 500L    # cluster-bootstrap replicates for the RIF decomposition
N_BOOT_OAXACA   <- 500L    # bootstrap replicates per imputation for Oaxaca-Blinder
N_BOOT_WEIGHTED <- 5000L   # bootstrap replicates for the weighted odds ratio in H1a

# The survey is stratified with sex-by-age quotas and no primary sampling unit
# identifier is released, so single-observation strata have to be handled
# explicitly rather than left to fail.
options(survey.lonely.psu = "adjust")

SECTIONS <- c("descriptive", "measurement", "imputation",
              "H1a", "H1b", "H2", "H3a", "H3b")

for (section in SECTIONS) {
  dir.create(file.path(OUT_ROOT, "tables",  section),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT_ROOT, "figures", section),
             recursive = TRUE, showWarnings = FALSE)
}

MANIFEST_PATH <- file.path(OUT_ROOT, "tables", "_manifest.html")

# Every table produced by the analysis is registered in this manifest, so that
# the appendix can be assembled from a complete inventory rather than from
# memory. The manifest is itself written out as HTML at the end of Block 9.
TABLE_MANIFEST <- tibble(section = character(), table = character(),
                         type = character(), path = character(),
                         caption = character())

# -----------------------------------------------------------------------------
# 0.2 Table and figure exporters
# -----------------------------------------------------------------------------
# Two exporters are used throughout. `save_table()` handles the data frames I
# build myself (descriptive summaries, pooled estimates, decompositions), and
# `save_model_table()` handles fitted model objects, for which stargazer
# produces its standard regression layout. Both write HTML into
# output/tables/<section>/ and register the result in the manifest.

register_table <- function(section, name, type, path, caption) {
  TABLE_MANIFEST <<- TABLE_MANIFEST %>%
    filter(!(section == !!section & table == !!name)) %>%
    add_row(section = section, table = name, type = type,
            path = path, caption = caption)
  message(sprintf("  saved  %-12s %-38s [%s]", section, name, type))
  invisible(path)
}

# For data-frame tables I write the HTML directly rather than through
# stargazer. stargazer's exporter parses every cell to decide how to align and
# round it, and that parser fails on some of the tables produced here with
# "missing value where TRUE/FALSE needed" even when the data contain no missing
# values at all. Writing the markup myself removes an unpredictable dependency
# from the step that produces the deliverable. The output keeps stargazer's
# appearance -- centred table, rules above and below the header, caption on top
# -- so the tables paste into the manuscript exactly as before. stargazer is
# still used for fitted models, where its regression layout is the point.

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

write_html_table <- function(df, caption, path) {
  n_col <- ncol(df)
  rule  <- sprintf('<tr><td colspan="%d" style="border-bottom: 1px solid black"></td></tr>',
                   n_col)

  cell_style <- 'style="text-align:left; padding: 3px 10px"'

  header <- paste0(
    "<tr>",
    paste0("<td ", cell_style, "><strong>", html_escape(names(df)),
           "</strong></td>", collapse = ""),
    "</tr>")

  body <- vapply(seq_len(nrow(df)), function(i) {
    cells <- html_escape(as.character(unlist(df[i, ], use.names = FALSE)))
    paste0("<tr>",
           paste0("<td ", cell_style, ">", cells, "</td>", collapse = ""),
           "</tr>")
  }, character(1))

  writeLines(c('<table style="text-align:center">',
               paste0("<caption><strong>", html_escape(caption),
                      "</strong></caption>"),
               rule, header, rule, body, rule,
               "</table>"),
             con = path, useBytes = TRUE)
  invisible(path)
}

save_table <- function(df, section, name, caption, digits = 3) {
  stopifnot(section %in% SECTIONS, is.data.frame(df), nrow(df) > 0)
  path <- file.path(OUT_ROOT, "tables", section, paste0(name, ".html"))

  # Every column is rendered to character before export, with missing values
  # written as an empty cell, which is how they should print in the manuscript
  # in any case. Whole numbers keep their integer appearance and the rest are
  # fixed to `digits` decimals, so that a column of counts is not printed as
  # 923.000.
  fmt_column <- function(x, column_name = "") {
    # A p-value rounded to a fixed number of decimals prints as 0 when it is
    # small, which reads as a probability of zero rather than as a very small
    # one. Any numeric column named `p` is therefore formatted as "<.001"
    # below that threshold.
    if (is.numeric(x) && column_name == "p") return(fmt_p(x))
    if (is.numeric(x)) {
      whole <- all(is.na(x) | x == round(x))
      out <- if (whole) {
        format(x, trim = TRUE, scientific = FALSE)
      } else {
        formatC(x, format = "f", digits = digits)
      }
    } else {
      out <- as.character(x)
    }
    ifelse(is.na(x), "", out)
  }

  out <- df %>%
    ungroup() %>%
    as.data.frame()
  out[] <- Map(fmt_column, out, names(out))

  new_names <- ifelse(names(out) %in% names(COLUMN_LABELS),
                      unname(COLUMN_LABELS[names(out)]), names(out))
  names(out) <- make.unique(new_names, sep = " ")

  write_html_table(out, caption, path)
  register_table(section, name, "data frame", path, caption)
}

save_model_table <- function(models, section, name, caption, ...) {
  stopifnot(section %in% SECTIONS)

  # A fitted model is itself a list, so `is.list()` cannot distinguish one
  # model from a list of models: testing it that way passes the model's own
  # components to stargazer as though each were a separate model. A plain list
  # has class "list" and a fitted object does not, which is the test used here.
  if (!identical(class(models), "list")) models <- list(models)

  path <- file.path(OUT_ROOT, "tables", section, paste0(name, ".html"))

  # stargazer covers the model classes used here (clogit and svyglm). Should a
  # class not be recognised in a future version of the package, the run must
  # not stop: the fallback writes the same information as a coefficient table.
  ok <- tryCatch({
    invisible(utils::capture.output(
      do.call(stargazer::stargazer,
              c(unname(models),
                list(type = "html", title = caption, out = path, ...)))
    ))
    TRUE
  }, error = function(e) {
    warning("stargazer could not render '", name, "' as a model table (",
            conditionMessage(e), "); falling back to a coefficient table.",
            call. = FALSE)
    FALSE
  })

  if (!ok) {
    labels <- if (is.null(names(models))) {
      paste("Model", seq_along(models))
    } else {
      names(models)
    }
    coefs <- map2_dfr(models, labels, function(m, label) {
      s <- summary(m)$coefficients
      tibble(model = label, term = rownames(s),
             estimate = s[, 1], se = s[, 2], p = s[, ncol(s)])
    })
    return(save_table(coefs, section, name, caption))
  }
  register_table(section, name, "model", path, caption)
}

save_figure <- function(plot_expr, section, name, caption,
                        width = 1000, height = 620, res = 110) {
  stopifnot(section %in% SECTIONS)
  path <- file.path(OUT_ROOT, "figures", section, paste0(name, ".png"))
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_expr)
  message(sprintf("  figure %-12s %-38s", section, name))
  invisible(path)
}

# -----------------------------------------------------------------------------
# 0.2b Presentation layer
# -----------------------------------------------------------------------------
# The analysis objects use the variable names of the dataset, which are not the
# names the tables carry in the manuscript. Rather than renaming by hand at
# every export, I keep a single lookup for column headers and another for model
# terms, plus a small set of formatters that assemble the composite cells the
# thesis tables use ("0.310 [0.160, 0.460]", "-0.690 (0.199)"). Everything
# downstream of this section works on the raw names; the translation happens
# only at the point of export, so no filter or join depends on it.

COLUMN_LABELS <- c(
  task = "Task", sex = "Sex", item = "Item", term = "Term", model = "Model",
  n = "N", se = "S.E.", p = "p", estimate = "Estimate", fmi = "FMI",
  ci_low = "95% CI lower", ci_high = "95% CI upper",
  odds_ratio = "Odds ratio", or_low = "OR 95% lower", or_high = "OR 95% upper",
  positive_hours = "Positive hours", zero_hours = "Zero hours",
  cant_say = "Can't say", refused = "Refused",
  mean = "Mean", median = "Median", max = "Max", p25 = "P25", p75 = "P75",
  p90 = "P90", variable = "Variable", level = "Level", block = "Block",
  unweighted_pct = "Unweighted %", weighted_pct = "Weighted %",
  pct_missing = "% missing", n_valid = "N valid", statistic = "Statistic",
  value = "Value", caption = "Caption", section = "Section", type = "Type",
  path = "File", check = "Check", observed = "Observed", expected = "Expected",
  pass = "Pass", target = "Target person", source = "Source", quantile = "Quantile",
  cell = "Cell", composition = "Composition", composition_se = "Composition SE",
  composition_p = "Composition p", structure = "Structure",
  structure_se = "Structure SE", structure_p = "Structure p",
  outcome = "Outcome", specification = "Specification",
  direction = "Direction", scenario = "Scenario", measure = "Measure",
  panel = "Panel", gap = "Gap (hours)", capped = "Capped p95",
  untruncated = "Untruncated", log_scale = "Log scale",
  median_reg = "Median regression", discordant = "Discordant pairs",
  mean_raw = "Mean (raw)", mean_winsorised = "Mean (winsorised)",
  cap_p95 = "Cap at p95", max_raw = "Max (raw)", pct_capped = "% capped",
  conflict_fw = "Family-to-work conflict", conflict_wf = "Work-to-family conflict",
  work_hours = "Paid-work hours", econ_strain = "Ease of making ends meet",
  ame_pp = "AME (percentage points)", se_pp = "SE",
  peak_task_load = "Peak declared load", peak_probability = "Peak probability",
  range_lo = "Observed range, lower", range_hi = "Observed range, upper",
  at_edge = "Peak at edge of range",
  n_persons = "N (persons)", n_pairs = "N (pairs)", n_conseq = "N (CONSEQ)",
  n_roles = "N (ROLES)", n_valid_roster = "N (roster valid)",
  conseq = "CONSEQ", roles = "ROLES", load_band = "Declared load band",
  item = "Item", group = "Group", set = "Set", delta = "Delta",
  income_quartile = "Income quartile", component = "Component",
  men_several_per_week = "Men, several times a week",
  women_several_per_week = "Women, several times a week",
  men_mostly_self = "Men, mostly themselves",
  women_mostly_self = "Women, mostly themselves",
  pct_someone_else = "Done by a third party (%)",
  does_much_more = "Does much more than their share",
  does_more = "Does more", fair_share = "Does roughly their share",
  does_less = "Does less", does_much_less = "Does much less",
  mean_load_when_fair = "Mean load when judged fair",
  median_load_when_fair = "Median load when judged fair",
  pct_fair = "Perceived as proportionate (%)")

# Model terms as they should read in the manuscript. Unknown terms fall through
# to a set of regular expressions and, failing those, are left untouched: a term
# I have forgotten must appear as it is rather than silently disappear.
TERM_LABELS <- c(
  "female"                 = "Woman",
  "work_hours"             = "Paid-work hours",
  "age"                    = "Age",
  "has_minor"              = "Minor in household",
  "has_elderly"            = "Person aged 75+ in household",
  "partnered"              = "Cohabiting partner",
  "econ_strain"            = "Ease of making ends meet",
  "modepaper"              = "Survey mode: paper",
  "conseq_index"           = "Discourse: consequences of maternal employment (CONSEQ)",
  "roles_index"            = "Discourse: normative prescription of roles (ROLES)",
  "task_load"              = "Declared task load",
  "I(task_load^2)"         = "Declared task load squared",
  "mental_load"            = "Mental load (V40)",
  "outsourced"             = "Task outsourced to a third party",
  "female:conseq_index"    = "Woman x CONSEQ",
  "female:roles_index"     = "Woman x ROLES",
  "female:task_load"       = "Woman x declared task load",
  "female:I(task_load^2)"  = "Woman x declared task load squared",
  "taskcare"               = "Care work (ref: domestic work)",
  "sourceproxy"            = "Proxy report (ref: self-report)")

relabel_terms <- function(x) {
  out <- unname(TERM_LABELS[x])
  out[is.na(out)] <- x[is.na(out)]
  # Regular patterns that would otherwise need one entry per level.
  out <- sub("^female:income_q(Q[1-4])$", "Woman x income \\1", out)
  out <- sub("^income_q(Q[1-4]):female$", "Woman x income \\1", out)
  out <- sub("^income_q(Q[1-4])$",        "Income quartile \\1", out)
  out <- sub("^educ4(\\d)$",              "Education level \\1", out)
  out
}

# Formatters for the composite cells used in the manuscript tables.
fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<.001",
                sub("^0", "", formatC(p, format = "f", digits = digits))))
}

fmt_est_ci <- function(est, lo, hi, digits = 3) {
  paste0(fmt_num(est, digits), " [", fmt_num(lo, digits), ", ",
         fmt_num(hi, digits), "]")
}

fmt_ci_dash <- function(lo, hi, digits = 2) {
  paste0(fmt_num(lo, digits), " - ", fmt_num(hi, digits))
}

fmt_b_se <- function(b, se, digits = 3) {
  paste0(fmt_num(b, digits), " (", fmt_num(se, digits), ")")
}

fmt_or_ci <- function(or, lo, hi, digits = 2) {
  paste0(fmt_num(or, digits), " [", fmt_num(lo, digits), ", ",
         fmt_num(hi, digits), "]")
}

# -----------------------------------------------------------------------------
# 0.3 Recoding helpers
# -----------------------------------------------------------------------------
# The CIS file uses distinct special codes for non-applicability, inability to
# answer and refusal. Collapsing them into a single missing category would
# destroy the information on which H1a rests, so each family of items has its
# own recoding function and the codes are handled explicitly.

# Set the listed codes to NA.
na_if_in <- function(x, codes) replace(x, x %in% codes, NA)

# Hours items (V34-V37). Three distinct outcomes must be kept apart:
#   996 not applicable | 998 "couldn't say" | 999 refused | 0-500 a figure.
# Code 998 is not missing data: it is the dependent variable of H1a.
hours_status <- function(x) {
  factor(case_when(x == 996         ~ "nap",
                   x == 998         ~ "cant_say",
                   x == 999         ~ "refused",
                   x >= 0 & x < 500 ~ "figure"),
         levels = c("figure", "cant_say", "refused", "nap"))
}

hours_value <- function(x) ifelse(x >= 0 & x < 500, as.numeric(x), NA_real_)

# Attitude items (V1-V6): 1 strongly agree ... 5 strongly disagree; 8 and 9 are
# missing. I recode them so that HIGH means stronger adherence to the
# naturalising discourse. V1 ("a working mother can establish just as warm a
# relationship") is the reverse-keyed item of the ISSP scale and therefore
# keeps its raw direction.
trad_recode <- function(x, reverse_keyed = FALSE) {
  x <- na_if_in(x, c(0, 8, 9))
  x <- ifelse(x >= 1 & x <= 5, x, NA_real_)
  if (reverse_keyed) x else 6 - x
}

# Task-division items (V39-V44): 1 always me ... 5 always partner, 6 "someone
# else", 0 not applicable. Recoded so that HIGH means the respondent does more.
# Code 6 is not a point on the scale and is therefore not rescaled.
task_recode <- function(x) {
  x <- ifelse(x >= 1 & x <= 5, x, NA_real_)
  6 - x
}

# Conflict items (V46-V49): 1 several times a week ... 4 never; 8 not
# applicable or not in work, 9 refused. Recoded so that HIGH means more
# frequent conflict.
conflict_recode <- function(x) {
  x <- ifelse(x >= 1 & x <= 4, x, NA_real_)
  5 - x
}

# Education in three levels. ESTUDIOS and SPESTUDIOS share a coding scheme,
# which is what makes the symmetric respondent/partner comparison of H1b
# possible. NIVELEDU is the CIS-derived four-level variable and has no partner
# counterpart.
recode_educ3 <- function(x) {
  x <- ifelse(x >= 1 & x <= 7, x, NA_real_)
  factor(case_when(x <= 2 ~ "low", x <= 5 ~ "mid", x <= 7 ~ "high"),
         levels = c("low", "mid", "high"))
}

# Household and personal income bands. Code 97 ("no income") is a substantive
# zero and is assigned to the lowest band; codes 0, 77 and 99 are missing.
income_band <- function(x) {
  case_when(x >= 1 & x <= 10 ~ as.numeric(x),
            x == 97          ~ 1,
            TRUE             ~ NA_real_)
}

# Row mean requiring a minimum number of available items. With a two-item scale
# the minimum has to be two: a single item is not a scale.
mean_min <- function(df, min_items) {
  ok  <- rowSums(!is.na(df)) >= min_items
  out <- rowMeans(df, na.rm = TRUE)
  ifelse(ok, out, NA_real_)
}

# Weighted quantiles. svyby() combined with svyquantile() no longer returns
# uncertainty estimates in current versions of the survey package, and quantile
# standard errors are not needed here: the quantiles are descriptive and the
# design-based standard error is reported for the mean, which is the statistic
# being compared.
wquantile <- function(x, w, probs) {
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which.max(cw >= p)], numeric(1))
}

pct    <- function(x) round(100 * x, 1)
sexlab <- function(f) factor(f, levels = 0:1, labels = c("Men", "Women"))

tidy_ci <- function(est, se, digits = 3) {
  tibble(estimate = round(est, digits),
         se       = round(se, digits),
         ci_low   = round(est - 1.96 * se, digits),
         ci_high  = round(est + 1.96 * se, digits),
         p        = round(2 * pnorm(-abs(est / se)), 4))
}

# -----------------------------------------------------------------------------
# 0.4 Reproducibility checks
# -----------------------------------------------------------------------------
# Each construction step is checked against the figures established during the
# exploratory phase. The checks are deliberately non-fatal: a failed check
# signals that the input file or a package default has changed and must be
# investigated, but it should not prevent the rest of the pipeline from running
# and reporting where exactly the divergence appears.

CHECKS <- list()

check <- function(label, observed, expected, tol = 0) {
  pass <- abs(observed - expected) <= tol
  CHECKS[[length(CHECKS) + 1]] <<- tibble(check = label, observed = observed,
                                          expected = expected, pass = pass)
  if (!pass) {
    warning(sprintf("CHECK FAILED — %s: observed %s, expected %s",
                    label, observed, expected), call. = FALSE)
  }
  invisible(pass)
}


# =============================================================================
# BLOCK 1. DATA PREPARATION
# =============================================================================
# In this block I build the respondent-level analysis file and the two
# long-format files required by H1. Two decisions are deliberately deferred to
# later blocks and are documented here so that the sequence is transparent.
#
# ON THE ATTITUDE SCALES. No attitude index is constructed here. This block
# only recodes the six items; their dimensionality is not a preparation
# decision but an empirical question, and it is tested in Block 2. Building an
# index at this point would presuppose the very structure that Block 2 is
# meant to establish, and the one-factor solution is in fact rejected.
#
# ON INCOME. No income-derived variable is constructed here either. The bands
# and the household composition are prepared as inputs; the equivalence scale,
# household income in euros, equivalised income and its quartiles are built
# passively inside each imputed dataset in Block 3 through derive_income(),
# which is defined at the end of this block.

message("\n=== BLOCK 1. DATA PREPARATION ===")
set.seed(SEED_PREP)

# -----------------------------------------------------------------------------
# 1.1 Loading the microdata
# -----------------------------------------------------------------------------

raw <- read_delim(RAW_PATH, delim = ";",
                  locale = locale(encoding = "UTF-8", decimal_mark = ","),
                  show_col_types = FALSE, na = character())

check("raw rows", nrow(raw), 1722)

# The three weights are distributed as comma-decimal strings and have to be
# converted before they can be used in any design object.
raw <- raw %>%
  mutate(across(c(PESO, PESODIS, PESOFIN),
                ~ as.numeric(str_replace(as.character(.x), ",", "."))))

# -----------------------------------------------------------------------------
# 1.2 Base filter and core identifiers
# -----------------------------------------------------------------------------
# The eight respondents who selected "prefer not to answer" for sex (SEXO = 7)
# are excluded. Sex is constitutive of every hypothesis and cannot be imputed
# here without determining the very quantity under analysis.
#
# PESOFIN is the final post-stratification weight and is the one used in all
# design-based estimates; PESO is the stratification weight and PESODIS the
# design weight, and neither is used in the analysis.

d <- raw %>%
  filter(SEXO %in% c(1, 2)) %>%
  mutate(id     = row_number(),
         female = as.integer(SEXO == 2),
         w      = PESOFIN,
         mode   = factor(MODE, levels = c(41, 34), labels = c("cawi", "paper")),
         age    = ifelse(EDAD >= 18 & EDAD < 120, EDAD, NA_real_))

check("working sample", nrow(d), 1714)
check("cohabiting partner (PARTLIV == 1)", sum(d$PARTLIV == 1, na.rm = TRUE), 1097)

# -----------------------------------------------------------------------------
# 1.3 Household roster
# -----------------------------------------------------------------------------
# The roster covers the OTHER members of the household; the respondent is not
# listed among them. Three parallel blocks of ten slots record the relationship
# to the respondent (HMP_SB_i, where 1 identifies the spouse or partner), the
# sex of the member (HMP_GD_i) and the year of birth (HMP_YR_i). I reshape them
# into a long file, clean the special codes and summarise each household back
# to one row per respondent.

slots <- 1:10
sb <- paste0("HMP_SB_", slots)
gd <- paste0("HMP_GD_", slots)
yr <- paste0("HMP_YR_", slots)

roster <- d %>%
  select(id, HOMPOP_TOT, all_of(c(sb, gd, yr))) %>%
  pivot_longer(-c(id, HOMPOP_TOT),
               names_to = c(".value", "slot"),
               names_pattern = "HMP_(SB|GD|YR)_(\\d+)") %>%
  mutate(slot       = as.integer(slot),
         rel        = na_if_in(SB, c(0, 9, 99)),
         sex        = na_if_in(GD, c(0, 9, 99)),
         birth_yr   = ifelse(YR > 1900 & YR < REF_YEAR + 1, YR, NA_real_),
         member_age = REF_YEAR - birth_yr)

roster_summary <- roster %>%
  group_by(id) %>%
  summarise(
    n_roster_valid = sum(!is.na(birth_yr)),
    n_minors       = sum(member_age <  18, na.rm = TRUE),
    n_under14      = sum(member_age <  14, na.rm = TRUE),
    n_adults14     = sum(member_age >= 14, na.rm = TRUE),
    n_elderly75    = sum(member_age >= 75, na.rm = TRUE),
    partner_sex    = first(sex[rel == 1 & !is.na(rel)], default = NA_real_),
    partner_age    = first(member_age[rel == 1 & !is.na(rel)], default = NA_real_),
    .groups = "drop")

d <- d %>%
  left_join(roster_summary, by = "id") %>%
  mutate(
    hh_size = na_if_in(HOMPOP_TOT, c(77, 99)),
    # TRUE when every co-resident has a usable year of birth, that is, when the
    # OECD equivalence scale can be built without recourse to imputation.
    roster_complete = !is.na(hh_size) &
                      (hh_size == 1 | n_roster_valid == hh_size - 1),
    has_minor   = as.integer(n_minors    > 0),
    has_elderly = as.integer(n_elderly75 > 0))

check("households with a minor",      sum(d$has_minor,   na.rm = TRUE), 366)
check("households with a person 75+", sum(d$has_elderly, na.rm = TRUE), 131)
check("hh_size missing",              sum(is.na(d$hh_size)),            101)

# -----------------------------------------------------------------------------
# 1.4 Couple structure
# -----------------------------------------------------------------------------
# H1b requires a partner who can be identified in the roster and whose sex
# differs from the respondent's, since the contrast compares self-reports with
# reports made by a partner of the other sex.

d <- d %>%
  mutate(
    partnered       = as.integer(PARTLIV == 1),
    partner_known   = as.integer(partnered == 1 & partner_sex %in% c(1, 2)),
    same_sex_couple = as.integer(partner_known == 1 & partner_sex == SEXO),
    diff_sex_couple = as.integer(partner_known == 1 & partner_sex != SEXO))

n_multi_partner <- roster %>%
  group_by(id) %>%
  summarise(k = sum(rel == 1, na.rm = TRUE)) %>%
  filter(k > 1) %>%
  nrow()

check("rosters listing more than one partner", n_multi_partner, 0)
check("partner identifiable in roster", sum(d$partner_known,   na.rm = TRUE), 974)
check("different-sex couples",          sum(d$diff_sex_couple, na.rm = TRUE), 952)
check("same-sex couples",               sum(d$same_sex_couple, na.rm = TRUE),  22)

# -----------------------------------------------------------------------------
# 1.5 Education, employment and material position
# -----------------------------------------------------------------------------

d <- d %>%
  mutate(
    educ3           = recode_educ3(ESTUDIOS),
    partner_educ3   = recode_educ3(SPESTUDIOS),
    educ4           = factor(na_if_in(NIVELEDU, c(0, 8, 9))),
    in_work         = as.integer(WORK == 1),
    partner_in_work = as.integer(SPWORK == 1),
    work_hours      = ifelse(WRKHRS >= 0 & WRKHRS <= 120, WRKHRS, NA_real_),
    econ_strain     = na_if_in(DIFICULT_ECO, c(0, 8, 9)),
    subj_status     = na_if_in(TOPBOT, c(98, 99)))

# -----------------------------------------------------------------------------
# 1.6 Hours of reproductive labour (V34-V37)
# -----------------------------------------------------------------------------
# The status variables preserve the distinction between an inability to
# quantify and a refusal, which is what allows the response process itself to
# be treated as an outcome. The dependent variables of H1a are built from that
# distinction: refusals are excluded rather than coded as zero, because a
# refusal is not a declaration of incommensurability.

d <- d %>%
  mutate(
    st_hours_dom_self     = hours_status(V34),
    st_hours_care_self    = hours_status(V35),
    st_hours_dom_partner  = hours_status(V36),
    st_hours_care_partner = hours_status(V37),
    hours_dom_self        = hours_value(V34),
    hours_care_self       = hours_value(V35),
    hours_dom_partner     = hours_value(V36),
    hours_care_partner    = hours_value(V37),
    noquant_dom  = case_when(st_hours_dom_self  == "cant_say" ~ 1L,
                             st_hours_dom_self  == "figure"   ~ 0L),
    noquant_care = case_when(st_hours_care_self == "cant_say" ~ 1L,
                             st_hours_care_self == "figure"   ~ 0L))

check("valid figures, own domestic hours", sum(!is.na(d$hours_dom_self)),  1344)
check("valid figures, own care hours",     sum(!is.na(d$hours_care_self)), 1149)
check("cant_say, own domestic hours", sum(d$st_hours_dom_self  == "cant_say"), 331)
check("cant_say, own care hours",     sum(d$st_hours_care_self == "cant_say"), 414)

# -----------------------------------------------------------------------------
# 1.7 Gender-role attitude items (recoding only)
# -----------------------------------------------------------------------------
# att_warm      V1  a working mother can have as warm a relationship (reverse-keyed)
# att_suffer    V2  a pre-school child suffers if the mother works
# att_famlife   V3  family life suffers when the woman works full time
# att_wantshome V4  what women really want is a home and children
# att_fulfil    V5  being a housewife is as fulfilling as paid work
# att_breadwin  V6  a man's job is to earn money, a woman's to look after the home
#
# All six items are carried forward, att_fulfil included: whether it belongs
# with the others is a question for the CFA in Block 2, not an assumption to be
# hard-coded here.

d <- d %>%
  mutate(
    att_warm      = trad_recode(V1, reverse_keyed = TRUE),
    att_suffer    = trad_recode(V2),
    att_famlife   = trad_recode(V3),
    att_wantshome = trad_recode(V4),
    att_fulfil    = trad_recode(V5),
    att_breadwin  = trad_recode(V6))

ATTITUDE_ITEMS <- c("att_warm", "att_suffer", "att_famlife",
                    "att_wantshome", "att_fulfil", "att_breadwin")
CORE_ITEMS     <- setdiff(ATTITUDE_ITEMS, "att_fulfil")

# Ordered-factor copies for the WLSMV estimator used in Block 2.
items_attitudes <- d %>%
  select(id, w, female, educ4, age, all_of(ATTITUDE_ITEMS)) %>%
  mutate(across(all_of(ATTITUDE_ITEMS), ~ ordered(.x, levels = 1:5)))

check("respondents with at least one attitude item",
      sum(rowSums(!is.na(d[, ATTITUDE_ITEMS])) > 0), 1700)
check("complete on the five core items",
      sum(complete.cases(d[, CORE_ITEMS])), 1518)

# -----------------------------------------------------------------------------
# 1.8 Division of household tasks (V39-V44)
# -----------------------------------------------------------------------------
# V39 laundry | V40 planning and organising family activities |
# V41 care of sick family members | V42 grocery shopping |
# V43 cleaning | V44 preparing meals
#
# The declared task load averages the four routine material tasks, matching the
# referent of V45 ("the household tasks"). V40 is kept apart as an indicator of
# mental load, because it captures cognitive and organisational responsibility
# rather than execution; V41 is excluded because it measures care, not
# housework.

d <- d %>%
  mutate(
    task_laundry   = task_recode(V39),
    task_planning  = task_recode(V40),
    task_sickcare  = task_recode(V41),
    task_groceries = task_recode(V42),
    task_cleaning  = task_recode(V43),
    task_meals     = task_recode(V44))

TASK_ITEMS <- c("task_laundry", "task_groceries", "task_cleaning", "task_meals")

d <- d %>%
  mutate(
    task_load   = mean_min(across(all_of(TASK_ITEMS)), min_items = 3),
    mental_load = task_planning,
    # Code 6 is not a point on the scale: it flags tasks performed by a third
    # party, which is the outsourcing indicator used as a control in H3.
    outsourced  = as.integer(rowSums(across(c(V39, V40, V41, V42, V43, V44),
                                            ~ .x == 6), na.rm = TRUE) > 0))

alpha_task <- psych::alpha(as.data.frame(d[, TASK_ITEMS]),
                           check.keys = FALSE, warnings = FALSE)
check("alpha, task load index", round(alpha_task$total$raw_alpha, 3),
      0.759, tol = 0.01)

# -----------------------------------------------------------------------------
# 1.9 Perceived proportionality of the split (V45)
# -----------------------------------------------------------------------------
# The scale runs from 1 ("much more than my fair share") through 3 ("roughly my
# fair share") to 5 ("much less"). The quantity described by H3a and H3b is the
# midpoint, so the dependent variable is binary rather than ordinal: an ordinal
# model would impose a monotone ordering on a category that is substantively
# the centre, not an extreme.

d <- d %>%
  mutate(
    fair_raw   = ifelse(V45 >= 1 & V45 <= 5, V45, NA_real_),
    fair_share = as.integer(fair_raw == 3),
    fair_dir   = factor(case_when(fair_raw %in% c(1, 2) ~ "does_more",
                                  fair_raw == 3         ~ "fair",
                                  fair_raw %in% c(4, 5) ~ "does_less"),
                        levels = c("does_more", "fair", "does_less")))

check("valid V45 among partnered",
      sum(d$partnered == 1 & !is.na(d$fair_raw), na.rm = TRUE), 1030)

# -----------------------------------------------------------------------------
# 1.10 Work-family interference (V46-V49)
# -----------------------------------------------------------------------------
# V46 and V47 capture the work-to-family direction; V48 and V49 the
# family-to-work direction, which is the one that operationalises the double
# presence and therefore the dependent variable of H2. Code 8 means "not
# applicable / not in work" and is not a zero.

d <- d %>%
  mutate(
    conf_tired_home  = conflict_recode(V46),
    conf_duties_home = conflict_recode(V47),
    conf_tired_work  = conflict_recode(V48),
    conf_concentrate = conflict_recode(V49),
    conflict_wf = mean_min(across(c(conf_tired_home, conf_duties_home)), 2),
    conflict_fw = mean_min(across(c(conf_tired_work, conf_concentrate)), 2))

alpha_fw <- psych::alpha(as.data.frame(d[, c("conf_tired_work", "conf_concentrate")]),
                         check.keys = FALSE, warnings = FALSE)
check("alpha, family-to-work conflict", round(alpha_fw$total$raw_alpha, 3),
      0.719, tol = 0.01)

check("H2 analytic sample",
      sum(d$in_work == 1 & !is.na(d$conflict_fw) & !is.na(d$work_hours),
          na.rm = TRUE), 888)

# -----------------------------------------------------------------------------
# 1.11 Income inputs
# -----------------------------------------------------------------------------
# NAT_INC records household income and NAT_RINC personal income, both in ten
# brackets, with 0 for non-applicability (single-person households skip the
# household question), 77 for "prefer not to answer", 97 for "no income" and 99
# for refusal. For single-person households the personal figure substitutes for
# the household one.

d <- d %>%
  mutate(
    hh_income_band  = ifelse(hh_size == 1, income_band(NAT_RINC),
                                           income_band(NAT_INC)),
    own_income_band = income_band(NAT_RINC),
    hh_income_src   = if_else(hh_size == 1, NAT_RINC, NAT_INC),
    no_income       = as.integer(hh_income_src == 97),
    # An explicit refusal (code 77) is informative behaviour, not inattention.
    # I flag it here so that the delta sensitivity analysis in Block 3 can be
    # documented against the composition of the non-response.
    income_refused  = as.integer(hh_income_src == 77))

check("household income band available", sum(!is.na(d$hh_income_band)), 1250,
      tol = 15)

# -----------------------------------------------------------------------------
# 1.12 Survey design
# -----------------------------------------------------------------------------
# The sample is stratified by autonomous community and municipality size with
# sex-by-age quotas, and administered in mixed mode. No primary sampling unit
# identifier is released, so no clusters can be declared and the design is
# specified with strata and weights only.

d <- d %>% mutate(stratum = interaction(CCAA, TAMUNI, drop = TRUE))

# From this point onwards the respondent-level analysis file is called `base`.
base <- d

# -----------------------------------------------------------------------------
# 1.13 Long file for H1a: person by task
# -----------------------------------------------------------------------------
# Each respondent contributes two rows and constitutes their own stratum, so
# every time-invariant individual characteristic is absorbed by construction.
# Only complete within-person pairs are retained, since an incomplete pair
# contributes nothing to a conditional likelihood.

h1a_long <- base %>%
  select(id, w, female, age, has_minor, has_elderly, educ4,
         noquant_dom, noquant_care) %>%
  pivot_longer(c(noquant_dom, noquant_care),
               names_to = "task", values_to = "cant_say") %>%
  mutate(task = factor(if_else(task == "noquant_dom", "domestic", "care"),
                       levels = c("domestic", "care"))) %>%
  drop_na(cant_say) %>%
  group_by(id) %>%
  filter(n() == 2) %>%
  ungroup()

check("H1a respondents", n_distinct(h1a_long$id), 1544)

discordant <- h1a_long %>%
  select(id, task, cant_say) %>%
  pivot_wider(names_from = task, values_from = cant_say) %>%
  summarise(dom_only  = sum(domestic == 1 & care == 0),
            care_only = sum(domestic == 0 & care == 1))

check("discordant: care only",     discordant$care_only, 190)
check("discordant: domestic only", discordant$dom_only,   57)

# -----------------------------------------------------------------------------
# 1.14 Long file for H1b: report level
# -----------------------------------------------------------------------------
# One row per combination of informant, target person and task. Each respondent
# supplies up to four rows, which is why the standard errors of the triple
# interaction are clustered on the respondent in Block 6. The covariates
# describe the TARGET of the report, not the informant, since the comparison is
# between two descriptions of comparable people rather than between two
# descriptions of the same person: only one member of each household is
# interviewed.

base_cols <- base %>%
  filter(diff_sex_couple == 1) %>%
  select(id, w, has_minor, mode,
         SEXO, age, educ3, in_work,
         partner_sex, partner_age, partner_educ3, partner_in_work,
         hours_dom_self, hours_care_self, hours_dom_partner, hours_care_partner)

self_rows <- base_cols %>%
  select(id, w, has_minor, mode, target_sex = SEXO, target_age = age,
         target_educ = educ3, target_works = in_work,
         domestic = hours_dom_self, care = hours_care_self) %>%
  mutate(source = "self")

proxy_rows <- base_cols %>%
  select(id, w, has_minor, mode, target_sex = partner_sex,
         target_age = partner_age, target_educ = partner_educ3,
         target_works = partner_in_work,
         domestic = hours_dom_partner, care = hours_care_partner) %>%
  mutate(source = "proxy")

h1b_long <- bind_rows(self_rows, proxy_rows) %>%
  pivot_longer(c(domestic, care), names_to = "task", values_to = "hours") %>%
  drop_na(hours) %>%
  mutate(source        = factor(source, levels = c("self", "proxy")),
         task          = factor(task,   levels = c("domestic", "care")),
         target_female = as.integer(target_sex == 2)) %>%
  # Winsorised within task, pooling both sources so that the cap is identical
  # across the arms being compared and cannot itself generate a gap.
  group_by(task) %>%
  mutate(hours_w   = pmin(hours, quantile(hours, 0.95, na.rm = TRUE)),
         log_hours = log1p(hours_w)) %>%
  ungroup()

check("H1b rows (all reports)", nrow(h1b_long),         2892)
check("H1b respondents",        n_distinct(h1b_long$id), 842)

# Rows surviving listwise deletion on the target covariates, that is, the
# sample the models of Block 6 are actually estimated on. Reported here for
# transparency rather than enforced at this stage.
h1b_complete <- h1b_long %>%
  drop_na(hours_w, target_age, target_educ, target_works)

check("H1b rows with complete covariates", nrow(h1b_complete), 2728)
check("H1b respondents with complete covariates",
      n_distinct(h1b_complete$id), 821)

# -----------------------------------------------------------------------------
# 1.15 derive_income(): passive imputation of the income chain
# -----------------------------------------------------------------------------
# This function is applied inside each completed dataset in Block 3 and must
# never be run on the observed data alone. Imputing a derived variable, or
# cutting quartiles before imputing, would understate the uncertainty and break
# the Rubin pooling: a quartile is a position in a distribution that does not
# yet exist at the point of imputation.
#
# Band 10 is open-ended (above 5,000 euros). The value of 6,500 euros comes
# from a Pareto fit to the upper tail and is the figure reported in the thesis.
#
# ON FRACTIONAL BANDS. The band is converted into euros by linear interpolation
# between adjacent midpoints rather than by rounding to the nearest integer
# band. For an observed or imputed band, which is always a whole number, the two
# procedures return exactly the same figure, so nothing in the baseline changes.
# They differ only for the fractional bands produced by the delta sensitivity
# analysis of section 3.5, and there the difference matters: R rounds halves to
# the nearest EVEN number, so adding 0.5 to a whole band moves the odd bands up
# and leaves the even ones exactly where they were. Rounding would therefore
# turn a half-band shift applied to every imputed case into a whole-band shift
# applied to whichever half of them happened to sit on an odd band, which is not
# the sensitivity analysis the thesis describes. Interpolation shifts every
# imputed case by half a band, as intended.

derive_income <- function(df, top_band_value = 6500) {
  mids <- c(600, 1025, 1350, 1750, 2150, 2450, 2900, 3600, 4500, 6500)
  mids[10] <- top_band_value

  band_to_euros <- function(b) {
    b  <- pmin(pmax(b, 1), 10)
    lo <- floor(b)
    hi <- ceiling(b)
    mids[lo] + (b - lo) * (mids[hi] - mids[lo])
  }

  dplyr::mutate(
    df,
    # Modified OECD scale, rebuilt from the imputed household composition.
    eq_scale      = 1 + 0.5 * pmax(n_adults14, 0) + 0.3 * pmax(n_under14, 0),
    hh_income_eur = band_to_euros(hh_income_band),
    eq_income     = hh_income_eur / eq_scale,
    income_q      = cut(eq_income,
                        breaks = stats::quantile(eq_income,
                                                 probs = seq(0, 1, 0.25),
                                                 na.rm = TRUE),
                        include.lowest = TRUE,
                        labels = paste0("Q", 1:4)))
}

# -----------------------------------------------------------------------------
# 1.16 Preparation report
# -----------------------------------------------------------------------------

check_report <- bind_rows(CHECKS)

cat("\n--- Block 1 benchmark checks ---\n")
print(as.data.frame(check_report), row.names = FALSE)
cat(sprintf("\n%d of %d checks passed.\n",
            sum(check_report$pass), nrow(check_report)))

cat("\n--- Winsorisation caps (95th percentile within task) ---\n")
h1b_long %>%
  group_by(task) %>%
  summarise(cap_hours  = max(hours_w),
            max_raw    = max(hours),
            pct_capped = round(100 * mean(hours > cap_hours), 1)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n--- Missing data entering the imputation block ---\n")
base %>%
  summarise(hh_income_band  = mean(is.na(hh_income_band)),
            hh_size         = mean(is.na(hh_size)),
            roster_complete = mean(!roster_complete),
            task_load       = mean(is.na(task_load))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "prop_missing") %>%
  mutate(prop_missing = round(prop_missing, 3)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

save_table(check_report, "descriptive", "D0_preparation_checks",
           paste("Reproducibility checks on the construction of the analysis file.",
                 "Each figure is compared with the benchmark established during the",
                 "exploratory phase; a failed check indicates that the input file or a",
                 "package default has changed."))


# =============================================================================
# BLOCK 2. MEASUREMENT OF THE GENDER-ATTITUDE BATTERY
# =============================================================================
# The thesis argues that survey instruments measure differently depending on
# the structural position of the respondent. H1a tests that claim on the clock;
# this block tests it on the attitude scale, and its result determines which
# comparisons are legitimate in H3a and H3b.
#
# I proceed through three questions, in order:
#   1. Dimensionality  — is the battery one construct or several?
#   2. Invariance      — do men and women use the scale in the same way?
#   3. Scoring         — what enters the substantive models?

message("\n=== BLOCK 2. MEASUREMENT ===")
set.seed(SEED_PREP)

items <- items_attitudes %>%
  mutate(sex = factor(female, levels = 0:1, labels = c("men", "women")))

CONSEQ_ITEMS <- c("att_warm", "att_suffer", "att_famlife")
ROLES_ITEMS  <- c("att_wantshome", "att_breadwin")
CORE_5       <- c(CONSEQ_ITEMS, ROLES_ITEMS)
ALL_6        <- c(CORE_5, "att_fulfil")

FIT_INDICES <- c("chisq.scaled", "df.scaled", "pvalue.scaled",
                 "cfi.scaled", "rmsea.scaled", "srmr")

fit_row <- function(fit, label) {
  m <- fitMeasures(fit, FIT_INDICES)
  tibble(model = label,
         chisq = round(unname(m[1]), 1), df = unname(m[2]),
         p     = round(unname(m[3]), 4),
         CFI   = round(unname(m[4]), 4),
         RMSEA = round(unname(m[5]), 4),
         SRMR  = round(unname(m[6]), 4))
}

# NOTE ON LAVAAN VERSIONS. Under parameterization = "theta" with ordered data,
# the scalar model must free the group-2 residual variances of all but one
# reference indicator per factor (Millsap and Yun-Tein, 2004). Versions 0.6-17
# and earlier did not do this and held every residual variance fixed, producing
# an over-constrained scalar model whose misfit looks like threshold
# non-invariance but is an artefact of the additional restriction. From 0.7-2
# the identification is correct and the scalar model has three fewer degrees of
# freedom here (21 rather than 24). Results obtained under 0.6-x are therefore
# not comparable with those reported in the thesis, and the degrees of freedom
# of the scalar model are printed below as a check.
#
# The appropriate chi-square difference method depends on which test statistic
# lavaan has produced. WLSMV yields "scaled.shifted", for which Satorra (2000)
# is correct; the Satorra-Bentler (2001) correction applies only to the
# Satorra-Bentler statistic and is rejected by lavaan from 0.6-18 onwards. I
# detect the statistic rather than hard-coding the method.
lrt_method <- function(fit) {
  tst <- unlist(lavInspect(fit, "options")$test)
  if (any(grepl("satorra.bentler|yuan.bentler", tst))) {
    "satorra.bentler.2001"
  } else {
    "satorra.2000"
  }
}

cfa_ord <- function(model, items_used = CORE_5, ...) {
  cfa(model, data = items, ordered = items_used,
      estimator = "WLSMV", parameterization = "theta", ...)
}

# -----------------------------------------------------------------------------
# 2.1 Dimensionality
# -----------------------------------------------------------------------------

cat("\n=== 2.1 Polychoric correlations ===\n")
poly <- lavCor(items[ALL_6], ordered = ALL_6)
print(round(poly, 2))

model_1f <- 'GEN =~ att_warm + att_suffer + att_famlife + att_wantshome + att_breadwin'
model_2f <- 'CONSEQ =~ att_warm + att_suffer + att_famlife
             ROLES  =~ att_wantshome + att_breadwin'
model_2f_v5 <- 'CONSEQ =~ att_warm + att_suffer + att_famlife
                ROLES  =~ att_wantshome + att_breadwin + att_fulfil'

fit_1f    <- cfa_ord(model_1f)
fit_2f    <- cfa_ord(model_2f)
fit_2f_v5 <- cfa_ord(model_2f_v5, items_used = ALL_6)

dim_tab <- bind_rows(
  fit_row(fit_1f,    "one factor, 5 items"),
  fit_row(fit_2f,    "two factors, 5 items"),
  fit_row(fit_2f_v5, "two factors + att_fulfil on ROLES"))

cat("\n=== 2.1 Dimensionality ===\n")
print(as.data.frame(dim_tab), row.names = FALSE)

# A second-order model over two first-order factors is NOT identified: with
# only two factors there is no way to estimate a general factor above them.
# This is the formal reason why no single composite index can be defended, and
# it is a matter of identification rather than of fit. I keep the specification
# here as a record of what was attempted:
#   model_2nd <- paste(model_2f, "GEN =~ CONSEQ + ROLES")

cat("\nFactor correlation (two-factor, 5 items): ",
    round(lavInspect(fit_2f, "std")$psi["ROLES", "CONSEQ"], 3), "\n")

cat("\nStandardised loadings:\n")
print(standardizedSolution(fit_2f) %>%
        filter(op == "=~") %>%
        transmute(factor = lhs, item = rhs, loading = round(est.std, 3)) %>%
        as.data.frame(), row.names = FALSE)

cat("\nInternal consistency:\n")
alpha_conseq <- psych::alpha(as.data.frame(base[, CONSEQ_ITEMS]), warnings = FALSE)
alpha_roles  <- psych::alpha(as.data.frame(base[, ROLES_ITEMS]),  warnings = FALSE)
cat(sprintf("  CONSEQ alpha = %.3f (3 items)\n  ROLES  alpha = %.3f (2 items)\n",
            alpha_conseq$total$raw_alpha, alpha_roles$total$raw_alpha))

# -----------------------------------------------------------------------------
# 2.2 Measurement invariance by sex
# -----------------------------------------------------------------------------
# I estimate the standard sequence: configural, then metric (equal loadings),
# then scalar (equal loadings and thresholds).
#
# Decision rule: the scaled chi-square difference test, complemented by the
# incremental criteria dCFI <= .010 and dRMSEA <= .015 (Chen, 2007). The
# incremental criteria take precedence, because chi-square is sensitive to
# sample size and here n = 1,518.

cfa_grp <- function(...) {
  cfa(model_2f, data = items, group = "sex", ordered = CORE_5,
      estimator = "WLSMV", parameterization = "theta", ...)
}

fit_configural <- cfa_grp()
fit_metric     <- cfa_grp(group.equal = "loadings")
fit_scalar     <- cfa_grp(group.equal = c("loadings", "thresholds"))

inv_tab <- bind_rows(fit_row(fit_configural, "configural"),
                     fit_row(fit_metric,     "metric"),
                     fit_row(fit_scalar,     "scalar")) %>%
  mutate(dCFI   = c(NA, round(diff(CFI), 4)),
         dRMSEA = c(NA, round(diff(RMSEA), 4)))

cat("\n=== 2.2 Invariance by sex ===\n")
print(as.data.frame(inv_tab), row.names = FALSE)
cat("\nScalar model df:", fitMeasures(fit_scalar, "df.scaled"),
    "- expect 21 with the correct categorical identification (24 signals",
    "an over-constrained model; check the lavaan version).\n")

LRT_METHOD <- lrt_method(fit_configural)
cat("\nTest statistic:",
    paste(unlist(lavInspect(fit_configural, "options")$test), collapse = ", "),
    "-> difference-test method:", LRT_METHOD, "\n")
cat("Scaled chi-square difference tests:\n")
print(lavTestLRT(fit_configural, fit_metric, fit_scalar, method = LRT_METHOD))

# -----------------------------------------------------------------------------
# 2.3 Locating differential item functioning
# -----------------------------------------------------------------------------
# The score test asks, for every equality constraint in the scalar model, how
# much the fit would improve were that single constraint released. Large values
# identify the parameters that behave differently across sexes: this is
# item-level differential functioning expressed in CFA rather than IRT terms.

score <- lavTestScore(fit_scalar, epc = FALSE)$uni

# From lavaan 0.7-2, lavTestScore() also returns scaled, adjusted and robust
# versions of the statistic when the estimator is robust, so the column names
# differ by version. I prefer the robust statistic where it is available and
# never assume a fixed name.
cat("\nColumns returned by lavTestScore():",
    paste(names(score), collapse = ", "), "\n")

pick <- function(nms, prefer) {
  hit <- prefer[prefer %in% nms]
  if (!length(hit)) stop("no usable column among: ", paste(nms, collapse = ", "))
  hit[1]
}

STAT_COL <- pick(names(score), c("X2.robust", "X2.scaled", "X2"))
P_COL    <- pick(names(score), c("p.value.robust", "p.value.scaled", "p.value"))
cat("Using:", STAT_COL, "/", P_COL, "\n")

ptab <- parTable(fit_scalar)
score$parameter <- vapply(score$lhs, function(lab) {
  r <- ptab[ptab$plabel == lab, ]
  paste0(r$lhs, r$op, r$rhs)
}, character(1))

dif_tab <- score %>%
  as_tibble() %>%
  transmute(parameter,
            X2 = round(.data[[STAT_COL]], 2),
            p  = round(.data[[P_COL]], 4)) %>%
  arrange(desc(X2))

cat("\n=== 2.3 Score test on the scalar constraints (top 10) ===\n")
print(as.data.frame(head(dif_tab, 10)), row.names = FALSE)

FREED <- dif_tab %>% filter(p < 0.01) %>% pull(parameter)

if (length(FREED) == 0) {
  cat("\nNo constraint is released: nothing reaches p < .01, so full scalar",
      "invariance is retained and no partial model is estimated.\n")
  fit_final   <- fit_scalar
  final_label <- "scalar (retained)"
} else {
  cat("\nConstraints released for the partial model (p < .01):\n  ",
      paste(FREED, collapse = ", "), "\n")
  fit_final   <- cfa_grp(group.equal = c("loadings", "thresholds"),
                         group.partial = FREED)
  final_label <- "partial scalar"
  inv_tab <- bind_rows(inv_tab, fit_row(fit_final, final_label)) %>%
    mutate(dCFI   = c(NA, round(diff(CFI), 4)),
           dRMSEA = c(NA, round(diff(RMSEA), 4)))
  cat("\n=== 2.3 Partial scalar invariance ===\n")
  print(as.data.frame(inv_tab), row.names = FALSE)
  cat("\nPartial scalar versus metric:\n")
  print(lavTestLRT(fit_metric, fit_final, method = LRT_METHOD))
}

# -----------------------------------------------------------------------------
# 2.4 Appendix table B7: CFA fit and measurement invariance
# -----------------------------------------------------------------------------
# The dimensionality models and the invariance models are reported in separate
# panels because they answer different questions: Panel A asks how many
# dimensions the battery contains, Panel B whether the retained two-factor
# structure is comparable across sexes.

appendix_fit_row <- function(fit, panel, model, specification,
                             reference_fit = NULL) {
  m <- fitMeasures(fit, FIT_INDICES)

  # nobs is a scalar for a single-group CFA and a vector for a multigroup one.
  n_fit <- sum(unlist(lavInspect(fit, "nobs")))

  # Delta fit indices are only meaningful for nested invariance models.
  if (is.null(reference_fit)) {
    d_cfi   <- NA_real_
    d_rmsea <- NA_real_
  } else {
    ref     <- fitMeasures(reference_fit, c("cfi.scaled", "rmsea.scaled"))
    d_cfi   <- unname(m["cfi.scaled"]   - ref["cfi.scaled"])
    d_rmsea <- unname(m["rmsea.scaled"] - ref["rmsea.scaled"])
  }

  tibble(panel         = panel,
         model         = model,
         specification = specification,
         N             = n_fit,
         chi2_scaled   = round(unname(m["chisq.scaled"]), 2),
         df            = unname(m["df.scaled"]),
         p             = round(unname(m["pvalue.scaled"]), 4),
         CFI           = round(unname(m["cfi.scaled"]), 3),
         RMSEA         = round(unname(m["rmsea.scaled"]), 3),
         SRMR          = round(unname(m["srmr"]), 3),
         delta_CFI     = round(d_cfi, 3),
         delta_RMSEA   = round(d_rmsea, 3))
}

cfa_dimensionality_table <- dim_tab %>%
  mutate(
    panel = "A. Dimensionality",
    specification = case_when(
      model == "one factor, 5 items"               ~ "GEN: V1-V4, V6",
      model == "two factors, 5 items"              ~ "CONSEQ: V1-V3 | ROLES: V4, V6",
      model == "two factors + att_fulfil on ROLES" ~ "CONSEQ: V1-V3 | ROLES: V4-V6",
      TRUE                                         ~ NA_character_),
    model = recode(model,
                   "one factor, 5 items"               = "One-factor model",
                   "two factors, 5 items"              = "Two-factor model",
                   "two factors + att_fulfil on ROLES" = "Two-factor model including V5"),
    N           = sum(unlist(lavInspect(fit_2f, "nobs"))),
    delta_CFI   = NA_real_,
    delta_RMSEA = NA_real_) %>%
  rename(chi2_scaled = chisq) %>%
  select(panel, model, specification, N,
         chi2_scaled, df, p, CFI, RMSEA, SRMR, delta_CFI, delta_RMSEA)

cfa_invariance_table <- bind_rows(
  appendix_fit_row(fit_configural,
                   panel = "B. Measurement invariance by sex",
                   model = "Configural",
                   specification = "Same two-factor structure"),
  appendix_fit_row(fit_metric,
                   panel = "B. Measurement invariance by sex",
                   model = "Metric",
                   specification = "Equal factor loadings",
                   reference_fit = fit_configural),
  appendix_fit_row(fit_scalar,
                   panel = "B. Measurement invariance by sex",
                   model = "Scalar",
                   specification = "Equal loadings + thresholds",
                   reference_fit = fit_metric))

# Should full scalar invariance fail and a partial model be estimated, it is
# added to the panel and compared against the metric model.
if (final_label == "partial scalar") {
  cfa_invariance_table <- bind_rows(
    cfa_invariance_table,
    appendix_fit_row(fit_final,
                     panel = "B. Measurement invariance by sex",
                     model = "Partial scalar",
                     specification = "Equal loadings + partially equal thresholds",
                     reference_fit = fit_metric))
}

cfa_appendix_table <- bind_rows(cfa_dimensionality_table, cfa_invariance_table)

cat("\n=== Appendix B7: CFA fit and measurement invariance ===\n")
print(as.data.frame(cfa_appendix_table), row.names = FALSE)

save_table(cfa_appendix_table, "measurement", "B7_cfa_fit_and_invariance",
           paste("Ordinal confirmatory factor analysis of the gender-attitude battery",
                 "(WLSMV). Panel A reports the dimensionality models, Panel B the",
                 "measurement invariance sequence by sex. Delta indices compare each",
                 "invariance model with the preceding, less restricted one."))

save_table(dif_tab, "measurement", "B7b_score_test_scalar_constraints",
           paste("Score test on each equality constraint of the scalar model, ordered by",
                 "the size of the statistic. Large values would identify item-level",
                 "differential functioning by sex; the threshold for releasing a",
                 "constraint was fixed at p < .01 before estimation."))

# -----------------------------------------------------------------------------
# 2.5 Scale scores
# -----------------------------------------------------------------------------
# Each dimension is measured by the simple mean of its items, requiring at
# least two available items. This is the measure used in the H3 models:
# transparent, reproducible, and it recovers almost all of the latent
# information. Scalar invariance having been retained in section 2.2, the
# comparison between men and women holding the same value on these scales is
# already licensed by the measurement model.

base <- base %>%
  mutate(conseq_index = mean_min(across(all_of(CONSEQ_ITEMS)), min_items = 2),
         roles_index  = mean_min(across(all_of(ROLES_ITEMS)),  min_items = 2))

cat("\n=== 2.5 Scale means by sex ===\n")
print(base %>%
        group_by(female) %>%
        summarise(conseq = round(mean(conseq_index, na.rm = TRUE), 3),
                  roles  = round(mean(roles_index,  na.rm = TRUE), 3),
                  n      = sum(!is.na(conseq_index))) %>%
        as.data.frame(), row.names = FALSE)

cat("\nMissingness on the new scales:\n")
cat(sprintf("  conseq_index %.1f%%   roles_index %.1f%%\n",
            100 * mean(is.na(base$conseq_index)),
            100 * mean(is.na(base$roles_index))))

cat("\n--- Environment ---\n")
cat("  lavaan", as.character(packageVersion("lavaan")),
    "| difference-test method:", LRT_METHOD,
    "| score statistic:", STAT_COL, "\n")

# The verdict is read off the estimated table rather than asserted in advance:
# the conclusion must follow the run, not the other way round.
invariance_ok <- function(i) {
  inv_tab$dCFI[i] >= -0.010 && inv_tab$dRMSEA[i] <= 0.015
}

metric_ok <- invariance_ok(2)
scalar_ok <- invariance_ok(3)

cat("\n--- Measurement decisions carried into the substantive blocks ---\n")
cat("  * the one-factor solution is rejected; no composite index is defensible\n")
cat("  * att_fulfil (V5) stays out of the scales\n")
cat("  * conseq_index and roles_index are the primary predictors in H3\n")
cat(if (metric_ok)
      "  * metric invariance holds -> within-group associations are comparable\n"
    else
      "  * METRIC INVARIANCE FAILS -> associations are not comparable by sex\n")
cat(if (scalar_ok)
      "  * scalar invariance holds -> observed means ARE comparable by sex\n"
    else
      "  * scalar invariance fails -> observed means are NOT comparable by sex\n")
cat("  * retained invariance model:", final_label, "\n")


# =============================================================================
# BLOCK 3. MULTIPLE IMPUTATION OF INCOME AND HOUSEHOLD COMPOSITION
# =============================================================================
# H2 requires equivalised household income, which does not exist in the
# questionnaire and has to be built in a chain:
#
#   income band          -> midpoint -> hh_income_eur
#   hh_size + age counts -> OECD     -> eq_scale
#   hh_income_eur / eq_scale         -> eq_income -> income_q
#
# Every link in that chain contains missing data, and listwise deletion would
# cost roughly a third of the sample, non-randomly.
#
# WHAT IS AND IS NOT IMPUTED. Only the raw inputs are imputed: the income band,
# household size and the count of co-residents under 14. The derived variables
# are recomputed deterministically inside each completed dataset by
# derive_income(), that is, by passive imputation. Imputing a derived variable
# directly would discard information known with certainty (a single-person
# household has eq_scale = 1 whatever its income), could produce combinations
# that no real household could exhibit, and in the case of quartiles would be
# circular, since a quartile is a position in a distribution that does not yet
# exist.

message("\n=== BLOCK 3. MULTIPLE IMPUTATION ===")
set.seed(SEED_PREP)

imp_base <- base %>%
  mutate(
    # Observed information that must not be discarded: where only SOME
    # co-residents lack a year of birth, the observed counts are a floor on the
    # true counts. They enter the imputation as predictors and, in section 3.4,
    # as bounds.
    n_under14_obs  = coalesce(n_under14,  0),
    n_adults14_obs = coalesce(n_adults14, 0),
    n_unknown_age  = pmax(hh_size - 1 - n_roster_valid, 0),
    # The counts themselves count as known only when the roster is complete.
    n_under14        = if_else(roster_complete, n_under14,  NA_real_),
    n_adults14       = if_else(roster_complete, n_adults14, NA_real_),
    hh_income_band_f = ordered(hh_income_band, levels = 1:10),
    no_income_answer = hh_income_src %in% c(77, 99))

# -----------------------------------------------------------------------------
# 3.1 Diagnosing the missingness
# -----------------------------------------------------------------------------
# Two things must be established before imputing: how much is missing, and
# whether it is plausibly missing at random. The MAR argument rests on
# respondents and non-respondents not differing on the outcome once the
# covariates are held constant. This is evidence for the assumption, not a test
# of it, and I treat it as such.

miss_tab <- imp_base %>%
  summarise(across(c(hh_income_band, hh_size, n_under14, n_adults14),
                   ~ round(100 * mean(is.na(.x)), 1))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_missing")

cat("\n=== 3.1 Missingness on the variables to impute ===\n")
print(as.data.frame(miss_tab), row.names = FALSE)

cat("\nMissingness is CLUSTERED, not scattered:\n")
print(imp_base %>%
        group_by(roster_complete) %>%
        summarise(n = n(),
                  income_missing = round(mean(is.na(hh_income_band)), 3),
                  .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("The same respondents withhold household and economic detail, so imputation\n",
    "is doing more work here than the marginal rates on their own suggest.\n")

cat("\n=== 3.1b Evidence for MAR: respondents versus non-respondents on income ===\n")
mar_tab <- imp_base %>%
  filter(!is.na(hh_income_src)) %>%
  mutate(group = if_else(no_income_answer, "no answer", "answered")) %>%
  group_by(group) %>%
  summarise(n           = n(),
            age         = round(mean(age, na.rm = TRUE), 1),
            female      = round(mean(female), 3),
            university  = round(mean(educ4 == 4, na.rm = TRUE), 3),
            in_work     = round(mean(in_work, na.rm = TRUE), 3),
            hours_dom   = round(mean(hours_dom_self,  na.rm = TRUE), 1),
            hours_care  = round(mean(hours_care_self, na.rm = TRUE), 1),
            conflict_fw = round(mean(conflict_fw, na.rm = TRUE), 3),
            .groups = "drop")
print(as.data.frame(mar_tab), row.names = FALSE)

# A formal check on the H2 outcome, since that is the association at stake.
mar_test <- t.test(conflict_fw ~ no_income_answer,
                   data = imp_base %>% filter(!is.na(hh_income_src)))
cat(sprintf("\nDifference in conflict_fw (H2 outcome): %.3f, p = %.3f\n",
            diff(rev(mar_test$estimate)), mar_test$p.value))
cat("Non-response is associated with education and employment, both of which are\n",
    "observed and both of which are used as predictors. The difference on the\n",
    "outcome is small but not negligible, so the MAR argument rests on the\n",
    "covariates carrying it rather than on response being unrelated to the\n",
    "outcome. Code 77 is moreover an EXPLICIT refusal, which is informative\n",
    "behaviour and makes MAR untestable in principle. The delta analysis in\n",
    "section 3.5 is therefore not optional here.\n")

save_table(miss_tab, "imputation", "I1_missingness_on_imputed_variables",
           paste("Percentage of missing values on the variables entering the imputation",
                 "model. Household income and household composition are the two inputs",
                 "of the equivalised-income chain required by H2."))

save_table(mar_tab, "imputation", "I2_mar_evidence",
           paste("Respondents who answered the household-income question compared with",
                 "those who did not, on the covariates and on the H2 outcome. Presented",
                 "as evidence bearing on the plausibility of the MAR assumption, not as",
                 "a test of it."))

# -----------------------------------------------------------------------------
# 3.2 Specifying the imputation model
# -----------------------------------------------------------------------------

TO_IMPUTE <- c("hh_income_band_f", "hh_size", "n_under14")

# STRUCTURALLY missing predictors. conflict_fw does not exist for people who
# are not in work; task_load and fair_share do not exist for people without a
# cohabiting partner. These are not unknown values but non-applicable ones, and
# imputing them would invent a housework split for someone who lives alone.
# Each therefore enters as a value-plus-indicator pair: the value with
# non-applicability coded to zero, and a flag marking that it was
# non-applicable.
STRUCTURAL <- c("conflict_fw", "conflict_wf", "work_hours", "fair_share",
                "task_load", "hours_dom_self", "hours_care_self")

# Predictors with ordinary, stochastic missingness. These ARE imputed, because
# a predictor that is never imputed blocks every row in which it is missing:
# that was the reason an earlier version of this model left 302 households
# uncompleted.
AUXILIARY <- c("age", "educ4", "in_work", "econ_strain", "subj_status",
               "marital_r", "conseq_index", "roles_index")

COMPLETE_PRED <- c("female", "has_minor", "has_elderly",
                   "n_under14_obs", "n_adults14_obs", "tamuni_r", "mode")

imp_data <- imp_base %>%
  mutate(marital_r = factor(na_if(MARITAL, 9)),
         tamuni_r  = factor(TAMUNI)) %>%
  select(id, all_of(TO_IMPUTE), all_of(AUXILIARY), all_of(COMPLETE_PRED),
         all_of(STRUCTURAL)) %>%
  mutate(across(all_of(STRUCTURAL),
                list(na = ~ as.integer(is.na(.x)),
                     z  = ~ coalesce(.x, 0)),
                .names = "{.col}_{.fn}")) %>%
  select(-all_of(STRUCTURAL))

STRUCT_COLS <- c(paste0(STRUCTURAL, "_na"), paste0(STRUCTURAL, "_z"))
PREDICTORS  <- c(AUXILIARY, COMPLETE_PRED, STRUCT_COLS)

cat("\n=== 3.2 Imputation model ===\n")
cat("Targets   :", paste(TO_IMPUTE, collapse = ", "), "\n")
cat("Auxiliary variables also imputed:", paste(AUXILIARY, collapse = ", "), "\n")
cat("Structural (value + indicator)  :", paste(STRUCTURAL, collapse = ", "), "\n")
cat("Predictors in total:", length(PREDICTORS), "\n")

# -----------------------------------------------------------------------------
# 3.3 Running mice
# -----------------------------------------------------------------------------
# METHOD FOR INCOME: a documented decision rather than a default. The band is
# ordinal with ten categories, so the proportional-odds model (polr) is the
# method that matches its measurement level, and it is the one used. It is
# slow, taking roughly three minutes here with twenty imputations and ten
# iterations. The number of iterations was reduced from twenty to ten because
# mice converges well before the tenth on this problem; the convergence trace
# produced below is the check on that decision.
#
# Predictive mean matching on the numeric band is the fallback should polr
# fail. PMM is standard for bracketed income, since it draws from observed
# donors and can therefore only return bands that exist, but it treats the
# bands as equally spaced, which polr does not. polr is preferred while it
# remains affordable.

ini  <- mice(imp_data, maxit = 0, printFlag = FALSE)

meth <- ini$method                     # mice's defaults for the auxiliaries
meth[STRUCT_COLS]        <- ""         # never impute a structural indicator
meth["id"]               <- ""
meth["hh_income_band_f"] <- "polr"
meth["hh_size"]          <- "pmm"
meth["n_under14"]        <- "pmm"

pred <- ini$predictorMatrix
pred[, "id"] <- 0
pred["id", ] <- 0
diag(pred)   <- 0

run_mice <- function(method) {
  mice(imp_data, m = N_IMPUTATIONS, maxit = N_ITERATIONS, method = method,
       predictorMatrix = pred, seed = SEED_PREP, printFlag = FALSE)
}

imp <- try(run_mice(meth), silent = TRUE)
INCOME_METHOD <- "polr"

if (inherits(imp, "try-error")) {
  cat("\npolr failed; falling back to pmm on the numeric band.\n")
  imp_data$hh_income_band_f <- as.numeric(as.character(imp_data$hh_income_band_f))
  meth["hh_income_band_f"]  <- "pmm"
  INCOME_METHOD <- "pmm"
  imp <- run_mice(meth)
}

cat("\n=== 3.3 Convergence ===\n")
cat("Income method:", INCOME_METHOD, "| m =", imp$m, "| maxit =", imp$iteration, "\n")

# mice returns NULL rather than an empty data frame when nothing was logged, so
# nrow() cannot be called on it directly.
logged   <- imp$loggedEvents
n_logged <- if (is.null(logged)) 0L else nrow(logged)

if (n_logged > 0) {
  cat("Logged events:", n_logged,
      "- predictors dropped for collinearity, repeated once per iteration\n",
      "  and per imputation. Distinct causes:\n")
  print(logged %>% count(dep, out, name = "times") %>% arrange(desc(times)) %>%
          as.data.frame(), row.names = FALSE)
} else {
  cat("No logged events.\n")
}

# The trace lines for the three targets should mix without trend across the ten
# iterations. Without the `y` argument mice plots every imputed variable and
# paginates, and png() would capture only the last page, which is an auxiliary
# rather than one of the variables whose convergence is at stake.
png(file.path(OUT_ROOT, "figures", "imputation", "I3_convergence.png"),
    width = 900, height = 700)
print(plot(imp, y = TO_IMPUTE, layout = c(2, 3)))
dev.off()
cat("Convergence trace written to output/figures/imputation/I3_convergence.png",
    "- inspect it before trusting the imputations.\n")

# The targets must come out fully completed. If they do not, a predictor is
# still blocking rows and the diagnosis above is where to look.
completeness <- complete(imp, 1) %>%
  summarise(across(all_of(TO_IMPUTE), ~ sum(is.na(.x))))
cat("\nRemaining NAs in the targets, first completed dataset:\n")
print(as.data.frame(completeness), row.names = FALSE)
stopifnot(all(completeness == 0))

# -----------------------------------------------------------------------------
# 3.4 Completing the datasets and deriving the income variables
# -----------------------------------------------------------------------------
# derive_income() is applied here, once per completed dataset. The quartile
# cut-points are recalculated within each dataset, because the distribution of
# equivalised income differs across them; that variation is precisely what
# Rubin's rules convert into a wider standard error.
#
# COHERENCE CONSTRAINT. Household size and the count of under-14s are imputed
# separately, so nothing guarantees that the counts add up. They are clamped to
# the feasible range implied by what was actually observed,
#     n_under14_obs <= n_under14 <= hh_size - 1,
# which is a documented approximation to a properly constrained imputation and
# never overrides an observed value.

harmonise_counts <- function(df) {
  df %>%
    mutate(hh_size    = pmax(round(hh_size), 1),
           n_under14  = pmin(pmax(round(n_under14), n_under14_obs),
                             pmax(hh_size - 1, 0)),
           n_adults14 = pmax(hh_size - 1 - n_under14, 0),
           hh_income_band = as.numeric(as.character(hh_income_band_f)))
}

completed <- complete(imp, action = "long", include = TRUE) %>%
  as_tibble() %>%
  harmonise_counts() %>%
  group_by(.imp) %>%
  group_modify(~ derive_income(.x)) %>%
  ungroup()

cat("\n=== 3.4 Completed data ===\n")
cat("Datasets:", n_distinct(completed$.imp) - 1, "plus the original\n")
cat("eq_income available after imputation: ",
    round(100 * mean(!is.na(completed$eq_income[completed$.imp > 0])), 1),
    "%\n", sep = "")

cat("\nEquivalised income by imputation (mean, first five datasets):\n")
print(completed %>%
        filter(.imp %in% 1:5) %>%
        group_by(.imp) %>%
        summarise(mean_eq = round(mean(eq_income, na.rm = TRUE)),
                  q1_cut  = round(quantile(eq_income, .25, na.rm = TRUE)),
                  q3_cut  = round(quantile(eq_income, .75, na.rm = TRUE)),
                  .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)
cat("Between-dataset variation in the cut-points is expected and is the point:\n",
    "it is what Rubin's rules convert into honest standard errors.\n")

# Sanity check: observed values must be untouched by the imputation.
obs_check <- completed %>%
  filter(.imp == 1) %>%
  left_join(base %>% select(id, obs_band = hh_income_band), by = "id") %>%
  filter(!is.na(obs_band)) %>%
  summarise(altered = sum(obs_band != hh_income_band))
stopifnot(obs_check$altered == 0)
cat("\nObserved income bands unchanged by imputation: OK\n")

# -----------------------------------------------------------------------------
# 3.5 Delta sensitivity for non-response on income
# -----------------------------------------------------------------------------
# Code 77 is "prefer not to answer": an explicit refusal rather than
# inattention. That is informative behaviour, so MAR cannot be verified and a
# sensitivity analysis is required rather than optional. The imputed bands are
# shifted by half a band in each direction and H2 is re-estimated on the
# shifted sets in Block 7.
#
# The shift is applied to EVERY imputed band, not only to the explicit
# refusals. Of the missing bands, the majority are code 77, a further group are
# code 99, and the remainder belong to households whose size was itself
# undeclared, so that the source variable could not even be identified.
# Restricting the shift to the explicit refusals would leave the least
# informative third of the non-response unperturbed, which is the opposite of a
# conservative check.
#
# THE FLAG IS BUILT ON THE PRE-IMPUTATION DATA. An earlier version derived it
# from is.na() inside the completed data, where nothing is missing any more:
# the flag was empty, the shifted datasets were identical copies of the
# original, and nothing failed loudly. Hence the guard below.
#
# THE SHIFT IS A GENUINE HALF BAND. The shifted band is passed to
# derive_income(), which interpolates between adjacent midpoints instead of
# rounding to the nearest whole band; see the note there. Under rounding, half
# of the imputed cases would not have moved at all, and which half would have
# been decided by the parity of the band rather than by anything substantive.

delta_flag <- imp_base %>%
  transmute(id,
            band_imputed     = is.na(hh_income_band),
            explicit_refusal = coalesce(income_refused == 1, FALSE))

completed <- completed %>% left_join(delta_flag, by = "id")

n_delta <- sum(completed$band_imputed[completed$.imp == 1])
cat("\n=== 3.5 Delta sensitivity ===\n")
cat("Bands imputed and therefore shifted:", n_delta, "\n")
cat("  of which explicit refusals (code 77):",
    sum(completed$explicit_refusal[completed$.imp == 1]), "\n")
stopifnot(n_delta > 300)

delta_shift <- function(delta) {
  completed %>%
    mutate(hh_income_band = if_else(.imp > 0 & band_imputed,
                                    pmin(pmax(hh_income_band + delta, 1), 10),
                                    hh_income_band)) %>%
    group_by(.imp) %>%
    group_modify(~ derive_income(.x)) %>%
    ungroup()
}

delta_sets <- list(minus = delta_shift(-0.5), plus = delta_shift(+0.5))

# The shifted sets must actually differ from the original; a silent no-op here
# was the original bug.
delta_check <- tibble(
  set = c("delta -0.5", "baseline", "delta +0.5"),
  mean_eq_income = c(mean(delta_sets$minus$eq_income[delta_sets$minus$.imp > 0]),
                     mean(completed$eq_income[completed$.imp > 0]),
                     mean(delta_sets$plus$eq_income[delta_sets$plus$.imp > 0])),
  q1_cut = c(quantile(delta_sets$minus$eq_income[delta_sets$minus$.imp > 0], .25),
             quantile(completed$eq_income[completed$.imp > 0], .25),
             quantile(delta_sets$plus$eq_income[delta_sets$plus$.imp > 0], .25))) %>%
  mutate(across(where(is.numeric), ~ round(.x)))

cat("\nShifted sets versus baseline:\n")
print(as.data.frame(delta_check), row.names = FALSE)
stopifnot(delta_check$mean_eq_income[1] < delta_check$mean_eq_income[2],
          delta_check$mean_eq_income[3] > delta_check$mean_eq_income[2])
cat("Both shifted sets differ from the baseline in the expected direction: OK\n")

save_table(delta_check, "imputation", "I4_delta_shift_check",
           paste("Mean equivalised income and lower quartile cut-point in the baseline",
                 "imputations and in the two shifted sets. Reported to document that the",
                 "delta sensitivity analysis actually perturbs the data in the intended",
                 "direction."))

cat("\n--- Carried into the model blocks ---\n")
cat("  * income method:", INCOME_METHOD, "| m =", N_IMPUTATIONS,
    "| seed =", SEED_PREP, "\n")
cat("  * models on eq_income / income_q are fitted per imputation and pooled\n")
cat("    with Rubin's rules, never on the average of the imputed datasets\n")
cat("  * the fraction of missing information is reported with the coefficients\n")


# =============================================================================
# BLOCK 4. DESCRIPTIVE ANALYSIS
# =============================================================================
# This block produces two kinds of output. The first is the set of sample-level
# tables that document the data and the cleaning decisions, most of which
# belong in the appendix. The second is the raw descriptive that precedes each
# hypothesis and that I ran before estimating any model. I retain the whole set
# here, including the exploratory tables that are not reproduced in the final
# text, because they record the sequence in which the evidence was examined and
# because several of them are the unadjusted form of results that the models
# later present in adjusted form.

message("\n=== BLOCK 4. DESCRIPTIVES ===")

design <- svydesign(ids = ~1, strata = ~stratum, weights = ~w,
                    data = base, nest = TRUE)

# -----------------------------------------------------------------------------
# D1. Response structure of the hours questions, by sex
# -----------------------------------------------------------------------------
# This table does double duty: it describes the sample and it is the raw form
# of H1a. It must be read VERTICALLY, since the comparison that matters is
# between tasks within a sex, not between sexes within a task.
#
# "Gave a figure" is split into positive and zero hours because they are
# different answers: a substantial share of those who quantify care report zero
# hours, which is not an inability to say but the absence of care work.

response_structure <- base %>%
  select(female, w,
         domestic = st_hours_dom_self, care = st_hours_care_self,
         h_dom = hours_dom_self, h_care = hours_care_self) %>%
  pivot_longer(c(domestic, care), names_to = "task", values_to = "status") %>%
  mutate(hours = if_else(task == "domestic", h_dom, h_care),
         category = case_when(status == "figure" & hours >  0 ~ "positive_hours",
                              status == "figure" & hours == 0 ~ "zero_hours",
                              status == "cant_say"            ~ "cant_say",
                              status == "refused"             ~ "refused"),
         task = factor(task, levels = c("domestic", "care"),
                       labels = c("Domestic work", "Care work")),
         sex  = sexlab(female)) %>%
  filter(!is.na(category)) %>%
  group_by(task, sex) %>%
  summarise(n = n(),
            positive_hours = pct(weighted.mean(category == "positive_hours", w)),
            zero_hours     = pct(weighted.mean(category == "zero_hours", w)),
            cant_say       = pct(weighted.mean(category == "cant_say", w)),
            refused        = pct(weighted.mean(category == "refused", w)),
            .groups = "drop") %>%
  select(task, sex, positive_hours, zero_hours, cant_say, refused, n)

save_table(response_structure, "descriptive", "TableA1_response_structure",
           paste("Response structure of the weekly-hours questions by task and sex.",
                 "Weighted percentages (PESOFIN), unweighted n.",
                 "Read down the task column rather than across sexes."))

# -----------------------------------------------------------------------------
# D2. Hours reported, by sex
# -----------------------------------------------------------------------------
# Both the mean and the median are reported, because they diverge sharply and
# the divergence is itself the point: a small group of respondents report more
# than one hundred weekly hours of care and some report the entire week. I use
# medians in the text and retain the means to document the skew.

hours_summary <- map_dfr(
  list(list(var = "hours_dom_self",  label = "Domestic work"),
       list(var = "hours_care_self", label = "Care work")),
  function(v) {
    f  <- as.formula(paste0("~", v$var))
    mn <- svyby(f, ~female, design, svymean, na.rm = TRUE)
    map_dfr(0:1, function(fem) {
      d <- base %>% filter(female == fem, !is.na(.data[[v$var]]))
      q <- wquantile(d[[v$var]], d$w, c(.25, .5, .75))
      tibble(task   = v$label,
             sex    = sexlab(fem),
             mean   = round(mn[[2]][mn$female == fem], 1),
             se     = round(mn[[3]][mn$female == fem], 2),
             p25    = q[1], median = q[2], p75 = q[3],
             max    = max(d[[v$var]]),
             n      = nrow(d))
    })
  })

save_table(hours_summary, "descriptive", "TableA2_hours_by_sex",
           paste("Weekly hours of reproductive labour among respondents who supplied a",
                 "figure. Design-based weighted means and weighted quantiles. Means and",
                 "medians diverge because of a long upper tail, so medians are reported",
                 "in the text."))

# -----------------------------------------------------------------------------
# D3. Sample composition, weighted and unweighted
# -----------------------------------------------------------------------------
# This table documents the quota design and shows what the weight does. Sex,
# age and education are calibration variables for PESOFIN, so their weighted
# margins are targets rather than findings.

composition <- bind_rows(
  base %>% count(level = sexlab(female)) %>% mutate(variable = "Sex"),
  base %>% mutate(level = cut(age, c(17, 29, 44, 59, 120),
                              labels = c("18-29", "30-44", "45-59", "60+"))) %>%
    count(level) %>% mutate(variable = "Age"),
  base %>% count(level = factor(educ4)) %>% mutate(variable = "Education (NIVELEDU)"),
  base %>% count(level = mode) %>% mutate(variable = "Mode")) %>%
  filter(!is.na(level)) %>%
  mutate(level = as.character(level)) %>%
  left_join(
    bind_rows(
      base %>% group_by(level = sexlab(female)) %>% summarise(wn = sum(w)),
      base %>% mutate(level = cut(age, c(17, 29, 44, 59, 120),
                                  labels = c("18-29", "30-44", "45-59", "60+"))) %>%
        group_by(level) %>% summarise(wn = sum(w)),
      base %>% group_by(level = factor(educ4)) %>% summarise(wn = sum(w)),
      base %>% group_by(level = mode) %>% summarise(wn = sum(w))) %>%
      mutate(level = as.character(level)) %>%
      filter(!is.na(level)),
    by = "level") %>%
  group_by(variable) %>%
  mutate(unweighted_pct = pct(n / sum(n)),
         weighted_pct   = pct(wn / sum(wn))) %>%
  ungroup() %>%
  select(variable, level, n, unweighted_pct, weighted_pct)

save_table(composition, "descriptive", "D3_sample_composition",
           paste("Sample composition, unweighted and weighted with PESOFIN. Sex, age and",
                 "education are calibration variables for the weight, so their weighted",
                 "margins are targets rather than findings."))

# -----------------------------------------------------------------------------
# D4. Audit of the special codes
# -----------------------------------------------------------------------------
# How many cases each special code actually affects. This makes the cleaning
# auditable: without it the analytic n cannot be reproduced from the raw file.

code_audit <- bind_rows(
  tibble(block = "Hours", variable = c("V34", "V35", "V36", "V37")) %>%
    mutate(not_applicable = map_int(variable, ~ sum(base[[.x]] == 996)),
           cant_say       = map_int(variable, ~ sum(base[[.x]] == 998)),
           refused        = map_int(variable, ~ sum(base[[.x]] == 999))),
  tibble(block = "Attitudes", variable = c("V1", "V2", "V3", "V4", "V5", "V6")) %>%
    mutate(not_applicable = 0L,
           cant_say       = map_int(variable, ~ sum(base[[.x]] == 8)),
           refused        = map_int(variable, ~ sum(base[[.x]] == 9))),
  tibble(block = "Task division",
         variable = c("V39", "V40", "V41", "V42", "V43", "V44")) %>%
    mutate(not_applicable = map_int(variable, ~ sum(base[[.x]] == 0)),
           cant_say       = map_int(variable, ~ sum(base[[.x]] == 6)),
           refused        = map_int(variable, ~ sum(base[[.x]] %in% c(8, 9)))),
  tibble(block = "Work-family conflict", variable = c("V46", "V47", "V48", "V49")) %>%
    mutate(not_applicable = map_int(variable, ~ sum(base[[.x]] == 8)),
           cant_say       = 0L,
           refused        = map_int(variable, ~ sum(base[[.x]] == 9))),
  tibble(block = "Fairness", variable = "V45") %>%
    mutate(not_applicable = sum(base$V45 == 0),
           cant_say       = sum(base$V45 == 8),
           refused        = sum(base$V45 == 9)),
  tibble(block = "Income", variable = c("NAT_INC", "NAT_RINC")) %>%
    mutate(not_applicable = map_int(variable, ~ sum(base[[.x]] == 0)),
           cant_say       = map_int(variable, ~ sum(base[[.x]] == 77)),
           refused        = map_int(variable, ~ sum(base[[.x]] == 99)))) %>%
  rename(code_A_not_applicable   = not_applicable,
         code_B_cantsay_or_other = cant_say,
         code_C_refused          = refused)

save_table(code_audit, "descriptive", "D4_special_codes_audit",
           paste("Number of cases affected by each special code, by block of variables.",
                 "Column B holds 'couldn't say' for hours and attitudes, 'someone else'",
                 "for the task-division items, and code 77 'prefer not to answer' for",
                 "income."))

# -----------------------------------------------------------------------------
# D5. Availability of the derived variables
# -----------------------------------------------------------------------------

derived_inventory <- tibble(
  variable = c("conseq_index", "roles_index", "task_load", "mental_load",
               "fair_share", "conflict_fw", "conflict_wf", "work_hours",
               "hh_income_band", "eq_scale_inputs_complete", "outsourced",
               "n_minors", "partner_known", "diff_sex_couple")) %>%
  mutate(
    n_valid = map_int(variable, function(v) {
      if (v == "eq_scale_inputs_complete") sum(base$roster_complete)
      else sum(!is.na(base[[v]]))
    }),
    pct_missing = round(100 * (nrow(base) - n_valid) / nrow(base), 1))

save_table(derived_inventory, "descriptive", "D5_derived_variables",
           paste("Availability of the derived analysis variables. Missingness on",
                 "task_load, mental_load and fair_share is structural: those items are",
                 "only asked of respondents with a cohabiting partner."))

# -----------------------------------------------------------------------------
# D8. Cells with fewer than thirty unweighted cases
# -----------------------------------------------------------------------------
# Identified in advance so that no percentage computed on a thin cell is read
# as a finding.

small_cells <- bind_rows(
  base %>% count(variable = "sex x education",
                 level = paste(sexlab(female), educ4)),
  base %>% filter(partnered == 1) %>%
    count(variable = "sex x fairness (V45)",
          level = paste(sexlab(female), fair_raw)),
  base %>% filter(in_work == 1) %>%
    count(variable = "sex x income band (1-10, pre-imputation)",
          level = paste(sexlab(female), hh_income_band))) %>%
  filter(n < 30, !grepl("NA", level))

save_table(small_cells, "descriptive", "D8_small_cells",
           paste("Cross-tabulation cells with fewer than thirty unweighted cases.",
                 "Percentages computed within these cells should not be interpreted."))

# -----------------------------------------------------------------------------
# H1a descriptive: the discordance table is itself the result
# -----------------------------------------------------------------------------
# Only the off-diagonal cells identify the conditional logit estimated in Block
# 5. The concordant respondents contribute nothing, because their individual
# propensity to quantify already accounts for both of their answers. With two
# observations per person and a binary predictor the conditional estimate has a
# closed form, which is what this table reports.

message("\n-- H1a descriptives --")

h1a_wide <- h1a_long %>%
  select(id, task, cant_say) %>%
  pivot_wider(names_from = task, values_from = cant_say)

tab_h1a <- table(domestic = ifelse(h1a_wide$domestic == 1, "cant_say", "figure"),
                 care     = ifelse(h1a_wide$care == 1, "cant_say", "figure"))

discordance <- as.data.frame.matrix(addmargins(tab_h1a)) %>%
  rownames_to_column("domestic") %>%
  rename(care_figure = figure, care_cant_say = cant_say, total = Sum)

b_cell  <- tab_h1a["figure", "cant_say"]   # quantifies domestic work, not care
c_cell  <- tab_h1a["cant_say", "figure"]   # the reverse
or_h1a  <- b_cell / c_cell
se_h1a  <- sqrt(1 / b_cell + 1 / c_cell)
mcnemar <- mcnemar.test(tab_h1a)

discordance_stats <- tibble(
  statistic = c("Discordant: domestic figure / care can't say",
                "Discordant: domestic can't say / care figure",
                "Total discordant pairs",
                "Concordant pairs (contribute nothing to the estimate)",
                "Conditional odds ratio (b/c)",
                "95% CI lower", "95% CI upper",
                "McNemar chi-squared", "McNemar p-value"),
  value = c(b_cell, c_cell, b_cell + c_cell, sum(tab_h1a) - b_cell - c_cell,
            round(or_h1a, 3),
            round(exp(log(or_h1a) - 1.96 * se_h1a), 3),
            round(exp(log(or_h1a) + 1.96 * se_h1a), 3),
            round(unname(mcnemar$statistic), 2),
            format.pval(mcnemar$p.value, digits = 3)))

save_table(discordance, "H1a", "TableB1_discordance_table",
           paste("Within-person cross-tabulation of the ability to quantify domestic and",
                 "care work in weekly hours. Only the off-diagonal cells identify the",
                 "effect; the concordant respondents contribute nothing."))

save_table(discordance_stats, "H1a", "H1a_discordance_statistics",
           paste("Conditional odds ratio computed from the discordant pairs, with the",
                 "McNemar test. With two observations per person and a binary predictor",
                 "the conditional logit has this closed form."))

# -----------------------------------------------------------------------------
# H1b descriptive: the four cells, unadjusted
# -----------------------------------------------------------------------------

message("\n-- H1b descriptives --")

cells_raw <- h1b_long %>%
  group_by(target = if_else(target_female == 1, "Woman", "Man"), task, source) %>%
  summarise(n = n(),
            mean_raw        = round(mean(hours), 1),
            mean_winsorised = round(mean(hours_w), 1),
            median          = median(hours),
            p90             = quantile(hours, .9),
            .groups = "drop") %>%
  arrange(target, task, source)

save_table(cells_raw, "H1b", "H1b_cells_unadjusted",
           paste("Reported weekly hours by target sex, task and source of the report,",
                 "unadjusted. The adjusted gaps and the RIF decomposition are estimated",
                 "in Block 6."))

winsorisation_caps <- h1b_long %>%
  group_by(task) %>%
  summarise(cap_p95    = max(hours_w),
            max_raw    = max(hours),
            pct_capped = round(100 * mean(hours > max(hours_w)), 1))

save_table(winsorisation_caps, "H1b", "H1b_winsorisation_caps",
           paste("Winsorisation caps (95th percentile within task, pooled across",
                 "sources) and the share of reports affected by them."))

# -----------------------------------------------------------------------------
# H2 descriptive: the conflict battery
# -----------------------------------------------------------------------------

message("\n-- H2 descriptives --")

workers <- base %>% filter(in_work == 1)

conflict_items <- tibble(
  variable = c("V46", "V47", "V48", "V49"),
  item = c("Too tired from work for household tasks",
           "Hard to meet family responsibilities because of work time",
           "Arrived at work too tired because of household tasks",
           "Hard to concentrate at work because of family responsibilities"),
  direction = c("work -> family", "work -> family",
                "family -> work", "family -> work")) %>%
  mutate(
    men_several_per_week = map_dbl(variable, ~ pct(weighted.mean(
      workers[[.x]][workers$female == 0 & workers[[.x]] %in% 1:4] == 1,
      workers$w[workers$female == 0 & workers[[.x]] %in% 1:4]))),
    women_several_per_week = map_dbl(variable, ~ pct(weighted.mean(
      workers[[.x]][workers$female == 1 & workers[[.x]] %in% 1:4] == 1,
      workers$w[workers$female == 1 & workers[[.x]] %in% 1:4]))),
    n = map_int(variable, ~ sum(workers[[.x]] %in% 1:4)))

save_table(conflict_items, "H2", "H2_conflict_items_by_sex",
           paste("Work-family conflict items among employed respondents: weighted",
                 "percentage answering 'several times a week'. The double-presence claim",
                 "concerns the family-to-work direction."))

conflict_indices <- workers %>%
  filter(!is.na(conflict_fw), !is.na(work_hours)) %>%
  group_by(sex = sexlab(female)) %>%
  summarise(n           = n(),
            conflict_fw = round(weighted.mean(conflict_fw, w), 3),
            conflict_wf = round(weighted.mean(conflict_wf, w, na.rm = TRUE), 3),
            work_hours  = round(weighted.mean(work_hours, w), 1),
            .groups = "drop")

save_table(conflict_indices, "H2", "H2_conflict_indices_by_sex",
           paste("Family-to-work and work-to-family conflict indices (1-4, higher means",
                 "more conflict) among employed respondents, with mean weekly hours of",
                 "paid employment."))

strain <- workers %>%
  filter(!is.na(econ_strain), !is.na(conflict_fw)) %>%
  group_by(econ_strain) %>%
  summarise(n = n(),
            conflict_fw = round(weighted.mean(conflict_fw, w), 3),
            .groups = "drop")

save_table(strain, "H2", "H2_conflict_by_economic_strain",
           paste("Family-to-work conflict by declared economic difficulty. Specified in",
                 "advance as the second measure of material position alongside",
                 "equivalised income; both are reported whatever the results."))

# -----------------------------------------------------------------------------
# H3a descriptive: discourse and the division of tasks
# -----------------------------------------------------------------------------

message("\n-- H3a descriptives --")

subscales <- base %>%
  group_by(sex = sexlab(female)) %>%
  summarise(n_conseq = sum(!is.na(conseq_index)),
            conseq   = round(weighted.mean(conseq_index, w, na.rm = TRUE), 3),
            n_roles  = sum(!is.na(roles_index)),
            roles    = round(weighted.mean(roles_index, w, na.rm = TRUE), 3),
            .groups = "drop")

save_table(subscales, "H3a", "H3a_discourse_subscales_by_sex",
           paste("The two dimensions of the naturalising discourse by sex (1-5, higher",
                 "means stronger adherence). Scalar invariance across sexes was",
                 "established in Block 2, so these means are comparable. A single index",
                 "masked these differences, since they run in opposite directions and",
                 "cancelled out."))

partnered_sample <- base %>% filter(partnered == 1)

TASKS <- c(task_laundry   = "Laundry",
           task_planning  = "Planning family activities",
           task_sickcare  = "Care of sick relatives",
           task_groceries = "Grocery shopping",
           task_cleaning  = "Cleaning",
           task_meals     = "Preparing meals")

task_division <- imap_dfr(TASKS, function(label, v) {
  d <- partnered_sample %>% filter(!is.na(.data[[v]]))
  tibble(item = label,
         men_mostly_self   = pct(weighted.mean(d[[v]][d$female == 0] >= 4,
                                               d$w[d$female == 0])),
         women_mostly_self = pct(weighted.mean(d[[v]][d$female == 1] >= 4,
                                               d$w[d$female == 1])),
         n = nrow(d))
})

save_table(task_division, "H3a", "H3a_task_division_by_sex",
           paste("Share reporting that they always or usually perform each task",
                 "themselves, among respondents with a cohabiting partner. Planning",
                 "family activities is the mental-load item and is available only in",
                 "this wave."))

# Code 6 has to be read off the RAW items. task_recode() keeps only the five
# points of the scale and sends 6 to NA, precisely because it is not a point on
# the scale, so the recoded variables can never equal 6 and this table would
# come out as a column of zeros if it were built from them. The denominator is
# the set of valid answers (1-6): the non-applicable, "don't know" and refusal
# codes are not couples in which the task is done by somebody else.
TASK_RAW <- c(task_laundry   = "V39",
              task_planning  = "V40",
              task_sickcare  = "V41",
              task_groceries = "V42",
              task_cleaning  = "V43",
              task_meals     = "V44")

outsourcing <- imap_dfr(TASKS, function(label, v) {
  answers <- partnered_sample[[TASK_RAW[[v]]]]
  answers <- answers[answers %in% 1:6]
  tibble(item             = label,
         pct_someone_else = pct(mean(answers == 6)),
         n                = length(answers))
})

save_table(outsourcing, "H3a", "H3a_outsourcing_by_task",
           paste("Share of couples in which the task is performed by a third party (code",
                 "6). This is not a point on the scale: it means that neither partner",
                 "does it, which is why cleaning shows the highest item missingness in",
                 "the declared task-load index."))

# -----------------------------------------------------------------------------
# H3b descriptive: perceived proportionality by sex and declared load
# -----------------------------------------------------------------------------

message("\n-- H3b descriptives --")

fairness_sample <- base %>% filter(partnered == 1, !is.na(fair_raw))

fairness_dist <- fairness_sample %>%
  group_by(sex = sexlab(female)) %>%
  summarise(n = n(),
            does_much_more = pct(weighted.mean(fair_raw == 1, w)),
            does_more      = pct(weighted.mean(fair_raw == 2, w)),
            fair_share     = pct(weighted.mean(fair_raw == 3, w)),
            does_less      = pct(weighted.mean(fair_raw == 4, w)),
            does_much_less = pct(weighted.mean(fair_raw == 5, w)),
            .groups = "drop")

fairness_dist <- fairness_dist %>%
  rename(Sex = sex,
         `Does much more than their share` = does_much_more,
         `Does more`                       = does_more,
         `Does roughly their share`        = fair_share,
         `Does less`                       = does_less,
         `Does much less`                  = does_much_less)

save_table(fairness_dist, "H3b", "Table5_fairness_distribution_by_sex",
           paste("Distribution of V45 by sex among respondents with a cohabiting",
                 "partner. The quantity described by the hypotheses is the middle",
                 "category, which is why the dependent variable is binary rather than",
                 "ordinal."))

fair_by_load <- fairness_sample %>%
  filter(!is.na(task_load)) %>%
  mutate(load_band = cut(task_load, c(0, 2, 2.75, 3.25, 4, 5),
                         labels = c("<=2", "2-2.75", "~3 (equal)",
                                    "3.25-4", ">4"))) %>%
  group_by(load_band, sex = sexlab(female)) %>%
  summarise(n = n(), pct_fair = pct(mean(fair_share)), .groups = "drop") %>%
  pivot_wider(names_from = sex, values_from = c(n, pct_fair))

save_table(fair_by_load, "H3b", "H3b_fairness_by_load_and_sex",
           paste("Percentage perceiving their own contribution as proportionate, by",
                 "declared task load and sex. This is the descriptive form of the central",
                 "finding: at the same declared load, men and women judge differently."))

equity_threshold <- fairness_sample %>%
  filter(fair_share == 1, !is.na(task_load)) %>%
  group_by(sex = sexlab(female)) %>%
  summarise(n = n(),
            mean_load_when_fair   = round(weighted.mean(task_load, w), 3),
            median_load_when_fair = median(task_load),
            .groups = "drop")

save_table(equity_threshold, "H3b", "H3b_equity_threshold",
           paste("Declared task load among respondents who say the split is what",
                 "corresponds to them, by sex, on a 1-5 scale in which 3 is an equal",
                 "split. Men locate the equity point below parity and women only above",
                 "it."))


# =============================================================================
# BLOCK 5. H1a — THE MEASURABILITY OF CARE
# =============================================================================
# H1a states that care is harder to quantify in units of time than domestic
# work. The design is within-person: each respondent answers both questions, so
# every time-invariant individual characteristic is absorbed by construction.

message("\n=== BLOCK 5. H1a ===")
set.seed(SEED_MODELS)

# -----------------------------------------------------------------------------
# M1. Conditional logit with person fixed effects
# -----------------------------------------------------------------------------
# Each respondent constitutes their own stratum. Only the discordant pairs
# identify the coefficient; the respondents who answered both questions in the
# same way contribute nothing, because their individual propensity to quantify
# already accounts for both answers.
#
# This design disarms the objection that V35 was asked without filtering on
# whether the respondent has anyone to care for. Whatever that circumstance may
# be, it is the same when the person answers about domestic work and about
# care, so it cannot generate a within-person difference.

m_clogit <- clogit(cant_say ~ task + strata(id), data = h1a_long)

h1a_main <- tidy_ci(coef(m_clogit)[1], sqrt(vcov(m_clogit)[1, 1])) %>%
  mutate(term       = "Care work (ref: domestic work)",
         odds_ratio = round(exp(estimate), 3),
         or_low     = round(exp(ci_low), 3),
         or_high    = round(exp(ci_high), 3),
         n_persons  = n_distinct(h1a_long$id)) %>%
  select(term, estimate, se, odds_ratio, or_low, or_high, p, n_persons)

save_table(h1a_main, "H1a", "H1a_M1_conditional_logit",
           paste("Conditional logit with person fixed effects: the odds of being unable",
                 "to quantify care relative to domestic work. Identified by the",
                 "discordant pairs alone; concordant respondents contribute nothing."))

# The same model in stargazer's regression layout, which reports the fitted
# object directly rather than a table I have assembled myself.
save_model_table(m_clogit, "H1a", "H1a_M1_conditional_logit_model",
                 paste("Conditional logit of the inability to quantify, on task, with the",
                       "respondent as stratum (survival::clogit)."),
                 dep.var.labels = "Cannot quantify weekly hours",
                 covariate.labels = "Care work (ref: domestic work)")

# -----------------------------------------------------------------------------
# M1b. Sensitivity to the survey weight
# -----------------------------------------------------------------------------
# clogit() offers method = "approximate" for weighted fits, but that
# approximation is designed for ties in survival models and attenuates
# matched-pair estimates severely: on these data it returns an odds ratio of
# 1.49 instead of 3.33 even with no weights at all, so it is the approximation
# and not the weighting that breaks down. I therefore do not use it.
#
# With two observations per person and a binary predictor the exact conditional
# logit has a closed form: the coefficient is the log of the ratio of the two
# discordant counts. The weighted analogue therefore sums the weights of the
# discordant respondents instead of counting them, so that each person enters
# according to what they represent in the population.
#
# The interval comes from a bootstrap over discordant respondents. Applying the
# McNemar formula to summed weights would treat the weighted total as a count
# of people and ignore the variability of the weights themselves, which yields
# a standard error roughly a third too small here.

disc_pairs <- h1a_long %>%
  select(id, w, task, cant_say) %>%
  pivot_wider(names_from = task, values_from = cant_say) %>%
  filter(domestic != care) %>%
  mutate(type = if_else(domestic == 0 & care == 1, "b", "c"))

w_b <- sum(disc_pairs$w[disc_pairs$type == "b"])
w_c <- sum(disc_pairs$w[disc_pairs$type == "c"])

boot_or <- replicate(N_BOOT_WEIGHTED, {
  s  <- disc_pairs[sample.int(nrow(disc_pairs), nrow(disc_pairs), replace = TRUE), ]
  bb <- sum(s$w[s$type == "b"])
  cc <- sum(s$w[s$type == "c"])
  if (cc == 0) NA_real_ else log(bb / cc)
})
ci_weighted <- exp(quantile(boot_or, c(.025, .975), na.rm = TRUE))

b_unweighted <- sum(disc_pairs$type == "b")
c_unweighted <- sum(disc_pairs$type == "c")
se_unweighted <- sqrt(vcov(m_clogit)[1, 1])

h1a_weighting <- tibble(
  Weighting = c("Unweighted", "PESOFIN"),
  `Quantifies domestic only` = c(b_unweighted, round(w_b, 1)),
  `Quantifies care only`     = c(c_unweighted, round(w_c, 1)),
  `Odds Ratio` = c(round(exp(coef(m_clogit)[1]), 3), round(w_b / w_c, 3)),
  `C.I. 95%` = c(
    fmt_ci_dash(exp(coef(m_clogit)[1] - 1.96 * se_unweighted),
                exp(coef(m_clogit)[1] + 1.96 * se_unweighted)),
    fmt_ci_dash(ci_weighted[1], ci_weighted[2])),
  Interval = c("Model-based",
               paste0("Bootstrap, ", format(N_BOOT_WEIGHTED, big.mark = ","),
                      " replicates")))

save_table(h1a_weighting, "H1a", "TableA3_weighting_sensitivity",
           paste("Capacity to quantify care work relative to domestic work.",
                 "'Quantifies domestic only' are respondents who marked 'Can't say' for",
                 "care but gave hours for domestic work; 'Quantifies care only' is the",
                 "reverse. The conditional logit does not admit survey weights, so the",
                 "weighted row is computed on the discordant pairs. Respondents who can",
                 "quantify care but not domestic work are under-represented in the",
                 "sample, which attenuates the estimate; the intervals overlap and the",
                 "conclusion is unchanged."))

# -----------------------------------------------------------------------------
# M2b. Conditional logit estimated separately by sex
# -----------------------------------------------------------------------------
# Estimated because the claim that the inability to quantify care is a property
# of women rather than of the task is common in the literature. Each estimate
# rests on roughly one hundred and twenty discordant pairs, so the comparison
# between the two is underpowered and is reported as description rather than as
# a formal contrast.

h1a_stratified <- map_dfr(c(0, 1), function(fem) {
  d <- h1a_long %>% filter(female == fem)
  m <- clogit(cant_say ~ task + strata(id), data = d)
  disc_sex <- d %>%
    select(id, task, cant_say) %>%
    pivot_wider(names_from = task, values_from = cant_say)
  tibble(Sex               = if_else(fem == 1, "Women", "Men"),
         `Odds Ratio`      = round(exp(coef(m)[1]), 3),
         `Standard Error`  = round(sqrt(vcov(m)[1, 1]), 3),
         p                 = fmt_p(summary(m)$coefficients[1, 5]),
         `Discordant pairs` = sum(disc_sex$domestic != disc_sex$care))
})

save_table(h1a_stratified, "H1a", "TableB2_conditional_logit_by_sex",
           paste("Conditional logit estimated separately for women and men. The last",
                 "column reports the number of discordant pairs on which each estimate",
                 "rests, not the size of the subsample. The two estimates are practically",
                 "identical and their intervals overlap completely: the irreducibility of",
                 "care to clock time is a property of the task rather than of who",
                 "performs it."))


# =============================================================================
# BLOCK 6. H1b — SELF-REPORTS VERSUS PROXY-INFORMANT REPORTS
# =============================================================================
# H1b states that the hours of reproductive labour attributed to women are
# lower when reported by a male proxy informant than when women report their
# own. Only one member of each household is interviewed, so this is not a
# within-household paired comparison: self-reports and proxy reports describe
# different people, and the contrast is adjusted for the characteristics of the
# target person rather than identified by matching.

message("\n=== BLOCK 6. H1b ===")
set.seed(SEED_MODELS)

H1B_COVARIATES <- c("target_age", "target_educ", "target_works", "has_minor")
h1b_model_data <- h1b_long %>% drop_na(all_of(H1B_COVARIATES), hours_w)

# -----------------------------------------------------------------------------
# M3. Adjusted gap in each of the four cells
# -----------------------------------------------------------------------------
# The two control cells are what identify the hypothesis. A generic informant
# deficit would deflate all four cells alike; it does not.
#
# ON WEIGHTING. These models are estimated unweighted, and the decision is a
# considered one rather than an omission. PESOFIN calibrates the sample of
# RESPONDENTS to the population margins of sex, age and education. Half of the
# rows of this file, however, do not describe a respondent: they describe the
# respondent's partner, a person who was never sampled and whose own
# probability of selection the weight says nothing about. Applying the
# respondent's calibration factor to a row describing somebody else would not
# make that row representative of anything. The quantity estimated here is
# moreover a conditional contrast between two ways of describing comparable
# people, not a population mean, and the covariates on which the comparison is
# adjusted are the ones the weight itself is built from. The weighted analogue
# is nonetheless reported below as a sensitivity check, exactly as in H1a.

cell_gap <- function(df, yvar, weighted = FALSE) {
  f <- reformulate(c("source", H1B_COVARIATES), response = yvar)
  m <- if (weighted) lm(f, data = df, weights = df$w) else lm(f, data = df)
  V <- sandwich::vcovHC(m, "HC1")
  b <- -unname(coef(m)["sourceproxy"])   # sign flipped: self minus proxy
  s <- sqrt(V["sourceproxy", "sourceproxy"])
  tibble(gap = b, se = s, p = 2 * pnorm(-abs(b / s)), n = nobs(m))
}

h1b_cells <- h1b_model_data %>%
  group_by(target = if_else(target_female == 1, "Woman", "Man"), task) %>%
  group_modify(~ cell_gap(.x, "hours_w")) %>%
  ungroup()

h1b_cells_display <- h1b_cells %>%
  transmute(`Target person` = target,
            Task            = recode(as.character(task),
                                     domestic = "Domestic work",
                                     care     = "Care work"),
            `Gap (hours)`   = sprintf("%+.2f", gap),
            `S.E.`          = round(se, 2),
            `95 % CI`       = fmt_ci_dash(gap - 1.96 * se, gap + 1.96 * se),
            p               = fmt_p(p),
            N               = n)

save_table(h1b_cells_display, "H1b", "TableA4_adjusted_gaps",
           paste("Adjusted self-minus-proxy gap in weekly hours of reproductive labour,",
                 "controlling for the target person's age, education and employment and",
                 "for the presence of minors in the household. Positive values indicate",
                 "that the target person reports more hours for themselves than their",
                 "cohabiting partner reports for them.",
                 "Heteroskedasticity-robust standard errors. Estimated unweighted;",
                 "the weighted analogue is reported as a sensitivity check."))

# -----------------------------------------------------------------------------
# M3c. Sensitivity to the survey weight
# -----------------------------------------------------------------------------
# The same four adjusted gaps refitted with PESOFIN. The weight belongs to the
# respondent, so in a proxy row it is being applied to a description of the
# partner; that is the reason it is not used in the primary estimates. The
# check establishes that the choice does not carry the result, which is what
# the reader needs to know.

h1b_cells_weighted <- h1b_model_data %>%
  group_by(target = if_else(target_female == 1, "Woman", "Man"), task) %>%
  group_modify(~ cell_gap(.x, "hours_w", weighted = TRUE)) %>%
  ungroup()

h1b_weighting <- h1b_cells %>%
  select(target, task, gap_u = gap, se_u = se, p_u = p, n) %>%
  left_join(h1b_cells_weighted %>%
              select(target, task, gap_w = gap, se_w = se, p_w = p),
            by = c("target", "task")) %>%
  transmute(`Target person`     = target,
            Task                = recode(as.character(task),
                                         domestic = "Domestic work",
                                         care     = "Care work"),
            `Unweighted gap`      = sprintf("%+.2f", gap_u),
            `95 % CI (unweighted)` = fmt_ci_dash(gap_u - 1.96 * se_u,
                                                 gap_u + 1.96 * se_u),
            `p (unweighted)`      = fmt_p(p_u),
            `PESOFIN gap`         = sprintf("%+.2f", gap_w),
            `95 % CI (PESOFIN)`   = fmt_ci_dash(gap_w - 1.96 * se_w,
                                                gap_w + 1.96 * se_w),
            `p (PESOFIN)`         = fmt_p(p_w),
            N                     = n)

save_table(h1b_weighting, "H1b", "TableB3b_weighting_sensitivity",
           paste("The four adjusted self-minus-proxy gaps estimated unweighted, as in the",
                 "main analysis, and refitted with the final post-stratification weight",
                 "(PESOFIN). The weight calibrates the sample of respondents, whereas half",
                 "of the rows of this file describe the respondent's partner, who was",
                 "never sampled; the unweighted estimate is therefore the primary one and",
                 "this table documents that the conclusion does not depend on that",
                 "choice."))

# -----------------------------------------------------------------------------
# M3b. Triple interaction, clustered by respondent
# -----------------------------------------------------------------------------
# H1b is not a claim about the independent effect of source, target sex or task
# but about their three-way interaction, which is what this model tests. Each
# respondent supplies up to four rows, so the standard errors are clustered on
# the respondent.

h1b_interaction <- map_dfr(c("hours_w", "hours", "log_hours"), function(y) {
  f  <- reformulate(c("source*target_female*task", H1B_COVARIATES), response = y)
  m  <- lm(f, data = h1b_model_data)
  ct <- lmtest::coeftest(m, vcov = sandwich::vcovCL, cluster = ~ id)
  r  <- ct["sourceproxy:target_female:taskcare", ]
  tibble(outcome = recode(y,
                          hours_w   = "Winsorised at p95",
                          hours     = "Untruncated",
                          log_hours = "log(1 + hours)"),
         estimate = round(r[1], 3), se = round(r[2], 3), p = round(r[4], 4))
})

save_table(h1b_interaction, "H1b", "TableB4_triple_interaction",
           paste("Triple interaction of source, target sex and task, with standard errors",
                 "clustered on the respondent. The estimate is directionally consistent",
                 "across all three specifications and significant in none: the four-cell",
                 "pattern carries the argument and the formal interaction is",
                 "underpowered."))

# -----------------------------------------------------------------------------
# M6. Sensitivity to the specification
# -----------------------------------------------------------------------------
# The same gap estimated four ways, to establish whether it depends on the
# treatment of the highly skewed distribution of reported hours.

h1b_sensitivity <- h1b_model_data %>%
  group_by(target = if_else(target_female == 1, "Woman", "Man"), task) %>%
  group_modify(function(d, k) {
    a <- cell_gap(d, "hours_w")
    b <- cell_gap(d, "hours")
    l <- cell_gap(d, "log_hours")
    qm <- rq(reformulate(c("source", H1B_COVARIATES), response = "hours"),
             tau = .5, data = d)
    qs <- summary(qm, se = "boot", R = 300)$coefficients
    tibble(capped      = paste0(round(a$gap, 2), " (p=", round(a$p, 3), ")"),
           untruncated = paste0(round(b$gap, 2), " (p=", round(b$p, 3), ")"),
           log_scale   = paste0(round(l$gap, 3), " (p=", round(l$p, 3), ")"),
           median_reg  = paste0(round(-qs["sourceproxy", "Value"], 2),
                                " (p=", round(qs["sourceproxy", "Pr(>|t|)"], 3), ")"),
           n = a$n)
  }) %>%
  ungroup()

h1b_sensitivity <- h1b_sensitivity %>%
  mutate(task = recode(as.character(task),
                       domestic = "Domestic", care = "Care")) %>%
  rename(`Target person` = target, Task = task,
         `Capped p95` = capped, Untruncated = untruncated,
         `Log scale` = log_scale, `Median reg` = median_reg, N = n)

save_table(h1b_sensitivity, "H1b", "TableB3_specification_sensitivity",
           paste("The same adjusted gap estimated four ways. The care gap for female",
                 "targets survives untruncated and grows; the domestic gap for male",
                 "targets is the one that holds at the median. The two discrepancies have",
                 "different distributional shapes."))

# -----------------------------------------------------------------------------
# M4. RIF decomposition by quantile
# -----------------------------------------------------------------------------
# The recentred influence function turns a quantile, which is not a mean and so
# cannot be placed on the left-hand side of a regression, into something that
# can be modelled linearly:
#
#     RIF(y; q_tau) = q_tau + (tau - 1{y <= q_tau}) / f(q_tau)
#
# Once the outcome has been transformed, the usual Oaxaca counterfactual
# applies. The reference coefficient vector is the pooled model (Neumark),
# which avoids having to designate one group as the norm.
#
# INTERPRETATION. "Composition" is the part of the gap explained by
# self-reporters and proxied partners being different kinds of people;
# "structure" is what remains once they are matched on the covariates. The
# hypothesis predicts structure. The structural term is a residual and must not
# be read as a causal proxy effect.

rif_quantile <- function(y, tau) {
  q <- as.numeric(quantile(y, tau, na.rm = TRUE))
  f <- density(y, from = q, to = q, n = 1, na.rm = TRUE)$y
  if (!is.finite(f) || f <= 0) return(rep(NA_real_, length(y)))
  q + (tau - as.numeric(y <= q)) / f
}

# Twofold Oaxaca decomposition with a pooled reference vector.
decompose <- function(df, tau) {
  df$rif <- rif_quantile(df$hours, tau)
  if (all(is.na(df$rif))) return(c(total = NA, composition = NA, structure = NA))
  f  <- reformulate(H1B_COVARIATES, response = "rif")
  a  <- df[df$source == "self", ]
  bp <- df[df$source == "proxy", ]
  if (nrow(a) < 30 || nrow(bp) < 30) {
    return(c(total = NA, composition = NA, structure = NA))
  }

  m_a <- lm(f, data = a)
  m_b <- lm(f, data = bp)
  m_p <- lm(f, data = df)
  x_a <- colMeans(model.matrix(m_a))
  x_b <- colMeans(model.matrix(m_b))
  nm  <- union(names(coef(m_p)), union(names(x_a), names(x_b)))
  fill <- function(v) {
    out <- setNames(rep(0, length(nm)), nm)
    out[names(v)] <- v
    out
  }
  b_a <- fill(coef(m_a)); b_b <- fill(coef(m_b)); b_p <- fill(coef(m_p))
  x_a <- fill(x_a);       x_b <- fill(x_b)

  comp <- sum((x_a - x_b) * b_p)
  strc <- sum(x_a * (b_a - b_p)) + sum(x_b * (b_p - b_b))
  c(total = comp + strc, composition = comp, structure = strc)
}

# Cluster bootstrap: respondents, not rows, are resampled. The row indices are
# pre-split by respondent so that a replicate is a single vector concatenation
# rather than one subset per respondent, and each replicate is reused across
# all five quantiles instead of being drawn afresh for each of them.
boot_cell <- function(df, taus, b_reps = N_BOOT_RIF) {
  idx <- split(seq_len(nrow(df)), df$id)
  k   <- length(idx)
  reps <- vapply(seq_len(b_reps), function(b) {
    rows <- unlist(idx[sample.int(k, k, replace = TRUE)], use.names = FALSE)
    d <- df[rows, ]
    vapply(taus, function(tau) {
      tryCatch(decompose(d, tau),
               error = function(e) c(NA_real_, NA_real_, NA_real_))
    }, numeric(3))
  }, matrix(0, nrow = 3, ncol = length(taus)))
  # reps has dimensions 3 x length(taus) x b_reps
  apply(reps, c(1, 2), sd, na.rm = TRUE)
}

RIF_QUANTILES <- c(.10, .25, .50, .75, .90)

message("  RIF decomposition with ", N_BOOT_RIF, " cluster-bootstrap replicates...")

rif_tab <- h1b_model_data %>%
  mutate(target = if_else(target_female == 1, "Woman", "Man")) %>%
  group_by(target, task) %>%
  group_split() %>%
  map_dfr(function(d) {
    key <- d %>% distinct(target, task)
    message("    ", key$target, " / ", key$task, "  (n = ", nrow(d), ")")
    est <- vapply(RIF_QUANTILES, function(tau) decompose(d, tau), numeric(3))
    se  <- boot_cell(d, RIF_QUANTILES)
    tibble(cell     = paste0(key$target, " - ", as.character(key$task)),
           quantile = RIF_QUANTILES,
           composition    = round(est[2, ], 2),
           composition_se = round(se[2, ], 2),
           composition_p  = round(2 * pnorm(-abs(est[2, ] / se[2, ])), 3),
           structure      = round(est[3, ], 2),
           structure_se   = round(se[3, ], 2),
           structure_p    = round(2 * pnorm(-abs(est[3, ] / se[3, ])), 3))
  })

rif_tab_display <- rif_tab %>%
  mutate(cell = str_replace(cell, " - domestic", " / Domestic"),
         cell = str_replace(cell, " - care",     " / Care"))

save_table(rif_tab_display, "H1b", "TableB5_rif_decomposition",
           paste("RIF decomposition of the self-minus-proxy gap at five quantiles.",
                 "'Composition' is the part of the gap explained by differences between",
                 "self-reporters and proxied partners in the target person's age,",
                 "education and employment; 'structure' is the residual reporting",
                 "difference. Pooled reference coefficients; standard errors from",
                 "cluster-bootstrap replicates resampling respondents; p-values are Wald",
                 "tests based on those standard errors."))

# -----------------------------------------------------------------------------
# Figure 1. The decomposition across the distribution
# -----------------------------------------------------------------------------
# The bands are plus or minus 1.96 standard errors, matching the Wald p-values
# reported in the table. In the upper tail of the woman-care cell the bootstrap
# distribution is right-skewed, so the symmetric band there is indicative
# rather than exact.

fig_rif <- rif_tab %>%
  pivot_longer(c(composition, structure),
               names_to = "component", values_to = "value") %>%
  mutate(se = if_else(component == "composition", composition_se, structure_se),
         lo = value - 1.96 * se,
         hi = value + 1.96 * se,
         component = recode(component,
                            composition = "Composition",
                            structure   = "Structure (reporting)"),
         cell = factor(cell, levels = c("Man - domestic", "Man - care",
                                        "Woman - domestic", "Woman - care"))) %>%
  ggplot(aes(quantile, value, colour = component, fill = component)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_line(linewidth = .8) +
  geom_point(size = 1.6) +
  facet_wrap(~ cell, scales = "free_y") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("Composition" = "grey45",
                                 "Structure (reporting)" = "#1f4e79")) +
  scale_fill_manual(values = c("Composition" = "grey45",
                               "Structure (reporting)" = "#1f4e79")) +
  labs(x = "Quantile of reported weekly hours",
       y = "Contribution to the gap (hours)", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

save_figure(fig_rif, "H1b", "Figure1_rif_by_quantile",
            paste("Decomposition of the self-minus-proxy gap across the distribution of",
                  "reported hours. The vertical scales differ across panels: the",
                  "woman-care cell reaches twenty-two hours while the other three lie",
                  "between zero and three, and a common scale would flatten them out of",
                  "view."))


# =============================================================================
# BLOCK 7. H2 — FAMILY-TO-WORK INTERFERENCE, SEX AND MATERIAL POSITION
# =============================================================================
# H2 states that at equal working hours women perceive greater interference of
# reproductive labour in paid work, and that this worsens with lower
# equivalised income. The models are estimated on each of the twenty imputed
# datasets and pooled with Rubin's rules, and the mean sex gap is additionally
# decomposed with a twofold Oaxaca-Blinder procedure.

message("\n=== BLOCK 7. H2 ===")
set.seed(SEED_MODELS)

# -----------------------------------------------------------------------------
# 7.1 Pooling machinery
# -----------------------------------------------------------------------------
# Rubin's rules are implemented directly rather than through a package, so that
# the fraction of missing information is available for reporting and so that
# the pooling itself is auditable:
#   Q  = mean of the m estimates
#   U  = mean of the m sampling variances       (within-imputation variance)
#   Bv = variance of the m estimates            (between-imputation variance)
#   T  = U + (1 + 1/m) Bv
# FMI is the share of the total variance contributed by the missing data.

pool_rubin <- function(est, var) {
  m  <- length(est)
  q  <- mean(est)
  u  <- mean(var)
  bv <- if (m > 1) var(est) else 0
  tv <- u + (1 + 1 / m) * bv
  r  <- (1 + 1 / m) * bv / u
  df <- (m - 1) * (1 + 1 / r)^2
  tibble(estimate = q, se = sqrt(tv),
         ci_low  = q - qt(.975, df) * sqrt(tv),
         ci_high = q + qt(.975, df) * sqrt(tv),
         p   = 2 * pt(-abs(q / sqrt(tv)), df),
         fmi = (r + 2 / (df + 3)) / (r + 1))
}

# Fit a model on every imputed dataset and pool its coefficients.
#
# The number of observations each model is actually estimated on is carried
# back as an attribute rather than as a column, because it belongs to the model
# and not to any one coefficient. It is the mean across the imputations: the
# analytic sample varies slightly between completed datasets, since employment
# status is itself one of the imputed variables and the filter that defines the
# sample is applied after imputation.
fit_pooled <- function(data_long, formula, design_fun) {
  imps <- sort(unique(data_long$.imp))
  imps <- imps[imps > 0]
  fits <- map(imps, function(i) design_fun(data_long %>% filter(.imp == i), formula))
  terms_all <- unique(unlist(map(fits, ~ names(coef(.x)))))
  out <- map_dfr(terms_all, function(tm) {
    est  <- map_dbl(fits, ~ unname(coef(.x)[tm]))
    vr   <- map_dbl(fits, ~ unname(diag(vcov(.x))[tm]))
    keep <- !is.na(est) & !is.na(vr)
    pool_rubin(est[keep], vr[keep]) %>% mutate(term = tm, .before = 1)
  }) %>%
    mutate(across(c(estimate, se, ci_low, ci_high), ~ round(.x, 4)),
           p = round(p, 4), fmi = round(fmi, 3))
  attr(out, "n_fitted")     <- mean(map_dbl(fits, nobs))
  attr(out, "n_imputations") <- length(fits)
  out
}

svy_lm <- function(d, f) {
  des <- svydesign(ids = ~1, strata = ~stratum, weights = ~w, data = d, nest = TRUE)
  svyglm(f, design = des)
}

# -----------------------------------------------------------------------------
# 7.2 Analytic sample
# -----------------------------------------------------------------------------
# The completed datasets carry only what the imputation model used: the
# structural variables were passed in as value-plus-indicator pairs and the
# survey-design columns were not needed there. The analysis variables and the
# design information are restored by joining on the respondent identifier.

ANALYSIS_COLS <- c("id", "conflict_fw", "conflict_wf", "work_hours",
                   "partnered", "stratum", "w")

attach_design <- function(d) {
  d %>%
    select(-any_of(setdiff(ANALYSIS_COLS, "id"))) %>%
    left_join(base %>% select(all_of(ANALYSIS_COLS)), by = "id")
}

completed  <- attach_design(completed)
delta_sets <- map(delta_sets, attach_design)

h2_data <- completed %>%
  filter(in_work == 1, !is.na(conflict_fw), !is.na(work_hours)) %>%
  mutate(educ4 = factor(educ4), income_q = factor(income_q))

message("  analytic n per imputation: ",
        round(mean(table(h2_data$.imp[h2_data$.imp > 0]))))

# Exploratory coding check run before estimating anything: higher values of
# econ_strain must mean MORE economic difficulty, and the raw association with
# the outcome must be inspected before the variable enters a model as a control.
cat("\n--- Exploratory check on the coding of economic strain ---\n")
h2_data %>%
  filter(.imp == 1) %>%
  distinct(id, econ_strain) %>%
  count(econ_strain, sort = FALSE) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\n--- Raw association between economic strain and the H2 outcome ---\n")
h2_data %>%
  filter(.imp == 1) %>%
  group_by(econ_strain) %>%
  summarise(n = n(),
            mean_conflict_unweighted = mean(conflict_fw, na.rm = TRUE),
            mean_conflict_weighted   = weighted.mean(conflict_fw, w, na.rm = TRUE),
            .groups = "drop") %>%
  as.data.frame() %>%
  print(row.names = FALSE)

# -----------------------------------------------------------------------------
# 7.3 M7. Pooled models
# -----------------------------------------------------------------------------
# The specifications are introduced sequentially: sex and paid-work hours
# first, then equivalised-income quartile with its interaction with sex and the
# demographic controls, then subjective economic strain as a second indicator
# of material position.

h2_f1 <- conflict_fw ~ female + work_hours
h2_f2 <- conflict_fw ~ female + work_hours + income_q + female:income_q +
  age + educ4 + has_minor + partnered + mode
h2_f3 <- update(h2_f2, . ~ . + econ_strain)

# The three pooled tables are kept as a named list rather than bound together
# straight away, because bind_rows() drops the attributes and the analytic n of
# each specification travels as an attribute of its own table.
h2_fitted <- list(
  "M1 sex and working hours" = fit_pooled(h2_data, h2_f1, svy_lm),
  "M2 + equivalised income"  = fit_pooled(h2_data, h2_f2, svy_lm),
  "M3 + economic strain"     = fit_pooled(h2_data, h2_f3, svy_lm))

H2_MODEL_N <- map_dbl(h2_fitted, ~ round(attr(.x, "n_fitted")))

cat("\n--- H2 analytic n by specification (mean across imputations) ---\n")
print(as.data.frame(H2_MODEL_N))

h2_models <- imap_dfr(h2_fitted, ~ mutate(.x, model = .y)) %>%
  filter(term != "(Intercept)") %>%
  select(model, term, estimate, se, ci_low, ci_high, p, fmi) %>%
  mutate(n = unname(H2_MODEL_N[model]))

save_table(h2_models, "H2", "H2_M7_pooled_models_full",
           paste("Family-to-work interference among employed respondents, complete",
                 "coefficient list. Design-based weighted OLS fitted on each of the",
                 "twenty imputed datasets and pooled with Rubin's rules. FMI is the",
                 "fraction of missing information, that is, the share of each",
                 "coefficient's variance contributed by the imputation. N is the mean",
                 "number of observations across the twenty completed datasets."))

# Panel A as it appears in the manuscript: one column of estimates and one of
# p-values per specification. The demographic controls are collapsed into a
# single Yes/No row, as in the text, because their individual coefficients are
# not part of the argument; the complete list remains available in the table
# exported immediately above.

H2_PANEL_TERMS <- c("female", "work_hours",
                    "female:income_qQ2", "female:income_qQ3", "female:income_qQ4",
                    "econ_strain", "has_minor", "age")

H2_MODEL_LABELS <- c("M1 sex and working hours" = "M1: Sex + working hours",
                     "M2 + equivalised income"  = "M2: + income & controls",
                     "M3 + economic strain"     = "M3: + economic conditions")

h2_panel_column <- function(model_key) {
  label <- unname(H2_MODEL_LABELS[model_key])
  h2_models %>%
    filter(model == model_key, term %in% H2_PANEL_TERMS) %>%
    transmute(term,
              !!label := fmt_est_ci(estimate, ci_low, ci_high),
              !!paste0("p (", sub(":.*", "", label), ")") := fmt_p(p))
}

h2_panel_a <- reduce(map(names(H2_MODEL_LABELS), h2_panel_column),
                     full_join, by = "term") %>%
  mutate(term = factor(term, levels = H2_PANEL_TERMS)) %>%
  arrange(term) %>%
  mutate(term = relabel_terms(as.character(term))) %>%
  mutate(across(everything(), ~ replace_na(.x, "—")))

# The controls row states which specifications carry the demographic block.
h2_controls_row <- tibble(term = "Education, partnership, survey mode")
for (key in names(H2_MODEL_LABELS)) {
  label <- unname(H2_MODEL_LABELS[key])
  present <- any(grepl("^educ4|^partnered$|^mode",
                       h2_models$term[h2_models$model == key]))
  h2_controls_row[[label]] <- if (present) "Yes" else "No"
  h2_controls_row[[paste0("p (", sub(":.*", "", label), ")")]] <- ""
}

# The sample size row. Employment status is one of the imputed variables, so
# the analytic sample is not identical in every completed dataset and the
# figure reported is the mean over the twenty of them.
h2_n_row <- tibble(term = "n")
for (key in names(H2_MODEL_LABELS)) {
  label <- unname(H2_MODEL_LABELS[key])
  h2_n_row[[label]] <- format(unname(H2_MODEL_N[key]), trim = TRUE)
  h2_n_row[[paste0("p (", sub(":.*", "", label), ")")]] <- ""
}

h2_panel_a <- bind_rows(h2_panel_a, h2_controls_row, h2_n_row) %>%
  rename(Predictor = term)

save_table(h2_panel_a, "H2", "PanelA_h2_pooled_models",
           paste("Pooled OLS models of family-to-work interference among employed",
                 "respondents. Cells report the pooled coefficient with its 95 per cent",
                 "confidence interval in brackets. In M2 and M3 the coefficient for Woman",
                 "is the female-male gap within the lowest income quartile, because sex",
                 "is interacted with the income quartile. n is the mean number of",
                 "observations across the twenty completed datasets, which differs",
                 "slightly between them because employment status is itself imputed."))

# -----------------------------------------------------------------------------
# 7.4 M7b. Female-male gap within each equivalised-income quartile
# -----------------------------------------------------------------------------
# In M2 and M3 the coefficient on `female` is the sex gap in Q1 alone, because
# sex is interacted with the income quartile. The quantities that directly test
# the second clause of H2 are the four linear combinations
#
#   Q1: female
#   Q2: female + female:income_qQ2
#   Q3: female + female:income_qQ3
#   Q4: female + female:income_qQ4
#
# Each contrast is computed inside every imputed dataset using the full
# variance-covariance matrix of the coefficients, and only then pooled.

income_gap_vector <- function(m, q, reference_q) {
  b <- coef(m)
  l <- setNames(rep(0, length(b)), names(b))

  if (!"female" %in% names(b)) {
    stop("`female` coefficient not found in the H2 model.")
  }
  l["female"] <- 1

  if (q != reference_q) {
    # Robust to either ordering of the interaction term's name.
    candidates <- c(paste0("female:income_q", q),
                    paste0("income_q", q, ":female"))
    hit <- intersect(candidates, names(b))
    if (length(hit) != 1) {
      stop("Could not uniquely identify the female-by-income interaction for ", q)
    }
    l[hit] <- 1
  }
  l
}

pool_income_sex_gaps <- function(data_long, formula, model_label) {
  imps <- sort(unique(data_long$.imp))
  imps <- imps[imps > 0]
  fits <- map(imps, function(i) svy_lm(data_long %>% filter(.imp == i), formula))

  quartiles   <- levels(data_long$income_q)
  reference_q <- quartiles[1]

  map_dfr(quartiles, function(q) {
    est <- map_dbl(fits, function(m) {
      l <- income_gap_vector(m, q, reference_q)
      sum(l * coef(m))
    })
    # Variance of the linear combination: Var(L'beta) = L' Var(beta) L.
    vr <- map_dbl(fits, function(m) {
      l <- income_gap_vector(m, q, reference_q)
      as.numeric(t(l) %*% vcov(m) %*% l)
    })
    pool_rubin(est, vr) %>%
      mutate(model = model_label, income_quartile = q, .before = 1)
  })
}

h2_income_gaps <- bind_rows(
  pool_income_sex_gaps(h2_data, h2_f2, "M2 + equivalised income"),
  pool_income_sex_gaps(h2_data, h2_f3, "M3 + economic strain")) %>%
  mutate(across(c(estimate, se, ci_low, ci_high), ~ round(.x, 4)),
         p = round(p, 4), fmi = round(fmi, 3))

save_table(h2_income_gaps, "H2", "H2_M7b_gender_gap_by_income_full",
           paste("Pooled female-minus-male contrasts in family-to-work interference",
                 "within each equivalised-income quartile, under both the income model",
                 "and the model that adds economic strain. Each contrast is computed from",
                 "the full coefficient covariance matrix in each imputed dataset and then",
                 "combined with Rubin's rules."))

# Panel B as it appears in the manuscript, reporting the fully adjusted
# specification (M3) alone.
h2_panel_b <- h2_income_gaps %>%
  filter(model == "M3 + economic strain") %>%
  transmute(`Income quartile` = recode(income_quartile,
                                       Q1 = "Q1 - lowest", Q4 = "Q4 - highest"),
            `Female - male gap` = round(estimate, 3),
            `95% CI`            = paste0("[", fmt_num(ci_low, 3), ", ",
                                         fmt_num(ci_high, 3), "]"),
            p                   = fmt_p(p),
            FMI                 = sub("^0", "", formatC(fmi, format = "f", digits = 3)))

save_table(h2_panel_b, "H2", "PanelB_h2_gender_gap_by_income",
           paste("Adjusted female-male gap in family-to-work interference by",
                 "equivalised-income quartile, from linear combinations of the sex and",
                 "sex-by-income coefficients using their full covariance matrix. Q1 is",
                 "the lowest income quartile. FMI is the fraction of missing",
                 "information."))

# -----------------------------------------------------------------------------
# 7.5 M8. Oaxaca-Blinder decomposition of the sex gap
# -----------------------------------------------------------------------------
# Twofold decomposition with a pooled reference coefficient vector (Neumark),
# which avoids having to designate one sex as the norm and sidesteps the choice
# of index that otherwise drives the detailed results.
#
# WHAT THE UNEXPLAINED PART IS NOT. It is not "the effect of sex" and it is not
# a measure of discrimination. It is the residual once the listed covariates
# are equalised: omitted variables, differences in how the question is answered
# and any unmeasured heterogeneity all land in it. I report it as the part not
# explained by the covariates included.

OB_COVARIATES <- c("work_hours", "age", "educ4", "has_minor",
                   "partnered", "income_q", "econ_strain")

weighted_colmeans <- function(x, w) {
  if (any(!is.finite(w)) || sum(w) <= 0) {
    stop("Invalid survey weights in the Oaxaca-Blinder calculation.")
  }
  colSums(x * w) / sum(w)
}

ob_once <- function(d) {
  f <- reformulate(OB_COVARIATES, response = "conflict_fw")

  is_woman <- d$female == 1
  is_man   <- d$female == 0
  if (!any(is_woman) || !any(is_man)) {
    stop("Both sexes must be represented in the Oaxaca-Blinder sample.")
  }

  # The model matrix is constructed once on the whole sample so that both
  # groups have exactly the same columns and the same factor coding.
  x  <- model.matrix(f, data = d)
  nm <- colnames(x)

  women <- d[is_woman, , drop = FALSE]
  men   <- d[is_man,   , drop = FALSE]

  m_women  <- lm(f, data = women, weights = women$w)
  m_men    <- lm(f, data = men,   weights = men$w)
  m_pooled <- lm(f, data = d,     weights = d$w)

  # A missing or aliased coefficient would leave the decomposition
  # unidentified for that imputation.
  if (anyNA(coef(m_women)) || anyNA(coef(m_men)) || anyNA(coef(m_pooled))) {
    stop("Aliased coefficient in the Oaxaca-Blinder regression; ",
         "inspect the support of the factors.")
  }

  fill <- function(v) {
    out <- setNames(rep(0, length(nm)), nm)
    hit <- intersect(nm, names(v))
    out[hit] <- v[hit]
    out
  }

  # The covariate means are weighted, since the regressions themselves are
  # weighted with PESOFIN; using unweighted means here would make the
  # decomposition internally inconsistent.
  x_women <- setNames(weighted_colmeans(x[is_woman, , drop = FALSE], d$w[is_woman]), nm)
  x_men   <- setNames(weighted_colmeans(x[is_man,   , drop = FALSE], d$w[is_man]),   nm)

  b_women  <- fill(coef(m_women))
  b_men    <- fill(coef(m_men))
  b_pooled <- fill(coef(m_pooled))

  explained   <- sum((x_women - x_men) * b_pooled)
  unexplained <- sum(x_women * (b_women - b_pooled)) +
                 sum(x_men   * (b_pooled - b_men))
  total       <- explained + unexplained

  # Audit check: with an intercept and coherent weighting, the Oaxaca total
  # must reproduce the direct weighted difference in means.
  direct_gap <- weighted.mean(women$conflict_fw, women$w) -
                weighted.mean(men$conflict_fw,   men$w)

  if (!isTRUE(all.equal(unname(total), unname(direct_gap), tolerance = 1e-7))) {
    warning("The Oaxaca decomposition does not reproduce the direct weighted ",
            "mean gap.", call. = FALSE)
  }

  c(total = total, explained = explained, unexplained = unexplained)
}

# STRATIFIED BOOTSTRAP. There is no primary sampling unit identifier in the CIS
# file, so respondents are resampled within the observed strata rather than
# from the complete sample treated as one undifferentiated pool.
resample_within_strata <- function(d) {
  idx  <- split(seq_len(nrow(d)), d$stratum, drop = TRUE)
  rows <- unlist(
    lapply(idx, function(ii) sample(ii, size = length(ii), replace = TRUE)),
    use.names = FALSE)
  d[rows, , drop = FALSE]
}

ob_boot_var <- function(d, b_reps = N_BOOT_OAXACA) {
  reps <- vapply(seq_len(b_reps), function(b) {
    tryCatch(ob_once(resample_within_strata(d)),
             error = function(e) c(total = NA_real_, explained = NA_real_,
                                   unexplained = NA_real_))
  }, numeric(3))
  apply(reps, 1, var, na.rm = TRUE)
}

message("  Oaxaca-Blinder over ", N_IMPUTATIONS, " imputations x ",
        N_BOOT_OAXACA, " bootstrap replicates...")

imps_h2 <- sort(unique(h2_data$.imp))
imps_h2 <- imps_h2[imps_h2 > 0]

ob_raw <- map(imps_h2, function(i) {
  d <- h2_data %>%
    filter(.imp == i) %>%
    drop_na(all_of(OB_COVARIATES), conflict_fw)
  list(est = ob_once(d), var = ob_boot_var(d))
})

h2_oaxaca <- map_dfr(c("total", "explained", "unexplained"), function(cmp) {
  pool_rubin(map_dbl(ob_raw, ~ .x$est[cmp]), map_dbl(ob_raw, ~ .x$var[cmp])) %>%
    mutate(component = recode(cmp,
                              total       = "Total gap (women - men)",
                              explained   = "Explained by covariates",
                              unexplained = "Not explained by the covariates included"),
           .before = 1)
}) %>%
  mutate(across(c(estimate, se, ci_low, ci_high), ~ round(.x, 4)),
         p = round(p, 4), fmi = round(fmi, 3))

h2_panel_c <- h2_oaxaca %>%
  mutate(`Approx. share` = paste0(
    round(100 * estimate / estimate[component == "Total gap (women - men)"]), "%")) %>%
  transmute(Component     = component,
            Estimate      = estimate,
            `95% CI`      = paste0("[", fmt_num(ci_low, 3), ", ",
                                   fmt_num(ci_high, 3), "]"),
            p             = fmt_p(p),
            `Approx. share`)

save_table(h2_panel_c, "H2", "PanelC_h2_oaxaca_blinder",
           paste("Twofold Oaxaca-Blinder decomposition of the female-male gap in",
                 "family-to-work interference, with pooled reference coefficients.",
                 "Estimated on each imputed dataset with bootstrap variances and pooled",
                 "with Rubin's rules. The unexplained term is a residual, not an estimate",
                 "of a sex effect."))

# -----------------------------------------------------------------------------
# 7.6 M10. Delta sensitivity of the income specification
# -----------------------------------------------------------------------------
# Code 77 is an explicit refusal, so MAR cannot be verified. The income clause
# of H2 can only be asserted if it holds in all three scenarios.

delta_rows <- map_dfr(names(delta_sets), function(nm) {
  d <- delta_sets[[nm]] %>%
    filter(in_work == 1, !is.na(conflict_fw), !is.na(work_hours)) %>%
    mutate(educ4 = factor(educ4), income_q = factor(income_q))
  fit_pooled(d, h2_f2, svy_lm) %>%
    filter(grepl("^female$|income_q", term)) %>%
    mutate(scenario = paste0("delta ", if_else(nm == "minus", "-0.5", "+0.5")))
})

h2_delta <- bind_rows(
  fit_pooled(h2_data, h2_f2, svy_lm) %>%
    filter(grepl("^female$|income_q", term)) %>%
    mutate(scenario = "baseline"),
  delta_rows) %>%
  select(scenario, term, estimate, se, ci_low, ci_high, p)

h2_delta_display <- h2_delta %>%
  filter(term == "female" | grepl("^female:income_q", term)) %>%
  mutate(term = if_else(term == "female", "Woman, Q1", relabel_terms(term)),
         cell = fmt_est_ci(estimate, ci_low, ci_high),
         scenario = recode(scenario,
                           baseline      = "Baseline",
                           `delta -0.5`  = "Income imputations -0.5 band",
                           `delta +0.5`  = "Income imputations +0.5 band")) %>%
  select(term, scenario, cell) %>%
  pivot_wider(names_from = scenario, values_from = cell) %>%
  select(Term = term, Baseline,
         `Income imputations -0.5 band`, `Income imputations +0.5 band`)

save_table(h2_delta_display, "H2", "TableB6_delta_sensitivity",
           paste("The income specification re-estimated with the imputed bands shifted",
                 "half a band in each direction. Code 77 is an explicit refusal, so MAR",
                 "cannot be verified; the income clause of H2 can only be asserted if it",
                 "holds in all three scenarios."))


# =============================================================================
# BLOCK 8. H3a AND H3b — DISCOURSE, SEX AND PERCEIVED PROPORTIONALITY
# =============================================================================
# H3a states that the naturalising discourse operates, in both sexes, as an
# incremental factor in perceiving one's own contribution as proportionate.
# H3b states that being male increases that perception independently of the
# discourse. Both are tested on the same dependent variable and the same
# population: respondents with a cohabiting partner who place their own
# contribution at the midpoint of V45.
#
# H3 uses no income variable, so it is estimated on the observed data with
# listwise deletion rather than on the imputed datasets.

message("\n=== BLOCK 8. H3a and H3b ===")
set.seed(SEED_MODELS)

h3 <- base %>%
  filter(partnered == 1, !is.na(fair_share)) %>%
  mutate(educ4 = droplevels(factor(educ4)))

des_h3 <- svydesign(ids = ~1, strata = ~stratum, weights = ~w, data = h3, nest = TRUE)

# The models are estimated by listwise deletion, so the eligible population and
# the analytic sample are not the same number and the difference has to be
# accounted for rather than left for the reader to notice. This table records
# where each case is lost, item by item, and is the source of the sentence in
# the text that explains the drop.
H3_MODEL_VARS <- c("conseq_index", "roles_index", "task_load",
                   "mental_load", "age", "educ4", "outsourced", "has_minor")

h3_attrition <- tibble(
  variable = H3_MODEL_VARS,
  n_missing = map_int(H3_MODEL_VARS, ~ sum(is.na(h3[[.x]]))))

h3_n_eligible <- nrow(h3)
h3_n_complete <- sum(complete.cases(h3[, H3_MODEL_VARS]))

cat("\n--- H3 analytic sample ---\n")
cat(sprintf("Partnered respondents with a valid answer to V45: %d\n", h3_n_eligible))
print(as.data.frame(h3_attrition), row.names = FALSE)
cat(sprintf("Complete on every model variable: %d (%.1f%% of the eligible cases)\n",
            h3_n_complete, 100 * h3_n_complete / h3_n_eligible))

h3_attrition_display <- h3_attrition %>%
  transmute(Variable  = variable,
            Missing   = n_missing,
            `N valid` = h3_n_eligible - n_missing) %>%
  add_row(Variable = "Eligible (partnered, valid V45)",
          Missing  = NA_integer_,
          `N valid` = h3_n_eligible, .before = 1) %>%
  add_row(Variable = "Complete on every model variable",
          Missing  = NA_integer_,
          `N valid` = h3_n_complete)

save_table(h3_attrition_display,
           "H3a", "H3_sample_attrition",
           paste("Item non-response on each variable entering the H3 models, among",
                 "partnered respondents who gave a valid answer to V45. The models are",
                 "estimated by listwise deletion, so this table accounts for the",
                 "difference between the eligible cases and the number of observations",
                 "the models report."))

# -----------------------------------------------------------------------------
# 8.1 M11. Series of design-weighted logits
# -----------------------------------------------------------------------------
# M1 enters the two discourse subscales, sex and the declared task load. M2
# adds age, education in four levels, the mental-load item, outsourcing and the
# presence of a minor. M3 interacts each subscale with sex, which is what tests
# the "in both sexes" clause of H3a. M4 adds a quadratic term in declared load,
# and M5 interacts both load terms with sex.
#
# The quadratic term is not decoration: estimated separately by sex, the linear
# load coefficient takes opposite signs, and these are the two arms of a single
# inverted U.

h3_g1 <- fair_share ~ female + conseq_index + roles_index + task_load
h3_g2 <- update(h3_g1, . ~ . + mental_load + age + educ4 + outsourced + has_minor)
h3_g3 <- update(h3_g2, . ~ . + female:conseq_index + female:roles_index)
h3_g4 <- update(h3_g2, . ~ . + I(task_load^2))
h3_g5 <- update(h3_g4, . ~ . + female:task_load + female:I(task_load^2))

fit_logit <- function(f) svyglm(f, design = des_h3, family = quasibinomial())

tidy_svy <- function(m, label) {
  s <- summary(m)$coefficients
  tibble(model = label, term = rownames(s),
         estimate   = round(s[, 1], 3),
         se         = round(s[, 2], 3),
         odds_ratio = round(exp(s[, 1]), 3),
         ci_low     = round(exp(s[, 1] - 1.96 * s[, 2]), 3),
         ci_high    = round(exp(s[, 1] + 1.96 * s[, 2]), 3),
         p          = round(s[, 4], 4)) %>%
    filter(term != "(Intercept)")
}

h3_specifications <- list("M1 discourse and load" = h3_g1,
                          "M2 + controls"         = h3_g2,
                          "M3 + discourse x sex"  = h3_g3,
                          "M4 + quadratic load"   = h3_g4,
                          "M5 + quadratic x sex"  = h3_g5)

h3_fits <- imap(h3_specifications, ~ fit_logit(.x))

h3_all <- imap_dfr(h3_fits, ~ tidy_svy(.x, .y)) %>%
  mutate(n = map_dbl(model, ~ nobs(h3_fits[[.x]])))

save_table(h3_all, "H3a", "H3_M11_logit_series_full",
           paste("Series of design-weighted logits on perceiving one's own contribution",
                 "to household work as proportionate, among respondents with a cohabiting",
                 "partner, complete coefficient list for all five specifications."))

# Panel A as it appears in the manuscript: the two specifications that carry
# the argument, M2 and M4, each with the coefficient and its standard error,
# the odds ratio with its interval, and the p-value. The remaining
# specifications are in the full table above and, for M5, in Panel B and
# Figure 2.

H3_PANEL_TERMS <- c("female", "conseq_index", "roles_index",
                    "task_load", "I(task_load^2)", "mental_load")

h3_panel_column <- function(model_key, prefix) {
  h3_all %>%
    filter(model == model_key, term %in% H3_PANEL_TERMS) %>%
    transmute(term,
              !!paste0(prefix, " b (SE)")      := fmt_b_se(estimate, se),
              !!paste0(prefix, " OR [95 % CI]") := fmt_or_ci(odds_ratio,
                                                             ci_low, ci_high),
              !!paste0(prefix, " p")           := fmt_p(p))
}

h3_panel_a <- full_join(h3_panel_column("M2 + controls", "M2"),
                        h3_panel_column("M4 + quadratic load", "M4"),
                        by = "term") %>%
  mutate(term = factor(term, levels = H3_PANEL_TERMS)) %>%
  arrange(term) %>%
  mutate(term = relabel_terms(as.character(term))) %>%
  mutate(across(everything(), ~ replace_na(.x, "—")))

h3_panel_a <- bind_rows(
  h3_panel_a,
  tibble(term = "n",
         `M2 b (SE)` = as.character(nobs(h3_fits[["M2 + controls"]])),
         `M2 OR [95 % CI]` = "", `M2 p` = "",
         `M4 b (SE)` = as.character(nobs(h3_fits[["M4 + quadratic load"]])),
         `M4 OR [95 % CI]` = "", `M4 p` = "")) %>%
  rename(Term = term)

save_table(h3_panel_a, "H3a", "PanelA_h3_logit_coefficients",
           paste("Design-weighted logits on perceiving one's own contribution as",
                 "proportionate. M2 adds age, education, the mental-load item,",
                 "outsourcing and the presence of a minor to the discourse subscales,",
                 "sex and declared task load; M4 adds a quadratic term in declared load,",
                 "whose linear and squared coefficients are read jointly through the",
                 "shape of the predicted relationship rather than separately."))

# The same series in stargazer's regression layout, which is the form used for
# the coefficient panel reported in the text.
save_model_table(h3_fits, "H3a", "H3_M11_logit_series_model",
                 paste("Design-weighted logistic regressions on perceiving one's own",
                       "contribution as proportionate (survey::svyglm, quasibinomial)."),
                 dep.var.labels = "Perceives own contribution as proportionate",
                 column.labels = names(h3_fits),
                 model.numbers = FALSE)

# -----------------------------------------------------------------------------
# 8.1b Appendix table B10: M5 in full, with the joint test of the interactions
# -----------------------------------------------------------------------------
# M5 is the specification plotted in Figure 2, and the claim the figure carries
# is that the shape of the relationship between declared task load and
# perceived proportionality differs by sex. That claim rests on the two
# interaction terms acting together, so the quantity that tests it is a joint
# Wald test on both, not the two separate p-values: each individual test asks
# whether one term could be dropped with the other retained, which is not the
# hypothesis. I report the full M5 coefficient vector and append the joint test
# as the last row, so that the figure can be checked against the model that
# produced it without reproducing all five specifications.

m5 <- h3_fits[["M5 + quadratic x sex"]]

m5_joint <- survey::regTermTest(m5, ~ female:task_load + female:I(task_load^2))

cat("\n=== Joint test of the sex-by-load interactions in M5 ===\n")
print(m5_joint)

h3_m5_table <- tidy_svy(m5, "M5") %>%
  transmute(Term           = relabel_terms(term),
            `b (SE)`       = fmt_b_se(estimate, se),
            `OR [95 % CI]` = fmt_or_ci(odds_ratio, ci_low, ci_high),
            p              = fmt_p(p))

h3_m5_table <- bind_rows(
  h3_m5_table,
  tibble(Term = "n", `b (SE)` = as.character(nobs(m5)),
         `OR [95 % CI]` = "", p = ""),
  # regTermTest returns an F statistic when the denominator degrees of freedom
  # are finite and a chi-square otherwise, so the label is read off the result
  # rather than assumed.
  {
    if (!is.null(m5_joint$Ftest)) {
      stat_label <- sprintf("Joint Wald test of the two sex-by-load interactions (F(%g, %g))",
                            m5_joint$df, m5_joint$ddf)
      stat_value <- as.numeric(m5_joint$Ftest)
    } else {
      stat_label <- sprintf("Joint Wald test of the two sex-by-load interactions (chi2(%g))",
                            m5_joint$df)
      stat_value <- as.numeric(m5_joint$chisq)
    }
    tibble(Term = stat_label,
           `b (SE)` = formatC(stat_value, format = "f", digits = 2),
           `OR [95 % CI]` = "",
           p = fmt_p(as.numeric(m5_joint$p)))
  })

save_table(h3_m5_table, "H3b", "TableB10_m5_coefficients",
           paste("Full coefficient vector of M5, the specification underlying Figure 2,",
                 "with a joint Wald test of the two sex-by-load interaction terms. The",
                 "joint test is the quantity that supports the figure: it asks whether",
                 "the shape of the relationship between declared task load and perceived",
                 "proportionality differs by sex, which neither interaction term answers",
                 "on its own. The average marginal effect of sex is uninformative in this",
                 "specification, because where the sex difference varies across the range",
                 "of load its average over that range is not the estimand of interest."))

h3_discourse <- h3_all %>%
  filter(grepl("conseq_index|roles_index", term)) %>%
  transmute(Model           = model,
            Term            = relabel_terms(term),
            `b (SE)`        = fmt_b_se(estimate, se),
            `OR [95 % CI]`  = fmt_or_ci(odds_ratio, ci_low, ci_high),
            p               = fmt_p(p))

save_table(h3_discourse, "H3a", "TableB8_discourse_coefficients",
           paste("The two dimensions of the naturalising discourse across the five",
                 "specifications, together with their interactions with sex in M3. The",
                 "interaction terms are what test the 'in both sexes' clause of H3a: if",
                 "they are not distinguishable from zero, the association between",
                 "discourse and perceived proportionality does not differ between women",
                 "and men, and the null result for H3a cannot be attributed to opposite",
                 "effects in the two groups cancelling out."))

# -----------------------------------------------------------------------------
# 8.2 Average marginal effect of sex
# -----------------------------------------------------------------------------
# The raw `female` coefficient is not comparable across these models. In M3 and
# M5 sex is interacted, so the main effect is the log-odds difference at a task
# load of zero and at zero discourse, a point outside the data entirely; in M5
# it comes out as an odds ratio in the hundreds of thousands, which is an
# extrapolation artefact rather than a result. The average marginal effect is
# the quantity that means the same thing in every specification: the average
# change in the predicted probability of perceiving the split as proportionate
# when sex is switched, holding everything else at each respondent's own
# values.

ame_female <- function(m, data, weights) {
  d0 <- d1 <- data
  d0$female <- 0
  d1$female <- 1
  p0 <- predict(m, newdata = d0, type = "response")
  p1 <- predict(m, newdata = d1, type = "response")
  ok <- !is.na(p0) & !is.na(p1)

  # The average marginal effect is a smooth function of the coefficients, so
  # its variance is obtained by simulating from their asymptotic distribution
  # rather than by the delta method.
  v     <- vcov(m)
  bhat  <- coef(m)
  draws <- MASS::mvrnorm(400, bhat, v)
  sims  <- apply(draws, 1, function(bb) {
    m2 <- m
    m2$coefficients <- bb
    q1 <- predict(m2, newdata = d1, type = "response")
    q0 <- predict(m2, newdata = d0, type = "response")
    weighted.mean((q1 - q0)[ok], weights[ok])
  })
  c(ame = weighted.mean((p1 - p0)[ok], weights[ok]), se = sd(sims))
}

h3_sex <- imap_dfr(h3_fits, function(m, label) {
  d    <- model.frame(m)
  rows <- as.integer(rownames(d))
  a    <- ame_female(m, h3[rows, ], h3$w[rows])
  tibble(model   = label,
         ame_pp  = round(100 * a["ame"], 1),
         se_pp   = round(100 * a["se"], 1),
         ci_low  = round(100 * (a["ame"] - 1.96 * a["se"]), 1),
         ci_high = round(100 * (a["ame"] + 1.96 * a["se"]), 1),
         p       = round(2 * pnorm(-abs(a["ame"] / a["se"])), 4),
         n       = nobs(m))
})

save_table(h3_sex, "H3b", "H3b_sex_average_marginal_effect_full",
           paste("Average marginal effect of being a woman on the probability of",
                 "perceiving one's own contribution as proportionate, in percentage",
                 "points, for all five specifications."))

# Panel B as it appears in the manuscript. M5 is included alongside M2 and M4:
# it is the specification underlying Figure 2, and reporting its marginal
# effect is what allows the figure to be read as evidence rather than as an
# illustration, since M5's own coefficients are not interpretable on their own.
h3_panel_b <- h3_sex %>%
  filter(model %in% c("M2 + controls", "M4 + quadratic load",
                      "M5 + quadratic x sex")) %>%
  transmute(Specification = model,
            `AME (percentage points)` = ame_pp,
            SE = se_pp,
            `95 % CI` = paste0("[", fmt_num(ci_low, 1), ", ",
                               fmt_num(ci_high, 1), "]"),
            p = fmt_p(p),
            n = .data$n)

save_table(h3_panel_b, "H3b", "PanelB_h3_average_marginal_effect",
           paste("Average marginal effect of being a woman on the probability of",
                 "perceiving one's own contribution as proportionate. The marginal effect",
                 "is used instead of the coefficient on sex because sex is interacted",
                 "with declared task load in M5, where that coefficient refers to a task",
                 "load of zero, a value the scale cannot take. The M5 row is reported for",
                 "completeness: in that specification the sex difference varies with",
                 "declared task load, and its average over the range of load is not the",
                 "estimand of interest. What M5 establishes is the dependence of the gap",
                 "on load, tested jointly in Table B10 and displayed in Figure 2."))

save_table(h3_all %>%
             filter(term == "female") %>%
             select(model, estimate, se, odds_ratio, p, n),
           "H3b", "H3b_sex_raw_coefficients",
           paste("Raw `female` coefficients, retained for completeness only. In M3 and M5",
                 "they are conditional on a task load of zero and are not interpretable;",
                 "the average marginal effects should be used instead."))

# -----------------------------------------------------------------------------
# 8.3 Models estimated separately by sex
# -----------------------------------------------------------------------------
# Estimated to verify that the pooled coefficients do not conceal opposite
# signs. The load coefficient does exactly that, which is the reason the pooled
# model requires the quadratic term.

h3_stratified <- map_dfr(c(0, 1), function(fem) {
  d  <- h3 %>% filter(female == fem)
  de <- svydesign(ids = ~1, strata = ~stratum, weights = ~w, data = d, nest = TRUE)
  m  <- svyglm(fair_share ~ conseq_index + roles_index + task_load + age + educ4,
               design = de, family = quasibinomial())
  tidy_svy(m, if_else(fem == 1, "Women", "Men")) %>% mutate(n = nobs(m))
})

save_table(h3_stratified, "H3a", "H3_M6_stratified_by_sex",
           paste("The H3 model estimated separately for women and men, as a check that",
                 "the pooled coefficients do not conceal opposite signs. The declared",
                 "task-load coefficient does precisely that, which is what motivates the",
                 "quadratic specification."))

# -----------------------------------------------------------------------------
# 8.4 Figure 2. Predicted probability by sex across declared load
# -----------------------------------------------------------------------------
# The quadratic must not be read outside the range in which each sex actually
# has data. Very few women report a load at or below two and very few men
# report four or above, so the curves are drawn between the fifth and
# ninety-fifth percentile within each sex and the peak - if there is any - is located inside that
# range.

support <- h3 %>%
  filter(!is.na(task_load)) %>%
  group_by(female) %>%
  summarise(lo = quantile(task_load, .05),
            hi = quantile(task_load, .95),
            .groups = "drop")

grid <- support %>%
  group_by(female) %>%
  reframe(task_load = seq(lo, hi, length.out = 60)) %>%
  mutate(conseq_index = weighted.mean(h3$conseq_index, h3$w, na.rm = TRUE),
         roles_index  = weighted.mean(h3$roles_index,  h3$w, na.rm = TRUE),
         mental_load  = weighted.mean(h3$mental_load,  h3$w, na.rm = TRUE),
         age          = weighted.mean(h3$age,          h3$w, na.rm = TRUE),
         educ4        = factor(levels(h3$educ4)[which.max(table(h3$educ4))],
                               levels = levels(h3$educ4)),
         outsourced   = 0,
         has_minor    = 0)

pred_h3 <- predict(h3_fits[["M5 + quadratic x sex"]], newdata = grid,
                   type = "link", se.fit = TRUE)

grid <- grid %>%
  mutate(fit   = as.numeric(pred_h3),
         se    = sqrt(as.numeric(attr(pred_h3, "var"))),
         p_hat = plogis(fit),
         lo    = plogis(fit - 1.96 * se),
         hi    = plogis(fit + 1.96 * se),
         Sex   = factor(female, 0:1, c("Men", "Women")))

fig_h3 <- ggplot(grid, aes(task_load, p_hat, colour = Sex, fill = Sex)) +
  geom_vline(xintercept = 3, linetype = "dashed", linewidth = .3, colour = "grey45") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_line(linewidth = .9) +
  annotate("text", x = 3, y = .02, label = "equal split", hjust = -0.08,
           size = 3, colour = "grey35") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_colour_manual(values = c(Men = "#8c5a2b", Women = "#1f4e79")) +
  scale_fill_manual(values = c(Men = "#8c5a2b", Women = "#1f4e79")) +
  labs(x = paste("Declared share of household tasks",
                 "(1 = partner does all, 5 = respondent does all)"),
       y = "Probability of perceiving the split as proportionate",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

save_figure(fig_h3, "H3b", "Figure2_predicted_probability",
            paste("Predicted probability of perceiving one's own contribution as",
                  "proportionate, by sex, across the declared task load, with the",
                  "covariates held at their weighted means. The horizontal distance",
                  "between the two peaks is the finding: men and women locate the equity",
                  "point differently."))

# Where each curve peaks, within the observed range of each sex.
ranges <- grid %>%
  group_by(Sex) %>%
  summarise(range_lo = round(min(task_load), 2),
            range_hi = round(max(task_load), 2),
            .groups = "drop")

peaks <- grid %>%
  group_by(Sex) %>%
  slice_max(p_hat, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Sex,
            peak_task_load   = round(task_load, 2),
            peak_probability = round(p_hat, 3)) %>%
  left_join(ranges, by = "Sex") %>%
  mutate(at_edge = peak_task_load <= range_lo + 0.05 |
                   peak_task_load >= range_hi - 0.05)

peaks_display <- peaks %>%
  transmute(Sex,
            `Peak declared load`     = peak_task_load,
            `Peak probability`       = peak_probability,
            `Observed range, lower`  = range_lo,
            `Observed range, upper`  = range_hi,
            `Peak at edge of range`  = if_else(at_edge, "Yes", "No"))

save_table(peaks_display, "H3b", "TableB9_equity_point_from_model",
           paste("Declared task load at which the predicted probability of perceiving",
                 "the split as proportionate is highest, by sex, from M5 and restricted",
                 "to the fifth to ninety-fifth percentile of load observed within each",
                 "sex. The last column is a diagnostic: 'Yes' means the maximum falls at",
                 "the boundary of the observed range, so the curve is still rising or",
                 "falling there and its position must not be read as an interior",
                 "optimum or compared as such across sexes."))


# =============================================================================
# BLOCK 9. TABLE MANIFEST AND SESSION INFORMATION
# =============================================================================
# Every table produced above is registered in the manifest, which is written
# out as HTML so that the appendix can be assembled from a complete inventory
# rather than from memory.

message("\n=== BLOCK 9. MANIFEST ===")

cat("\n--- Tables written ---\n")
print(as.data.frame(TABLE_MANIFEST %>% select(section, table, type)),
      row.names = FALSE)

write_html_table(
  as.data.frame(TABLE_MANIFEST),
  paste("Inventory of every table produced by the analysis, with the section it",
        "belongs to and the caption used in the manuscript."),
  MANIFEST_PATH)

cat("\n", nrow(TABLE_MANIFEST), "tables written under",
    file.path(OUT_ROOT, "tables"), "\n")
cat("Manifest:", MANIFEST_PATH, "\n")

# -----------------------------------------------------------------------------
# 9.1 Single file for transfer into the manuscript
# -----------------------------------------------------------------------------
# Every table is also concatenated into one document, in the order of the
# manifest and with its own heading. 

ALL_TABLES_PATH <- file.path(OUT_ROOT, "tables", "_all_tables.html")

combine_tables_html <- function(manifest = TABLE_MANIFEST,
                                path = ALL_TABLES_PATH) {
  m <- manifest %>%
    mutate(section = factor(section, levels = SECTIONS)) %>%
    arrange(section, table)

  parts <- unlist(lapply(seq_len(nrow(m)), function(i) {
    if (!file.exists(m$path[i])) return(character(0))
    c(paste0("<h3>", html_escape(as.character(m$section[i])), " &mdash; ",
             html_escape(m$table[i]), "</h3>"),
      readLines(m$path[i], warn = FALSE),
      "<p>&nbsp;</p>")
  }))

  writeLines(c("<html><head><meta charset='utf-8'></head>",
               "<body style=\"font-family: Times New Roman, serif\">",
               parts, "</body></html>"),
             con = path, useBytes = TRUE)
  invisible(path)
}

combine_tables_html()
cat("All tables in one file:", ALL_TABLES_PATH, "\n")

sessionInfo()


######################################################## END OF SCRIPT

