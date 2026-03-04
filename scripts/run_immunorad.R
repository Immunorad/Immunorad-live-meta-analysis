###############################################################################
# Immunorad live meta-analysis
# - Reads Excel
# - Produces forest plots (PDF + PNG) into repo/figures
# - Writes pooled per-tumour results into repo/results
# - Rewrites repo/index.html to display the latest outputs
###############################################################################

# ---- Install missing packages (optional convenience) ----
pkgs <- c("readxl", "metafor", "pdftools")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(metafor)
library(pdftools)

# =============================================================================
# 0) USER SETTINGS (only edit these 2 paths)
# =============================================================================

# 0a) Path to your Excel (unchanged)
file_path <- "C:/Users/P087325/OneDrive - Amsterdam UMC/Documenten/Immunorad review/Analyses/Data excel + analyses R/Data_immunorad_FSC_PRO_24032025_IO_inclGBM_clean.xlsx"

# 0b) Path to the LOCAL folder where your GitHub repo is on your computer
#     Example: "C:/Users/<you>/Documents/GitHub/Immunorad-live-meta-analysis"
repo_dir  <- "C:/PATH/TO/YOUR/LOCAL/Immunorad-live-meta-analysis"

# Derived output folders inside the repo
fig_dir   <- file.path(repo_dir, "figures")
res_dir   <- file.path(repo_dir, "results")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1) Fixed tumour type order (lowercase)
# =============================================================================

tumor_order <- c(
  "glioblastoma",
  "head and neck",
  "nasopharyngeal",
  "non-small-cell lung",
  "small-cell lung",
  "esophageal",
  "pancreatic",
  "cervical",
  "prostate",
  "sarcoma",
  "cutaneous squamous-cell carcinoma"
)

tumor_pretty <- c(
  "glioblastoma" = "Glioblastoma",
  "head and neck" = "Head and neck",
  "nasopharyngeal" = "Nasopharyngeal",
  "non-small-cell lung" = "NSCLC",
  "small-cell lung" = "SCLC",
  "esophageal" = "Esophageal",
  "pancreatic" = "Pancreatic",
  "cervical" = "Cervical",
  "prostate" = "Prostate",
  "sarcoma" = "Sarcoma",
  "cutaneous squamous-cell carcinoma" = "Cutaneous squamous-cell carcinoma"
)

# =============================================================================
# 2) Load data
# =============================================================================

dat <- read_excel(file_path)
dat$Tumortype_clean <- tolower(trimws(dat$Tumortype))
dat$Tumortype_clean <- factor(dat$Tumortype_clean, levels = tumor_order)

# =============================================================================
# Formatting helpers
# =============================================================================

format_p <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p < 0.001) return("p<0.001")
  paste0("p=", formatC(p, digits = 3, format = "f"))
}

format_ci <- function(lo, hi, digits = 2) {
  paste0("[", round(lo, digits), "-", round(hi, digits), "]")
}

# =============================================================================
# Sample size: STRICTLY from TotalParticipants (robust parsing)
# =============================================================================

get_total_participants <- function(d, colname = "TotalParticipants") {
  nm_low <- tolower(names(d))
  target <- tolower(colname)

  idx <- which(nm_low == target)
  idx <- as.integer(unlist(idx))

  if (length(idx) != 1) {
    stop(
      "Column 'TotalParticipants' not found (or not unique). ",
      "Please ensure the dataset contains exactly one column named 'TotalParticipants'."
    )
  }

  n_raw <- as.character(d[[idx]])
  n_raw <- gsub("\\s+", "", n_raw)
  n_raw <- gsub("\\.", "", n_raw)
  n_raw <- gsub(",", ".", n_raw)
  n <- suppressWarnings(as.numeric(n_raw))

  n[!is.finite(n) | n <= 0] <- NA_real_
  n
}

# =============================================================================
# Helpers: endpoint prep, RE model, heterogeneity (I²)
# =============================================================================

prepare_endpoint_data <- function(d, hr_col, lo_col, hi_col) {
  d <- d[!is.na(d[[hr_col]]) & !is.na(d[[lo_col]]) & !is.na(d[[hi_col]]), ]
  if (nrow(d) == 0) return(NULL)

  d$N_total <- get_total_participants(d, "TotalParticipants")

  d$HR  <- as.numeric(d[[hr_col]])
  d$L   <- as.numeric(d[[lo_col]])
  d$U   <- as.numeric(d[[hi_col]])
  d$yi  <- log(d$HR)
  d$sei <- (log(d$U) - log(d$L)) / (2 * 1.96)

  author <- if ("First Author" %in% names(d)) as.character(d$`First Author`) else NA_character_
  acr    <- if ("Acronym" %in% names(d)) as.character(d$Acronym) else NA_character_

  author[is.na(author) | author == ""] <- "Study"
  acr[is.na(acr) | acr == ""] <- NA_character_

  d$study_label <- ifelse(is.na(acr), author, paste0(author, " (", acr, ")"))
  d
}

re_model <- function(d) rma(yi, sei, data = d, method = "REML")

compute_I2_manual <- function(yi, sei) {
  wi <- 1 / (sei^2)
  mu_hat <- sum(wi * yi) / sum(wi)
  Q <- sum(wi * (yi - mu_hat)^2)
  df <- length(yi) - 1

  tau2 <- max(0, (Q - df) / (sum(wi) - sum(wi^2) / sum(wi)))
  I2 <- ifelse(Q > 0, max(0, ((Q - df) / Q)) * 100, 0)

  list(Q = Q, df = df, tau2 = tau2, I2 = I2)
}

# =============================================================================
# Plot sizing helpers
# =============================================================================

calc_plot_height <- function(n_studies, n_groups, extra_rows = 10,
                             inches_per_row = 0.28, min_h = 10, max_h = 18) {
  h <- (n_studies + n_groups * 2 + extra_rows) * inches_per_row
  h <- max(min_h, min(max_h, h))
  h
}

diamond_height <- function(total_n, ref_n, base_step,
                           min_frac = 0.20, max_frac = 1.60) {
  if (!is.finite(total_n) || total_n <= 0 || !is.finite(ref_n) || ref_n <= 0) {
    return(base_step * 0.60)
  }

  s <- log1p(total_n) / log1p(ref_n)
  s <- pmax(0, pmin(1, s))

  frac <- min_frac + (max_frac - min_frac) * s
  base_step * frac
}

box_psize <- function(n_vec, min_psize = 0.7, max_psize = 1.8) {
  n <- suppressWarnings(as.numeric(n_vec))
  n[!is.finite(n) | n <= 0] <- NA_real_

  if (all(is.na(n))) return(rep(1.0, length(n)))

  max_n <- max(n, na.rm = TRUE)
  if (!is.finite(max_n) || max_n <= 0) return(rep(1.0, length(n)))

  s <- sqrt(n / max_n)
  s[!is.finite(s)] <- 0
  s <- pmax(0, pmin(1, s))

  min_psize + (max_psize - min_psize) * s
}

# =============================================================================
# Forest plotting function
# =============================================================================

forest_with_two_bottom_pools <- function(dat, tumor_order, tumor_pretty,
                                         hr_col, lo_col, hi_col,
                                         main_title,
                                         xlim, at,
                                         base_step = 1.15,
                                         gap_after_pooled = 2.2,
                                         label_offset = 0.85,
                                         pval_x = NULL,
                                         show_excl_gbm = TRUE,
                                         study_col_shift = 0.00) {

  dat$Tumortype_clean <- factor(tolower(trimws(dat$Tumortype)), levels = tumor_order)
  dat <- dat[!is.na(dat$Tumortype_clean), ]

  d <- prepare_endpoint_data(dat, hr_col, lo_col, hi_col)
  if (is.null(d) || nrow(d) == 0) stop("No complete cases for endpoint: ", hr_col)

  d$Tumortype_clean <- factor(as.character(d$Tumortype_clean), levels = tumor_order)
  d <- d[order(d$Tumortype_clean, d$HR), ]

  tumor_types_present <- tumor_order[tumor_order %in% as.character(unique(d$Tumortype_clean))]

  y_current <- 200
  row_pos <- rep(NA_real_, nrow(d))

  pooled_rows <- setNames(rep(NA_real_, length(tumor_types_present)), tumor_types_present)
  has_pooled  <- setNames(rep(FALSE, length(tumor_types_present)), tumor_types_present)

  label_pos <- list()
  group_n_studies <- list()

  for (tt in tumor_types_present) {
    idx <- which(as.character(d$Tumortype_clean) == tt)
    if (length(idx) == 0) next

    n_tt <- length(idx)
    group_n_studies[[tt]] <- n_tt

    rows <- y_current - base_step * (0:(n_tt - 1))
    row_pos[idx] <- rows
    label_pos[[tt]] <- max(rows) + label_offset

    if (n_tt >= 2) {
      pooled_row <- min(rows) - base_step
      pooled_rows[tt] <- pooled_row
      has_pooled[tt]  <- TRUE
      y_current <- pooled_row - gap_after_pooled
    } else {
      y_current <- min(rows) - gap_after_pooled
    }
  }

  d$row_pos <- row_pos

  min_y <- min(d$row_pos, na.rm = TRUE)
  pooled_used <- pooled_rows[has_pooled]
  if (length(pooled_used) > 0) min_y <- min(min_y, pooled_used, na.rm = TRUE)

  bottom_all <- min_y - 2.2

  ylim_upper <- max(unlist(label_pos), na.rm = TRUE) + 2.5

  if (isTRUE(show_excl_gbm)) {
    bottom_nogbm <- bottom_all - 1.8
    ylim_lower <- bottom_nogbm - 2.5
  } else {
    bottom_nogbm <- NA_real_
    ylim_lower <- bottom_all - 2.5
  }

  par(mar = c(4, 2, 2, 4))
  par(xpd = NA)

  textpos_left <- xlim[1] + study_col_shift * diff(xlim)
  textpos_vec  <- c(textpos_left, xlim[2])

  psize_vec <- box_psize(d$N_total)

  forest(
    x = d$yi, sei = d$sei,
    slab = d$study_label,
    xlim = xlim,
    at = at,
    atransf = exp,
    xlab = "",
    rows = d$row_pos,
    refline = log(1),
    ylim = c(ylim_lower, ylim_upper),
    cex = 0.9,
    header = c("Study", "Hazard ratio [95% CI]"),
    textpos = textpos_vec,
    psize = psize_vec
  )

  x_left  <- textpos_left
  x_right <- par("usr")[2]
  if (is.null(pval_x)) pval_x <- x_right

  for (tt in tumor_types_present) {
    n_tt <- group_n_studies[[tt]]
    if (!is.null(n_tt) && n_tt > 0) {
      lab <- if (!is.na(tumor_pretty[tt])) tumor_pretty[tt] else tt
      text(
        x_left, label_pos[[tt]],
        labels = paste0(lab, " (n=", n_tt, ")"),
        pos = 4, font = 2, cex = 0.95
      )
    }
  }

  x_left2  <- par("usr")[1]
  x_right2 <- par("usr")[2]
  y_title  <- par("usr")[4] - 0.90
  title_shift <- 0.06 * (x_right2 - x_left2)

  text(
    x = (x_left2 + x_right2) / 2 + title_shift,
    y = y_title,
    labels = main_title,
    font = 2,
    cex = 1.0
  )

  group_totalN <- setNames(rep(NA_real_, length(tumor_types_present)), tumor_types_present)
  for (tt in tumor_types_present) {
    sub_idx <- which(as.character(d$Tumortype_clean) == tt)
    if (length(sub_idx) > 0) {
      group_totalN[tt] <- sum(d$N_total[sub_idx], na.rm = TRUE)
    }
  }

  totalN_all <- sum(d$N_total, na.rm = TRUE)

  nogbm <- d[as.character(d$Tumortype_clean) != "glioblastoma", ]
  totalN_nogbm <- sum(nogbm$N_total, na.rm = TRUE)

  refN <- max(group_totalN[is.finite(group_totalN)], na.rm = TRUE)
  if (!is.finite(refN) || refN <= 0) refN <- 1

  cat("\n--- Diamond scaling debug (", hr_col, ") ---\n", sep = "")
  print(data.frame(Tumor = names(group_totalN), SumN = as.numeric(group_totalN)))
  cat("refN (max tumour SumN) = ", refN, "\n", sep = "")

  per_tumor_tab <- data.frame(
    Endpoint = character(0),
    Tumor = character(0),
    n = integer(0),
    HR = numeric(0),
    CI_low = numeric(0),
    CI_high = numeric(0),
    p_value = numeric(0),
    p_text = character(0),
    I2 = numeric(0),
    stringsAsFactors = FALSE
  )

  for (tt in tumor_types_present) {
    if (!isTRUE(has_pooled[tt])) next
    sub <- d[as.character(d$Tumortype_clean) == tt, ]
    if (nrow(sub) >= 2) {
      m <- tryCatch(re_model(sub), error = function(e) NULL)
      if (!is.null(m)) {

        thisN  <- group_totalN[tt]
        h_diam <- diamond_height(thisN, refN, base_step)

        cat("Tumor: ", tt, " | SumN=", thisN, " | height=", h_diam, "\n", sep = "")

        addpoly(m, row = pooled_rows[tt], atransf = exp, cex = 0.9, mlab = "", height = h_diam)
        text(x_left, pooled_rows[tt], "Pooled HR", pos = 4, cex = 0.9)

        hr   <- as.numeric(exp(m$b))
        ci_l <- as.numeric(exp(m$ci.lb))
        ci_u <- as.numeric(exp(m$ci.ub))
        pval <- as.numeric(m$pval)
        i2   <- compute_I2_manual(sub$yi, sub$sei)$I2

        per_tumor_tab <- rbind(
          per_tumor_tab,
          data.frame(
            Endpoint = hr_col,
            Tumor = as.character(tt),
            n = nrow(sub),
            HR = hr,
            CI_low = ci_l,
            CI_high = ci_u,
            p_value = pval,
            p_text = format_p(pval),
            I2 = i2,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  m_all <- re_model(d)
  h_all <- diamond_height(totalN_all, refN, base_step)
  cat("ALL studies | SumN=", totalN_all, " | height=", h_all, "\n", sep = "")
  addpoly(m_all, row = bottom_all, atransf = exp, cex = 0.9, mlab = "", height = h_all)

  text(
    x_left, bottom_all,
    paste0("Pooled HR (all studies; n=", nrow(d), ")"),
    pos = 4, font = 2, cex = 0.9
  )
  text(
    pval_x, bottom_all - 0.9,
    labels = format_p(m_all$pval),
    pos = 2, cex = 0.85
  )

  i2_all <- compute_I2_manual(d$yi, d$sei)$I2
  if (!is.na(i2_all)) {
    text(
      x_left, bottom_all - 0.9,
      labels = paste0("I", "\u00B2", "=", formatC(i2_all, digits = 2, format = "f"), "%"),
      pos = 4, cex = 0.85, xpd = TRUE
    )
  }

  if (isTRUE(show_excl_gbm) && nrow(nogbm) >= 2) {
    m_nogbm <- re_model(nogbm)
    h_nogbm <- diamond_height(totalN_nogbm, refN, base_step)
    cat("Excl GBM | SumN=", totalN_nogbm, " | height=", h_nogbm, "\n", sep = "")
    addpoly(m_nogbm, row = bottom_nogbm, atransf = exp, cex = 0.9, mlab = "", height = h_nogbm)

    text(
      x_left, bottom_nogbm,
      paste0("Pooled HR (excluding glioblastoma; n=", nrow(nogbm), ")"),
      pos = 4, font = 2, cex = 0.9
    )
    text(
      pval_x, bottom_nogbm - 0.9,
      labels = format_p(m_nogbm$pval),
      pos = 2, cex = 0.85
    )

    i2_nogbm <- compute_I2_manual(nogbm$yi, nogbm$sei)$I2
    if (!is.na(i2_nogbm)) {
      text(
        x_left, bottom_nogbm - 0.9,
        labels = paste0("I", "\u00B2", "=", formatC(i2_nogbm, digits = 2, format = "f"), "%"),
        pos = 4, cex = 0.85, xpd = TRUE
      )
    }
  }

  x_left2  <- par("usr")[1]
  x_right2 <- par("usr")[2]
  y_favor  <- par("usr")[3] - 1.8
  favor_shift <- 0.17 * (x_right2 - x_left2)

  text(
    x = (x_left2 + x_right2) / 2 + favor_shift,
    y = y_favor,
    labels = expression(bold("Favors RIT") ~ "\u2190" ~ "   HR   " ~ "\u2192" ~ bold("Favors RT")),
    cex = 0.95,
    xpd = TRUE
  )

  invisible(list(
    n_studies = nrow(d),
    per_tumor = per_tumor_tab
  ))
}

# =============================================================================
# Console printing of per-tumour pooled results
# =============================================================================

print_per_tumor_results <- function(res, endpoint_name, tumor_pretty) {
  tab <- res$per_tumor
  if (is.null(tab) || nrow(tab) == 0) {
    cat("\n", "Per-tumour pooled results (", endpoint_name, "): none (no tumour groups with n>=2)\n", sep = "")
    return(invisible(NULL))
  }

  tab$Tumor_pretty <- ifelse(!is.na(tumor_pretty[tab$Tumor]), tumor_pretty[tab$Tumor], tab$Tumor)
  tab$HR_CI <- paste0(
    formatC(tab$HR, digits = 2, format = "f"),
    " ",
    format_ci(tab$CI_low, tab$CI_high, digits = 2)
  )
  tab$I2_text <- ifelse(is.na(tab$I2), "I²=NA", paste0("I²=", formatC(tab$I2, digits = 2, format = "f"), "%"))

  out <- tab[, c("Endpoint","Tumor_pretty", "n", "HR_CI", "p_text", "I2_text")]
  colnames(out) <- c("Endpoint","Tumour type", "n", "Pooled HR [95% CI]", "p-value", "I²")

  cat("\n", "Per-tumour pooled results (", endpoint_name, "):\n", sep = "")
  print(out, row.names = FALSE)
  invisible(out)
}

# =============================================================================
# 3) Run endpoints and save to repo/figures and repo/results
# =============================================================================

all_tables <- list()

# ---- OS ----
tmp_os <- dat
tmp_os <- tmp_os[!is.na(tmp_os$Tumortype_clean) &
                   !is.na(tmp_os$HR_OS) & !is.na(tmp_os$Lower_CI_OS) & !is.na(tmp_os$Upper_CI_OS), ]
h_os <- calc_plot_height(
  n_studies = nrow(tmp_os),
  n_groups  = length(tumor_order[tumor_order %in% as.character(unique(tmp_os$Tumortype_clean))])
)

pdf(file.path(fig_dir, "Forest_OS_by_tumor_type.pdf"), width = 8.5, height = h_os)
res_os <- forest_with_two_bottom_pools(
  dat = dat, tumor_order = tumor_order, tumor_pretty = tumor_pretty,
  hr_col = "HR_OS", lo_col = "Lower_CI_OS", hi_col = "Upper_CI_OS",
  main_title = "Overall survival by tumor type",
  xlim = c(log(0.01), log(10)),
  at   = log(c(0.1, 0.5, 1, 2, 3)),
  study_col_shift = 0.13
)
dev.off()
tab_os <- print_per_tumor_results(res_os, endpoint_name = "OS", tumor_pretty = tumor_pretty)
all_tables[["OS"]] <- tab_os

# ---- PFS ----
tmp_pfs <- dat
tmp_pfs <- tmp_pfs[!is.na(tmp_pfs$Tumortype_clean) &
                     !is.na(tmp_pfs$HR_PFS) & !is.na(tmp_pfs$lower_CI_PFS) & !is.na(tmp_pfs$upper_CI_PFS), ]
h_pfs <- calc_plot_height(
  n_studies = nrow(tmp_pfs),
  n_groups  = length(tumor_order[tumor_order %in% as.character(unique(tmp_pfs$Tumortype_clean))])
)

pdf(file.path(fig_dir, "Forest_PFS_by_tumor_type.pdf"), width = 8.5, height = h_pfs)
res_pfs <- forest_with_two_bottom_pools(
  dat = dat, tumor_order = tumor_order, tumor_pretty = tumor_pretty,
  hr_col = "HR_PFS", lo_col = "lower_CI_PFS", hi_col = "upper_CI_PFS",
  main_title = "Progression-free survival by tumor type",
  xlim = c(log(0.01), log(10)),
  at   = log(c(0.1, 0.5, 1, 2, 3)),
  study_col_shift = 0.13
)
dev.off()
tab_pfs <- print_per_tumor_results(res_pfs, endpoint_name = "PFS", tumor_pretty = tumor_pretty)
all_tables[["PFS"]] <- tab_pfs

# ---- EFS ----
tmp_efs <- dat
tmp_efs <- tmp_efs[!is.na(tmp_efs$Tumortype_clean) &
                     !is.na(tmp_efs$HR_EFS) & !is.na(tmp_efs$Lower_CI_EFS) & !is.na(tmp_efs$Upper_CI_EFS), ]
h_efs <- calc_plot_height(
  n_studies = nrow(tmp_efs),
  n_groups  = length(tumor_order[tumor_order %in% as.character(unique(tmp_efs$Tumortype_clean))])
)
h_efs <- max(7, h_efs * 0.25)

pdf(file.path(fig_dir, "Forest_EFS_by_tumor_type.pdf"), width = 8.5, height = h_efs)
res_efs <- forest_with_two_bottom_pools(
  dat = dat, tumor_order = tumor_order, tumor_pretty = tumor_pretty,
  hr_col = "HR_EFS", lo_col = "Lower_CI_EFS", hi_col = "Upper_CI_EFS",
  main_title = "Event-free survival by tumor type",
  xlim = c(log(0.01), log(10)),
  at   = log(c(0.1, 0.5, 1, 2, 3)),
  show_excl_gbm = FALSE,
  study_col_shift = 0.13
)
dev.off()
tab_efs <- print_per_tumor_results(res_efs, endpoint_name = "EFS", tumor_pretty = tumor_pretty)
all_tables[["EFS"]] <- tab_efs

# PDF -> PNG (into figures/)
pdf_convert(
  pdf = file.path(fig_dir, "Forest_OS_by_tumor_type.pdf"),
  format = "png", dpi = 300,
  filenames = file.path(fig_dir, "Forest_OS_by_tumor_type.png")
)
pdf_convert(
  pdf = file.path(fig_dir, "Forest_PFS_by_tumor_type.pdf"),
  format = "png", dpi = 300,
  filenames = file.path(fig_dir, "Forest_PFS_by_tumor_type.png")
)
pdf_convert(
  pdf = file.path(fig_dir, "Forest_EFS_by_tumor_type.pdf"),
  format = "png", dpi = 300,
  filenames = file.path(fig_dir, "Forest_EFS_by_tumor_type.png")
)

# =============================================================================
# 4) Write pooled tables to CSV
# =============================================================================

# Combine available tables (OS/PFS/EFS)
combined <- do.call(rbind, all_tables)
if (!is.null(combined) && nrow(combined) > 0) {
  write.csv(combined, file = file.path(res_dir, "per_tumour_pooled_results.csv"), row.names = FALSE)
}

# =============================================================================
# 5) Rewrite index.html to show the latest outputs
# =============================================================================

updated <- format(Sys.time(), "%Y-%m-%d %H:%M")

index_html <- paste0(
'<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Immunorad Living Meta-analysis</title>
</head>
<body style="font-family: Arial; max-width: 980px; margin: auto; line-height: 1.4;">
  <h1>Immunorad Living Meta-analysis</h1>
  <p><b>Last updated:</b> ', updated, '</p>

  <h2>Overall Survival</h2>
  <p><a href="figures/Forest_OS_by_tumor_type.pdf">Download PDF</a></p>
  <img src="figures/Forest_OS_by_tumor_type.png" style="width: 100%; max-width: 900px;">

  <h2>Progression-Free Survival</h2>
  <p><a href="figures/Forest_PFS_by_tumor_type.pdf">Download PDF</a></p>
  <img src="figures/Forest_PFS_by_tumor_type.png" style="width: 100%; max-width: 900px;">

  <h2>Event-Free Survival</h2>
  <p><a href="figures/Forest_EFS_by_tumor_type.pdf">Download PDF</a></p>
  <img src="figures/Forest_EFS_by_tumor_type.png" style="width: 100%; max-width: 900px;">

  <h2>Pooled results (per tumour)</h2>
  <p><a href="results/per_tumour_pooled_results.csv">Download CSV</a></p>
</body>
</html>'
)

writeLines(index_html, con = file.path(repo_dir, "index.html"))

cat("\nDone. Outputs written to:\n")
cat(" - ", fig_dir, "\n", sep = "")
cat(" - ", res_dir, "\n", sep = "")
cat(" - ", file.path(repo_dir, "index.html"), "\n", sep = "")
