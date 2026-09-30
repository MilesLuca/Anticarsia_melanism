# Pool-seq population-statistic plots

**Paper:** Figure 2B; population-statistic panels of S8–S9  

Visualize saved grenedalf FST, between-pool divergence DXY and combined-pool nucleotide diversity. Retain the historical windowed Fisher-statistic distribution that is visibly present in Figure 2B.

## Inputs

- PoolSeq/Popgen_stats/Grendalf_stats/sizes.genome
- PoolSeq/Popgen_stats/Grendalf_stats/PlDkfst.csv
- PoolSeq/Popgen_stats/Grendalf_stats/PlDkpi-between.csv
- PoolSeq/Popgen_stats/Grendalf_stats/PlDkpi-total.csv
- PoolSeq/Popgen_stats/Grendalf_stats/PlainDark_10k_1k.filt.fet
- PoolSeq/Popgen_stats/Grendalf_stats/NC_134752.1.gff

## Outputs

- GWA-versus-background density plots and summary table.
- Individual genome-wide and chromosome-8 FST/DXY/Pi-total panels.

## Software

R: ggplot2, dplyr, readr, cowplot and tidyr. grenedalf generated the inputs. [Analysis 20](20_grenedalf_population_statistics.md) provides a new calculation workflow, supported by the Hudson identity in the saved FST/pi tables.

## Execution and interpretation

1. Windows are 10 kb with a 1-kb slide. The GWA group contains only windows fully inside chr8:2,012,000–2,613,000: start >= 2,012,000 and end <= 2,613,000. Windows crossing either boundary remain in the rest-of-genome group. The recovered selection code is unchanged.
2. Negative FST values were removed. The archive now filters FST to values >= 0 immediately after loading, before all FST plots, distributions and summaries. Zero values are retained.
3. Export calls only save the selected population-statistic panels.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-01"></a>

### 1. Read, summarize and plot the published statistic series

Curated extraction; source ranges and all edits are recorded in the source manifest.

**Source:** `PoolSeq/Popgen_stats/plot_stacked_gwa_popgen.R`  
**Save as:** `plot_stacked_gwa_popgen.R`  
**Archive treatment:** Assemble source lines 9–37, 45–173, 188–218, 225–462, 473–483, 526–536, 577–670. Omit Pi-within from distribution rows/factor levels and repair the list comma. Omit superseded GWAS scatterplot calls, window-significance cutoff, mixed stacks and diagnostics. Add six explicitly labeled PDF export calls for the selected existing population-statistic plot objects. Exclude of negative FST values before plotting and summaries; retain zero values. No estimator or normalization change. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples.

```r
setwd("/path/to/project/PoolSeq/Popgen_stats/Grendalf_stats")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(cowplot)
  library(tidyr)
})

## ----------------------------------------------------------------------------
## 0. User settings
## ----------------------------------------------------------------------------

target_scaffold <- "NC_134752.1"   # chr8
zoom_start      <- 1500000L
zoom_end        <- 3500000L

window_size <- 10000L
step_size   <- 1000L

# Reconstruct 10 kb GWAS window bounds from the .fet midpoint-like position
fet_left_half  <- 4999L
fet_right_half <- 5000L

# Fixed, user-established GWA interval
gwa_start <- 2012000L
gwa_end   <- 2613000L


# Output directory
outdir <- "stacked_gwa_popgen_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ----------------------------------------------------------------------------
## 1. Helpers
## ----------------------------------------------------------------------------

fmt_bp <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

axis_bp_labels <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

save_pdf_plot <- function(plot_obj, filename, width, height) {
  pdf(filename, width = width, height = height, onefile = TRUE)
  print(plot_obj)
  dev.off()
}

## ----------------------------------------------------------------------------
## 2. Read genome sizes and build cumulative genome coordinates
## ----------------------------------------------------------------------------

sizes <- read_tsv(
  "sizes.genome",
  col_names = FALSE,
  show_col_types = FALSE
)
colnames(sizes) <- c("scaffold", "length")

sizes <- sizes %>%
  mutate(type = ifelse(grepl("^NC_", scaffold), "chrom", "unplaced"))

# NC_ scaffolds each become one chromosome in file order
chrom_scaffolds <- sizes %>%
  filter(type == "chrom") %>%
  mutate(
    chrom_label = as.character(row_number()),
    offset_within_chr = 0L
  )

# All NW_ scaffolds are concatenated into a single "Unplaced" pseudochromosome
unplaced_scaffolds <- sizes %>%
  filter(type == "unplaced") %>%
  arrange(scaffold) %>%
  mutate(
    chrom_label = "Unplaced",
    offset_within_chr = lag(cumsum(length), default = 0)
  )

scaffold_map <- bind_rows(chrom_scaffolds, unplaced_scaffolds)

n_chrom <- sum(sizes$type == "chrom")
chrom_levels <- c(as.character(seq_len(n_chrom)), "Unplaced")

chrom_lengths <- scaffold_map %>%
  group_by(chrom_label) %>%
  summarise(chr_len = sum(length), .groups = "drop") %>%
  mutate(chrom_label = factor(chrom_label, levels = chrom_levels)) %>%
  arrange(chrom_label) %>%
  mutate(
    cum_start = lag(cumsum(chr_len), default = 0),
    center    = cum_start + chr_len / 2,
    parity    = factor(seq_len(n()) %% 2)
  )

scaffold_map <- scaffold_map %>%
  left_join(
    chrom_lengths %>% select(chrom_label, cum_start, parity),
    by = "chrom_label"
  )

background_rects <- chrom_lengths %>%
  transmute(
    chrom_label,
    xmin = cum_start,
    xmax = cum_start + chr_len,
    center,
    parity
  )

add_genome_coords <- function(df) {
  df %>%
    mutate(
      start    = as.integer(start),
      end      = as.integer(end),
      midpoint = as.integer(floor((start + end) / 2))
    ) %>%
    left_join(
      scaffold_map %>%
        select(scaffold, chrom_label, offset_within_chr, cum_start, parity),
      by = "scaffold"
    ) %>%
    mutate(
      pos_chr = midpoint + offset_within_chr,
      pos_cum = pos_chr + cum_start
    )
}

## ----------------------------------------------------------------------------
## 3. Read GWAS .fet
## ----------------------------------------------------------------------------

fet <- read_tsv(
  "PlainDark_10k_1k.filt.fet",
  col_names = FALSE,
  show_col_types = FALSE
)

colnames(fet) <- c(
  "scaffold",
  "pos",
  "n_sites",
  "freq_diff",
  "avg_cov",
  "logp"
)

fet <- fet %>%
  mutate(
    pos      = as.integer(pos),
    n_sites  = as.integer(n_sites),
    start    = pmax(1L, pos - fet_left_half),
    end      = pos + fet_right_half
  ) %>%
  add_genome_coords()

read_metric_csv <- function(file, dataset_label) {
  df <- read_csv(file, show_col_types = FALSE)

  metric_col <- names(df)[ncol(df)]

  out <- df %>%
    transmute(
      scaffold = .data$chrom,
      start    = as.integer(.data$start),
      end      = as.integer(.data$end),
      metric   = as.numeric(.data[[metric_col]])
    ) %>%
    filter(!is.na(metric)) %>%
    add_genome_coords() %>%
    mutate(dataset = dataset_label)

  observed_widths <- unique(out$end - out$start + 1L)

  if (!all(observed_widths == window_size)) {
    warning(
      file, ": unexpected window sizes detected: ",
      paste(sort(unique(observed_widths)), collapse = ", ")
    )
  }

  out
}

# Discard negative FST estimates before plots and summaries.
# Keep zero estimates; do not replace negative estimates with zero.
fst_df       <- read_metric_csv("PlDkfst.csv",        "Fst") %>%
  filter(metric >= 0)
dxy_df       <- read_metric_csv("PlDkpi-between.csv", "Dxy")
pi_total_df  <- read_metric_csv("PlDkpi-total.csv",   "Pi total")

if (is.null(gff_file)) {
  gff_candidates <- list.files(
    pattern = "\\.(gff|gff3)$",
    ignore.case = TRUE
  )

  if ("NC_134752.1.gff" %in% gff_candidates) {
    gff_file <- "NC_134752.1.gff"
  } else if (length(gff_candidates) == 1) {
    gff_file <- gff_candidates[1]
  } else {
    stop(
      "Could not uniquely determine the GFF file. Candidates found: ",
      paste(gff_candidates, collapse = ", "),
      "\nSet gff_file <- 'your_file.gff' manually near the top of the script."
    )
  }
}

gff <- read_tsv(
  gff_file,
  comment = "#",
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(gff) < 9) {
  stop("GFF file does not appear to have at least 9 columns: ", gff_file)
}

gff <- gff[, 1:9]
colnames(gff) <- c(
  "seqid", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

gff <- gff %>%
  mutate(
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  filter(
    seqid == target_scaffold,
    end >= zoom_start,
    start <= zoom_end
  )

gene_feature_type <- if ("gene" %in% gff$type) {
  "gene"
} else if ("mRNA" %in% gff$type) {
  "mRNA"
} else {
  NA_character_
}

if (is.na(gene_feature_type)) {
  warning(
    "No gene or mRNA features found in ", gff_file,
    " for the zoom region; annotation panel may be sparse."
  )
}

gff_gene <- gff %>%
  filter(type == gene_feature_type) %>%
  transmute(
    xmin = start,
    xmax = end,
    y    = ifelse(strand == "+", 0.78, 0.22)
  )

gff_exon <- gff %>%
  filter(type %in% c("exon", "CDS")) %>%
  transmute(
    xmin = start,
    xmax = end,
    ymin = ifelse(strand == "+", 0.60, 0.10),
    ymax = ifelse(strand == "+", 0.90, 0.40)
  )

p_gff <- ggplot() +
  geom_rect(
    data = data.frame(xmin = gwa_start, xmax = gwa_end),
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "goldenrod",
    alpha = 0.14
  ) +
  geom_rect(
    data = gff_exon,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey20"
  ) +
  geom_segment(
    data = gff_gene,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    linewidth = 0.35
  ) +
  scale_x_continuous(
    limits = c(zoom_start, zoom_end),
    labels = axis_bp_labels,
    expand = c(0, 0)
  ) +
  labs(y = NULL) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

## ----------------------------------------------------------------------------
## 6. Plot helpers
## ----------------------------------------------------------------------------

make_genomewide_plot <- function(df, y_col, y_lab,
                                 sig_df = NULL,
                                 cutoff = NULL,
                                 point_size = 0.22,
                                 point_alpha = 0.75) {

  p <- ggplot() +
    geom_rect(
      data = background_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = parity),
      inherit.aes = FALSE,
      alpha = 0.22
    ) +
    geom_point(
      data = df,
      aes(x = pos_cum, y = .data[[y_col]], color = parity),
      size = point_size,
      alpha = point_alpha
    ) +
    scale_fill_manual(values = c("0" = "grey90", "1" = "grey80")) +
    scale_color_manual(values = c("0" = "grey45", "1" = "black")) +
    scale_x_continuous(
      breaks = chrom_lengths$center,
      labels = chrom_lengths$chrom_label,
      expand = c(0.001, 0.001)
    ) +
    labs(x = "Chromosome", y = y_lab) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )

  if (!is.null(sig_df) && !is.null(cutoff)) {
    p <- p +
      geom_point(
        data = sig_df,
        aes(x = pos_cum, y = .data[[y_col]]),
        inherit.aes = FALSE,
        color = "red",
        size = point_size * 1.8,
        alpha = 0.9
      ) +
      geom_hline(
        yintercept = cutoff,
        linetype = "dashed",
        color = "red",
        linewidth = 0.35
      )
  }

  p
}

make_zoom_plot <- function(df, y_col, y_lab,
                           sig_df = NULL,
                           cutoff = NULL,
                           line_color = "grey35",
                           point_color = "grey20") {

  zoom_df <- df %>%
    filter(
      scaffold == target_scaffold,
      end >= zoom_start,
      start <= zoom_end
    ) %>%
    arrange(midpoint)

  p <- ggplot(zoom_df, aes(x = midpoint, y = .data[[y_col]])) +
    geom_rect(
      data = data.frame(xmin = gwa_start, xmax = gwa_end),
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "goldenrod",
      alpha = 0.14
    ) +
    geom_line(color = line_color, linewidth = 0.35) +
    geom_point(color = point_color, size = 0.35, alpha = 0.75) +
    scale_x_continuous(
      limits = c(zoom_start, zoom_end),
      labels = axis_bp_labels,
      expand = c(0, 0)
    ) +
    labs(x = paste0(target_scaffold, " position (bp)"), y = y_lab) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 5.5, 0, 5.5)
    )

  if (!is.null(sig_df) && !is.null(cutoff)) {
    sig_zoom <- sig_df %>%
      filter(
        scaffold == target_scaffold,
        end >= zoom_start,
        start <= zoom_end
      )

    p <- p +
      geom_point(
        data = sig_zoom,
        aes(x = midpoint, y = .data[[y_col]]),
        inherit.aes = FALSE,
        color = "red",
        size = 0.6,
        alpha = 0.9
      ) +
      geom_hline(
        yintercept = cutoff,
        linetype = "dashed",
        color = "red",
        linewidth = 0.35
      )
  }

  p
}

p_fst_genome <- make_genomewide_plot(
  fst_df, "metric", "Fst"
)

p_dxy_genome <- make_genomewide_plot(
  dxy_df, "metric", "Dxy"
)

p_pitotal_genome <- make_genomewide_plot(
  pi_total_df, "metric", "Pi total"
)

p_fst_zoom <- make_zoom_plot(
  fst_df, "metric", "Fst"
)

p_dxy_zoom <- make_zoom_plot(
  dxy_df, "metric", "Dxy"
)

p_pitotal_zoom <- make_zoom_plot(
  pi_total_df, "metric", "Pi total"
)

make_dist_df <- function(df, value_col, dataset_label) {
  df %>%
    transmute(
      scaffold,
      start,
      end,
      value = .data[[value_col]]
    ) %>%
    filter(!is.na(value)) %>%
    mutate(
      label = ifelse(
        scaffold == target_scaffold & start >= gwa_start & end <= gwa_end,
        "GWA interval",
        "rest of genome"
      ),
      dataset = dataset_label
    )
}

dist_long <- bind_rows(
  make_dist_df(fet,          "logp",   "GWAS (-log10 P)"),
  make_dist_df(fst_df,       "metric", "Fst"),
  make_dist_df(dxy_df,       "metric", "Dxy"),
  make_dist_df(pi_total_df,  "metric", "Pi total")
) %>%
  mutate(
    dataset = factor(
      dataset,
      levels = c("GWAS (-log10 P)", "Fst", "Dxy", "Pi total")
    ),
    label = factor(label, levels = c("GWA interval", "rest of genome"))
  )

dist_summary <- dist_long %>%
  group_by(dataset, label) %>%
  summarise(
    n_windows = n(),
    mean      = mean(value, na.rm = TRUE),
    median    = median(value, na.rm = TRUE),
    sd        = sd(value, na.rm = TRUE),
    min       = min(value, na.rm = TRUE),
    max       = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  dist_summary,
  file.path(outdir, "distribution_summary_GWA_interval_vs_rest.tsv")
)

p_dist <- ggplot(
  dist_long,
  aes(x = value, fill = label)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    bins = 80,
    alpha = 0.55
  ) +
  facet_wrap(~ dataset, scales = "free", ncol = 1) +
  scale_fill_manual(
    values = c(
      "GWA interval"   = "#E6A19A",
      "rest of genome" = "#6CCBCD"
    )
  ) +
  labs(
    x = "Metric value",
    y = "density",
    fill = "label"
  ) +
  theme_gray() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

save_pdf_plot(
  p_dist,
  file.path(outdir, "distribution_histograms_GWA_interval_vs_rest.pdf"),
  width = 12,
  height = 16
)

ggsave(
  filename = file.path(outdir, "distribution_histograms_GWA_interval_vs_rest.png"),
  plot = p_dist,
  width = 12,
  height = 16,
  units = "in",
  dpi = 300
)


# Archive-added export calls: save only the selected population-statistic panels.
# These filenames are new archive exports, not names of the historical figure files.
archive_panels <- list(
  Fst_genome = p_fst_genome, Dxy_genome = p_dxy_genome,
  Pi_total_genome = p_pitotal_genome,
  Fst_chr8 = p_fst_zoom, Dxy_chr8 = p_dxy_zoom,
  Pi_total_chr8 = p_pitotal_zoom
)
for (archive_panel_name in names(archive_panels)) {
  save_pdf_plot(archive_panels[[archive_panel_name]],
                file.path(outdir, paste0("archive_", archive_panel_name, ".pdf")),
                width = 12, height = 4)
}
```
