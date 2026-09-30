# Merian-element chromosome painting

**Paper:** Figure S4A–C  

Assign complete BUSCO orthologues to ancestral Merian elements, plot chromosome-wide assignments and highlight differences from the dominant assignment on each chromosome.

## Inputs

- An uncompressed A. gemmatalis assembly FASTA for the optional Compleasm preparation step. Use the assembly and sequence names underlying S4 when reproducing the paper.
- full_table.tsv in BUSCO format, with Complete/Duplicated status labels. Section 0 supplies it from a new Compleasm run, or use the retained historical table when available.
- Anticarsia_gemmatalis_primary_pseudohap_scaffold.fasta.fai, corresponding to the assembly and sequence names used for that Compleasm run.
- Merian_elements_full_table.tsv from the pinned lep_busco_painter revision documented below.

## Outputs

- Compleasm assessment outputs and a painting folder containing full_table.tsv and the matching FASTA index.
- buscopainter_summary.tsv; buscopainter_complete_location.tsv; buscopainter_duplicated_location.tsv.
- buscopainter_complete_location.tsv_buscopainter.png and .pdf for the selected plotting mode; preserve each mode’s files separately.
- The summary table supplies Merian-assignment counts and fractions for S4C.

## Software

Compleasm 0.2.6, verified in the mamba environment named by the supplied workflow; Python 3 (the painting notes specify a Python 3.9 environment), lep_busco_painter, SAMtools; R with optparse, tidyverse and scales. Workflow notes: based on Charlotte Wright’s lep_busco_painter.

## Execution and interpretation

1. The supplied local notes contain an explicit ilAntGemm2 plotting command using the A. gemmatalis FASTA index. This establishes an organism-specific painting workflow.
2. Compleasm generated the ortholog assignments used for Figure S4. Section 0 adapts the supplied example to a user-supplied A. gemmatalis assembly, with Lepidoptera and 32 threads from that example. It explicitly uses Compleasm 0.2.6 and its busco mode with Lepidoptera odb10. The workflow notes activate `compleasm_0.2.6_JDA`; its package metadata and executable both confirm Compleasm 0.2.6 (Bioconda build `pyh7cba7a3_0`). The environment history records installation on 9 May 2025. The environment does not identify the original A. gemmatalis assembly, output table or lineage-download snapshot. The adapted example is not the original run script.
3. The installed buscopainter.py and Merian reference table match upstream commit 904ec735c9e78b3747dfe5dbed2c680c593f3d2a byte for byte. Download those exact files via the pinned raw links below; the notes’ wget /blob/ links would retrieve GitHub pages rather than raw source.
4. The custom R script in the shared scripts directory is byte-identical to the installed copy called by the workflow. Both are SHA-256 ba956c7818bb7330c4be6fddfd1cde409e2fbc57517d3888a1fa2e282823f3bd. The notes reference lab commit 8d869d052bd958bb1bc5fb072fb6a859e782d887, but that private remote revision could not be compared directly.
5. For S4A use the recorded -m TRUE invocation. For S4B rerun the same plotting command with -d TRUE added, as specified in the notes. Both modes write the same filenames: preserve/rename the S4A outputs before the S4B run.
6. The plotter requires at least three complete BUSCOs per displayed contig by default (-n 3), filters chromosome labels containing a colon, and uses FASTA-index lengths. The supplied S4 displays 31 chromosome scaffolds; this count has not been regenerated without the exact input table/index.
7. The custom script supports other generic modes, but only the two Merian-mode invocations are relevant here. The original R implementation is preserved, including its repeated final save calls.
8. Keep attribution to Charlotte Wright for the workflow and R edits, and Charlotte Wright for lep_busco_painter.
9. After section 0, run the section 1 painting command from the new output's painting directory. Set lep_scripts_path to the folder containing the pinned Python tool, Merian reference table and the supplied R script. The FASTA index and query table must use the same sequence names.

Save each code block under its stated filename and supply the listed inputs. Replace the example paths with locations on your own computer, install the listed tools, and run the steps in order. Bash scripts run directly with `bash`; any required sample index is explained at the top of the script.

## Code

<a id="step-00"></a>

### 0. Generate the BUSCO table and matching assembly index

**Save as:** `run_compleasm_agem.sh`

Install Compleasm 0.2.6 and SAMtools on your PATH. For example, use `conda create -n compleasm_026 -c conda-forge -c bioconda compleasm=0.2.6 samtools`, then `conda activate compleasm_026`. The package version is pinned explicitly here; an environment name alone does not pin a package version.

Run `bash run_compleasm_agem.sh /path/to/assembly.fasta /path/to/new_output_directory`. Set `THREADS` to your available CPU count if needed; the default is 32. The script downloads/caches Lepidoptera lineage data under the output directory, records versions/settings, and creates the index for the same assembly that Compleasm assesses. Use paths without whitespace for Compleasm 0.2.6.

The retained [Compleasm 0.2.6 interface](https://github.com/huangnengCSU/compleasm/blob/v0.2.6/README.md) supports the example's options. Its [implementation](https://github.com/huangnengCSU/compleasm/blob/v0.2.6/compleasm.py) uses Lepidoptera odb10 and writes both a native table and `full_table_busco_format.tsv`. The adapted code copies the latter to `painting/full_table.tsv`: the pinned painting script selects `Complete` records, whereas the native table uses `Single`. Other status labels and coordinates are preserved; no manual relabelling of the native table is needed.

The generated index uses the historical filename expected by section 1. This is an interface filename; the assembly actually analysed is the input supplied to this script and is recorded in `run_settings.tsv`. The original S4 assembly and exact lineage snapshot could not be determined from the mamba environment; use the assembly supplied to this portable example. The template and output-format connection have been checked, but a biological completeness analysis was not rerun for this archive.

```bash
#!/usr/bin/env bash
# Adapted from Jasmine D. Alqassar's 2024 example; original A. gemmatalis run unconfirmed.
# Usage: bash run_compleasm_agem.sh ASSEMBLY.fasta NEW_OUTPUT_DIRECTORY
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 ASSEMBLY.fasta NEW_OUTPUT_DIRECTORY" >&2
    exit 1
fi
THREADS="${THREADS:-32}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }
[[ -s "$1" ]] || { echo "Assembly FASTA not found: $1" >&2; exit 1; }
[[ ! -e "$2" ]] || { echo "Use a new output directory: $2" >&2; exit 1; }
for tool in compleasm samtools; do
    command -v "$tool" >/dev/null || { echo "Required tool not on PATH: $tool" >&2; exit 1; }
done
ASSEMBLY="$(cd -- "$(dirname -- "$1")" && pwd)/$(basename -- "$1")"
VERSION="$(compleasm --version)"
[[ "$VERSION" =~ [[:space:]]0\.2\.6$ ]] || {
    echo "This template expects Compleasm 0.2.6 (Lepidoptera odb10); found: $VERSION" >&2
    exit 1
}
# Compleasm 0.2.6 constructs some subprocess commands without quoting paths.
[[ ! "$ASSEMBLY" =~ [[:space:]] && ! "$2" =~ [[:space:]] && ! "$PWD" =~ [[:space:]] ]] || {
    echo "Use assembly, output and working-directory paths without whitespace" >&2
    exit 1
}
mkdir -p -- "$2"
OUTROOT="$(cd -- "$2" && pwd)"
LINEAGE_DIR="$OUTROOT/mb_downloads"
PAINTDIR="$OUTROOT/painting"
mkdir -p "$PAINTDIR" "$LINEAGE_DIR"
printf '%s\n' "$VERSION" > "$OUTROOT/compleasm_version.txt"
samtools --version > "$OUTROOT/samtools_version.txt"
printf 'assembly\t%s\nlineage\tlepidoptera_odb10\nmode\tbusco\nthreads\t%s\n' \
    "$ASSEMBLY" "$THREADS" > "$OUTROOT/run_settings.tsv"

# Lepidoptera and 32 threads follow the supplied example; busco is the 0.2.6 default.
compleasm run \
    -a "$ASSEMBLY" \
    -o "$OUTROOT/compleasm" \
    -t "$THREADS" \
    -l lepidoptera \
    -L "$LINEAGE_DIR" \
    -m busco \
    > "$OUTROOT/compleasm.log" 2>&1

# Use the BUSCO-format output: the painter selects status == "Complete".
TABLE="$OUTROOT/compleasm/lepidoptera_odb10/full_table_busco_format.tsv"
[[ -s "$TABLE" ]] || { echo "Expected BUSCO-format table not found: $TABLE" >&2; exit 1; }
cp "$TABLE" "$PAINTDIR/full_table.tsv"
samtools faidx \
    --fai-idx "$PAINTDIR/Anticarsia_gemmatalis_primary_pseudohap_scaffold.fasta.fai" \
    "$ASSEMBLY"
if [[ -f "$LINEAGE_DIR/lepidoptera_odb10/dataset.cfg" ]]; then
    cp "$LINEAGE_DIR/lepidoptera_odb10/dataset.cfg" "$OUTROOT/lineage_dataset.cfg"
fi
echo "Run the Merian-painting command from: $PAINTDIR"
```

<a id="step-01"></a>

### 1. Assign Merian elements and plot ilAntGemm2

Verbatim A. gemmatalis command block from the supplied notes. Adapt lep_scripts_path in a working copy to the directory containing the pinned tool, reference table and supplied custom R script.

**Source:** `Merian_element_painting/Merian_element_painting_JDA.md`  
**Save as:** `paint_ilAntGemm2.sh`  
**Archive treatment:** Verbatim fenced A. gemmatalis command block from the supplied Markdown; installation examples and the other-species Compleasm run are omitted. Replace installation-specific data paths with editable /path/to/project or /path/to/tools examples. Use tools available on PATH instead of local modules, environment activation and executable installation paths.

```bash

lep_scripts_path="/path/to/tools/lep_busco_painter"


${lep_scripts_path}/buscopainter.py -r ${lep_scripts_path}/Merian_elements_full_table.tsv -q full_table.tsv #full_table.tsv is your compleasm output


${lep_scripts_path}/plot_buscopainter_editedJA.R \
    -f buscopainter_complete_location.tsv -p "ilAntGemm2" -i Anticarsia_gemmatalis_primary_pseudohap_scaffold.fasta.fai -m TRUE

```

<a id="step-02"></a>

### 2. Custom chromosome-painting implementation

Jasmine Alqassar’s modified plotter, based on Charlotte Wright’s code. The two recovered copies were identical.

**Source:** `lep_busco_painter/plot_buscopainter_editedJA.R`  
**Save as:** `plot_buscopainter_editedJA.R`  
**Archive treatment:** Complete source, verbatim.

```r
#!/usr/bin/env Rscript


### Load packages
library(optparse)
suppressMessages(library(tidyverse))
suppressMessages(library(scales))


### Funcions for making busco paints in R ####

prepare_data <- function(args1){
  locations <- read_tsv(args1, col_types = cols())
  locations <- locations %>% filter(!grepl(':', query_chr)) # format location data 
  locations <- locations %>% filter(!grepl(':', assigned_chr)) # format location data 
  locations <- locations %>% group_by(query_chr) %>% mutate(length = max(position)) %>% ungroup()
  locations$start <- 0
  return(locations)
}

prepare_data_with_index <- function(args1, args2){
  locations <- read_tsv(args1, col_types=cols())
  contig_lengths <- read_tsv(args2, col_names=FALSE, col_types = cols())
  colnames(contig_lengths) <- c('Seq', 'length', 'offset', 'linebases', 'linewidth')
  locations <- locations %>% filter(!grepl(':', query_chr)) # format location data 
  locations <- locations %>% filter(!grepl(':', assigned_chr)) # format location data 
  locations <- merge(locations, contig_lengths, by.x="query_chr", by.y="Seq")
  locations$start <- 0
  return(locations)
}

filter_buscos <- function(locations, minimum){ # minimum of buscos to be present
  locations_filt <- locations  %>% 
    group_by(query_chr) %>%   # filter df to only keep query_chr with >=3 buscos to remove shrapnel
    mutate(n_busco = n()) %>% # make a new column reporting number buscos per query_chr
    ungroup() %>%
    filter(n_busco >= minimum)
  return(locations_filt)
}
set_merian_colour_mapping <- function(location_set){ # Set mapping of Merian element to colour when only plot 
  merian_order = c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
colour_palette <- c("#999999","#D50062","#F10059","#FF104E","#FF4F44","#FF7839","#FF9E2F",
  "#FF9700","#CF8F00","#9C8600","#6A7A00","#336C00","#0A8500","#009F00",
  "#00B846","#00D27C","#00EBB3","#00D5BF","#00BEC9","#00A6D1","#008FD5",
  "#0079D4","#008BFC","#009BFF","#00A8FF","#00B3FF","#87BCFF","#B198FF",
  "#C971FF","#D546E1","#D500B0","#CC007E", "grey")
  status_merians <- unique(location_set$status)
  subset_merians <- subset(colour_palette, merian_order %in% status_merians)
  return(subset_merians)
}
busco_paint_theme <- theme(legend.position="right",
                           strip.text.x = element_text(margin = margin(0,0,0,0, "cm")), 
                           panel.background = element_rect(fill = "white", colour = "white"), 
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           axis.line.x = element_line(color="black", size = 0.5),
                           axis.text.x = element_text(size=15),
                           axis.title.x = element_text(size=15),
                           strip.text.y = element_text(angle=0),
                           strip.background = element_blank(),
                           plot.title = element_text(hjust = 0.5, face="italic", size=20),
                           plot.subtitle = element_text(hjust = 0.5, size=20))

busco_paint_no_facet_labels_theme <- theme(legend.position="right",
                           strip.text.x = element_blank(), 
                           panel.background = element_rect(fill = "white", colour = "white"), 
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           axis.line.x = element_line(color="black", size = 0.5),
                           axis.text.x = element_text(size=15),
                           axis.title.x = element_text(size=15),
                           strip.text.y = element_text(angle=0),
                           strip.background = element_blank(),
                           plot.title = element_text(hjust = 0.5, face="italic", size=20),
                           plot.subtitle = element_text(hjust = 0.5, size=20))
# plot only buscos that have moved - paint by Merians
paint_merians_differences_only <- function(spp_df, subset_merians, num_col, title, karyotype){
  merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  spp_df$status_f =factor(spp_df$status, levels=merian_order)
  chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(length, decreasing=TRUE) 
  chr_levels <- chr_levels$query_chr
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=rev(chr_levels)) # set chr order as order for plotting, JA added rev to reverse the order
  sub_title <- paste("n contigs =", karyotype) 
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=subset_merians, aesthetics=c("colour", "fill")) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=status_f)) +
    facet_wrap(query_chr_f ~., ncol=1, strip.position="right")  + guides(scale="none") + 
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + 
    ggtitle(label=title, subtitle= sub_title)  +
    guides(fill=guide_legend("Merian element"), color = "none")  +
    busco_paint_theme
  return(the_plot)
}


# plot only buscos that have moved - paint by species
paint_species_differences_only <- function(spp_df, num_col, title, karyotype){
  chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(length, decreasing=TRUE)
  chr_levels <- chr_levels$query_chr
  chr_levels = chr_levels [! chr_levels %in% "self"]
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=chr_levels) # set chr order as order for plotting query chr
  legend_levels <- unique(spp_df$status)
  legend_levels <- legend_levels[legend_levels != 'self'] # remove 'self' from list
  legend_levels <- c('self',legend_levels) # then put 'self' back in to have it in first position as want 'self' to always be painted grey.
  num_colours <- length(legend_levels)
  col_palette <- hue_pal()(num_colours) 
  col_palette[1] <- 'grey'
  spp_df$status_f = factor(spp_df$status, levels=legend_levels) # set chr order as order for plotting
  
  sub_title <- paste("n contigs =", karyotype) 
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=col_palette, aesthetics=c("fill"), breaks=legend_levels) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=status_f)) + 
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + 
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + 
    ggtitle(label=title, subtitle= sub_title)  +
    guides(fill=guide_legend("Query chromosome"), color = "none") +
    busco_paint_theme
  return(the_plot)
}

paint_merians_all <- function(spp_df, num_col, title, karyotype){
colour_palette <- c("#999999","#D50062","#F10059","#FF104E","#FF4F44","#FF7839","#FF9E2F",
  "#FF9700","#CF8F00","#9C8600","#6A7A00","#336C00","#0A8500","#009F00",
  "#00B846","#00D27C","#00EBB3","#00D5BF","#00BEC9","#00A6D1","#008FD5",
  "#0079D4","#008BFC","#009BFF","#00A8FF","#00B3FF","#87BCFF","#B198FF",
  "#C971FF","#D546E1","#D500B0","#CC007E", "grey")
  merian_order <- c('MZ', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20','M21', 'M22', 'M23', 'M24', 'M25', 'M26', 'M27', 'M28', 'M29', 'M30', 'M31', 'self')
  spp_df$assigned_chr_f =factor(spp_df$assigned_chr, levels=merian_order)
  chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(length, decreasing=TRUE) 
  chr_levels <- chr_levels$query_chr
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=rev(chr_levels)) # set chr order as order for plotting
  sub_title <- paste("n contigs =", karyotype) 
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=colour_palette, aesthetics=c("colour", "fill")) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=assigned_chr_f)) +
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + 
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + 
    ggtitle(label=title, subtitle= sub_title)  +
    guides(fill=guide_legend("Merian element"), color = "none") +
    busco_paint_theme
  return(the_plot)
}

# paint all buscos by species
paint_species_all <- function(spp_df, num_col, title, karyotype){
  chr_levels <- subset(spp_df, select = c(query_chr, length)) %>% unique() %>% arrange(length, decreasing=TRUE)
  chr_levels <- chr_levels$query_chr
  chr_levels = chr_levels [! chr_levels %in% "self"]
  spp_df$query_chr_f =factor(spp_df$query_chr, levels=chr_levels) # set chr order as order for plotting
  legend_levels <- subset(spp_df, select = c(assigned_chr)) %>% unique() 
  legend_levels <- legend_levels$assigned_chr
  num_colours <- length(legend_levels)
  col_palette <- hue_pal()(num_colours) 
  spp_df$assigned_chr_f = factor(spp_df$assigned_chr, levels=legend_levels) # set chr order as order for plotting
  
  sub_title <- paste("n contigs =", karyotype) 
  the_plot <- ggplot(data = spp_df) +
    scale_colour_manual(values=col_palette, aesthetics=c("fill"), breaks=legend_levels) +
    geom_rect(aes(xmin=start, xmax=length, ymax=0, ymin =12), colour="black", fill="white") + 
    geom_rect(aes(xmin=position-2e4, xmax=position+2e4, ymax=0, ymin =12, fill=assigned_chr_f)) + 
    facet_wrap(query_chr_f ~., ncol=num_col, strip.position="right") + guides(scale="none") + 
    xlab("Position (Mb)") +
    scale_x_continuous(labels=function(x)x/1e6, expand=c(0.005,1)) +
    scale_y_continuous(breaks=NULL) + 
    ggtitle(label=title, subtitle= sub_title)  +
    guides(fill=guide_legend("Query chromosome"), color = "none") +
    busco_paint_theme
  return(the_plot)
}

### get args
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
        help="location.tsv file", metavar="character"),
    make_option(c("-p", "--prefix"), type="character", default="Query species",
        help="prefix for plot title",  metavar="character"),
    make_option(c("-i", "--index"), type="character", default="False", 
        help="genome index file", metavar="character"),
    make_option(c("-m", "--merians"), type="character", default="False", 
        help="use this flag if you are comparing a genome to Merian elements", metavar="character"),
    make_option(c("-d", "--differences"), type="character", default="False",
        help="only colour buscos that have moved from the dominant chromosome", metavar="character"),
    make_option(c("-n", "--minimum"), type="integer", default=3,
        help="minimum number of buscos ", metavar="number")
        );
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

locations <- opt$file
prefix <- opt$prefix

index <- opt$index
merians <- opt$merians
differences_only <- opt$differences
minimum <- opt$minimum

if (index == "False"){ # if no index supplied
    location_set <- prepare_data(locations)
    locations_filt <- filter_buscos(location_set, minimum)
} else { # if index supplied
    location_set <- prepare_data_with_index(locations, index) 
    locations_filt <- filter_buscos(location_set, minimum)
}

total_contigs <- length(unique(location_set$query_chr))# total number of query_chr before filtering
num_contigs <- as.character(length(unique(locations_filt$query_chr))) # number of query_chr after filtering
num_removed_contigs <- length(unique(location_set$query_chr)) - length(unique(locations_filt$query_chr)) 
print(paste('Number of contigs before filtering by number of BUSCOs:', total_contigs))
print(paste('Number of contigs removed by filtering :', num_removed_contigs))
print(paste('Number of contigs post-filtering:', num_contigs))

if (merians != "False"){ # if Merian elements are being used as the comparator
    subset_merians <- set_merian_colour_mapping(locations_filt)
    }

# generate the plot - four possible options based on given arguments to script
# plot only buscos that have moved - paint by Merians
if (merians == "False"){ # if comparing two species
    if (differences_only == "False"){ # if colouring all orthologs
        p <- paint_species_all(locations_filt, 1, prefix, num_contigs)
    } else { # if only colouring orthologs that have moved
        p <- paint_species_differences_only(locations_filt, 1, prefix, num_contigs)
    }

} else { # comparing one species to Merian elements 
    if (differences_only == "False"){ # if colouring all orthologs
        p <- paint_merians_all(locations_filt, 1, prefix, num_contigs)
    } else { # if only colouring orthologs that have moved
    if (length(locations_filt$query_chr) < 100){
      p <- paint_merians_differences_only(locations_filt, subset_merians, 1, prefix, num_contigs) 
      p <- p + busco_paint_theme
    } else {
      p <- paint_merians_differences_only(locations_filt, subset_merians, 3, prefix, num_contigs) 
      #p <- p + busco_paint_theme
      p <- p + busco_paint_no_facet_labels_theme
    }
    }    
}


ggsave(paste(as.character(opt$file), "_buscopainter.png", sep = ""), plot = p, width = 15, height = 30, units = "cm", device = "png")
pdf(NULL)
ggsave(paste(as.character(opt$file), "_buscopainter.pdf", sep = ""), plot = p, width = 15, height = 30, units = "cm", device = "pdf")


ggsave(paste(as.character(opt$file), "_buscopainter.png", sep = ""), plot = p, width = 15, height = 30, units = "cm", device = "png")
pdf(NULL)
ggsave(paste(as.character(opt$file), "_buscopainter.pdf", sep = ""), plot = p, width = 15, height = 30, units = "cm", device = "pdf")

```

## Pinned dependencies and attribution


The upstream Python program and Merian reference table were compared with the installed files on 14 September 2026. Both match the revision named in Jasmine Alqassar’s workflow notes. Save the linked raw files under their original basenames. Ensure scripts invoked directly have executable permissions, or use Python/Rscript in a working copy.


| Dependency | SHA-256 |
|---|---|

| [buscopainter.py](https://raw.githubusercontent.com/charlottewright/lep_busco_painter/904ec735c9e78b3747dfe5dbed2c680c593f3d2a/buscopainter.py) | `4fa3379c5cf0fc2d630956ae1f8c720f0775cbbacad1a55a7b3de67647863010` |

| [Merian_elements_full_table.tsv](https://raw.githubusercontent.com/charlottewright/lep_busco_painter/904ec735c9e78b3747dfe5dbed2c680c593f3d2a/reference_data/Merian_elements_full_table.tsv) | `0646b1794446ca3570c378f9a94d8c7eec608d14aade215e6127374c0d765e73` |


[Charlotte Wright’s upstream usage documentation](https://github.com/charlottewright/lep_busco_painter/blob/904ec735c9e78b3747dfe5dbed2c680c593f3d2a/README.md). Lab workflow source: `Merian_element_painting/Merian_element_painting_JDA.md`, originally linked from the Martin-Lab-Bioinformatics repository.
