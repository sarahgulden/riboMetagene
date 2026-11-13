
<!-- README.md is generated from README.Rmd. Please edit that file -->

# riboMetagene

riboMetagene turns Ribo-seq BAM files into P-site offsets and
P-site-corrected metagene profiles around start and stop codons. It also
checks the read length distribution and the 3-nt periodicity of your
data. It is a work in progress and is developed and tested on mouse
tissue data (GENCODE annotation).

What it does:

- reads each BAM once and keeps a small table of read 5’ ends, which is
  easy to copy from a cluster to a laptop
- shows the read length distribution (data frame, line plot and heatmap)
  and suggests read lengths
- selects one transcript per gene from your annotation, using a rule you
  choose
- works in transcript coordinates, so introns never break a window
  around a codon
- finds the P-site offset for every read length, checks the periodicity
  in the CDS, and keeps the read lengths that pass both checks
- returns corrected profiles, one row per sample, read length, codon
  type and position, ready for your own plots

## Installation

``` r
# install.packages("remotes")
remotes::install_github("sarahgulden/riboMetagene")
```

The package needs these Bioconductor packages, which are installed
automatically if your setup allows it: GenomicAlignments, GenomicRanges,
GenomeInfoDb, BiocGenerics, IRanges, Rsamtools, S4Vectors and
GenomicFeatures. If the installation stops with a “package not
available” error, install them first:

``` r
BiocManager::install(c("GenomicAlignments", "GenomicRanges", "GenomeInfoDb", "BiocGenerics",
                       "IRanges", "Rsamtools", "S4Vectors", "GenomicFeatures"))
```

The plot functions also need `ggplot2` and `patchwork`. They are
optional, and the package tells you if they are missing.

On a cluster, use `remotes` (not `devtools`) and restart R after
installing, because an already loaded version is not replaced.

## What you need

| Input | Details |
|----|----|
| BAM files | Genomic alignments, for example from STAR. Only the extraction step reads them. rRNA removal and the handling of multi-mappers happen before this package |
| Annotation | A GFF3 or GTF of the same genome, imported with `rtracklayer::import()`. It needs a `type` column with exon, start_codon and stop_codon rows, and the columns `gene_id` and `transcript_id`. A `tag` column is optional. Chromosome names must match the BAMs |
| Experiment information | Optional: columns such as tissue, timepoint or replicate in your sample table (see below) |

You do not need chromosome sizes, bedgraph files or transcriptome BAMs.

## Workflow at a glance

    BAM files (cluster)
       list_bam_samples()  ->  samples
       extract_5prime()    ->  counts: one small 5' end table per sample
                                  |-- readlength_distribution() -> qc -> plots, suggest_readlengths()
    annotation (GFF3/GTF)         |
       prepare_annotation() ->  annot: one transcript per gene (exons + lengths)
                                  |
       metagene_pipeline(counts, annot, samples)
            1 map 5' ends to transcripts
            2 profiles around start / stop codons
            3 P-site offset per read length (the samples vote)
            4 periodicity in the CDS
            5 select read lengths, apply the correction

Stage 1 runs on the cluster, where the BAM files are. Everything after
it runs on small files and can run on a laptop.

## Quick start

``` r
library(riboMetagene)

## Stage 1, on the cluster: BAM files -> 5' end tables
samples <- list_bam_samples("path/to/bams")
counts  <- extract_5prime(samples, out_dir = "5prime_tables")   # saves one .rds per sample

## Copy the folder "5prime_tables" to your laptop, then:
counts  <- read_5prime_tables("5prime_tables")

## Your experiment information (see the next section)
samples <- data.frame(sample_name = names(counts))
samples$tissue    <- ifelse(grepl("KI", samples$sample_name), "Kidney", "Liver")
samples$timepoint <- gsub(".*_(KI|LIV)_(\\d+).*", "\\1_\\2", samples$sample_name)

## Read length QC
qc <- readlength_distribution(counts, samples)
plot_readlength_qc(qc, color_by = "timepoint", highlight = suggest_readlengths(qc))

## Annotation, once per annotation: one transcript per gene
annotation <- rtracklayer::import("annotation.gff3.gz")
annot <- prepare_annotation(annotation,
                            prefer = transcripts_with_tag(annotation, "Ensembl_canonical"))

## The pipeline
res <- metagene_pipeline(counts, annot, samples = samples, readlengths = 26:32)
```

## Adding your experiment information

The package knows nothing about your experiment. You describe it in the
`samples` table by adding columns, and the functions that take
`samples =` attach these columns to their results.

``` r
# from the sample names
samples$replicate <- gsub("_Ribo.*", "", samples$sample_name)

# or from a sample sheet (a CSV with a column sample_name)
samples <- merge(samples, read.csv("sample_sheet.csv"), by = "sample_name")
```

- Every name in `names(counts)` has to appear in `samples$sample_name`.
  Extra rows are ignored.
- Avoid these column names, because they are used in the results:
  `sample_name`, `readlength`, `feature`, `position`, `signal`, `n`,
  `pct`, `frame`, `fraction`, `score`, `peak_position`,
  `position_corrected`.
- A short `label` column makes the heatmap rows readable. A factor
  column with the levels in the order you want gives that order in the
  plots.

## Output directory of the extraction

`extract_5prime()` always returns the tables in memory. With `out_dir`
it also saves one `<sample_name>.rds` per sample. A rerun then loads the
finished samples instead of reading the BAM again, so an interrupted run
resumes. Load the folder later with `read_5prime_tables()`. Without
`out_dir` nothing is saved.

## Choosing one transcript per gene

The choice happens in `prepare_annotation()`, before the pipeline. Only
transcripts with exons, a start codon and a stop codon are candidates.

``` r
# a tag in your annotation (GENCODE: "Ensembl_canonical"; human: also "MANE_Select")
prepare_annotation(annotation, prefer = transcripts_with_tag(annotation, "Ensembl_canonical"))

# your own list of transcript IDs
prepare_annotation(annotation, prefer = my_transcript_ids)

# no preference: the transcript with the longest CDS per gene
prepare_annotation(annotation)

# random, as a robustness check (same seed = same result)
prepare_annotation(annotation, fallback = "random", seed = 1)
```

For each gene, the preferred transcripts come first, and the fallback
rule (longest CDS, or random) decides among the rest. The result `annot`
is a list with `exons_tx` (the exons in 5’ to 3’ order) and `tx_info`
(transcript length and the lengths of the 5’ UTR, CDS and 3’ UTR). Do
not copy `annot` between machines as an `.rds` file. Rebuild it from the
annotation file instead, which takes about a minute.

## Settings of the pipeline

``` r
metagene_pipeline(counts, annot, samples = NULL,
  readlengths = 26:32,                               # read lengths to consider
  select = TRUE, min_agree = 0.9, min_frame0 = 0.8,  # automatic selection, see below
  features = c("start", "stop"),                     # "start" is always needed for the offsets
  upstream = 25, downstream = 25,                    # nt before and after the codon
  threshold = 0.95,                                  # maximum share of zero positions per site
  common_set = FALSE,                                # TRUE: same genes for start and stop
  skip_codons = 5,                                   # codons skipped at both ends of the CDS
  trim_edges = TRUE,                                 # keep positions that all selected lengths cover
  keep_mapped = FALSE)                               # TRUE: also return the mapped tables (large)
```

- **`threshold`:** a site (one transcript’s window) is kept if its share
  of positions without reads is below this value. At 0.95 a site needs
  reads at 3 of 53 positions. Kept sites are normalized to 0-1 and
  averaged per position.
- **`select`:** with `TRUE`, a read length is kept only if at least
  `min_agree` of the samples have their start-codon peak at the
  consensus offset and frame 0 is the dominant CDS frame in at least
  `min_frame0` of the samples. With `FALSE`, all of `readlengths` are
  used. The table with the reasons is always printed.
- **`counts`:** you can pass tables that were already mapped with
  `map_to_transcripts()` (or `res$mapped`). The slowest step is then
  skipped, which helps when you try other settings.

## What you get back

`metagene_pipeline()` returns a list:

| Element | Content |
|----|----|
| `corrected` | The main result: profiles of the selected read lengths with `position_corrected`, one row per sample, read length, codon type and position, plus your experiment columns |
| `offsets` | The P-site offset of every considered read length. The votes of the samples are in `attr(res$offsets, "peaks")` |
| `readlengths` | The read lengths that were used |
| `periodicity` | Per sample and read length: number of P-sites, share in frame 0, dominant frame and a score |
| `frames` | P-sites per frame (0, 1, 2) for every sample and read length |
| `signal` | Profiles of all considered read lengths before the correction |
| `mapped` | Only with `keep_mapped = TRUE` |

In `signal` and `corrected`, position 0 is the first nucleotide of the
codon and negative positions are towards the 5’ end of the transcript.
After the correction, position 0 means that the P-site is at the first
nucleotide of the codon. Save the whole list with
`saveRDS(res, "res.rds")`. It only contains plain data frames.

## Plots

Compare timepoints (average over read lengths and replicates, one panel
per tissue and codon type):

``` r
library(ggplot2)

res$corrected |>
  dplyr::group_by(tissue, timepoint, feature, position_corrected) |>
  dplyr::summarise(signal = mean(signal), .groups = "drop") |>
  ggplot(aes(x = position_corrected, y = signal, colour = timepoint, group = timepoint)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_vline(xintercept = c(0, 3), linetype = "dashed", linewidth = 0.4) +
  facet_grid(tissue ~ feature, scales = "free_x") +
  scale_x_continuous(breaks = -30:40, labels = function(x) ifelse(x %% 3 == 0, x, "")) +
  labs(x = "Corrected position relative to codon (nt)", y = "Mean normalised signal") +
  theme_classic(base_size = 12)
```

Periodicity: the fraction of P-sites per frame. Without periodicity
every frame has 1/3 (dashed line).

``` r
ggplot(res$frames, aes(x = factor(readlength), y = fraction, fill = factor(frame))) +
  geom_boxplot(outlier.size = 0.8) +
  geom_hline(yintercept = 1/3, linetype = "dashed") +
  labs(x = "Read length (nt)", y = "Fraction of P-sites in the CDS", fill = "Frame") +
  theme_classic(base_size = 12)
```

The read length QC has its own plot functions (`plot_readlength_line()`,
`plot_readlength_heatmap()` and `plot_readlength_qc()`). They return
ggplot objects that you can change further.

## Function reference

### Stage 1: BAM files

| Function | Purpose and arguments |
|----|----|
| `list_bam_samples(base_dir, pattern = "\\.bam$")` | Finds the BAM files (subfolders too). Returns `filepath`, `filename`, `sample_name` |
| `bam_to_5prime_counts(bam_file, unique_only = FALSE, check_complete = TRUE)` | One BAM to a table of 5’ ends (chr, pos, strand, readlength, count). Keeps all read lengths. Reports multi-mappers (NH \> 1), drops them with `unique_only`, stops on a truncated BAM |
| `extract_5prime(samples, unique_only = FALSE, out_dir = NULL)` | Runs it for all samples. Failed samples are skipped with a message |
| `read_5prime_tables(dir)` | Loads the saved tables into a named list |

### Read length QC

| Function | Purpose and arguments |
|----|----|
| `readlength_distribution(counts, samples = NULL)` | Data frame with `n` and `pct` per sample and read length, plus your columns |
| `suggest_readlengths(qc, method = "cumulative", cutoff = 0.9, n = 5, contiguous = TRUE)` | Suggests read lengths: the most common ones until they cover `cutoff` of the reads (`method = "top_n"`: the `n` most common). `contiguous` fills gaps. Only a suggestion |
| `plot_readlength_line(qc, color_by = "sample_name", highlight = NULL, length_range = NULL)` | Line plot of % of reads per read length |
| `plot_readlength_heatmap(qc, label_by = "sample_name", order_by = NULL, highlight = NULL, length_range = NULL, high_colour = "#17848E")` | Heatmap with one row per sample. `label_by` must identify every sample uniquely |
| `plot_readlength_qc(qc, color_by, label_by, order_by, highlight, length_range, title, subtitle)` | Both plots in one figure |

`highlight` marks a range of read lengths with dashed lines, for example
the output of `suggest_readlengths()`. `length_range = c(24, 40)` shows
only part of the read lengths.

### Annotation

| Function | Purpose and arguments |
|----|----|
| `transcripts_with_tag(annotation, tag)` | IDs of the transcripts that carry a tag |
| `prepare_annotation(annotation, prefer = NULL, fallback = "longest_cds", seed = 1)` | One transcript per gene. Returns `exons_tx` and `tx_info`. `fallback = "random"` uses the seed |
| `build_transcript_table(annotation)`, `select_transcripts(tx_info, prefer, fallback, seed)` | The two parts of `prepare_annotation()` |

### Steps of the pipeline (each can be run alone)

| Function | Purpose and arguments |
|----|----|
| `map_to_transcripts(counts, exons_tx)` | 5’ ends to transcript positions: transcript_id, tx_pos, readlength, count |
| `build_signal_df_tx(mapped, tx_info, readlengths, features, upstream, downstream, threshold, common_set, samples)` | Profiles around the codons: sample_name, feature, readlength, position, signal |
| `normalize_mat(mat, threshold, label)` | Filters and normalizes the sites (called inside the function above) |
| `find_psite_offsets(signal_df, group_cols = "sample_name")` | Offset per read length: the highest peak before the start codon in each sample, then the most common position. Prints the votes |
| `apply_psite_correction(signal_df, offsets)` | Adds `peak_position` and `position_corrected` (position minus offset) |
| `psite_frames(mapped, tx_info, offsets, skip_codons = 5, samples = NULL)` | P-sites per frame of the CDS |
| `periodicity_score(frames)` | Per sample and read length: number of P-sites, share in frame 0, dominant frame and a score from 0 (no periodicity) to 1 (all in one frame) |
| `select_readlengths(offsets, score, min_agree = 0.9, min_frame0 = 0.8)` | Applies the two criteria and returns the read lengths that pass |
| `metagene_pipeline(...)` | Runs everything above (see the settings section) |

## How it works

- **5’ ends:** for every read only the first nucleotide in transcription
  direction is kept, together with the read length.
- **Transcript coordinates:** every 5’ end is converted to a position
  along its transcript, with introns removed. The 5’ UTR and CDS lengths
  tell where the start and stop codons are, so a window around a codon
  is always one continuous stretch of the transcript. A window is used
  only if it fits inside the transcript: for a start codon the 5’ UTR
  must be at least `upstream` and the CDS at least 3 + `downstream`. For
  a stop codon the CDS must be at least `upstream` + 3 and the 3’ UTR at
  least `downstream`.
- **Offset:** a footprint’s 5’ end lies a fixed distance upstream of the
  ribosome’s P-site, and the distance depends on the read length. For
  each sample and read length, the highest peak before the start codon
  gives the offset (for example -12). The samples vote, and the most
  common position is used for all samples.
- **Periodicity:** a P-site is the 5’ end shifted by the offset. In the
  CDS, true P-sites fall mostly on the first nucleotide of a codon
  (frame 0). The share of P-sites per frame therefore checks the offsets
  from a second angle.

## What good results look like

In our mouse kidney and liver data (28 samples), the main footprint
lengths of 27-29 nt mapped to the selected transcripts for about 94-95%
of their reads, and their offsets were -11 or -12 with at least 93% of
the samples agreeing. Frame 0 was the dominant frame for 27, 28 and 29
nt. Read lengths outside this range had unclear offsets or little
periodicity, and the automatic selection dropped them. Your data may
look different, and the selection table printed by the pipeline shows
why a read length was kept or dropped.

## Notes and limitations

- One transcript per gene is an approximation. A read from another
  isoform that is spliced differently near the start codon can get a
  slightly wrong distance.
- Start and stop codons are taken from `start_codon` and `stop_codon`
  rows of the annotation, so annotations without these rows are not
  supported yet.
- The offset comes from the start codon and is the same for all samples.
  Offsets per group are not implemented. If the votes of your samples
  disagree in a systematic way, check the vote table first.
- Soft-clipped bases are part of the read length. In our data this
  affected about 1.6% of the reads.
- The offset is shared per read length. A read length with scattered
  votes or no frame 0 dominance has no reliable offset and is left out
  by the selection.

## Troubleshooting

- **“BAM file looks truncated”:** the file is incomplete (for example a
  failed download). `extract_5prime()` skips it and reports it.
  Re-download the file.
- **“None of the chromosome names … are in the annotation”:** the BAM
  and the annotation use different names (for example `chr1` and `1`).
- **“These samples are missing in ‘samples’”:** a name in
  `names(counts)` has no row in your sample table.
- **An error about a package called `Seqinfo`:** an annotation object
  was saved on a machine with a newer Bioconductor. Import the GFF3/GTF
  file again on the machine where you use it.

## Status

This package is under development. Older functions based on genomic
windows are being retired in favour of the transcript-coordinate
functions described here.
