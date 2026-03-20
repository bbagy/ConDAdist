# ConDA-dist

`ConDA-dist` is a microbiome differential abundance framework that supports both:

- single-method DA runs with standardized outputs
- multi-method consensus scoring with beta-diversity-guided evidence integration

The main entrypoint is `Go_ConDaDist()`.

## What ConDA-dist does

`Go_ConDaDist()` can be used in two modes.

- Single mode:
  Run one DA method and keep its native behavior as much as possible.
- Consensus mode:
  Run multiple DA methods, summarize agreement across methods, optionally add beta-diversity contribution evidence, and export a final ranked table.

Supported DA methods:

- `deseq2`
- `aldex2`
- `ancombc2`
- `maaslin2`
- `corncob`

Supported distance metrics:

- `bray`
- `jaccard`
- `aitchison`

Why distances matter:

- DA methods detect feature-level abundance shifts
- distance-based beta evidence asks whether a feature also contributes to
  between-group community separation
- this gives `ConDA-dist` a second axis of evidence beyond DA p-values alone

`distances = NULL` turns off beta-diversity and runs DA-only mode.

## Design Principles

- Keep single-method runs close to the older Go DA family workflow
- Use one shared output schema for consensus integration
- Export all important results as CSV files
- Keep plotting decoupled through bridge files and `Go_volcanoPlot`
- Use a directory-first workflow rather than returning large in-memory objects

## Installation

This repository is currently structured like a development package.

```r
devtools::load_all("/path/to/ConDA-dist")
```

Or source files directly in an interactive workflow.

```r
source("ConDA-dist/R/Go_utils.R")
source("ConDA-dist/R/Go_DA_adapters.R")
source("ConDA-dist/R/Go_Consensus.R")
source("ConDA-dist/R/Go_IO.R")
source("ConDA-dist/R/Go_ConDaPlot.R")
source("ConDA-dist/R/Go_ConDaDist.R")
```

## Input

The main input is a standard `phyloseq` object.

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj"
)
```

Key arguments:

- `group_var`: metadata column defining the comparison groups
- `group_1`: baseline/reference group
- `group_2`: one or more target groups
- `project`: project name used for dated output folders
- `covariates`: optional fixed-effect metadata covariates
- `random_effects`: optional random-effect metadata variables
- `methods`: DA methods to run
- `distances`: up to 3 distance metrics per run

## Single-Method Mode

Single mode is triggered when:

- only one DA method is supplied
- `distances = NULL`

Example: `DESeq2`

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj",
  methods = "deseq2",
  distances = NULL
)
```

Example: native `MaAsLin2`

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = c("GLP-2", "D7", "D14"),
  project = "DemoProj",
  methods = "maaslin2",
  distances = NULL,
  covariates = c("Age", "Sex"),
  random_effects = c("SubjectID")
)
```

Notes:

- `maaslin2` single mode uses the native `MaAsLin2` workflow
- `ancombc2` supports mixed-effects through `random_effects`
- single mode returns the output directory path invisibly

## Consensus Mode

Consensus mode runs when:

- more than one DA method is provided
- and/or distance metrics are provided

Example:

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj",
  methods = c("deseq2", "aldex2", "ancombc2", "maaslin2", "corncob"),
  distances = c("bray", "jaccard", "aitchison")
)
```

Consensus combines:

- standardized DA outputs from all requested methods
- Fisher-style combined evidence
- direction and effect consistency
- beta-diversity contribution evidence when distances are enabled

## Method Signature Rules

Output directories use a method signature.

- Single method:
  full method name
- Multi-method:
  fixed-order initials in `D A N M C`

Examples:

- `deseq2`
- `aldex2`
- `DMC`
- `AC`
- `DANMC`

This avoids overwriting outputs from different method combinations.

## Output Layout

Outputs are organized under a dated project directory:

```text
project_YYMMDD/
├── table/
│   ├── ConDaDist/
│   │   └── <signature>/
│   │       └── <group1.vs.group2>/
│   └── ConDaDist_volcano/
└── pdf/
    ├── DA_plot/
    └── ConDa_plot/
```

In practice, the two most important directories are:

- `table/ConDaDist/`
  main analysis tables
- `table/ConDaDist_volcano/`
  bridge CSV files for `Go_volcanoPlot`

For a typical consensus comparison:

```text
DemoProj_260319/
├── table/
│   ├── ConDaDist/
│   │   └── DANMC/
│   │       └── Control.vs.GLP-2/
│   │           ├── DemoProj.filtered_feature_table.csv
│   │           ├── DemoProj.all_methods_standardized.csv
│   │           ├── DemoProj.feature_consensus_summary.csv
│   │           ├── DemoProj.beta_summary.csv
│   │           ├── DemoProj.beta_feature_contribution.csv
│   │           └── DemoProj.final_consensus_scores.csv
│   └── ConDaDist_volcano/
│       ├── condadist.DANMC.(Control.vs.GLP-2.DemoProj).volcano_bridge.csv
│       └── deseq2.(Control.vs.GLP-2.DemoProj).volcano_bridge.csv
└── pdf/
    ├── DA_plot/
    └── ConDa_plot/
```

For a single-method run, the signature directory uses the full method name:

```text
DemoProj_260319/
└── table/
    └── ConDaDist/
        └── deseq2/
            └── Control.vs.GLP-2/
```

## Typical Workflow

```r
# 1. Run analysis
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj",
  methods = c("deseq2", "aldex2", "ancombc2"),
  distances = c("bray", "jaccard", "aitchison"),
  volcano_plot = TRUE
)

# 2. Inspect exported tables
list.files(res_dir)

# 3. Review volcano output
project_dir <- normalizePath(file.path(res_dir, "..", "..", "..", ".."))
list.files(file.path(project_dir, "pdf", "ConDa_plot"))
```

This directory-first workflow is intentional: `ConDA-dist` writes its main
results to CSV files and returns the output path invisibly.

## Volcano Plot Integration

`ConDA-dist` exports bridge CSV files that can be read by the older Gotools `Go_volcanoPlot`.

Enable this during a run:

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj",
  methods = c("deseq2", "aldex2"),
  distances = NULL,
  volcano_plot = TRUE
)
```

Behavior:

- single DA method bridges are rendered into `pdf/DA_plot/`
- consensus bridges are rendered into `pdf/ConDa_plot/`

This keeps `Go_volcanoPlot` compatible with `ConDA-dist` outputs.

## QC Plot Integration

QC plots can be created automatically:

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj",
  qc_plot = TRUE
)
```

QC plots are intended as diagnostic summaries, not publication-ready final figures.

Current QC outputs include:

- consensus volcano
- DA vs beta scatter
- method overlap summary
- top taxa summary

## Distance Guidance

Use distances when you want consensus ranking to reflect not only differential
abundance, but also contribution to overall community structure differences.

At most 3 distances are allowed per run.

Recommended combinations:

- tree-free data:
  `c("bray", "jaccard", "aitchison")`
- when phylogenetic distances are needed:
  choose a smaller representative set per run

Why the limit exists:

- beta contribution is the most computationally expensive step
- feature-level leave-one-feature-out scoring scales quickly with both feature count and number of distances

When to set `distances = NULL`:

- when you only want a single-method DA result
- when you want a fast DA-only consensus run
- when beta-diversity contribution is not part of the current question

## Return Value

`Go_ConDaDist()` returns an output directory path invisibly.

- for one comparison:
  the comparison directory
- for multiple comparisons:
  the method-signature directory containing all comparison subdirectories

Example:

```r
res_dir <- Go_ConDaDist(...)
res_dir
```

## Current Status

Single-mode runs have been exercised for:

- `deseq2`
- `aldex2`
- `ancombc2`
- `maaslin2`
- `corncob`

Consensus smoke tests have been exercised across all non-single combinations of:

- `D`
- `A`
- `N`
- `M`
- `C`

## Recommended Workflow

1. Start with a single-method run to inspect native behavior.
2. Run a consensus configuration that matches your biological question.
3. Review exported tables first.
4. Use `Go_volcanoPlot` bridge outputs for standardized volcano plots.
5. Use QC plots only as diagnostic summaries.

## Development Notes

`ConDA-dist` is still under active development. The current implementation focuses on:

- stable single-mode behavior
- consensus orchestration
- backward-compatible output contracts
- Git-friendly exported tables and plotting bridges

The most important interface to keep stable is:

- `Go_ConDaDist()`
- exported CSV tables
- volcano bridge CSV files consumed by `Go_volcanoPlot`
