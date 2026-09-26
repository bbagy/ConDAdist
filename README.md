# ConDA-dist

`ConDA-dist` is a microbiome differential abundance framework for single-method
analysis and multi-method consensus scoring, with optional distance-guided
evidence integration.

- single-method DA runs with standardized outputs
- multi-method consensus scoring with beta-diversity-guided evidence integration

The main entrypoint is `Go_ConDaDist()`.

## V5 configuration presets

The package remains `ConDAdist`; the current internal consensus engine is V5.
The public API exposes one frozen preset and one user-defined mode:

- `preset = "broad_panel"` (default): five DA tests grouped into four method
  families, no distance contribution, and all-but-one family partial conjunction.
- `preset = "custom"`: user-selected methods, distances, weights, or combiner;
  custom combinations are not claimed to be independently calibrated.

The Stage-X exhaustive-search candidate is retained only as a reported negative
benchmark result and is not exposed as a named preset. Explicitly supplying
`methods`, `distances`, `weights`, or `p_combine` without a
`preset` automatically selects `custom`, preserving older calls.

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
- `corncob_wald`
- `corncob_lrt`

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

### Covariate adjustment and restricted permutations

When distances are enabled, `covariates` are used in both parts of the
workflow. DA adapters receive the requested fixed effects, and the distance
layer estimates the marginal group effect after those covariates. Its
leave-one-feature-out score is the change in the covariate-adjusted group
partial R-squared, rather than an unadjusted separation score.

Use `strata` for paired or clustered permutation designs. It names one
metadata column that defines exchangeability blocks for PERMANOVA and the
feature-level beta permutation test. `strata` restricts label permutations;
it does not replace a fixed covariate or a mixed-effects model.

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "Case",
  project = "AdjustedDemo",
  covariates = c("Age", "Sex"),
  strata = "SubjectID",
  distances = "bray"
)
```

## Design Principles

- Keep single-method runs close to the older Go DA family workflow
- Use one shared output schema for consensus integration
- Export all important results as CSV files
- Keep plotting decoupled through bridge files and `Go_volcanoPlot`
- Use a directory-first workflow rather than returning large in-memory objects

## Installation

### 1. Install Rust (required for ANCOMBC)

ANCOMBC depends on `CVXR`, which depends on `clarabel` — a package that must
be compiled from Rust source. Without Rust, ANCOMBC installation will fail.

Run this once in your **Terminal** (not in R):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After the installer finishes, reload your shell environment:

```bash
source ~/.cargo/env
```

Then **restart R** so it picks up the updated PATH.

If Rust is already installed but R cannot find `cargo`, add this to your
`~/.Renviron` or run before installing:

```r
Sys.setenv(PATH = paste0(path.expand("~/.cargo/bin"), ":", Sys.getenv("PATH")))
```

### 2. Install libomp (required for terra, macOS only)

Some Bioconductor packages pull in `terra` as a transitive dependency.
`terra` requires OpenMP to compile on macOS, which is not included with
Apple's default Xcode toolchain.

Run this once in your **Terminal**:

```bash
brew install libomp
```

Then configure R to find it by adding these lines to `~/.R/Makevars`
(create the file if it does not exist):

```
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp
```

You can do this from R:

```r
dir.create("~/.R", showWarnings = FALSE)
write(c(
  "LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp",
  "CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp"
), file = "~/.R/Makevars", append = TRUE)
```

Then **restart R**.

### 3. Install R dependencies

Load ConDA-dist, then run the dependency installer once:

```r
source("ConDA-dist/R/Go_utils.R")
# ... source remaining files ...
condadist_dependency()
```

This installs all required Bioconductor and CRAN packages automatically,
including ANCOMBC, DESeq2, ALDEx2, Maaslin2, corncob, phyloseq, and vegan.

### 3. Load ConDA-dist

```r
devtools::load_all("/path/to/ConDA-dist")
```

Or source files directly:

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
- `orders`: optional ordered levels used for full pairwise runs
- `pairwise_all`: if `TRUE`, run all pairwise contrasts in `orders`
- `filter_scope`: feature filtering scope for DA runs.
  - `"pairwise"`: subset to the current `group_1` vs `group_2` comparison first, then calculate prevalence/abundance filtering on that subset.
  - `"global"`: calculate prevalence/abundance filtering on the full aligned dataset before running pairwise DA.

## Mode Summary

`ConDA-dist` behaves differently depending on the combination of `methods` and
`distances`.

- `single + NULL`
  pure single-method mode
- `single + dist`
  one-method DA result plus beta-diversity evidence
- `multi + NULL`
  DA-only multi-method consensus
- `multi + dist`
  full ConDA mode

### Consensus inference default (V5)

The conservative `broad_panel` default uses five tests grouped into four
method families:

```r
methods = c(
  "deseq2", "aldex2", "ancombc2",
  "corncob_wald", "corncob_lrt"
)
p_combine = "family_partial_conjunction"
```

Within any V5 panel, related corncob Wald and LRT tests are first collapsed
into one Bonferroni family p-value. For `m` planned families, CDD tests the
all-but-one partial-conjunction hypothesis with `h = max(1, m - 1)` and
`min(1, (m - h + 1) * p_(h))`, then applies BH across features. Missing planned
tests occupy their slots with `p = 1`; one extreme method therefore cannot
determine a multi-family consensus result by itself.

`p_combine = "adaptive_cauchy"` remains available only for explicit legacy
or exploratory reproduction. It is not the production consensus default.
Custom multi-method subsets may use the same generalized rule, but are labeled
as not independently calibrated rather than inheriting the Stage-X claim.

Pure single-method and one-method-plus-distance runs retain their existing V2
behavior automatically because no cross-method p-value consensus is possible.

The `single + dist` case is especially useful when a user strongly prefers one
DA method but still wants to reinterpret that result with an additional
community-structure axis. In that setting, the preferred DA tool is kept, while
beta-diversity contribution is added as complementary evidence.

Comparison behavior:

- `pairwise_all = FALSE`
  run `group_1` against each value in `group_2`
- `pairwise_all = TRUE`
  ignore `group_1` and `group_2`, and run all pairwise contrasts in `orders`

Filtering behavior:

- `filter_scope = "pairwise"`
  default. Recommended for baseline-vs-target microbiome DA because feature filtering is evaluated within the current pairwise subset.
- `filter_scope = "global"`
  optional. Keeps one shared feature universe across comparisons, but may retain sparse features that are not well-supported within a given pairwise contrast.

Why this matters:

- In sparse microbiome data, global filtering can keep features that pass prevalence in the full cohort but are nearly absent within a specific pairwise contrast.
- Those features can inflate instability in DA methods, especially effect-size-heavy outputs such as `DESeq2` volcano plots.
- `pairwise` filtering keeps the tested feature universe closer to the actual comparison question.

### When to use `filter_scope = "global"`

For most analyses `"pairwise"` is the correct choice. Use `"global"` only when
one of the following applies:

**1. Multi-group `pairwise_all = TRUE` runs where feature universe consistency is required**

When running all pairwise contrasts with `pairwise_all = TRUE` (e.g., A vs B,
A vs C, B vs C), each pairwise subset may produce a different feature set
under `"pairwise"` filtering. If downstream comparisons of effect sizes or
rankings across contrasts must be based on an identical feature universe,
`"global"` ensures a single consistent set is used throughout.

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  orders = c("Control", "Low", "High"),
  pairwise_all = TRUE,
  filter_scope = "global"   # same features in all three contrasts
)
```

**2. Very small pairwise sample sizes**

When group sizes drop to 3–5 samples per arm, `"pairwise"` filtering can be
overly aggressive and discard biologically relevant features that happen to
appear at low prevalence in the small subset. `"global"` uses the full cohort
prevalence, which is more stable at small N.

**3. Reviewer or publication requirements for a consistent feature universe**

Some journals or analysis platforms require that the same features are reported
across all subgroup comparisons. In that case, `"global"` ensures that every
comparison table shares the same row set.

**Practical summary:**

| Scenario | Recommended `filter_scope` |
|---|---|
| Standard two-group comparison | `"pairwise"` (default) |
| All-pairwise multi-group run | `"global"` for cross-contrast consistency |
| Very small n (< 5 per group) | `"global"` to avoid over-filtering |
| Publication requiring identical feature lists | `"global"` |

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
  distances = NULL,
  filter_scope = "pairwise"
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
- if one DA method is supplied together with `distances`, the run becomes a
  one-method-plus-beta mode rather than pure single mode

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
  distances = c("bray", "jaccard", "aitchison"),
  filter_scope = "pairwise"
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
│   └── ConDaDist_plot_Tab/
└── pdf/
    ├── DA_plot/
    └── ConDa_plot/
```

In practice, the two most important directories are:

- `table/ConDaDist/`
  main analysis tables
- `table/ConDaDist_plot_Tab/`
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
│   └── ConDaDist_plot_Tab/
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
  distances = c("bray", "jaccard", "aitchison")
)

# 2. Inspect project outputs
list.files(res_dir)

# 3. Review volcano output
list.files(file.path(res_dir, "pdf", "ConDa_plot"))
```

This directory-first workflow is intentional: `ConDA-dist` writes its results
to CSV/PDF/HTML files beneath a dated project directory and returns that root
path invisibly.

## Volcano Plot Integration

`ConDA-dist` exports bridge CSV files automatically so they can be read by the
older Gotools `Go_volcanoPlot`.

```r
res_dir <- Go_ConDaDist(
  psIN = ps,
  group_var = "TreatmentGroup",
  group_1 = "Control",
  group_2 = "GLP-2",
  project = "DemoProj",
  methods = c("deseq2", "aldex2"),
  distances = NULL,
  filter_scope = "pairwise"
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

## Method Settings

`ConDA-dist` keeps one shared interface across DA methods, but each method still
uses its own default controls. These defaults are chosen to stay reasonably
close to microbiome use cases while avoiding overly aggressive transformation or
filtering inside the adapter itself.

### Native Adapter Thresholds

Before running a native method, `ConDA-dist` checks that there are enough
samples and features to support it. If the check fails, that method is skipped
and the reason is reported (see Failure Behavior below).

| Method    | Min samples per group | Min features |
|-----------|-----------------------|--------------|
| `ancombc2`| 3                     | 10           |
| `deseq2`  | 3                     | 5            |
| `corncob` | 3                     | 10           |
| `aldex2`  | 2                     | 2            |
| `maaslin2`| 2                     | 2            |

Note: group size is counted from the actual pairwise metadata after sample
alignment. If `group_var` is stored as a factor with additional levels not
part of the current comparison, those extra levels are ignored correctly.

### `ancombc2`

Default controls:

- `prv_cut = 0.1`
- `lib_cut = 1000`
- `p_adj_method = "BH"`
- `pseudo = 0`
- `pseudo_sens = FALSE`
- `struc_zero = TRUE`
- `neg_lb = TRUE`
- `global = TRUE`
- `pairwise = FALSE`
- `dunnet = FALSE`
- `trend = FALSE`

Interpretation:

- Uses ANCOM-BC2 as the primary microbiome-native DA method.
- Handles compositional bias through log-linear bias correction.
- Structural zeros are enabled by default, which is important for sparse data.
- Pairwise contrasts in `ConDA-dist` are handled by the outer comparison loop,
  so the adapter itself keeps `pairwise = FALSE`.
- Supports covariates and random effects via `fix_formula` and `rand_formula`.
- If the first run fails, the adapter automatically retries with incremental
  legacy relative-mean filtering (0.001 to 0.01) until convergence.

### `aldex2`

Default controls:

- `mc_samples = 128`
- `denom = "iqlr"`
- `use_mc = FALSE`
- `paired_test = FALSE`
- `zero_replace = FALSE`
- `zero_replace_value = 0.5`
- `seed = 1`

Interpretation:

- CLR/Dirichlet-based compositional DA method.
- Uses `iqlr` as the default denominator for robustness to asymmetric
  differential abundance, which is common in microbiome data.
- In simple two-group runs with no covariates, the adapter uses the native
  ALDEx2 Wilcoxon t-test path (`wi.ep`, `wi.eBH`).
- When covariates are present, the adapter automatically switches to the
  ALDEx2 GLM path.
- Seed is fixed at `1` by default for reproducibility of the Monte Carlo draws.

### `maaslin2`

Default controls:

- `min_abundance = 0`
- `min_prevalence = 0`
- `normalization = "TSS"`
- `transform = "LOG"`
- `analysis_method = "LM"`
- `max_significance = 1`
- `standardize = FALSE`

Interpretation:

- Linear mixed model DA with TSS+LOG normalization.
- Filtering (`min_abundance = 0`, `min_prevalence = 0`) is intentionally
  set to zero so that feature selection is controlled by the upstream
  `Go_FilterFeatures` step, not repeated inside the adapter.
- Supports covariates as fixed effects and `random_effects` for repeated
  measures or nested designs.
- Single-mode MaAsLin2 (`methods = "maaslin2"`, `distances = NULL`) uses a
  special native workflow that supports multi-level `group_2` comparisons.

### `corncob`

Default controls:

- `phi_formula = "auto"`
- `phi_null_formula = "auto"`
- `boot = FALSE`
- `fdr_cutoff = 1`
- `filter_discriminant = TRUE`
- `min_prevalence = 0.05`
- `min_total_count = 10`
- `legacy_filter_cutoff = 0`

Retry controls (applied automatically on first failure):

- `retry_min_prevalence = 0.1`
- `retry_min_total_count = 20`
- `retry_legacy_filter_cutoff = 0.001`
- `retry_phi_formula = "~ 1"`
- `retry_phi_null_formula = "~ 1"`

Interpretation:

- Beta-binomial DA model that handles overdispersion directly.
- Most stability-sensitive method in the set. Convergence can fail on very
  sparse or unbalanced data.
- `phi_formula = "auto"` uses the same formula as the mean model by default.
- The adapter applies a two-attempt retry strategy: if the first Wald test
  fails, it retries with stricter feature filtering and a simplified
  dispersion formula (`~ 1`), which improves convergence in difficult cases.
- In single-method mode, additional stability-focused defaults are applied
  automatically (stricter prevalence and count thresholds).

### `deseq2`

Default controls:

- `min_count = 0`
- `zero_replace = FALSE`
- `zero_replace_value = 1`
- `size_factors_type = "poscounts"`
- `lfc_shrink = TRUE`
- `lfc_shrink_type = "ashr"`

Size factor fallback chain:

`ConDA-dist` uses a robust DESeq2 wrapper (`Go_RunDESeqRobust`) that attempts
size factor estimation in order until one succeeds:

1. `poscounts` (default) — geometric mean of non-zero counts; more stable
   for sparse microbiome tables than the standard ratio method
2. `iterate` — iterative estimation; used when `poscounts` fails

This fallback chain is recorded in the `notes` column of the standardized
output (`size_factor_type=<method_used>`).

Interpretation:

- Uses `poscounts` size-factor estimation by default, which is more stable
  than the basic ratio method for sparse microbiome count tables.
- Applies log-fold-change shrinkage (`ashr`) by default to reduce extreme
  LFC artifacts caused by low-count or sparse features.
- LFC shrinkage is applied after `DESeq2::results()` using
  `DESeq2::lfcShrink()` with `type = "ashr"`. If `lfcShrink` fails, the
  unshrunken LFC is used and a warning is printed.
- Still best interpreted as a supportive count-model DA method. It is not
  microbiome-native but is widely used as a reference in benchmarks.

## Failure Behavior

When a native DA method cannot run — due to insufficient sample size, too few
features, or a runtime error — `ConDA-dist` skips that method instead of
silently substituting a different statistical test.

Skipped-method reporting:

- `method` retains the requested method name
- p-values, q-values and effect sizes are `NA`
- `notes` column records which method failed and why
- A `[ConDA] WARNING` message is printed to the console

Example console output:

```
[ConDA] WARNING: DESeq2 could not run and was skipped. Reason: Native deseq2 skipped: smallest group has 3 sample(s).
```

Skipped methods retain their planned consensus slot but contribute no evidence;
their missing family-level p-value is handled conservatively as `p = 1`.
No substitute Wilcoxon CSV or volcano plot is created.

This design ensures that:

- A successful-looking substitute analysis cannot hide native method failure.
- Results remain explicit about which requested methods actually ran.
- Consensus scoring is not polluted by evidence from a different test.

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
