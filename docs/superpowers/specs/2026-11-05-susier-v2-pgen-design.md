# SuSiE 2.0: PGEN inputs and gene batches

Status: approved for implementation on 2026-11-06. Default to 50 genes per batch,
with configurable batch size, PLINK multi-region extraction, sequential gene
fits and separate gene files. The PGEN reader has passed tests and review.

## Intent and agreed constraints

- Target Terra managed Cromwell with WDL 1.0.
- Input PGEN hardcalls rather than Tabix dosage tables.
- Typical cohort: 9,000–10,000 samples and about 30 million variants genome-wide.
- Gene counts vary by molecular trait.
- Run several independent gene analyses in each task, then aggregate summaries.
- Keep each gene's RDS and supporting files separate.
- Preserve the existing statistical settings, compatibility fallback, and
  covariate-explained variant/credible-set diagnostics.
- Use GitHub Actions for image smoke tests. Do not build Docker locally or
  submit Terra jobs without the user's approval.

## Implementation choices

The user approved the following simple batch structure. Implementation details
are engineering choices unless listed as user constraints above.

Use chromosome PGEN/PVAR/PSAM sets as storage inputs. Scatter batches of nearby
genes within each chromosome. Sort genes by
position and split them into batches of 50. Make that size configurable and keep
all outcomes for a gene together. No overlap optimizer or variant counting is
needed to choose these groups.

For each batch, supply every complete cis interval to PLINK 2 and extract their
union into one compact task-local PGEN/PVAR/PSAM set. For inclusive coordinates,
use `--extract bed1 batch.windows.tsv --make-pgen`. Keep all samples and do not
add MAF, missingness, allele-conversion or imputation options. The existing
per-phenotype R code remains responsible for statistical filtering.

Use one R process and fit genes in sequence within a batch. Open the PGEN
handle and parse variant/sample metadata once. Read only the current gene's
complete cis window and release its large working matrices after saving its
files. The first version need not cache genotype blocks across genes. A small
bounded raw-block cache can be added if the pilot shows a useful benefit.
Do not decode a chromosome or a growing batch union into a dense matrix.

## Simple batch planning and extraction

1. Separate analysis cohorts, PGEN inputs and chromosomes. Validate gene/outcome
   coordinates. Sort genes by genomic position, with gene ID as a stable tie
   breaker. Keep all selected outcomes for each gene together.
2. Split each sorted list into consecutive batches of at most 50 genes by
   default. Leave the last smaller batch intact. Expose the batch size as an
   input. The number counts genes, not phenotype fits.
3. Write a task-local interval list containing each selected outcome's complete
   cis window. Use 1-based inclusive bounds and clip the lower bound at 1 for
   the PLINK query. Preserve original requested coordinates in reports. Do not
   replace disjoint windows with one min-to-max span.
4. Use PLINK 2's multi-interval extraction to write the union of matching variants
   to a batch PGEN. Shared variants are retained once per source record, while
   distinct records at the same position remain distinct. Keep the complete
   sample cohort and preserve source variant/allele/sample order. Validate these
   properties with real hardcall round trips before claiming numerical parity.
5. Open that batch set once in R. Run genes in sequence. Decode only the current
   gene/outcome's complete window, apply the existing filtering and fit, write
   its separate files, then release its large matrices.
6. Select memory from resident inputs, the maximum active gene working set and
   a reserve, not the union of batch windows. A compact batch PGEN on disk is
   not a dense matrix in R. Never truncate a gene's window to meet a budget.
7. Assign each gene to exactly one batch. Gene order, batch size and batch ID must
   not change its statistical inputs, seed or output identity.

Map each requested window to global row indices of the currently open PGEN in
the reader. Preserve duplicate positions and original record order. This
mapping is required for correct reads, not for choosing batch membership.

Fixed gene counts do not guarantee equal task times: the number of selected
outcomes can differ. Measure batch duration and localization costs before
adding more scheduling logic. No overlap threshold, raw cross-gene cache or
dynamic programming is required for the initial version.

Any optional genotype cache shares raw hardcalls only. Each gene/outcome keeps its own
sample alignment, filtering, imputation, centering and residualization. In
particular, do not substitute residualized columns from a larger shared union
into the compatibility fallback: different multiplication dimensions can change
rounding. Other caches need exact input/design/sample keys and separate parity
tests before use.

Overlap reduces repeated reads, not the number of SuSiE fits. The measured
single-gene PGEN run spent 2.076 seconds reading and 32.391 seconds fitting,
out of 48.884 seconds in the entry point. Even perfect removal of that read
stage alone would save about 4.2% of that local run. Sharing process/input
startup may add benefit but has not been measured for a batch. Fifty comparable
fits still represent about 27 minutes of fit time. A larger dense-region batch
is appropriate when its total work stays within the chosen budget.

The existing local benchmark took 48.884 seconds and 6.670 GiB for one gene
with 10,000 samples and 20,000 variants. It is a starting point for testing,
not a measured batch/cloud runtime or a memory guarantee. A dense matrix for
10,000 samples and 30 million variants would use 2.4 TB before temporary copies.

## Alternatives

1. One task per chromosome: minimizes repeated input transfer, but limits
   parallel work and produces long, unequal task times. Gene checkpoints help
   retries, but do not improve this load balance.
2. Smaller gene batches with chromosome inputs: recommended first version.
   More parallel work and smaller retries; chromosome input files are localized
   for each batch.
3. Physical regional PGEN files: add only if repeated chromosome-file
   localization dominates. Assign each gene to exactly one core region and
   include enough flanking variants for every assigned gene's full cis window.
   Regional boundaries must never truncate a gene's variant window.

## Workflow structure

1. Sort gene work items by chromosome and position and split into batches of
   50 genes, configurable. Grouping does not require PVAR/PGEN reads. Create
   each batch's interval list inside a task, then run PLINK multi-region
   extraction before the R gene loop. Input files remain typed WDL Files.
2. Scatter batch tasks. Load common phenotype/covariate data once where practical,
   align sample IDs, and run the unchanged per-gene statistical procedure.
3. Save each gene locally before advancing. Return its artifacts as separate
   typed File outputs; Cromwell collects them when the batch succeeds. The first
   version uses native batch retries/call caching. A failed VM reruns that batch;
   task-local gene files are not durable checkpoints. Durable per-gene recovery
   is a possible extension, not an implicit promise of local file writing.
   Keep scientific identity and RNG seeds separate from batch identity.
4. Produce batch summary tables and a complete gene/outcome status index.
5. Aggregate summary tables and indexes. Reject missing/duplicate expected
   outcomes and unreadable input files. No-credible-set outcomes remain in the
   index. Do not aggregate RDS objects across genes.

## Gene artifacts

Each gene has its own directory with a stable safe key and its original ID in
the index. The index records status, original IDs, artifact basenames and hashes;
Cromwell's typed File outputs provide the remote output locations. Keep the
current gene RDS contract: a named list of outcome fits,
including multiple outcomes where applicable. Keep per-gene credible-set,
full-SuSiE and log-Bayes-factor Parquets, the flagged-variant TSV, variant-position
TSV, status/provenance manifest, and logs. The final gene RDS contains only that
gene's outcomes.

Carry `has_covariate_explained_variant` into both per-gene and aggregated
credible-set summaries. Keep diagnostic lists specific to the phenotype's final
sample mask and covariate design.

## Numerical and input parity

- Use the current legacy fitting/preparation path, not the separate checkpointed
  runner's QR/filtering kernel. Reuse storage/control patterns through adapters.
- Preserve cis coordinates, variant/allele order, sample order, phenotype
  selection, MAF/list filters, imputation population, and model settings.
- Read missing hardcalls without imputation. Align each outcome's samples before
  calculating imputation means and covariate adjustments.
- Require an explicit legacy missing-code setting, `NA` or `-1`, for numerical
  parity testing. Map PGEN missing hardcalls to that code before the existing
  variant/MAF filters. Test both settings; synthetic examples use `NA`. PGEN
  cannot recover which missing cells originally used which code in a mixed
  source. Exact mixed-source parity needs the original missing-code information;
  do not infer it or promise that parity from a PGEN alone.
- Validate allele/ploidy support. The tested reader covers diploid biallelic
  hardcalls; multiallelic ALT counts must not be silently collapsed.
- Give each gene/outcome a stable seed independent of task order and batch size.
- Preserve all selected outcomes in each gene RDS. Do not turn task batching into
  a joint model across genes.

## Terra file and image contract

Declare PGEN/PVAR/PSAM, phenotype/covariate files, worklists, and aggregation
inputs as typed File values, including inside structs. Keep File types until
command rendering. Pass named CLI arguments; create any newline file lists in
the task command. Workflow-scope file-writing functions remain forbidden.
Output manifests are actual structured metadata, not JSON CLI argument wrappers.

Read local files with the R PGEN reader. A regional read limits decoded RAM;
it does not avoid localization of the chromosome file. Check input readability
and report unresolved cloud URIs as localization failures.

Use a separate pinned micromamba-based v2 image with pinned dependencies,
including the tested pgenlibr/susieR versions. Log gene start, input dimensions,
variant counts, fallback use, completion, checkpoint/recovery, and failures.
Build/test the image in GitHub Actions only when its contained source changes.

## Validation before choosing defaults

- Compare identical genes run singly and in different batch sizes/orders.
- Compare PGEN hardcalls and all previous numerical fit fields with the current
  reference, with fixed gene/outcome seeds and explicit missing-value policy.
- Test sample reordering, per-outcome missingness, multiple outcomes per gene,
  chromosome/region boundaries, and the compatibility fallback/diagnostics.
- Test failed-batch reruns, corrupt or incomplete outputs, and gene ownership.
- Record all relevant source/image/input identities in provenance. Validate
  successful batch outputs before aggregation; do not aggregate failed attempts.
- Test typed cloud-to-local input rendering and statically reject workflow-level
  file writers. A local syntax check is not Terra execution validation.
- Measure dense and sparse regions, multi-outcome genes, large windows, and
  several work/time budgets. Include batches with more than 50 genes. Record
  decoded variants, cache hits, per-gene fit time and peak memory before fixing
  resource defaults. The complete v2 workflow has not been tested on Terra.

Open sizing detail: actual chromosome PGEN sizes are not yet known. Measure
localization cost before deciding whether physical PGEN regions are necessary.
