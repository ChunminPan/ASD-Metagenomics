# Running the code

The repository is readable without installing anything. Running sequencing analysis requires Linux or WSL, suitable storage and memory, external bioinformatics tools, databases, and a private sample sheet. No actual study data are distributed here.

## Prepare local inputs

From the repository root:

```bash
mkdir -p private
cp templates/sample_metadata.csv private/sample_metadata.csv
cp config/config.example.json config/local.json
```

Edit the private CSV to provide the selected public run accessions and groups. Its required columns are `Run` and `group`. Each accession must appear exactly once. This project-specific workflow expects 18 rows with `ASD_baseline` and 20 rows with `TD_control`. The empty template contains only its header. Cohort selection must be based on source-study metadata, not inferred from accession order.

`private/` and `config/local.json` are ignored by Git. The pipeline requires its `work_dir` to be outside the repository, so generated data and results stay outside the code tree. All relative paths in the JSON are resolved from the repository root, including when the command is launched from another directory.

In `config/local.json`, set the host-removal, MetaPhlAn, ChocoPhlAn, UniRef protein, KO mapping, and EC mapping paths to resources installed on your machine. The values beginning `/path/to/` are placeholders. Thread counts are editable local resource settings, not an assertion about the original final run.

## Software and databases

Python 3.10+ is sufficient for the project utilities and their tests. The full pipeline requires SRA Toolkit (`prefetch`, `fasterq-dump`), `pigz`, KneadData and its dependencies, MetaPhlAn and its merging utility, and HUMAnN with its table utilities.

Follow the official installation instructions linked in [PROVENANCE.md](PROVENANCE.md). This repository does not silently install packages or download annotation databases. The `download` stage explicitly downloads selected sequencing runs from SRA.

If the tools are in different Conda environments, put the existing environment names in the `conda_environments` object. The keys are `sra`, `qc`, `taxonomy`, `function`, and `regroup`. Leave a value empty to use executables already on `PATH`. A nonempty value runs the corresponding tool through `conda run --no-capture-output -n NAME`.

The taxonomic command uses the recovered MetaPhlAn 4.2-style `--db_dir` and `--mapout` flags. The historical index in the example is `mpa_vOct22_CHOCOPhlAnSGB_202403`. Ensure it is installed and appropriate for the chosen environment. The new functional stage uses a documented HUMAnN CLI with an existing taxonomic profile, ChocoPhlAn and UniRef databases. All three need compatible versions. The final historical HUMAnN version was not established, so no environment lockfile is claimed to reproduce it.

The EC mapping name in the example follows the recovered local script. Its mapping members must match the actual gene-family identifiers. A filename match alone does not establish compatibility. The portable regroup stage records mapping-file hashes, while the original mapping snapshots remain unavailable.

## Check and preview

After filling in local paths and the sample sheet, run:

```bash
bash run_pipeline.sh check --config config/local.json
```

This checks the sample sheet, required database paths, executable availability, and key CLI flags. It prints available tool versions. It does not validate full database contents or run profiling. The final line should start with `CHECK_COMPLETE`.

Preview all commands without executing the external tools:

```bash
bash run_pipeline.sh all --config config/local.json --dry-run
```

The preview requires a valid local sample sheet. It prints sample-specific commands and creates no analysis outputs. A successful preview is not evidence that the real bioinformatics tools will complete.

## Execute by stage

```bash
bash run_pipeline.sh download --config config/local.json
bash run_pipeline.sh dehost --config config/local.json
bash run_pipeline.sh taxonomy --config config/local.json
bash run_pipeline.sh function --config config/local.json
bash run_pipeline.sh merge --config config/local.json
bash run_pipeline.sh normalize --config config/local.json
bash run_pipeline.sh regroup --config config/local.json
bash run_pipeline.sh export --config config/local.json
```

Alternatively, `bash run_pipeline.sh all --config config/local.json` runs those eight stages in order. Samples are processed sequentially within each stage; each external tool uses its configured threads. Sequencing downloads and profiling can be expensive in disk space and compute time. Plan resources for the source dataset before starting.

## Files created locally

| Stage directory under `work_dir` | Main outputs |
| --- | --- |
| `01_download` | SRA files and paired compressed FASTQ files |
| `02_dehost` | KneadData output and cleaned paired reads |
| `03_taxonomy` | Per-sample MetaPhlAn profiles and mapping output |
| `04_function` | Per-sample HUMAnN gene and pathway profiles |
| `05_merge` | Combined species, gene-family and pathway tables |
| `06_normalize` | Gene CPM and pathway relative-abundance tables, split by stratification |
| `07_regroup` | KO and EC tables plus mapping hashes |
| `08_export/ready` | Aligned MicrobiomeAnalyst input tables and a validation report |

Each executed stage writes a command log and a config/sample-sheet snapshot. A `COMPLETE.json` marker is written only after that stage succeeds. Downstream stages require matching completion records. The fingerprint covers config, sample-sheet and portable-code contents; it is not a hash of all reads, executable binaries or database files.

Existing stage directories are never automatically overwritten or skipped. If a stage fails, inspect its `commands.log`, preserve or move that incomplete stage directory, and rerun the stage. If code, metadata, or config changes, use a new `work_dir` rather than mixing configurations. Reusing already generated files without these completion markers is supported only by calling the standalone formatter below, not by bypassing upstream checks.

The HUMAnN stage concatenates cleaned paired reads into a temporary FASTQ. The new concatenated copy is removed after both main HUMAnN output files are present. Original SRA, raw FASTQ, and cleaned reads are retained. The portable workflow uses cleaned paired reads; it does not add singleton files. This is a documented release implementation choice, not a recovered statement of every original read-handling detail.

The workflow recognizes both unnumbered HUMAnN gene/pathway filenames and the documented numbered variants (`_2_genefamilies.tsv` and `_4_pathabundance.tsv`). It stops if both versions exist for the same sample instead of choosing one. Recognizing filenames does not establish scientific compatibility between HUMAnN versions.

## Format existing local matrices

To use previously generated data without rerunning sequencing analysis:

```bash
python3 scripts/prepare_inputs.py \
  --metadata private/sample_metadata.csv \
  --pathway /absolute/path/to/pathabundance_relab_unstratified.tsv \
  --species /absolute/path/to/metaphlan_merged_with_sample_headers.tsv \
  --ko /absolute/path/to/ko_raw.tsv \
  --ec /absolute/path/to/ec_raw.tsv \
  --outdir /absolute/path/to/new_prepared_directory
```

Species/KO/EC inputs are optional for this standalone utility. An optional `--gene` argument checks the gene matrix's sample header only. All supplied abundance tables must have explicit sample identifiers. The formatter validates and aligns columns; it does not apply additional normalization, estimate missing values, run hypothesis tests, or reconstruct missing sample labels. The output directory must be empty.

## Tests and downstream analysis

```bash
python3 -m unittest discover -s tests -v
```

The tests create small synthetic fixtures in temporary directories. They do not contain the study's true sample accessions, assignments, or results, and do not run third-party profiling tools.

Subsequent MicrobiomeAnalyst analyses occur outside this code. Record the web-session filtering, transformations, statistical settings and software date separately. No missing historical session settings or significance tests are invented in this repository.
