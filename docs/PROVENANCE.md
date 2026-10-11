# Code provenance

This release contains recovered project files, commands retrieved from historical conversations, and new code that organizes them into a portable workflow. These have different levels of historical evidence.

## Stage-by-stage record

| Stage | Evidence recovered | Publication treatment |
| --- | --- | --- |
| Download and FASTQ conversion | A 4 May 2026 pilot block with `prefetch`, `fasterq-dump --split-files --outdir ... --temp ... -e 4`, and `pigz -p 4`; later records confirm local SRA files | Adapted to a caller-supplied sample sheet and per-sample directories; not the recovered final 38-sample launcher |
| Quality control and host removal | June local and July cloud KneadData commands, paired `-i1/-i2`, `-db`, and `SLIDINGWINDOW:4:20 MINLEN:50` | Retains those analysis arguments; adds an explicit output prefix and configurable thread count |
| Taxonomic profiling | July corrected MetaPhlAn invocation with `--db_dir`, `-x`, `--mapout`, and paired-file input | Retains the corrected CLI form and fixed database index; paths and threads are configurable |
| Functional profiling | Workflow records, database setup records, and saved HUMAnN matrix identities; final full production invocation unavailable | Reconstructed from documented HUMAnN interfaces; concatenates cleaned mates, reuses the per-sample taxonomic profile, and retains nucleotide and translated search |
| Merging | Historical reference to `merge_metaphlan_tables.py` and partial `humann_join_tables` calls | Full staging/merge code newly organized for release |
| Normalization and splitting | Saved matrix names identify gene CPM, pathway relative abundance, and unstratified outputs; exact commands unavailable | Explicit reconstructed `community` normalization with special rows included, then stratification splitting; historical parameter equivalence is not asserted |
| KO/EC regrouping | Complete local Shell script and historical completion discussion | Preserves the original script; portable version requires explicit custom mapping files and records their hashes |
| MicrobiomeAnalyst formatting | Complete Python helper and cloud Shell wrapper | Preserves originals and supplies an adapted converter with stricter validation |

## Preserved files

The following are byte-for-byte copies of recovered saved files, with the latest local KO/EC filename simplified for the repository:

| Repository filename | Recovered source filename |
| --- | --- |
| `archive/original_scripts/prepare_microbiomeanalyst_inputs.py` | `prepare_microbiomeanalyst_inputs.py` |
| `archive/original_scripts/run_prepare_microbiomeanalyst_inputs.sh` | `run_prepare_microbiomeanalyst_inputs.sh` |
| `archive/original_scripts/run_local_KO_EC.sh` | `run_local_KO_EC(3).sh` |

File hashes are in `archive/original_scripts/SHA256SUMS.txt`. These files contain the original local/cloud paths, output naming, and software-discovery behavior. They are historical reference copies. The portable entry point uses the separate code under `scripts/`.

The old wrappers may search local storage, download mapping resources, create an environment if explicitly enabled, and overwrite their named archives. The new orchestrator uses explicit inputs and fresh stage directories. It does not copy those machine-management behaviors.

## Recovered command fragments

These are transcriptions of retrieved command fragments with filesystem locations represented by variables. They are not complete historical batch scripts.

```bash
# Download/conversion pilot; May 2026
prefetch --option-file "$RUN_LIST" -O "$SRA_DIR" --max-size 100g
fasterq-dump "$SRA_DIR/${srr}/${srr}.sra" --split-files --outdir "$FASTQ_DIR" --temp "$TMP_DIR" -e 4
pigz -p 4 "$FASTQ_DIR"/*.fastq

# Cloud quality control; July 2026
kneaddata -i1 "fastq/${SRR}_1.fastq.gz" -i2 "fastq/${SRR}_2.fastq.gz" \
  -db "$HOST_DB" -o "kneaddata_out/${SRR}" -t 12 \
  --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50"

# Corrected MetaPhlAn call; July 2026
metaphlan "${CLEAN1},${CLEAN2}" --input_type fastq --nproc "$META_THREADS" \
  --db_dir "$MPA_DB" -x "$MPA_INDEX" \
  --mapout "metaphlan_out/${SRR}.bowtie2out.bz2" \
  -o "metaphlan_out/${SRR}_metaphlan_profile.txt"
```

The MetaPhlAn correction matters: earlier conversation examples used `--bowtie2out`; the later recorded installation required `--mapout`. The portable preflight checks the available CLI rather than silently switching interfaces. The names `knead1_master` and `knead2_master` in old audit summaries are insufficient evidence of two sequential dehosting passes, so no second pass has been invented.

The release download stage uses one accession at a time, with the `prefetch -O` destination ending in that accession as specified by current SRA Toolkit documentation. This is an adaptation of the pilot fragment above.

## Software versions

June local setup records mention HUMAnN `4.0.0.alpha.1` and MetaPhlAn `4.1.1`. July cloud setup records mention HUMAnN `3.9` and MetaPhlAn `4.2.4`, including an environment named `humann4`. An environment name is not a version. The exact final HUMAnN executable and all database snapshots remain unresolved.

The historical MetaPhlAn index `mpa_vOct22_CHOCOPhlAnSGB_202403` is included in the configuration example. The ChocoPhlAn, protein and annotation mapping paths must be supplied locally and checked for compatibility with the chosen HUMAnN installation. The local KO/EC script sought `map_ko_uniref90.txt*` and `map_level4ec_uniclust90.txt*`; the example preserves those names without assuming that a current download is identical to the original mapping.

## Scope of new code

`scripts/pipeline.py`, `run_pipeline.sh`, configuration templates, and tests were prepared for this publication. `scripts/prepare_inputs.py` adapts the recovered helper, retaining accession alignment and feature-formatting operations while rejecting invalid numeric values and stratified functional input. Execution logs, config snapshots, and completion markers are generated locally during a new run. None are bundled as study results.

The rebuilt HUMAnN and normalization steps are usable workflow code based on official interfaces, but they are not evidence of the exact original invocation or its output. Full raw-read reproduction has not been performed for this release.

## Technical references for reconstruction

- [HUMAnN manual, including paired reads, custom taxonomic profiles, and table utilities](https://github.com/biobakery/humann)
- [HUMAnN normalization options](https://github.com/biobakery/humann/blob/master/humann/tools/renorm_table.py)
- [SRA Toolkit download and conversion](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
- [KneadData paired-end interface](https://github.com/biobakery/kneaddata)
- [EasyMetagenome workflow](https://github.com/YongxinLiu/EasyMetagenome/blob/master/1Pipeline.sh)
