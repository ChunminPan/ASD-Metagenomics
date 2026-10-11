# Software and analysis settings

This workflow calls external bioinformatics tools. Install compatible tools and databases before running it, and supply their locations in `config/local.json`. The repository contains project scripts, not copies of those tools or databases.

## Tools by stage

| Stage | Software | Main settings or behavior |
| --- | --- | --- |
| Download | SRA Toolkit and pigz | Downloads each accession to its own directory, converts to paired FASTQ files, and compresses them. |
| Quality control and host removal | KneadData | Uses paired reads, a supplied human reference database, and `SLIDINGWINDOW:4:20 MINLEN:50` trimming settings. |
| Taxonomic profiling | MetaPhlAn | Uses paired cleaned reads, a specified database index, and the `--db_dir` and `--mapout` interface. |
| Functional profiling | HUMAnN | Concatenates cleaned mates, supplies the corresponding taxonomic profile, and runs nucleotide and translated searches with the configured databases. |
| Merging | MetaPhlAn and HUMAnN utilities | Combines per-sample profiles into common tables. |
| Normalization | HUMAnN table utilities | Produces gene-family CPM and pathway relative-abundance tables, then splits stratified and unstratified features. |
| Functional regrouping | HUMAnN regrouping utility | Uses explicit KO and EC mapping files and records their hashes. |
| Input preparation | Project Python utility | Validates sample labels, aligns columns, and formats tables for downstream use. |

## Version and database compatibility

The configuration example specifies the MetaPhlAn index `mpa_vOct22_CHOCOPhlAnSGB_202403`. The installed tool must support the selected index and the required command-line flags. Run the `check` command described in [RUNNING.md](RUNNING.md) to inspect executable availability, key flags, and database paths.

HUMAnN requires compatible taxonomic profiles, ChocoPhlAn nucleotide references, and UniRef protein references. The workflow recognizes unnumbered gene/pathway filenames and the numbered variants `_2_genefamilies.tsv` and `_4_pathabundance.tsv`. Matching filenames alone does not establish database or scientific compatibility. Record tool versions and database versions for each analysis.

The example mapping paths use `map_ko_uniref90.txt.gz` and `map_level4ec_uniclust90.txt.gz`. Choose mapping resources whose members match the identifiers in the gene-family table. Their filenames alone cannot establish compatibility.

## Normalization and input preparation

The example configuration uses HUMAnN `community` normalization with special features included (`special: y`). Gene-family values are converted to copies per million (`cpm`), and pathway values to relative abundance (`relab`). These options are visible in the local configuration and command logs.

The formatter removes special summary rows from analysis inputs, checks for unstratified functional features, and aligns sample columns. It preserves the supplied numerical values without applying another normalization step. The optional gene-family input is used to check sample headers only.

## Archived scripts

`archive/original_scripts/` contains the earlier standalone input-preparation and KO/EC scripts. Their checksums are stored alongside them. These copies retain machine-specific paths and software-discovery behavior. Some archived wrappers can download mapping files, create environments when enabled, or overwrite named output archives; review them before standalone use. The main entry point calls the code in `scripts/`.

## Official documentation

- [SRA Toolkit download and conversion](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
- [KneadData](https://github.com/biobakery/kneaddata)
- [MetaPhlAn](https://github.com/biobakery/MetaPhlAn)
- [HUMAnN and its table utilities](https://github.com/biobakery/humann)
- [HUMAnN normalization options](https://github.com/biobakery/humann/blob/master/humann/tools/renorm_table.py)
- [EasyMetagenome workflow guide](https://github.com/YongxinLiu/EasyMetagenome)
