# From 16S to shotgun metagenomics

### Exploring the gut microbiome in autism and gastrointestinal conditions

**Chunmin Pan | Independent research | Python, Bash, and Linux**

## Why I began this research

While working with autistic musicians, I heard parents describe their children's gastrointestinal discomfort and uncertainty about dietary choices. These conversations led me to study the gut microbiome, the community of microorganisms living in the digestive tract. I wanted to investigate whether public biological data could help me examine the questions raised by those everyday experiences.

My earlier project analyzed 325 publicly available 16S microbiome samples across four groups: healthy controls, individuals with autism, individuals with intestinal disorders, and individuals with both conditions. I compared the composition of their gut microbial communities and examined its association with whole-grain intake. The samples came from an existing public study; I carried out the secondary analysis.

## Why continue beyond 16S?

16S sequencing uses a marker gene to help identify bacteria. My first analysis focused on bacterial genera and how their abundances varied across groups. It helped me study which microorganisms were present, but it did not directly measure the metabolic genes carried by those communities. Discussions of functions such as fiber fermentation therefore depended on what previous research had reported about particular bacteria.

That limitation shaped my next question: **what metabolic capabilities are represented in the microbial DNA, and how do those capabilities vary between groups?**

I extended the research to shotgun metagenomics, which samples DNA across the microbial community. This approach can provide more detailed information about microbial species and the genes involved in metabolic processes. I examined an independent public dataset of 38 gut metagenomes, comprising 18 baseline samples from autistic children with gastrointestinal symptoms and 20 typically developing controls. I excluded post-treatment and donor samples from the selected comparison.

My areas of interest included short-chain fatty acid production, which is linked to microbial fermentation; the breakdown and use of dietary fiber; and sulfur metabolism. These provided biological questions for interpreting functional profiles. The extension uses a different dataset from the 16S study and does not directly test the earlier whole-grain association.

## What I did

I taught myself Linux command-line tools and worked through the practical steps needed to analyze sequencing data. Using EasyMetagenome as a workflow guide, I configured bioinformatics software and reference databases, processed public sequencing files, and worked through errors in file formats, tool arguments, and software compatibility. I used MetaPhlAn to examine microbial composition and HUMAnN to examine gene-family and pathway profiles, then prepared the tables for downstream analysis in MicrobiomeAnalyst.

This repository presents the computational workflow behind that research. It contains code and documentation only. Sequencing files, actual sample assignments, abundance tables, figures, and research results are not included.

## What the pipeline does

| Stage | In plain language | Main tool |
| --- | --- | --- |
| 1. Download | Retrieve the selected public sequencing files and convert them into readable sequence files | SRA Toolkit |
| 2. Quality control and host removal | Filter low-quality reads and remove reads matching the human reference | KneadData |
| 3. Species profiling | Estimate which microbes are represented in each sample | MetaPhlAn |
| 4. Functional profiling | Estimate the abundance of microbial gene families and metabolic pathways | HUMAnN |
| 5. Merge | Bring individual sample profiles into shared tables | MetaPhlAn and HUMAnN utilities |
| 6. Normalize | Put gene and pathway abundances on comparable scales across samples | HUMAnN utilities |
| 7. Group functional annotations | Organize gene families into KEGG Orthology and Enzyme Commission categories | HUMAnN regrouping |
| 8. Prepare analysis inputs | Align sample identifiers and format tables for downstream analysis | Project Python code |

KEGG Orthology groups genes by shared functions. Enzyme Commission categories describe enzyme-catalyzed reactions. Metabolic pathways describe connected biochemical processes. These annotations offer ways to examine microbial functional potential; they do not directly measure metabolite concentrations or establish effects on human health.

## About the code

The repository preserves three recovered project scripts in `archive/original_scripts/`. The portable workflow in `scripts/` organizes recovered commands and documented publication edits into a common interface. Some production command blocks, including the final HUMAnN invocation and normalization commands, could not be recovered verbatim. Their replacements are explicitly identified as reconstructions in [PROVENANCE.md](docs/PROVENANCE.md). This is not a claim that every file here was used unchanged in the original analysis.

For technical use, start with [RUNNING.md](docs/RUNNING.md). A local sample sheet and installed reference databases are required; the included sample-sheet template is empty. The pipeline can print its planned commands before executing them:

```bash
bash run_pipeline.sh all --config config/local.json --dry-run
```

See [VALIDATION.md](docs/VALIDATION.md) for the checks performed and their limits. The original MicrobiomeAnalyst sessions are not recreated by this code.

## Data sources and acknowledgments

The earlier 16S work used GMrepo abundance data and metadata associated with [PRJEB11419](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB11419). The shotgun extension used [PRJNA782533](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA782533). I did not collect the original biological samples or generate the sequencing data.

The workflow relies on [EasyMetagenome](https://github.com/YongxinLiu/EasyMetagenome), [SRA Toolkit](https://github.com/ncbi/sra-tools), [KneadData](https://github.com/biobakery/kneaddata), [MetaPhlAn](https://github.com/biobakery/MetaPhlAn), [HUMAnN](https://github.com/biobakery/humann), and [MicrobiomeAnalyst](https://www.microbiomeanalyst.ca/). Project scripts were developed and refined with AI coding assistance. The underlying bioinformatics methods and tools belong to their respective developers.

Project code and documentation are available under the [MIT License](LICENSE). Third-party tools and databases retain their own licenses.
