# Code checks

The release was checked with Python syntax parsing, Bash syntax checks, and ten offline unit tests. The tests generate artificial sample labels and small synthetic inputs in temporary directories. No study measurements, actual sample assignments, or research results are bundled as test fixtures.

The tests cover:

- Correct sample alignment and preservation of supplied abundance values.
- Rejection of mismatched sample sets, nonfinite numbers, and stratified functional rows.
- Planning all eight stages without writing analysis outputs, with both HUMAnN search stages retained and the recovered MetaPhlAn CLI flags.
- Accession-specific SRA download paths.
- Recognition of documented HUMAnN output filenames and rejection of ambiguous duplicates.
- Rejection of changed configuration across stages and prevention of overwriting an existing stage.
- Keeping the analysis work directory outside the repository.
- Concatenation of plain and compressed reads while preserving input files.

The three archived scripts were compared byte-for-byte with their recovered source files. Their SHA-256 hashes are included beside them. The distribution was inspected for actual study accessions, sample tables, sequencing files, result matrices, figures, and local execution logs; only an empty CSV header is included as a data-entry template.

## Re-run the offline checks

From the repository root, with Python 3.10 or newer:

```bash
python3 -m unittest discover -s tests -v
```

## What these checks do not establish

SRA downloads, KneadData, MetaPhlAn, HUMAnN, and their reference databases were not executed end-to-end for this release. Tool interfaces were checked against recovered records and official documentation; command planning tests do not prove runtime compatibility. The exact historical HUMAnN invocation, normalization settings, executable version, and every reference snapshot were not recovered.

No claim is made that the reconstructed workflow reproduces the historical numerical results. A real run requires locally supplied sample assignments and compatible software/databases, followed by inspection of logs and outputs. The original MicrobiomeAnalyst web sessions are outside the automated workflow.
