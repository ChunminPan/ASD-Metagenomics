#!/usr/bin/env python3
"""Publication adaptation of the recovered July 2026 input-preparation scripts.

Convert saved, already profiled abundance tables. This does not run sequencing
QC, taxonomic profiling, HUMAnN profiling, or statistical hypothesis tests.
Only the Python standard library is required. See docs/PROVENANCE.md.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

SRR = re.compile(r"(?<![A-Za-z0-9])SRR\d+(?!\d)")
SPECIAL = {"UNMAPPED", "UNINTEGRATED", "UNGROUPED"}
FILENAMES = {
    "pathway": "MA_pathway_abundance_MetaCyc_unstratified.txt",
    "species": "MA_species_abundance_MetaPhlAn.txt",
    "ko": "MA_KO_abundance_from_HUMAnN.txt",
    "ec": "MA_EC_abundance_from_HUMAnN.txt",
}


def open_text(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8-sig", newline="")


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_metadata(path):
    with open_text(path) as f:
        first = f.readline()
        f.seek(0)
        reader = csv.DictReader(f, delimiter="\t" if "\t" in first else ",")
        fields = reader.fieldnames or []
        run_key = next((k for k in ("Run", "#NAME", "NAME") if k in fields), None)
        group_key = next((k for k in ("group", "Group") if k in fields), None)
        if not run_key or not group_key:
            raise ValueError("Metadata requires Run/#NAME and group columns")
        result = {}
        for row in reader:
            run = (row.get(run_key) or "").strip()
            group = (row.get(group_key) or "").strip()
            if not re.fullmatch(r"SRR\d+", run) or not group:
                raise ValueError(f"Invalid sample/group in metadata: {run!r}")
            if run in result:
                raise ValueError(f"Duplicate metadata sample: {run}")
            result[run] = group
    expected = Counter({"ASD_baseline": 18, "TD_control": 20})
    if Counter(result.values()) != expected:
        raise ValueError(f"This project requires 18 ASD_baseline and 20 TD_control samples; got {dict(Counter(result.values()))}")
    return result


def sample_header(row):
    samples = []
    for value in row[1:]:
        matches = SRR.findall(value)
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one SRR accession in header: {value!r}")
        samples.append(matches[0])
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("Missing or duplicate sample columns")
    return samples


def table_rows(path):
    """Skip leading single-field comments; never infer missing sample order."""
    with open_text(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        for row in reader:
            if not row or not any(row):
                continue
            if len(row) == 1 and row[0].startswith("#"):
                continue
            header = row
            break
        if header is None:
            raise ValueError(f"Empty table: {path}")
        samples = sample_header(header)
        yield samples
        for row in reader:
            if not row or not any(row):
                continue
            if len(row) != len(header):
                raise ValueError(f"Wrong field count in {path}, line {reader.line_num}")
            for value in row[1:]:
                try:
                    x = float(value)
                except ValueError as exc:
                    raise ValueError(f"Non-numeric abundance in {path}, line {reader.line_num}") from exc
                if not math.isfinite(x) or x < 0:
                    raise ValueError(f"Abundance must be finite and nonnegative in {path}, line {reader.line_num}")
            yield row


def clean_id(raw, kind):
    raw = raw.strip()
    if kind != "species" and raw.split("|", 1)[0] in SPECIAL:
        return None, "special"
    if kind != "species" and "|" in raw:
        raise ValueError(f"Expected an unstratified {kind} table; found {raw!r}")
    if kind == "pathway":
        feature = raw.split(":", 1)[0].strip()
    elif kind == "ko":
        feature = raw.split(":", 1)[0].strip()
        if not re.fullmatch(r"K\d{5}", feature):
            raise ValueError(f"Invalid KO ID: {raw!r}")
    elif kind == "ec":
        feature = re.sub(r"^EC:", "", raw).split(":", 1)[0].strip()
        if not re.fullmatch(r"\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)", feature):
            raise ValueError(f"Invalid EC ID: {raw!r}")
    else:
        feature = raw.replace("|", ";")
        ranks = feature.split(";")
        if not any(r.startswith("s__") for r in ranks) or any(r.startswith("t__") for r in ranks):
            return None, "other_rank"
    if not feature:
        raise ValueError("Empty feature ID")
    return feature, None


def convert(path, kind, samples, metadata, output):
    rows = table_rows(path)
    observed = next(rows)
    if set(observed) != set(metadata):
        raise ValueError(f"Sample set mismatch in {kind}")
    indices = [observed.index(s) + 1 for s in samples]
    seen, skipped, mapping = set(), Counter(), []
    with output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["#NAME", *samples])
        for row in rows:
            feature, reason = clean_id(row[0], kind)
            if feature is None:
                skipped[reason] += 1
                continue
            if feature in seen:
                raise ValueError(f"Duplicate {kind} ID after cleanup: {feature}")
            seen.add(feature)
            writer.writerow([feature, *(row[i] for i in indices)])
            if kind == "pathway":
                mapping.append((feature, row[0].strip()))
    if not seen:
        raise ValueError(f"No {kind} features retained")
    return {"features": len(seen), "skipped": dict(skipped), "input_sha256": digest(path)}, mapping


def run(args):
    metadata = read_metadata(args.metadata)
    samples = next(table_rows(args.pathway))
    if set(samples) != set(metadata):
        raise ValueError("Pathway sample set does not match metadata")
    if args.gene:
        # As in the recovered script, gene input is used only for sample alignment.
        if set(next(table_rows(args.gene))) != set(samples):
            raise ValueError("Gene sample set does not match pathway samples")
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Output directory must be empty to avoid mixing runs")
    outdir.mkdir(parents=True, exist_ok=True)
    report = {"status": "IN_PROGRESS", "samples": len(samples),
              "groups": dict(Counter(metadata.values())),
              "metadata_sha256": digest(args.metadata),
              "gene_header_checked": bool(args.gene), "tables": {}}
    try:
        with (outdir / "MA_metadata_ASD_TD.txt").open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["#NAME", "group"])
            writer.writerows((s, metadata[s]) for s in samples)
        for kind, name in FILENAMES.items():
            path = getattr(args, kind)
            if not path:
                continue
            stats, mapping = convert(path, kind, samples, metadata, outdir / name)
            report["tables"][kind] = stats
            if mapping:
                with (outdir / "MetaCyc_feature_name_map.tsv").open("w", encoding="utf-8", newline="") as f:
                    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                    writer.writerow(["MetaCyc_ID", "Original_HUMAnN_feature"])
                    writer.writerows(mapping)
        report["status"] = "PASS"
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = str(exc)
        raise
    finally:
        (outdir / "validation_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    files = sorted(p for p in outdir.iterdir() if p.is_file())
    (outdir / "SHA256SUMS.txt").write_text("".join(f"{digest(p)}  {p.name}\n" for p in files), encoding="utf-8")
    print(f"PASS: {len(samples)} samples; " + ", ".join(f"{k}={v['features']}" for k,v in report["tables"].items()))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", required=True)
    p.add_argument("--pathway", required=True)
    p.add_argument("--species")
    p.add_argument("--ko")
    p.add_argument("--ec")
    p.add_argument("--gene", help="Optional gene matrix; checks its sample header only")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    try:
        run(args)
    except (ValueError, OSError) as exc:
        p.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
