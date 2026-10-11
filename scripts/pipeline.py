#!/usr/bin/env python3
"""Code-only publication workflow. Historical recovery is documented in PROVENANCE.md.

Download, QC/host removal, taxonomy, function, merge, normalization, regrouping,
and table export. This orchestrator is new release code. In particular HUMAnN
invocation and normalization are reconstructions, not recovered production files.
No study data or result values are embedded. Commands execute without a shell.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from prepare_inputs import read_metadata

ROOT = Path(__file__).resolve().parents[1]
STAGES = ("download", "dehost", "taxonomy", "function", "merge", "normalize", "regroup", "export")
ENV_FOR = {
    "prefetch": "sra", "fasterq-dump": "sra", "pigz": "sra",
    "kneaddata": "qc", "metaphlan": "taxonomy", "merge_metaphlan_tables.py": "taxonomy",
    "humann": "function", "humann_join_tables": "function",
    "humann_renorm_table": "function", "humann_split_stratified_table": "function",
    "humann_regroup_table": "regroup",
}


def local_path(value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty file: {path}")


def concatenate_reads(paths, destination):
    """Concatenate paired reads as independent records for HUMAnN, per its manual."""
    with destination.open("xb") as output:
        for path in paths:
            opener = gzip.open if str(path).endswith(".gz") else open
            with opener(path, "rb") as source:
                shutil.copyfileobj(source, output)


class Pipeline:
    def __init__(self, config_path, dry_run=False):
        self.config_path = config_path.resolve()
        self.config = json.loads(config_path.read_text(encoding="utf-8"))
        self.dry_run = dry_run
        self.work = local_path(self.config["work_dir"])
        if self.work == ROOT or ROOT.is_relative_to(self.work) or self.work.is_relative_to(ROOT):
            raise ValueError("work_dir must be outside the repository and must not be its ancestor")
        self.metadata = local_path(self.config["metadata"])
        self.groups = read_metadata(self.metadata)
        self.samples = list(self.groups)
        self.threads = self.config["threads"]
        for name in ("download", "qc", "taxonomy", "function"):
            if type(self.threads.get(name)) is not int or self.threads[name] < 1:
                raise ValueError(f"threads.{name} must be a positive integer")
        self.db = self.config["databases"]
        self.log = None
        self.fingerprint = hashlib.sha256(
            config_path.read_bytes() + self.metadata.read_bytes() + Path(__file__).read_bytes()
            + (ROOT / "scripts/prepare_inputs.py").read_bytes()
        ).hexdigest()

    def directory(self, stage):
        return self.work / f"{STAGES.index(stage)+1:02d}_{stage}"

    def prefix(self, tool):
        environment = self.config.get("conda_environments", {}).get(ENV_FOR.get(tool, ""), "")
        if environment:
            return ["conda", "run", "--no-capture-output", "-n", environment, tool]
        return [tool]

    def command(self, tool, *arguments):
        return self.prefix(tool) + [str(a) for a in arguments]

    def run(self, cmd):
        print(shlex.join([str(x) for x in cmd]), flush=True)
        if not self.dry_run:
            self.log.write(shlex.join([str(x) for x in cmd]) + "\n")
            self.log.flush()
            subprocess.run(cmd, check=True, stdout=self.log, stderr=subprocess.STDOUT)

    def database(self, key, is_file=False):
        value = self.db.get(key, "")
        if not value:
            raise ValueError(f"Set databases.{key} in the local config")
        path = local_path(value)
        if not self.dry_run:
            if is_file:
                require_file(path)
            elif not path.is_dir():
                raise ValueError(f"Database directory does not exist: {path}")
        return path

    def completed(self, stage):
        if self.dry_run:
            return
        marker = self.directory(stage) / "COMPLETE.json"
        require_file(marker)
        record = json.loads(marker.read_text())
        if record.get("fingerprint") != self.fingerprint:
            raise ValueError(f"Config/manifest/code changed after {stage}; choose a fresh work_dir")

    def clean_pair(self, sample):
        folder = self.directory("dehost") / sample
        paths = []
        for mate in (1, 2):
            base = folder / f"{sample}_kneaddata_paired_{mate}.fastq"
            if self.dry_run:
                paths.append(base)
                continue
            candidates = [p for p in (base, Path(str(base) + ".gz")) if p.is_file()]
            if len(candidates) != 1:
                raise ValueError(f"Expected one cleaned mate {mate} for {sample}, got {candidates}")
            require_file(candidates[0])
            paths.append(candidates[0])
        return paths

    def functional_table(self, sample, kind):
        """Resolve documented HUMAnN output names without guessing between files."""
        numbers = {"genefamilies": "2", "pathabundance": "4"}
        folder = self.directory("function") / sample / "profile"
        plain = folder / f"{sample}_{kind}.tsv"
        if self.dry_run:
            return plain
        numbered = folder / f"{sample}_{numbers[kind]}_{kind}.tsv"
        candidates = [p for p in (plain, numbered) if p.is_file()]
        if len(candidates) != 1:
            raise ValueError(f"Expected one HUMAnN {kind} table for {sample}, got {candidates}")
        require_file(candidates[0])
        return candidates[0]

    def check(self):
        """Read-only tool/CLI checks; no installation, database download or analysis."""
        print(f"MANIFEST_OK: {len(self.samples)} samples, 18 ASD_baseline and 20 TD_control")
        self.database("host")
        self.database("metaphlan")
        self.database("chocophlan")
        self.database("uniref")
        self.database("ko_map", is_file=True)
        self.database("ec_map", is_file=True)
        checks = {
            "prefetch": ["--max-size"], "fasterq-dump": ["--split-files"], "pigz": [],
            "kneaddata": ["--input1", "--input2", "--output-prefix"],
            "metaphlan": ["--mapout", "--db_dir"],
            "merge_metaphlan_tables.py": [],
            "humann": ["--taxonomic-profile", "--nucleotide-database", "--protein-database"],
            "humann_join_tables": ["--file_name"],
            "humann_renorm_table": ["--units", "--mode", "--special"],
            "humann_split_stratified_table": ["--input"],
            "humann_regroup_table": ["--custom"],
        }
        for tool, flags in checks.items():
            cmd = self.command(tool, "--help")
            if self.dry_run:
                print(shlex.join(cmd))
                continue
            p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90)
            if p.returncode:
                raise ValueError(f"Cannot run {tool} --help:\n{p.stdout[-2000:]}")
            for flag in flags:
                if flag not in p.stdout:
                    raise ValueError(f"Installed {tool} does not advertise {flag}; see docs/RUNNING.md")
            print(f"TOOL_OK: {tool}")
            if tool in ("prefetch", "fasterq-dump", "pigz", "kneaddata", "metaphlan", "humann"):
                version = subprocess.run(self.command(tool, "--version"), text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90, check=True)
                print(version.stdout.strip())
        print("CHECK_COMPLETE: CLI availability does not verify database compatibility or scientific results.")

    def execute(self, stage):
        folder = self.directory(stage)
        if not self.dry_run:
            if folder.exists():
                raise ValueError(f"Stage directory already exists: {folder}. No files will be overwritten.")
            folder.mkdir(parents=True)
            self.log = (folder / "commands.log").open("x", encoding="utf-8")
            (folder / "config_snapshot.json").write_text(json.dumps(self.config, indent=2) + "\n")
            shutil.copyfile(self.metadata, folder / "sample_metadata.csv")
        print(f"\nSTAGE {stage}", flush=True)
        try:
            getattr(self, stage)(folder)
            if not self.dry_run:
                record = {"stage":stage, "fingerprint":self.fingerprint,
                          "completed_utc":datetime.now(timezone.utc).isoformat(), "samples":len(self.samples)}
                (folder / "COMPLETE.json").write_text(json.dumps(record, indent=2) + "\n")
        finally:
            if self.log:
                self.log.close()
                self.log = None

    def download(self, out):
        # Adapted from recovered pilot prefetch/fasterq-dump commands.
        sra = out / "sra"
        for sample in self.samples:
            fastq, temp = out / "fastq" / sample, out / "tmp" / sample
            if not self.dry_run:
                fastq.mkdir(parents=True)
                temp.mkdir(parents=True)
            self.run(self.command("prefetch", sample, "-O", sra / sample, "--max-size", self.config["sra_max_size"]))
            self.run(self.command("fasterq-dump", sra / sample / f"{sample}.sra",
                "--split-files", "--outdir", fastq, "--temp", temp, "-e", self.threads["download"]))
            for mate in (1, 2):
                fq = fastq / f"{sample}_{mate}.fastq"
                if not self.dry_run:
                    require_file(fq)
                self.run(self.command("pigz", "-p", self.threads["download"], fq))
                if not self.dry_run:
                    require_file(Path(str(fq) + ".gz"))

    def dehost(self, out):
        # Recovered cloud call, with a new explicit output prefix for stable paths.
        self.completed("download")
        host = self.database("host")
        for sample in self.samples:
            source = self.directory("download") / "fastq" / sample
            r1, r2 = [source / f"{sample}_{m}.fastq.gz" for m in (1,2)]
            if not self.dry_run:
                require_file(r1)
                require_file(r2)
            self.run(self.command("kneaddata", "-i1", r1, "-i2", r2, "-db", host,
                "-o", out / sample, "-t", self.threads["qc"],
                "--output-prefix", f"{sample}_kneaddata",
                "--trimmomatic-options", "SLIDINGWINDOW:4:20 MINLEN:50"))
            if not self.dry_run:
                self.clean_pair(sample)

    def taxonomy(self, out):
        # Recovered corrected --mapout/--db_dir invocation, used in July 2026.
        self.completed("dehost")
        database = self.database("metaphlan")
        index = self.config["metaphlan_index"]
        if not index:
            raise ValueError("metaphlan_index must identify a fixed installed database")
        for sample in self.samples:
            reads = self.clean_pair(sample)
            profile = out / f"{sample}_metaphlan_profile.txt"
            self.run(self.command("metaphlan", ",".join(map(str, reads)), "--input_type", "fastq",
                "--nproc", self.threads["taxonomy"], "--db_dir", database, "-x", index,
                "--mapout", out / f"{sample}.bowtie2out.bz2", "-o", profile))
            if not self.dry_run:
                require_file(profile)

    def function(self, out):
        # RECONSTRUCTED for release: exact production HUMAnN command not recovered.
        # Both nucleotide and translated search are retained. No bypass flags.
        self.completed("taxonomy")
        self.completed("dehost")
        nucleotide, protein = self.database("chocophlan"), self.database("uniref")
        for sample in self.samples:
            folder = out / sample
            combined = folder / f"{sample}.fastq"
            reads = self.clean_pair(sample)
            profile = self.directory("taxonomy") / f"{sample}_metaphlan_profile.txt"
            if self.dry_run:
                print(f"[Python concatenate] {reads[0]} + {reads[1]} -> {combined}")
            else:
                require_file(profile)
                folder.mkdir()
                concatenate_reads(reads, combined)
            self.run(self.command("humann", "--input", combined, "--output", folder / "profile",
                "--threads", self.threads["function"], "--output-basename", sample,
                "--taxonomic-profile", profile, "--nucleotide-database", nucleotide,
                "--protein-database", protein))
            if not self.dry_run:
                for kind in ("genefamilies", "pathabundance"):
                    self.functional_table(sample, kind)
                # Only this newly created concatenated copy is removed after success.
                combined.unlink()

    def merge(self, out):
        # Reconstructed full calls around the recovered merge utility names.
        self.completed("function")
        self.completed("taxonomy")
        profiles = [self.directory("taxonomy") / f"{s}_metaphlan_profile.txt" for s in self.samples]
        merged_species = out / "metaphlan_merged_all_taxa.tsv"
        self.run(self.command("merge_metaphlan_tables.py", *profiles, "-o", merged_species))
        if not self.dry_run:
            require_file(merged_species)
        for kind in ("genefamilies", "pathabundance"):
            inputs = out / f"input_{kind}"
            if not self.dry_run:
                inputs.mkdir()
                for sample in self.samples:
                    source = self.functional_table(sample, kind)
                    (inputs / source.name).symlink_to(source)
            target = out / f"{kind}.tsv"
            self.run(self.command("humann_join_tables", "--input", inputs, "--file_name", kind, "--output", target))
            if not self.dry_run:
                require_file(target)

    def normalize(self, out):
        # RECONSTRUCTED explicit settings; see PROVENANCE.md before reuse.
        self.completed("merge")
        options = self.config["normalization"]
        if options["mode"] not in ("community", "levelwise") or options["special"] not in ("y", "n"):
            raise ValueError("Invalid normalization mode or special-feature setting")
        for kind, units in (("genefamilies", "cpm"), ("pathabundance", "relab")):
            source = self.directory("merge") / f"{kind}.tsv"
            target = out / f"{kind}_{units}.tsv"
            self.run(self.command("humann_renorm_table", "--input", source, "--output", target,
                "--units", units, "--mode", options["mode"], "--special", options["special"]))
            self.run(self.command("humann_split_stratified_table", "--input", target, "--output", out))
            if not self.dry_run:
                require_file(out / f"{kind}_{units}_unstratified.tsv")

    def regroup(self, out):
        # Adaptation of the recovered KO/EC script; requires explicit mapping files.
        self.completed("normalize")
        source = self.directory("normalize") / "genefamilies_cpm_unstratified.tsv"
        for kind in ("ko", "ec"):
            mapping = self.database(f"{kind}_map", is_file=True)
            if not self.dry_run:
                h = hashlib.sha256()
                with mapping.open("rb") as f:
                    for block in iter(lambda:f.read(1024*1024), b""):
                        h.update(block)
                (out / f"{kind}_mapping.sha256").write_text(f"{h.hexdigest()}  {mapping.name}\n")
            target = out / f"{kind}_raw.tsv"
            self.run(self.command("humann_regroup_table", "--input", source, "--custom", mapping, "--output", target))
            if not self.dry_run:
                require_file(target)

    def export(self, out):
        self.completed("regroup")
        self.completed("normalize")
        self.completed("merge")
        normalized = self.directory("normalize")
        self.run([sys.executable, str(ROOT / "scripts/prepare_inputs.py"),
            "--metadata", str(self.metadata),
            "--gene", str(normalized / "genefamilies_cpm_unstratified.tsv"),
            "--pathway", str(normalized / "pathabundance_relab_unstratified.tsv"),
            "--species", str(self.directory("merge") / "metaphlan_merged_all_taxa.tsv"),
            "--ko", str(self.directory("regroup") / "ko_raw.tsv"),
            "--ec", str(self.directory("regroup") / "ec_raw.tsv"),
            "--outdir", str(out / "ready")])
        if not self.dry_run:
            require_file(out / "ready/validation_report.json")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("check", "all", *STAGES))
    p.add_argument("--config", type=Path, default=ROOT / "config/local.json")
    p.add_argument("--dry-run", action="store_true", help="Print planned commands without analysis/downloads or output writes")
    args = p.parse_args()
    try:
        pipeline = Pipeline(args.config, args.dry_run)
        if args.stage == "check":
            pipeline.check()
        else:
            for stage in STAGES if args.stage == "all" else (args.stage,):
                pipeline.execute(stage)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        p.exit(1, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
