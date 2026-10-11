"""Offline checks using synthetic fixtures generated in temporary directories.

No real study accessions, assignments, measurements, or output tables are stored
here. Third-party tools are not executed by these tests.
"""
import argparse
import contextlib
import csv
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_inputs as prep
from pipeline import Pipeline, STAGES, concatenate_reads


def write_table(path, samples, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["#NAME", *samples])
        writer.writerows(rows)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        # Artificial accession-shaped labels. These are not selected study runs.
        self.samples = [f"SRR{900000000+i}" for i in range(38)]
        self.meta = self.directory / "metadata.csv"
        with self.meta.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Run", "group"])
            writer.writerows((s, "ASD_baseline" if i<18 else "TD_control") for i,s in enumerate(self.samples))
        self.pathway = self.directory / "pathway.tsv"
        write_table(self.pathway, self.samples,
            [["UNMAPPED", *(["0.1"]*38)], ["PWY-TEST: synthetic label", *(["0.3"]*38)]])
        self.ko = self.directory / "ko.tsv"
        write_table(self.ko, self.samples[::-1], [["K00001", *(str(i) for i in reversed(range(38)))]])
        self.ec = self.directory / "ec.tsv"
        write_table(self.ec, self.samples, [["1.1.1.1", *(["2.5"]*38)]])
        self.species = self.directory / "species.tsv"
        write_table(self.species, self.samples,
            [["k__Bacteria|g__Synthetic|s__Synthetic_example", *(["5"]*38)]])

    def arguments(self):
        return argparse.Namespace(metadata=self.meta,pathway=self.pathway,ko=self.ko,ec=self.ec,
            species=self.species,gene=None,outdir=self.directory/"prepared")

    def configuration(self):
        config = json.loads((ROOT/"config/config.example.json").read_text())
        config["metadata"] = str(self.meta)
        config["work_dir"] = str(self.directory/"analysis-work")
        path = self.directory/"local.json"
        path.write_text(json.dumps(config))
        return path

    def test_table_alignment_and_no_renormalization(self):
        args = self.arguments()
        report = prep.run(args)
        with (args.outdir/prep.FILENAMES["ko"]).open() as f:
            rows=list(csv.reader(f,delimiter="\t"))
        self.assertEqual(rows[0][1:],self.samples)
        self.assertEqual(rows[1][1:],[str(i) for i in range(38)])
        with (args.outdir/prep.FILENAMES["pathway"]).open() as f:
            rows=list(csv.reader(f,delimiter="\t"))
        self.assertEqual(rows[1],["PWY-TEST", *(["0.3"]*38)])
        self.assertEqual(report["tables"]["pathway"]["skipped"],{"special":1})

    def test_sample_mismatch_fails(self):
        changed=self.samples.copy()
        changed[0]="SRR999999999"
        write_table(self.ko,changed,[["K00001",*(["1"]*38)]])
        with self.assertRaisesRegex(ValueError,"Sample set mismatch"):
            prep.run(self.arguments())

    def test_nonfinite_values_fail(self):
        write_table(self.ko,self.samples,[["K00001", "nan", *(["1"]*37)]])
        with self.assertRaisesRegex(ValueError,"finite and nonnegative"):
            prep.run(self.arguments())

    def test_stratification_is_not_double_counted(self):
        write_table(self.pathway,self.samples,[["PWY-TEST|g__Synthetic.s__example",*(["1"]*38)]])
        with self.assertRaisesRegex(ValueError,"unstratified"):
            prep.run(self.arguments())

    def test_all_stage_dry_run_has_no_writes_or_bypass(self):
        pipeline=Pipeline(self.configuration(),dry_run=True)
        commands=[]
        with contextlib.redirect_stdout(io.StringIO()), patch.object(pipeline,"run",side_effect=commands.append):
            for stage in STAGES:
                pipeline.execute(stage)
        self.assertFalse(pipeline.work.exists())
        downloads=[c for c in commands if c[0]=="prefetch"]
        self.assertEqual(len(downloads),38)
        self.assertTrue(all(Path(c[c.index("-O")+1]).name==c[1] for c in downloads))
        humann=[c for c in commands if c[0]=="humann"]
        self.assertEqual(len(humann),38)
        self.assertTrue(all("--taxonomic-profile" in c for c in humann))
        self.assertFalse(any(any(str(arg).startswith("--bypass") for arg in c) for c in commands))
        taxonomy=[c for c in commands if c[0]=="metaphlan"]
        self.assertTrue(all("--mapout" in c and "--db_dir" in c for c in taxonomy))
        normalized=[c for c in commands if c[0]=="humann_renorm_table"]
        self.assertEqual({c[c.index("--units")+1] for c in normalized},{"cpm","relab"})

    def test_changed_config_blocks_downstream_stage(self):
        config_path=self.configuration()
        initial=Pipeline(config_path)
        stage=initial.directory("download")
        stage.mkdir(parents=True)
        (stage/"COMPLETE.json").write_text(json.dumps({"fingerprint":initial.fingerprint}))
        config=json.loads(config_path.read_text())
        config["threads"]["download"]+=1
        config_path.write_text(json.dumps(config))
        changed=Pipeline(config_path)
        with self.assertRaisesRegex(ValueError,"changed after download"):
            changed.completed("download")

    def test_humann_output_names_and_ambiguity(self):
        pipeline=Pipeline(self.configuration())
        sample=self.samples[0]
        folder=pipeline.directory("function")/sample/"profile"
        folder.mkdir(parents=True)
        for kind,number in (("genefamilies","2"),("pathabundance","4")):
            plain=folder/f"{sample}_{kind}.tsv"
            numbered=folder/f"{sample}_{number}_{kind}.tsv"
            plain.write_text("synthetic fixture\n")
            self.assertEqual(pipeline.functional_table(sample,kind),plain)
            plain.rename(numbered)
            self.assertEqual(pipeline.functional_table(sample,kind),numbered)
            plain.write_text("synthetic fixture\n")
            with self.assertRaisesRegex(ValueError,"Expected one HUMAnN"):
                pipeline.functional_table(sample,kind)

    def test_existing_stage_is_not_overwritten(self):
        pipeline=Pipeline(self.configuration())
        directory=pipeline.directory("download")
        directory.mkdir(parents=True)
        marker=directory/"keep.txt"
        marker.write_text("keep")
        with self.assertRaisesRegex(ValueError,"already exists"):
            pipeline.execute("download")
        self.assertEqual(marker.read_text(),"keep")

    def test_work_directory_cannot_contain_repository_data(self):
        config_path=self.configuration()
        config=json.loads(config_path.read_text())
        config["work_dir"]=str(ROOT/"results")
        config_path.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError,"outside the repository"):
            Pipeline(config_path)

    def test_concatenate_plain_and_gzip_reads(self):
        first=self.directory/"r1.fastq"
        second=self.directory/"r2.fastq.gz"
        left,right=b"@synthetic/1\nAC\n+\nII\n", b"@synthetic/2\nGT\n+\nII\n"
        first.write_bytes(left)
        with gzip.open(second,"wb") as f:
            f.write(right)
        output=self.directory/"both.fastq"
        concatenate_reads([first,second],output)
        self.assertEqual(output.read_bytes(),left+right)
        self.assertEqual(first.read_bytes(),left)
        self.assertTrue(second.exists())


if __name__ == "__main__":
    unittest.main()
