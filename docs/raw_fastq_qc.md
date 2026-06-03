# Optional raw FASTQ QC, contamination screening, and MultiQC

This scRNA-seq workflow normally starts from processed count matrices: `.h5ad`, 10X `.h5`, or 10X MEX directories. If raw FASTQs are also available, the workflow can run upstream QC before the count-matrix analysis starts.

## Enable raw FASTQ QC

Add `raw_fastq_r1` and `raw_fastq_r2` columns to `config/samples.tsv`. Leave the cells blank for samples without raw FASTQs.

```tsv
sample	path	batch	condition	description	raw_fastq_r1	raw_fastq_r2
pbmc_ctrl	data/pbmc_ctrl_raw.h5ad	batch1	control	PBMC donor 1	data/fastq/pbmc_ctrl_R1.fastq.gz	data/fastq/pbmc_ctrl_R2.fastq.gz
```

Configure the feature in `config/config.yaml`:

```yaml
raw_fastq_qc:
  enabled: true
  fastq_screen_conf: "config/fastq_screen.conf"
  fastq_screen_subset: 100000
```

## Configure FastQ Screen

Copy the template and replace every `/path/to/...` entry with a real Bowtie2 index prefix:

```bash
cp config/fastq_screen.conf.example config/fastq_screen.conf
```

Keep `config/fastq_screen.conf` local; it is ignored by Git because it usually contains site-specific reference paths.

## Outputs

For each sample with raw FASTQs configured, the workflow writes:

- `results/<sample>/raw_fastq_qc/fastqc/` - raw FASTQ FastQC HTML and ZIP outputs
- `results/<sample>/raw_fastq_qc/fastq_screen/` - FastQ Screen text and HTML reports

The workflow also writes:

- `results/multiqc/multiqc_report.html` - consolidated MultiQC report

If no raw FASTQs are configured, MultiQC writes a small placeholder HTML file so the default workflow target still completes.
