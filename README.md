# Short-Read-Genome-Assembler

Batch pipeline for Illumina paired‑end reads:
FastQC → fastp → SPAdes → Racon → seqtk → CheckM2.
Resumable, per‑sample, and skips completed steps.

## Requirements

- Bash + Conda
- Tools:
  - ShortReadAssembler.yml: assembly-stats, fastqc, fastp, spades.py, minimap2, racon, seqtk
  - checkm2.yml: checkm2 (with a valid CheckM2 database)
- Sufficient CPU/RAM (defaults are very high; adjust before running)

## Input

- Place paired‑end reads in a directory named `00_reads/`
- Filenames must end with `_R1.fastq.gz` and `_R2.fastq.gz` for each sample
  - Example: `00_reads/SAMPLE_A_R1.fastq.gz` and `00_reads/SAMPLE_A_R2.fastq.gz`

The sample ID is derived from the basename before `_R1.fastq.gz`.

## Configure

Edit these variables at the top of the script as needed:

```bash
THREADS=120
MEMORY_GB=800
LENGTH_FILTER=500
RACON_ROUNDS=2
checkm2_database=/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

## Usage

```bash
# make executable 
chmod +x ShortReadAssembler.sh

# run
./ShortReadAssembler.sh | tee run.log
```

The script activates the `ShortReadAssembler` environment first and later switches to `checkm2` for CheckM2.
It uses `eval "$(conda shell.bash hook)"`, so ensure conda is installed and accessible.

## Workflow

For each sample in `00_reads/`:

1. FastQC on raw reads → `01_fastqc/<ID>/`
2. fastp trim + reports → `02_trimmed_reads/<ID>/`
3. SPAdes assembly (`--isolate`) → `03_assemblies/<ID>/<ID>.fasta`
4. Racon polishing (minimap2 alignments, N rounds) → `04_racon_polish/<ID>/`
5. seqtk length filter (min length=LENGTH_FILTER) → copy to final
6. CheckM2 on all final assemblies → `05_final_assemblies/checkm2_results/`

Skips steps if outputs already exist; failed SPAdes samples are recorded in `failed.txt` and the pipeline continues (`set +e` with `pipefail`).

## Outputs

- `01_fastqc/<ID>/` — FastQC HTML reports for raw reads
- `02_trimmed_reads/<ID>/`
  - `<ID>_R1.trimmed.fastq.gz`, `<ID>_R2.trimmed.fastq.gz`
  - `<ID>.fastp.html`, `<ID>.fastp.json`
- `03_assemblies/<ID>/`
  - SPAdes output directory and `<ID>.fasta` (copied contigs)
- `04_racon_polish/<ID>/`
  - `<ID>_polished_roundN.fasta` per round, temporary PAFs removed
  - `<ID>_filtered.fasta` after seqtk length filter
- `05_final_assemblies/`
  - `<ID>.fasta` — final filtered assembly per sample
  - `checkm2_results/` — CheckM2 quality metrics for all assemblies
- `failed.txt` — sample IDs that failed SPAdes (or missing polished assembly)

## Notes

- The script is resumable: re‑running will skip completed outputs.
- Combined reads for polishing are concatenated gz files (valid multi‑member gzip).
- Minimises disk usage by cleaning temporary alignment files.
