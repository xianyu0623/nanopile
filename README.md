# Nanopile

**Nanopile** is a specialized tool for analyzing nanopore sequencing data. It functions similarly to `samtools mpileup` but is specifically designed to parse and utilize **move table values** from nanopore BAM files.

> [!NOTE]
> This tool is currently under active development. Bug reports and feature requests are welcome!

## Installation

### Prerequisites

You need to have **Rust** and **Cargo** installed to build this tool. You can install them from the [official Rust website](https://www.rust-lang.org/tools/install).

### Build from Source

```bash
# Clone the repository
git clone https://github.com/yourusername/nanopile.git
cd nanopile

# Build the project
cargo build --release

# The binary will be available at ./target/release/nanopile
./target/release/nanopile --help
```

## Usage

```bash
nanopile [OPTIONS] --bam_fp <BAM_FP> --region <REGION>...
# OR
nanopile [OPTIONS] --bam_fp <BAM_FP> --bed_fp <BED_FP>
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --bam_fp` | Input BAM file (must be sorted and indexed) | **Required** |
| `-r, --region` | Target region (1-based, inclusive, e.g., `chr1:100-200`). Can be specified multiple times. | Required if no BED |
| `-l, --bed_fp` | Input BED file (0-based, half-open). Mutually exclusive with `--region`. | Required if no Region |
| `-f, --ref_fp` | Reference FASTA file | Optional |
| `--buffer_size` | Buffer size for reading BAM file | `10000` |
| `--margin` | Margin for reading BAM file | `500` |
| `-q, --min_mapq` | Minimum mapping quality | `0` |
| `-Q, --min_baseq` | Minimum base quality | `13` |
| `--flag_filter` | SAM flag filter | `0` |

### Output Flags

Use these flags to include additional information in the output:

| Flag | Description |
|------|-------------|
| `--output_mv` | Output Move Table values |
| `--output_bq` | Output Base Quality scores |
| `--output_mapq` | Output Mapping Quality scores |
| `--output_read_name` | Output Read Names |

## Python API

The core pileup engine is also available from Python behind an optional feature flag.

1. Install [maturin](https://www.maturin.rs/) once (e.g. `pip install maturin`).
2. From the repository root run:

   ```bash
   maturin develop --release --features python
   ```

   This builds and installs the `nanopile` extension module into your current Python environment.

3. Call the binding from Python:

   ```python
   import nanopile

   pileup = nanopile.run_nanopile(
       bam_fp="reads.bam",
       ref_fp="reference.fa",
       regions=["chr1:1000-1100"],
       output_bq=True,
       output_read_name=True,
   )

   for pos in pileup:
       print(pos.chrom, pos.pos + 1, pos.depth, pos.bases)
   ```

`run_nanopile` mirrors the CLI flags: you must provide either `bed_fp` or `regions`, and you can toggle the optional outputs with the same boolean parameters. The function returns a Python `list` of `PyPileupPos` objects, so every position can be iterated over and its attributes accessed directly (`bases`, `read_names`, `map_qualities`, `quality_scores`, `mv_values`).

## Help

To see the full list of options, run:

```bash
nanopile --help
```