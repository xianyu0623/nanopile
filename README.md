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

## Help

To see the full list of options, run:

```bash
nanopile --help
```