use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;

mod nanopileup;
mod region;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[clap(
        short = 'i',
        long = "bam_fp",
        help = "Input BAM file, must be sorted and indexed"
    )]
    bam_fp: PathBuf,

    #[clap(short = 'f', long = "ref_fp", help = "Input reference FASTA file")]
    ref_fp: Option<PathBuf>,

    #[clap(
        short = 'l',
        long = "bed_fp",
        help = "Input BED file (zero-based and half-open interval)",
        conflicts_with = "region",
        required_unless_present = "region"
    )]
    bed_fp: Option<PathBuf>,

    #[clap(
        short = 'r',
        long = "region",
        help = "Input region (1-based and inclusive at both ends, e.g. chr1:100-200). Can be specified multiple times.",
        conflicts_with = "bed_fp",
        required_unless_present = "bed_fp"
    )]
    region: Option<Vec<String>>,

    #[clap(
        long = "buffer_size",
        default_value_t = 10000,
        help = "Buffer size used when reading BAM file"
    )]
    buffer_size: usize,

    #[clap(
        long = "margin",
        default_value_t = 500,
        help = "Margin used when reading BAM file"
    )]
    margin: usize,

    #[clap(
        short = 'q',
        long = "min_mapq",
        default_value_t = 0,
        help = "Minimum mapping quality"
    )]
    min_mapq: u8,

    #[clap(
        short = 'Q',
        long = "min_baseq",
        default_value_t = 13,
        help = "Minimum base quality"
    )]
    min_baseq: u8,

    #[clap(long = "flag_filter", default_value_t = 0, help = "Flag filter")]
    flag_filter: u32,

    #[clap(
        long = "output_mv",
        default_value_t = false,
        help = "Output Move table values"
    )]
    output_mv: bool,

    #[clap(
        long = "output_bq",
        default_value_t = false,
        help = "Output Base Quality scores"
    )]
    output_bq: bool,

    #[clap(
        long = "output_mapq",
        default_value_t = false,
        help = "Output Mapping Quality scores"
    )]
    output_mapq: bool,

    #[clap(
        long = "output_read_name",
        default_value_t = false,
        help = "Output Read Names"
    )]
    output_read_name: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Parse regions
    let regions = if let Some(bed_fp) = args.bed_fp {
        region::parse_bed_file(&bed_fp)?
    } else if let Some(region_strs) = args.region {
        region_strs
            .iter()
            .map(|s| s.parse())
            .collect::<Result<Vec<region::Region>>>()?
    } else {
        unreachable!("Either region or bed_fp must be provided");
    };

    for region in regions {
        // println!("Region: {:?}", region);
        let results = nanopileup::nanopileup(
            &args.bam_fp,
            &region,
            args.ref_fp.as_ref(),
            args.min_mapq,
            args.min_baseq,
            args.flag_filter,
            args.buffer_size,
            args.margin,
            args.output_bq,
            args.output_mapq,
            args.output_read_name,
            args.output_mv,
        )?;

        for p in results {
            let bases_str = p.bases.join("");
            // let quals_str = ".".repeat(bases_str.len()); // Placeholder
            let mut output = format!(
                "{}\t{}\t{}\t{}\t{}",
                p.chrom,
                p.pos + 1, // 1-based output
                p.ref_base,
                p.depth,
                bases_str,
            );

            if let Some(rn) = p.read_names {
                output.push('\t');
                output.push_str(&rn.join(","));
            }
            if let Some(mq) = p.map_qualities {
                output.push('\t');
                output.push_str(
                    &mq.iter()
                        .map(|q| q.to_string())
                        .collect::<Vec<_>>()
                        .join(""),
                );
            }
            if let Some(qs) = p.quality_scores {
                output.push('\t');
                output.push_str(
                    &qs.iter()
                        .map(|q| ((*q + 33) as char).to_string())
                        .collect::<Vec<_>>()
                        .join(""),
                );
            }
            if let Some(mvs) = p.mv_values {
                output.push('\t');
                output.push_str(&mvs.join(";"));
            }

            println!("{}", output);
        }
    }

    Ok(())
}
