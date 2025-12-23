use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub chromosome: String,
    pub start: usize, // 0-based, inclusive
    pub end: usize,   // 0-based, exclusive
}

impl Region {
    pub fn new(chromosome: String, start: usize, end: usize) -> Self {
        Self {
            chromosome,
            start,
            end,
        }
    }
}

impl FromStr for Region {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        // Expected format: chrom:start-end (1-based)
        let s = s.replace(',', "");
        let parts: Vec<&str> = s.split(':').collect();
        if parts.len() != 2 {
            return Err(anyhow::anyhow!("Invalid region format: {}", s));
        }
        let chrom = parts[0].to_string();
        let range_parts: Vec<&str> = parts[1].split('-').collect();
        if range_parts.len() != 2 {
            return Err(anyhow::anyhow!("Invalid region range format: {}", parts[1]));
        }

        let start_1based: usize = range_parts[0].parse().context("Invalid start coordinate")?;
        let end_1based: usize = range_parts[1].parse().context("Invalid end coordinate")?;

        if start_1based == 0 {
            return Err(anyhow::anyhow!(
                "Start coordinate must be > 0 for region string"
            ));
        }

        // Convert to 0-based
        Ok(Region::new(chrom, start_1based - 1, end_1based))
    }
}

pub fn parse_bed_file<P: AsRef<Path>>(path: P) -> Result<Vec<Region>> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref).with_context(|| {
        format!(
            "Failed to open BED file located at '{}'",
            path_ref.display()
        )
    })?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();

    for (line_idx, line_res) in reader.lines().enumerate() {
        let line_no = line_idx + 1;
        let line = line_res.with_context(|| {
            format!(
                "Failed to read line {} from BED file '{}'",
                line_no,
                path_ref.display()
            )
        })?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            eprintln!(
                "Warning: skipping malformed BED line {} in '{}': {}",
                line_no,
                path_ref.display(),
                line
            );
            continue;
        }

        let chrom = fields[0].to_string();
        let start: usize = fields[1].parse().with_context(|| {
            format!(
                "Invalid BED start coordinate at line {} in '{}'",
                line_no,
                path_ref.display()
            )
        })?;
        let end: usize = fields[2].parse().with_context(|| {
            format!(
                "Invalid BED end coordinate at line {} in '{}'",
                line_no,
                path_ref.display()
            )
        })?;

        regions.push(Region::new(chrom, start, end));
    }
    Ok(regions)
}
