use crate::nanopileup::PileupPos;
use crate::{nanopileup, region};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyModule;
use std::path::PathBuf;

const DEFAULT_BUFFER_SIZE: usize = 10_000;
const DEFAULT_MARGIN: usize = 500;
const DEFAULT_MIN_MAPQ: u8 = 0;
const DEFAULT_MIN_BASEQ: u8 = 13;
const DEFAULT_FLAG_FILTER: u32 = 0;

fn runtime_error(err: anyhow::Error) -> PyErr {
    PyRuntimeError::new_err(err.to_string())
}

fn collect_regions(
    bed_path: Option<&PathBuf>,
    region_strings: Option<Vec<String>>,
) -> PyResult<Vec<region::Region>> {
    match (bed_path, region_strings) {
        (Some(_), Some(_)) => Err(PyValueError::new_err(
            "Provide either `bed_fp` or `regions`, not both.",
        )),
        (None, None) => Err(PyValueError::new_err(
            "You must set `bed_fp` or supply at least one region string.",
        )),
        (Some(path), None) => region::parse_bed_file(path).map_err(runtime_error),
        (None, Some(region_list)) => region_list
            .into_iter()
            .map(|s| {
                s.parse::<region::Region>()
                    .map_err(|e| PyValueError::new_err(e.to_string()))
            })
            .collect(),
    }
}

#[pyclass]
pub struct PyPileupPos {
    #[pyo3(get)]
    chrom: String,
    #[pyo3(get)]
    pos: usize,
    #[pyo3(get)]
    ref_base: char,
    #[pyo3(get)]
    depth: usize,
    #[pyo3(get)]
    bases: Vec<String>,
    #[pyo3(get)]
    read_names: Option<Vec<String>>,
    #[pyo3(get)]
    map_qualities: Option<Vec<u8>>,
    #[pyo3(get)]
    quality_scores: Option<Vec<u8>>,
    #[pyo3(get)]
    mv_values: Option<Vec<String>>,
}

impl From<PileupPos> for PyPileupPos {
    fn from(pos: PileupPos) -> Self {
        Self {
            chrom: pos.chrom,
            pos: pos.pos,
            ref_base: pos.ref_base,
            depth: pos.depth,
            bases: pos.bases,
            read_names: pos.read_names,
            map_qualities: pos.map_qualities,
            quality_scores: pos.quality_scores,
            mv_values: pos.mv_values,
        }
    }
}

#[pyfunction]
#[pyo3(signature = (
    bam_fp,
    ref_fp=None,
    bed_fp=None,
    regions=None,
    buffer_size=DEFAULT_BUFFER_SIZE,
    margin=DEFAULT_MARGIN,
    min_mapq=DEFAULT_MIN_MAPQ,
    min_baseq=DEFAULT_MIN_BASEQ,
    flag_filter=DEFAULT_FLAG_FILTER,
    output_bq=false,
    output_mapq=false,
    output_read_name=false,
    output_mv=false,
))]
pub fn run_nanopile(
    bam_fp: &str,
    ref_fp: Option<&str>,
    bed_fp: Option<&str>,
    regions: Option<Vec<String>>,
    buffer_size: usize,
    margin: usize,
    min_mapq: u8,
    min_baseq: u8,
    flag_filter: u32,
    output_bq: bool,
    output_mapq: bool,
    output_read_name: bool,
    output_mv: bool,
) -> PyResult<Vec<PyPileupPos>> {
    let bam_path = PathBuf::from(bam_fp);
    let reference_path = ref_fp.map(PathBuf::from);
    let bed_path = bed_fp.map(PathBuf::from);

    let regions_to_process = collect_regions(bed_path.as_ref(), regions)?;

    let mut aggregated = Vec::new();
    for region in regions_to_process {
        let result = nanopileup::nanopileup(
            &bam_path,
            &region,
            reference_path.as_ref(),
            min_mapq,
            min_baseq,
            flag_filter,
            buffer_size,
            margin,
            output_bq,
            output_mapq,
            output_read_name,
            output_mv,
        )
        .map_err(runtime_error)?;

        aggregated.extend(result.into_iter().map(PyPileupPos::from));
    }

    Ok(aggregated)
}

#[pymodule]
fn nanopile(_py: Python, m: Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPileupPos>()?;
    m.add_function(wrap_pyfunction!(run_nanopile, m.clone())?)?;
    Ok(())
}

