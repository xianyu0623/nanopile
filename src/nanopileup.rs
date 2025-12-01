use crate::region;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct BaseInfo {
    pub base: char,
    pub qual: u8,
    pub is_reverse: bool,
    pub insertion: Option<String>,
    pub deletion_len: Option<u32>,
    pub is_head: bool,
    pub is_tail: bool,
    pub mapq: u8,
    pub mv_value: Option<Vec<i32>>,
}

#[derive(Debug)]
pub struct CachedRead {
    pub _read_id: String,
    pub ref_start: i64,
    pub ref_end: i64,
    pub mapq: u8,
    pub seq_data: Vec<Option<BaseInfo>>,
}

impl CachedRead {
    pub fn new(record: &bam::Record, output_mv: bool) -> Result<Self> {
        //check if read seq is in the record if no skip this read
        if record.seq().len() == 0 {
            return Ok(Self {
                _read_id: String::from_utf8_lossy(record.qname()).to_string(),
                ref_start: record.pos(),
                ref_end: record.cigar().end_pos(),
                mapq: record.mapq(),
                seq_data: vec![],
            });
        }
        let read_id = String::from_utf8_lossy(record.qname()).to_string();
        let ref_start = record.pos();
        let ref_end = record.cigar().end_pos();
        let mapq = record.mapq();
        let is_reverse = record.is_reverse();

        // Initialize Vec with None
        let len = (ref_end - ref_start).max(0) as usize;
        let mut seq_data = vec![None; len];

        let qseq = record.seq();
        let qual = record.qual();
        let cigar = record.cigar();
        let mut ref_pos = ref_start;
        let mut query_pos = 0;
        let mut mv_per_query_base: Option<Vec<i32>> = None;
        if output_mv {
            if let Ok(mv_tag) = record.aux(b"mv") {
                let raw_mv_values: Vec<u8> = match mv_tag {
                    bam::record::Aux::ArrayU8(val) => val.iter().skip(1).map(|x| x as u8).collect(),
                    bam::record::Aux::ArrayI8(val) => val.iter().skip(1).map(|x| x as u8).collect(),
                    _ => Vec::new(),
                };

                if !raw_mv_values.is_empty() {
                    let qlen = qseq.len();
                    let mut counts = vec![0; qlen];
                    let mut base_idx: i32 = -1;

                    for move_val in raw_mv_values {
                        if move_val == 1 {
                            base_idx += 1;
                        }
                        if base_idx >= 0 && (base_idx as usize) < qlen {
                            counts[base_idx as usize] += 1;
                        }
                    }

                    if is_reverse {
                        counts.reverse();
                    }
                    mv_per_query_base = Some(counts);
                }
            }
        }
        // println!("Seq data mv: {:?}", mv_per_query_base);
        for cigar_entry in cigar.iter() {
            match cigar_entry {
                bam::record::Cigar::Match(len)
                | bam::record::Cigar::Equal(len)
                | bam::record::Cigar::Diff(len) => {
                    for _ in 0..*len {
                        let base_char =
                            b"=ACMGRSVTWYHKDBN"[qseq.encoded_base(query_pos) as usize] as char;
                        let q = qual[query_pos];

                        let idx = (ref_pos - ref_start) as usize;
                        if idx < seq_data.len() {
                            seq_data[idx] = Some(BaseInfo {
                                base: base_char,
                                qual: q,
                                is_reverse,
                                insertion: None,
                                deletion_len: None,
                                is_head: false,
                                is_tail: false,
                                mapq,
                                mv_value: mv_per_query_base
                                    .as_ref()
                                    .and_then(|v| v.get(query_pos).map(|&x| vec![x])),
                            });
                        }

                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                bam::record::Cigar::Ins(len) => {
                    if let Some(last_pos) = ref_pos.checked_sub(1) {
                        let idx = (last_pos - ref_start) as usize;
                        if idx < seq_data.len() {
                            if let Some(info) = seq_data[idx].as_mut() {
                                let mut ins_seq = String::new();
                                for _ in 0..*len {
                                    let base_char = b"=ACMGRSVTWYHKDBN"
                                        [qseq.encoded_base(query_pos) as usize]
                                        as char;
                                    ins_seq.push(base_char);

                                    if let Some(mvs) = &mv_per_query_base {
                                        if let Some(mv_vec) = info.mv_value.as_mut() {
                                            if let Some(&val) = mvs.get(query_pos) {
                                                mv_vec.push(val);
                                            }
                                        }
                                    }

                                    query_pos += 1;
                                }
                                info.insertion = Some(ins_seq);
                            } else {
                                query_pos += *len as usize;
                            }
                        } else {
                            query_pos += *len as usize;
                        }
                    } else {
                        query_pos += *len as usize;
                    }
                }
                bam::record::Cigar::Del(len) => {
                    if let Some(last_pos) = ref_pos.checked_sub(1) {
                        let idx = (last_pos - ref_start) as usize;
                        if idx < seq_data.len() {
                            if let Some(info) = seq_data[idx].as_mut() {
                                info.deletion_len = Some(*len);
                                if info.mv_value.is_some() {
                                    info.mv_value.as_mut().unwrap().push(0);
                                }
                            }
                        }
                    }
                    ref_pos += *len as i64;
                }
                bam::record::Cigar::RefSkip(len) => {
                    ref_pos += *len as i64;
                }
                bam::record::Cigar::SoftClip(len) => {
                    query_pos += *len as usize;
                }
                bam::record::Cigar::HardClip(_) => {}
                bam::record::Cigar::Pad(_) => {}
            }
        }

        // Mark head and tail
        // Find first non-None
        if let Some(info) = seq_data.iter_mut().find_map(|x| x.as_mut()) {
            info.is_head = true;
        }
        // Find last non-None
        if let Some(info) = seq_data.iter_mut().rev().find_map(|x| x.as_mut()) {
            info.is_tail = true;
        }

        Ok(Self {
            _read_id: read_id,
            ref_start,
            ref_end,
            mapq,
            seq_data,
        })
    }
}

pub struct ReadCache {
    pub reads: HashMap<String, CachedRead>,
}

impl ReadCache {
    pub fn new() -> Self {
        Self {
            reads: HashMap::new(),
        }
    }

    pub fn prune(&mut self, min_ref_pos: i64) {
        self.reads.retain(|_, read| read.ref_end > min_ref_pos);
    }
}

#[derive(Debug)]
pub struct PileupPos {
    pub chrom: String,
    pub pos: usize, // 0-based
    pub ref_base: char,
    pub depth: usize,
    pub bases: Vec<String>,
    pub read_names: Option<Vec<String>>,
    pub map_qualities: Option<Vec<u8>>,
    pub quality_scores: Option<Vec<u8>>,
    pub mv_values: Option<Vec<String>>,
}

impl PileupPos {
    pub fn new(
        chrom: String,
        pos: usize,
        ref_base: char,
        output_bq: bool,
        output_mapq: bool,
        output_read_name: bool,
        output_mv: bool,
    ) -> Self {
        Self {
            chrom,
            pos,
            ref_base,
            depth: 0,
            bases: Vec::new(),
            read_names: if output_read_name {
                Some(Vec::new())
            } else {
                None
            },
            map_qualities: if output_mapq { Some(Vec::new()) } else { None },
            quality_scores: if output_bq { Some(Vec::new()) } else { None },
            mv_values: if output_mv { Some(Vec::new()) } else { None },
        }
    }
}

pub fn nanopileup(
    bam_path: &PathBuf,
    region: &region::Region,
    ref_fp: Option<&PathBuf>,
    min_mapq: u8,
    min_baseq: u8,
    flag_filter: u32,
    buffer_size: usize,
    margin: usize,
    output_bq: bool,
    output_mapq: bool,
    output_read_name: bool,
    output_mv: bool,
) -> Result<Vec<PileupPos>> {
    let mut bam = bam::IndexedReader::from_path(bam_path)?;
    // let _header = bam.header().clone(); // Clone needed?

    // Load reference sequence for the region
    let ref_seq = if let Some(path) = ref_fp {
        if path.exists() {
            let fa_reader = faidx::Reader::from_path(path)?;
            // fetch_seq returns Ok(Vec<u8>)
            Some(fa_reader.fetch_seq_string(
                &region.chromosome,
                region.start as usize,
                region.end as usize - 1,
            )?)
        } else {
            None
        }
    } else {
        None
    };

    let mut cache = ReadCache::new();
    let mut results = Vec::new();

    let start = region.start;
    let end = region.end;

    for window_start in (start..end).step_by(buffer_size) {
        // println!("Window start: {}", window_start);
        let window_end = (window_start + buffer_size).min(end);

        let fetch_start = window_start.saturating_sub(margin);
        let fetch_end = window_end + margin;

        // Fetch reads
        bam.fetch((
            region.chromosome.as_bytes(),
            fetch_start as i64,
            fetch_end as i64,
        ))?;

        for result in bam.records() {
            // println!("Record: {:?}", result);
            let record = result?;
            // Skip if already in cache
            let read_id = String::from_utf8_lossy(record.qname()).to_string();
            if cache.reads.contains_key(&read_id) {
                continue;
            }

            // Filter by overlap with current window (strict)
            // fetch() gives loose overlap, we might want to be sure
            if record.pos() >= (window_end + margin) as i64
                || record.cigar().end_pos() <= (window_start.saturating_sub(margin)) as i64
            {
                continue;
            }

            if record.mapq() < min_mapq {
                continue;
            }
            if (record.flags() as u32) & flag_filter != 0 {
                continue;
            }

            let cached_read = CachedRead::new(&record, output_mv)?;
            // println!("Cached read: {:?}", cached_read);
            cache.reads.insert(read_id, cached_read);
        }

        // Prune cache: remove reads that end before this window starts
        // actually we need to keep reads that overlap the window.
        // If read.ref_end <= window_start, it's done.
        cache.prune(window_start as i64);

        // Generate pileup for [window_start, window_end)
        for pos in window_start..window_end {
            // Get ref base
            let ref_base = if let Some(seq) = &ref_seq {
                let offset = (pos - start) as usize;
                if offset < seq.len() {
                    seq.as_bytes()[offset] as char
                } else {
                    'N'
                }
            } else {
                'N'
            };
            let mut p = PileupPos::new(
                region.chromosome.clone(),
                pos as usize,
                ref_base,
                output_bq,
                output_mapq,
                output_read_name,
                output_mv,
            );

            // Collect active reads for this position
            // We need deterministic order.
            // Sorting all reads every time is slow.
            // Better: get all keys, sort them, then iterate.
            let mut active_read_ids: Vec<&String> = cache.reads.keys().collect();
            active_read_ids.sort(); // Deterministic

            for read_id in active_read_ids {
                let read = &cache.reads[read_id];
                if pos as i64 >= read.ref_start && (pos as i64) < read.ref_end {
                    // Check if we have base info
                    // Calculate index
                    let idx = (pos as i64 - read.ref_start) as usize;
                    if let Some(Some(info)) = read.seq_data.get(idx) {
                        if info.qual < min_baseq {
                            continue;
                        }
                        // Construct base string
                        let mut base_str = String::new();

                        // Start marker
                        if info.is_head {
                            base_str.push('^');
                            // mapq as char, usually +33
                            base_str.push((info.mapq + 33) as char);
                        }
                        // println!("Info base: {}", info.base);
                        // The base itself
                        let b = if info.is_reverse {
                            match info.base {
                                '=' => ',', // match reverse
                                _ => info.base.to_ascii_lowercase(),
                            }
                        } else {
                            match info.base {
                                '=' => '.', // match forward
                                _ => info.base,
                            }
                        };
                        base_str.push(b);

                        // Insertion
                        if let Some(ins) = &info.insertion {
                            base_str.push('+');
                            base_str.push_str(&ins.len().to_string());
                            //positive strand uppercase
                            if info.is_reverse {
                                base_str.push_str(&ins.to_ascii_lowercase());
                            } else {
                                base_str.push_str(&ins.to_ascii_uppercase());
                            }
                        }

                        // Deletion
                        if let Some(del_len) = info.deletion_len {
                            base_str.push('-');
                            base_str.push_str(&del_len.to_string());
                            if let Some(seq) = &ref_seq {
                                let start_del = (pos + 1 - start) as usize;
                                let end_del = start_del + del_len as usize;
                                if end_del <= seq.len() {
                                    base_str.push_str(&seq[start_del..end_del]);
                                }
                            }
                        }

                        // End marker
                        if info.is_tail {
                            base_str.push('$');
                        }

                        p.bases.push(base_str);
                        p.depth += 1;

                        if output_read_name {
                            if let Some(rn) = p.read_names.as_mut() {
                                rn.push(read._read_id.clone());
                            }
                        }
                        if output_mapq {
                            if let Some(mq) = p.map_qualities.as_mut() {
                                mq.push(info.mapq);
                            }
                        }
                        if output_bq {
                            if let Some(qs) = p.quality_scores.as_mut() {
                                qs.push(info.qual);
                            }
                        }
                        if output_mv {
                            if let Some(mvs) = p.mv_values.as_mut() {
                                if let Some(vals) = &info.mv_value {
                                    let mut s = String::new();
                                    for (i, v) in vals.iter().enumerate() {
                                        if i > 0 {
                                            s.push_str(",+");
                                        }
                                        s.push_str(&v.to_string());
                                    }
                                    mvs.push(s);
                                } else {
                                    mvs.push("0".to_string());
                                }
                            }
                        }
                    }
                }
            }
            results.push(p);
        }
    }

    Ok(results)
}
