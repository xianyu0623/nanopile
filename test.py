import nanopile

pileup = nanopile.run_nanopile(
    bam_fp="/autofs/bal19/xyu/ont_open_data/HG001/hac_calls/downsample/HG001_10x.bam",
    ref_fp="/autofs/bal19/zxzheng/testData/ont/data/GRCh38_no_alt_analysis_set.fasta",
    regions=["chr20:63000-63100"],
    output_bq=True,
    output_read_name=True,
    output_mv=True,
    output_mapq=True,
)

for pos in pileup:
    print(pos.chrom, pos.pos + 1, pos.depth, pos.bases,  pos.mv_values)