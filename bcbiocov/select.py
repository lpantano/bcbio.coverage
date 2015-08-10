import os
import pandas as pd
import six

import pybedtools

from bcbio.provenance import do
from bcbio.utils import rbind, file_exists
from bcbio.distributed.transaction import file_transaction

def _calc_regional_coverage(in_bam, chrom, start, end, samplename, work_dir):
    """
    given a BAM and a region, calculate the coverage for each base in that
    region. returns a pandas dataframe of the format:

    chrom position coverage name

    where the samplename column is the coverage at chrom:position
    """
    region_bt = pybedtools.BedTool("%s\t%s\t%s\n" % (chrom, start, end), from_string=True).saveas()
    region_file = region_bt.fn
    coords = "%s:%s-%s" % (chrom, start, end)
    tx_tmp_file = os.path.join(work_dir, "coverage-%s-%s.txt" % (samplename, coords.replace(":", "_")))
    cmd = ("samtools view -b {in_bam} {coords} | "
           "bedtools coverage -abam - -b {region_file} -d > {tx_tmp_file}")
    do.run(cmd.format(**locals()), "Plotting coverage for %s %s" % (samplename, coords))
    names = ["chom", "start", "end", "offset", "coverage"]
    df = pd.io.parsers.read_table(tx_tmp_file, sep="\t", header=None,
                                  names=names)
    os.remove(tx_tmp_file)
    df["sample"] = samplename
    df["chrom"] = chrom
    df["position"] = df["start"] + df["offset"] - 1
    return df[["chrom", "position", "coverage", "sample"]]

def _combine_regional_coverage(in_bams, samplenames, chrom, start, end, work_dir):
    """
    given a list of bam files, sample names and a region, calculate the
    coverage in the region for each of the samples and return a tidy pandas
    dataframe of the format:

    chrom position coverage name
    """
    dfs = [_calc_regional_coverage(bam, chrom, start, end, sample, work_dir) for bam, sample
           in zip(in_bams, samplenames)]
    return rbind(dfs)

def save_multiple_regions_coverage(samples, out_file, region_bed=None):
    """
    given a list of bcbio samples and a bed file or BedTool of regions,
    makes a plot of the coverage in the regions for the set of samples

    if given a bed file or BedTool of locations in stem_bed with a label,
    plots lollipops at those locations
    Adapted form Rory Kirchner (@roryk)
    """
    PAD = 100
    if file_exists(out_file):
        return out_file
    in_bams = samples
    samplenames = [os.path.splitext(os.path.basename(fn))[0] for fn in samples]
    if isinstance(region_bed, six.string_types):
        region_bed = pybedtools.BedTool(region_bed)
    with file_transaction(out_file) as tx_out_file:
        for line in region_bed:
            chrom = line.chrom
            start = max(line.start - PAD, 0)
            end = line.end + PAD
            df = _combine_regional_coverage(in_bams, samplenames, chrom,
                                            start, end, os.path.dirname(tx_out_file))
            df.to_csv(tx_out_file, mode='a', index=False, header=None)
    return out_file

