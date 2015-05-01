"""
calculate coverage across a list of regions
"""
import os

import six
from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib
# import seaborn as sns
import pandas as pd
import numpy as np
# from collections import Counter, defaultdict

import pybedtools

from bcbio import bam
from bcbio.utils import rbind, file_exists
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction


def _calc_regional_coverage(in_bam, samplename, chrom, start, end, work_dir):
    """
    given a BAM and a region, calculate the coverage for each base in that
    region. returns a pandas dataframe of the format:

    chrom position coverage name

    where the samplename column is the coverage at chrom:position
    """
    # bam.index(in_bam, {'algorithm': {}})
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
    mean = np.mean(np.array(df["coverage"]))
    std = np.std(np.array(df["coverage"]))
    # print df["coverage"]
    # print mean
    ntup = sum(np.array(df["coverage"]) > mean + std * 3)
    ntdown = sum(np.array(df["coverage"]) < mean - std * 3)
    dfs = {}
    dfs["mean"] = mean
    dfs["std"] = std
    dfs["ntup"] = ntup
    dfs["ntdow"] = ntdown
    dfs["size"] = len(df["coverage"])
    dfs["sample"] = samplename
    dfs["region"] = coords
    df = pd.DataFrame(dfs, index=['1'])
    # print df
    return df

def _sample_bed(bed_file, n=1000, seed=None):
    out_file = os.path.splitext(bed_file)[0] + "-sample.bed"
    if seed:
        seed = " -seed " + str(seed)
    cmd = ("bedtools sample -n {n} -i {bed_file} {seed} > {tx_out} ")
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out:
        do.run(cmd.format(**locals()), "Sampling %s from %s" % (n, bed_file))
    return out_file

def calculate_bias_over_multiple_regions(in_bam, args):
    """
    given a list of bcbio samples and a bed file or BedTool of regions,
    makes a plot of the coverage in the regions for the set of samples

    if given a bed file or BedTool of locations in stem_bed with a label,
    plots lollipops at those locations
    Adapted form Rory Kirchner (@roryk)
    """
    PAD = 0

    samplename = os.path.splitext(os.path.basename(in_bam))[0]
    out_file = os.path.splitext(os.path.basename(args.out))[0] + "-%s.tsv" % samplename
    if file_exists(out_file):
        return out_file
    if isinstance(args.region, six.string_types):
        region_bed = _sample_bed(args.region, args.n_sample, args.seed)
        region_bed = pybedtools.BedTool(region_bed)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, 'w') as in_handle:
            print >> in_handle, "mean,ntdow,ntup,region,sample,size,std"

        for line in region_bed:
            chrom = line.chrom
            start = max(line.start - PAD, 0)
            end = line.end + PAD
            df = _calc_regional_coverage(in_bam, samplename, chrom,
                                            start, end, os.path.dirname(tx_out_file))
            df.to_csv(tx_out_file, mode='a', index=False, header=None)
    return out_file

