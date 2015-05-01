"""
calculate coverage across a list of regions
"""
import os

# import six
from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib
# import seaborn as sns
from ichwrapper import cluster, arguments
# import pandas as pd
# from collections import defaultdict

# import pybedtools

from bcbio.utils import file_exists
# from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

from ecov.total import _calc_total_exome_coverage
from ecov.bias import calculate_bias_over_multiple_regions
from ecov.variants import calc_variants_stats
from ecov.select import save_multiple_regions_coverage
from ecov.basic import calculate_bam, calculate_tstv

def calculate_genes_per_vcf(args):
    """
    count number of genes with variants
    """

def calculate_cg_depth_coverage(args):
    resources = {'name': 'vcf_stats', 'mem': 1, 'cores': 1}
    cluster.send_job(calc_variants_stats, args.bams, args, resources)

def bias_exome_coverage(args):
    resources = {'name': 'bias', 'mem': 1, 'cores': 1}
    cluster.send_job(calculate_bias_over_multiple_regions, args.bams, args, resources)

def average_exome_coverage(args):
    out_file = args.out
    if file_exists(out_file):
        return out_file
    with file_transaction(out_file) as tx_out_file:
        # dfs = [_calc_total_exome_coverage(bam, bed_file) for bam in in_bams]
        resources = {'name': 'bedtools', 'mem': 12, 'cores': 1}
        cluster.send_job(_calc_total_exome_coverage, args.bams, args, resources)
        # df = rbind(dfs)
        # df.to_csv(tx_out_file, mode='a', index=False, header=["r10", "r25", "r50", "region", "size", "sample"])
    return out_file


if __name__ == "__main__":
    parser = ArgumentParser(description="Create file with coverage of a region")
    parser = arguments.myargs(parser)
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("--reference", help="genome fasta file.")
    parser.add_argument("--out", required=True, help="output file.")
    parser.add_argument("bams", nargs="*", help="Bam files.")
    parser.add_argument("--run", required=True, help="type of analysis", choices=['basic-bam', 'basic-vcf', 'stats', 'tstv', 'bias', 'plot'])
    # parser.add_argument("--basic-bam", action="store_true", help="Calculate bam stats")
    # parser.add_argument("--basic-vcf", action="store_true", help="Calculate vcf stats")
    # parser.add_argument("--stats", action="store_true", help="Calculate stats for all regions.")
    # parser.add_argument("--tstv", action="store_true", help="Calculate ts/tv ratio.")
    # parser.add_argument("--bias", action="store_true", help="Calculate bias for all regions.")
    # parser.add_argument("--plot", action="store_true", help="Calculate nt coverage for given regions.")
    parser.add_argument("--n_sample", default=1000, help="sample bed files with this number of lines")
    parser.add_argument("--seed", help="replication of sampling")
    args = parser.parse_args()

    if os.path.exists(args.out):
        os.remove(args.out)

    if args.run == "stats":
        average_exome_coverage(args)
    elif args.run == "bias":
        bias_exome_coverage(args)
    elif args.run == "tstv":
        calculate_tstv(args)
    elif args.run == "basic-bam":
        calculate_bam(args)
    elif args.run == "basic-vcf":
        calculate_cg_depth_coverage(args)
    elif args.run == "plot":
        save_multiple_regions_coverage(args.bams, args.out, args.region)
