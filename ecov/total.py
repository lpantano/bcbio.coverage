import os
import os.path as op
import pandas as pd
import subprocess
from collections import Counter

import pysam
import pybedtools

from bcbio.utils import file_exists, tmpfile, chdir
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger


class cov_class:

    def __init__(self, size, name, sample):
        self.size = size
        self.name = name
        self.sample = sample
        self.cov = {'10': 0, '25': 0, '50': 0}
        self.total = Counter()

    def update(self, size):
        self.size += size

    def save(self, cov, pt):
        if cov > 10:
            self.cov['10'] += pt
        if cov > 25:
            self.cov['25'] += pt
        if cov > 50:
            self.cov['50'] += pt

    def save_coverage(self, cov, nt):
        if cov > 100:
            cov = 100
        # self.size += size
        self.total[cov] += nt

    def write_coverage(self, out_file):
        # names = ["region", "size", "sample", "10", "25", "50"]
        df = pd.DataFrame({'depth': self.total.keys(), 'nt': self.total.values()})
        df["size"] = self.size
        df["sample"] = self.sample
        df.to_csv(out_file, mode='a', header=False, index=False, sep="\t")

    def write_completeness(self, out_file):
        # names = ["region", "size", "sample", "10", "25", "50"]
        df = pd.DataFrame(self.cov, index=["1"])
        df["region"] = self.name
        df["size"] = self.size
        df["sample"] = self.sample
        df.to_csv(out_file, mode='a', header=False, index=False, sep="\t")
        #return df

def _get_exome_coverage_stats(fn, sample, out_file, total_cov):
    tmp_region = ""
    stats = ""
    with open(fn) as in_handle:
        for line in in_handle:
            if line.startswith("all"):
                continue
            cols = line.strip().split()
            cur_region = "_".join(cols[0:2])
            if cur_region != tmp_region:
                if tmp_region != "":
                    stats.write_completeness(out_file)
                stats = cov_class(cols[-2], cur_region, sample)
            stats.save(int(cols[-4]), float(cols[-1]))
            total_cov.save_coverage(int(cols[-4]), int(cols[-3]))
            tmp_region = cur_region
        total_cov.update(int(cols[-2]))
        stats.write_completeness(out_file)
    return total_cov

def _silence_run(cmd):
    do._do_run(cmd, False)

def coverage(data):
    with chdir("coverage"):
        in_bam = data['bam']['ready']
        bed_file = data['region']
        sample = op.splitext(os.path.basename(in_bam))[0]
        region_bed = pybedtools.BedTool(bed_file)
        parse_file = op.join(sample + "_cov.tsv")
        parse_total_file = op.join(sample + "_cov_total.tsv")
        if not file_exists(parse_file):
            total_cov = cov_class(0, None, sample)
            bam_api = pysam.AlignmentFile(in_bam)
            with file_transaction(parse_file) as out_tx:
                with open(out_tx, 'w') as out_handle:
                    print >>out_handle, "q10\tq25\tq50\tregion\tsize\tsample"
                with tmpfile() as tx_tmp_file:
                    # tx_tmp_file = "tmpintersect"
                    for line in region_bed:
                        chrom = line.chrom
                        start = max(line.start, 0)
                        end = line.end
                        region_file = pybedtools.BedTool("%s\t%s\t%s\n" % (chrom, start, end), from_string=True).saveas().fn
                        coords = "%s:%s-%s" % (chrom, start, end)
                        cmd = ("samtools view -b {in_bam} {coords} | "
                               "bedtools coverage -a {region_file} -b - -hist > {tx_tmp_file}")
                        _silence_run(cmd.format(**locals()))
                        total_cov = _get_exome_coverage_stats(op.abspath(tx_tmp_file), sample, out_tx, total_cov)
            total_cov.write_coverage(parse_total_file)

        return [[data]]

