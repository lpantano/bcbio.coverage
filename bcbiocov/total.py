import os
import os.path as op
import pandas as pd
import subprocess
from collections import Counter

import numpy as np
import math
import pysam
import pybedtools

from bcbio.utils import file_exists, tmpfile, chdir
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger


class cov_class:

    def __init__(self, size, name, sample):
        self.size = int(size)
        self.name = name
        self.sample = sample
        self.cov = {'4': 0, '10': 0, '20': 0, '50': 0}
        self.total = Counter()
        self.raw = 0

    def update(self, size):
        self.size += size

    def save(self, cov, pt):
        self.raw += cov
        self.total[cov] = pt
        for cut in [4, 10, 20, 50]:
            if cov > cut:
                self.cov[str(cut)] += pt

    def save_coverage(self, cov, nt):
        if cov > 100:
            cov = 100
        elif cov > 10:
            cov =  int(math.ceil(cov / 10.0)) * 10
        # self.size += size
        self.total[cov] += nt

    def write_coverage(self, out_file):
        # names = ["region", "size", "sample", "10", "25", "50"]
        df = pd.DataFrame({'depth': self.total.keys(), 'nt': self.total.values()})
        df["size"] = self.size
        df["sample"] = self.sample
        df.to_csv(out_file, mode='a', header=False, index=False, sep="\t")

    def _noise(self):
        m = np.average(map(int, self.total.keys()), weights=self.total.values())
        x = []
        [x.extend([k] * int(float(v) * self.size)) for k, v in self.total.items()]
        sd = np.std(x)
        return m, sd

    def write_completeness(self, out_file):
        # names = ["region", "size", "sample", "10", "25", "50"]
        m, sd = self._noise()
        df = pd.DataFrame(self.cov, index=["1"])
        df['mean'] = m
        df['sdt'] = sd
        df['total'] = self.raw
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
                    print >>out_handle, "q10\tq20\tq4\tq50\tmean\tsdt\ttotal\tregion\tsize\tsample"
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

