import os
import pandas as pd
# from collections import Counter

from bcbio.utils import rbind, file_exists
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger


class cov_class:

    def __init__(self, size, name, sample):
        self.size = size
        self.name = name
        self.sample = sample
        self.cov = {'10': 0, '25': 0, '50': 0}

    def save(self, cov, pt):
        if cov > 10:
            self.cov['10'] += pt
        if cov > 25:
            self.cov['25'] += pt
        if cov > 50:
            self.cov['50'] += pt

    def dataframe(self, out_file):
        # names = ["region", "size", "sample", "10", "25", "50"]
        df = pd.DataFrame(self.cov, index=["1"])
        df["region"] = self.name
        df["size"] = self.size
        df["sample"] = self.sample
        df.to_csv(out_file, mode='a', index=False, header=None)
        #return df

def _get_exome_coverage_stats(fn, sample, out_file):
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
                    stats.dataframe(out_file)
                stats = cov_class(cols[-2], cur_region, sample)
            stats.save(int(cols[-4]), float(cols[-1]))
            tmp_region = cur_region
        stats.dataframe(out_file)
    # return [st.dataframe(out_file) for st in stats.values()]

def _calc_total_exome_coverage(in_bam, args):
    bed_file = args.region
    sample = os.path.splitext(os.path.basename(in_bam))[0]
    cov_file = sample + ".dat"
    parse_file = sample + "_cov.csv"
    if not file_exists(cov_file):
        with file_transaction(cov_file) as cov_tx:
            cmd = ("bedtools coverage -abam {in_bam} -b {bed_file} -hist > {cov_tx}")
            do.run(cmd.format(**locals()), "exome coverage for %s" % in_bam)
    if not file_exists(parse_file):
        with file_transaction(parse_file) as out_tx:
            logger.info('parsing coverage: %s' % sample)
            _get_exome_coverage_stats(cov_file, sample, out_tx)
    # return df
    return cov_file

