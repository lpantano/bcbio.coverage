import os
import os.path as op

from collections import defaultdict
import pandas as pd

from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, safe_makedir


def calculate_bam(data):
    """
    samtools flagstat output
    """
    bam = data[0]["bam"]["ready"]
    sample = os.path.basename(bam).split("-")[0].replace(".bam", "")
    out_sample = op.join(sample + ".flagstat")
    if not file_exists(out_sample):
        with file_transaction(out_sample) as tx_out:
            cmd = ("samtools flagstat {bam} > {tx_out}")
            do.run(cmd.format(**locals()), "bam stats for %s" % bam)
            # print stats[sample]
    data[0]["flagstat"] = op.abspath(out_sample)
    return [data]

def calculate_tstv(args):
    """
    get tstv from bcftools stat for all, known and new variants
    """
    tstv = defaultdict(list)
    for in_vcf in args.bams:
        out_file = os.path.splitext(in_vcf)[0] + ".stats"
        known_file = os.path.splitext(in_vcf)[0] + ".known.stats"
        new_file = os.path.splitext(in_vcf)[0] + ".new.stats"
        sample = os.path.basename(in_vcf).split("-")[0]
        if not file_exists(out_file):
            with file_transaction(out_file) as tx_out:
                cmd = ("bcftools stats {in_vcf} > {tx_out}")
                do.run(cmd.format(**locals()), "ts/tv ratio for %s" % in_vcf)
        if not file_exists(new_file):
            with file_transaction(new_file) as tx_new:
                cmd = ("bcftools filter -i DB=0 {in_vcf} | bcftools stats /dev/stdin > {tx_new}")
                do.run(cmd.format(**locals()), "ts/tv ratio for %s" % in_vcf)
        if not file_exists(known_file):
            with file_transaction(known_file) as tx_known:
                cmd = ("bcftools filter -i DB=1 {in_vcf} | bcftools stats /dev/stdin > {tx_known}")
                do.run(cmd.format(**locals()), "ts/tv ratio for %s" % in_vcf)
        for fn, name in zip([out_file, known_file, new_file], ['all', 'known', 'new']):
            with open(fn) as in_handle:
                for line in in_handle:
                    if line.startswith("TSTV"):
                        tstv[sample].append(line.split()[4])
                        break
    df = pd.DataFrame(tstv, index=['all', 'known', 'new'])
    df.to_csv(args.out)

