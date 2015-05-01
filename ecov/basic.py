import os

from collections import defaultdict
import pandas as pd

from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists


def calculate_bam(args):
    """
    samtools flagstat output
    """
    stats = defaultdict(list)
    samples = []
    for bam in args.bams:
        # sample = os.path.splitext(bam)[0].split("-")[0]
        sample = os.path.basename(bam).split("-")[0]
        out_sample = sample + ".flagstat"
        if not file_exists(out_sample):
            with file_transaction(out_sample) as tx_out:
                cmd = ("samtools flagstat {bam} > {tx_out}")
                do.run(cmd.format(**locals()), "bam stats for %s" % bam)
        with open(out_sample) as in_handle:
            for line in in_handle:
                if line.find("mapQ") == -1:
                    stats[line.strip().split(" + 0 ")[1].split("(")[0].strip()].append(line.strip().split(" + ")[0])
                # print stats[sample]
        samples.append(sample)
    with open(args.out, 'w') as out_handle:
        out_handle.write("\t".join(['measure'] + samples) + '\n')
        for feature in stats:
            out_handle.write("\t".join([feature] + stats[feature]) + "\n")

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

