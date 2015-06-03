import os
import os.path as op
import pandas as pd
# from collections import Counter

from bcbio.utils import rbind, file_exists, splitext_plus
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio import broad

def calc_variants_stats(data, args):
    in_vcf = data['vcf']
    ref_file = args.reference
    # gatk_jar = '/groups/bcbio/bcbio/toolplus/gatk/3.2-2-gec30cee/GenomeAnalysisTK.jar'
    jvm_opts = broad.get_gatk_framework_opts(data['config'])
    gatk_jar = config_utils.get_program("gatk", data['config'], "dir")
    bed_file = args.region
    sample = splitext_plus(op.basename(in_vcf))[0]
    in_bam = data['bam']
    cg_file = op.join(args.out, sample + "_with-gc.vcf.gz")
    parse_file = op.join(args.out, sample + "_cg-depth-parse.tsv")
    if not file_exists(cg_file):
        with file_transaction(cg_file) as tx_out:
            cmd = ("java -jar {gatk_jar}/GenomeAnalysisTK.jar -T VariantAnnotator -R {ref_file} "
                   "-L {bed_file} -I {in_bam} "
                   "-A GCContent --variant {in_vcf} --out {tx_out}")
            do.run(cmd.format(**locals()), " cg for %s" % in_vcf)

    if not file_exists(parse_file):
        with file_transaction(parse_file) as out_tx:
            with open(out_tx, 'w') as out_handle:
                print >>out_handle, "CG\tdepth\tsample"
            cmd = ("bcftools query -f '[%GC][\\t%DP][\\t%SAMPLE]\\n' -R  {bed_file} {cg_file} >> {out_tx}")
            do.run(cmd.format(**locals()), " query for %s" % in_vcf)
            logger.info('parsing coverage: %s' % sample)
    # return df
    return parse_file

