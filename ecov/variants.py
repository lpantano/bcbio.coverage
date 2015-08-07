import os
import os.path as op
import pandas as pd

from bcbio.utils import rbind, file_exists, splitext_plus, chdir
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.log import logger
from bcbio.pipeline import config_utils
from bcbio import broad


def variants(data):
    with chdir("cg"):
        if not data['vcf']:
            return [[]]
        in_vcf = data['vcf'].values()[0]
        ref_file = data['reference']
        assert ref_file, "Need the reference genome fasta file."
        jvm_opts = broad.get_gatk_framework_opts(data['config'])
        gatk_jar = config_utils.get_program("gatk", data['config'], "dir")
        bed_file = data['region']
        sample = splitext_plus(op.basename(in_vcf))[0]
        in_bam = data['bam']['ready']
        cg_file = op.join(sample + "_with-gc.vcf.gz")
        parse_file = op.join(sample + "_cg-depth-parse.tsv")
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
        return [[data]]

