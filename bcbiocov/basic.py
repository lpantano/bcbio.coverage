import os
import os.path as op

from bcbio.utils import file_exists, tmpfile, chdir
from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction

def bams(data):
    """
    samtools flagstat output
    """
    with chdir("flagstat"):
        bam = data["bam"]["ready"]
        sample = os.path.basename(bam).split("-")[0].replace(".bam", "")
        out_sample = op.join(sample + ".flagstat")
        if not file_exists(out_sample):
            with file_transaction(out_sample) as tx_out:
                cmd = ("samtools flagstat {bam} > {tx_out}")
                do.run(cmd.format(**locals()), "bam stats for %s" % bam)
                # print stats[sample]
        data["flagstat"] = op.abspath(out_sample)
        return [[data]]

