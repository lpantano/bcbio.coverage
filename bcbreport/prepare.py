from collections import defaultdict
import string
import shutil
import glob

import os.path as op
# from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists
#, splitext_plus, safe_makedir, rbind, chdir


def _get_template(temp):
    template =  op.normpath(op.join(op.dirname(op.realpath(__file__)), temp + ".Rmd"))
    return open(template).read()

def report(out_dir):
    """
    create rmd template
    """
    content = _get_template("header")
    out_content = string.Template(content).safe_substitute({'path_results': op.abspath(out_dir)})
    out_file = op.join(out_dir, "report-ready.Rmd")
    with open(out_file, 'w') as out_handle:
        print >>out_handle, out_content

        if file_exists(op.join(out_dir, "qsignature.ma")):
            print >>out_handle, _get_template("qsignature")

        print >>out_handle, _get_template("qc")

        if glob.glob(op.join(out_dir, "coverage")):
            print >>out_handle, _get_template("regions")
        if glob.glob(op.join(out_dir, "variants")):
            print >>out_handle, _get_template("variants")

