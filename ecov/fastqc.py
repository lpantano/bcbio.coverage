import os.path as op

import pandas as pd
from fadapa import Fadapa

from bcbio.distributed.transaction import file_transaction
from bcbio.utils import rbind, file_exists, safe_makedir


def _get_module(fastq_list, module, wide=True):
    dt = []
    dt_together = []
    for sample in fastq_list:
        itern = fastq_list[sample].clean_data(module)
        header = itern[0]
        for data in itern[1:]:
            if data[0].startswith("#"):
                header = data
                continue
            if wide:
                if data[0].find("-") > -1:
                    f, s = map(int, data[0].split("-"))
                    for pos in range(f, s):
                        dt.append([str(pos)] + data[1:])
                else:
                    dt.append(data)
        dt = pd.DataFrame(dt)
        dt.columns = [h.replace(" ", "_") for h in header]
        dt['sample'] = sample
        dt_together.append(dt)
    dt_together = rbind(dt_together)
    return dt_together


def merge_fastq(data, args):
    """
    merge all fastqc samples into one by module
    """
    out_dir = safe_makedir(args.out)
    fastqc_list = {}
    for s in data:
        name = s[0]['name']
        fn = s[0]['fastqc']
        fastqc_list[name] = Fadapa(fn)

    module = [m[1] for m in fastqc_list[name].summary()][2:9]
    for m in module:
        out_file = op.join(out_dir, m.replace(" ", "_") + ".tsv")
        dt = _get_module(fastqc_list, m)
        dt.to_csv(out_file, sep="\t", index=False)
