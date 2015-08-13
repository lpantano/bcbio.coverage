from collections import defaultdict
import pandas as pd
import string
import shutil
import glob

import os.path as op
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, splitext_plus, safe_makedir, rbind, chdir
from bcbio.log import logger
from bcbio.provenance import profile

from bcbiocov.manage import process_coverage, process_bams,process_variants
from bcbiocov.bias import calculate_bias_over_multiple_regions
from bcbiocov.select import save_multiple_regions_coverage
from bcbiocov import basic
from bcbiocov.fastqc import merge_fastq

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

def _bcbio_metrics(yaml_data):
    """
    parse project.yaml file to get metrics for each bam
    """
    project = yaml_data
    out_file = op.join("metrics", "metrics.tsv")
    dt_together = []
    with file_transaction(out_file) as out_tx:
        for s in project['samples']:
            m = s['summary']['metrics']
            for me in m:
                if isinstance(m[me], list):
                    m[me] = ":".join(m[me])
            dt = pd.DataFrame(m, index=['1'])
            # dt = pd.DataFrame.from_dict(m)
            dt.columns = [k.replace(" ", "_").replace("(", "").replace(")", "") for k in dt.columns]
            dt['sample'] = s['description']
            dt_together.append(dt)
        dt_together = rbind(dt_together)
        dt_together.to_csv(out_tx, index=False, sep="\t")

def bcbio_complete(run_parallel, samples, summary, dirs, qsignature=None, region=None):
    """
    samples is a list for each sample. sample[0] is a dict with the following keys:
        [name] = name
        [bam][ready] = bam file
        [vcf][caller] = vcf file
        [qc][fastqc] = data.txt from fastqc
        [region] = region in bed file
        [reference] = reference in fasta format
    """

    logger.info("copy qsignature")
    if qsignature:
        if file_exists(qsignature) and not file_exists("qsignature.ma"):
            shutil.copy(qsignature, "qsignature.ma")

    with profile.report("basic bam metrics", dirs):
        out_dir = safe_makedir("flagstat")
        samples = run_parallel(process_bams, samples)

        with chdir(out_dir):
            out_file = op.join("flagstat.tsv")
            stats = defaultdict(list)
            samples_name = []
            for sample in samples:
                samples_name.append(sample[0]['name'])
                with open(sample[0]["flagstat"]) as in_handle:
                    for line in in_handle:
                        if line.find("mapQ") == -1:
                            stats[line.strip().split(" + 0 ")[1].split("(")[0].strip()].append(line.strip().split(" + ")[0])
            with open(out_file, 'w') as out_handle:
                out_handle.write("\t".join(['measure'] + samples_name) + '\n')
                for feature in stats:
                    out_handle.write("\t".join([feature] + stats[feature]) + "\n")

    with profile.report("bcbio bam metrics", dirs):
        _bcbio_metrics(summary)

    with profile.report("bcbio fastq metrics", dirs):
        out_dir = safe_makedir("fastq")
        with chdir(out_dir):
            merge_fastq(samples)

    if region:
        with profile.report("doing coverage regions", dirs):
            out_dir = safe_makedir("coverage")
            samples = run_parallel(process_coverage, samples)

    with profile.report("doing cg-depth in vcf files", dirs):
        out_dir = safe_makedir("cg")
        run_parallel(process_variants, samples)

    # print "doing bias-coverage"
    # new_args = ['--run', 'bias-coverage', '--out', 'bias', '--region', args.region, '--n_sample', str(args.n_sample)] + galaxy + bam + cluster
    # new_args = params().parse_args(new_args)
    # data = _prepare_samples(new_args)
    # bias_exome_coverage(data, new_args)


    # print "doing report"
    report("report")
    return samples
