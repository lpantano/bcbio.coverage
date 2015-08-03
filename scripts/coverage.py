"""
calculate coverage across a list of regions
"""
import os
from collections import namedtuple, defaultdict
import pandas as pd
import string
import glob
import shutil

# import six
from argparse import ArgumentParser
import yaml
import os.path as op
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib
# import seaborn as sns
from ichwrapper import cluster, arguments
# import pandas as pd
# from collections import defaultdict

# import pybedtools

from bcbio.utils import file_exists, splitext_plus, safe_makedir, rbind
# from bcbio.provenance import do
from bcbio.distributed.transaction import file_transaction
from bcbio.install import _get_data_dir

from ecov.total import _calc_total_exome_coverage
from ecov.bias import calculate_bias_over_multiple_regions
from ecov.variants import calc_variants_stats
from ecov.select import save_multiple_regions_coverage
from ecov.basic import calculate_bam, calculate_tstv
from ecov.fastqc import merge_fastq

def _find_bam(bam_files, sample):
    """
    Find the most similar file name
    """
    score = 0
    candidate = None
    for fn in bam_files:
        sc = sum(a == b for a, b in zip(op.basename(sample), op.basename(fn)))
        if sc > score:
            score = sc
            candidate = fn
    return candidate

def _find_fastqc(fq_files, sample):
    """
    Find the most similar file name
    """
    score = 0
    candidate = None
    for fn in fq_files:
        sc = sum(a == b for a, b in zip(op.dirname(sample).split("/")[-1], op.dirname(fn.replace("qc/fastqc/fastqc_data.txt", "")).split("/")[-1]))
        if sc > score:
            score = sc
            candidate = fn
    return candidate

def _update_algorithm(data, resources):
    """
    Update algorithm dict with new cores set
    """
    new_data = []
    for sample in data:
        sample[0]['config']['algorithm'] = resources
        new_data.append(sample)
    return new_data

def _config(args):
    if args.galaxy:
        system_config = args.galaxy
    else:
        system_config = op.join(_get_data_dir(), "galaxy", "bcbio_system.yaml")
    config = yaml.load(open(system_config))
    config['algorithm'] = {}
    return config

def _prepare_samples(args):
    """
    create dict for each sample having all information
    """
    config = _config(args)
    data = []
    vcf_files = [fn for fn in args.bams if fn.endswith('vcf.gz')]
    bam_files = [fn for fn in args.bams if fn.endswith('bam')]
    fastqc_files = [fn for fn in args.bams if fn.endswith('data.txt')]
    assert bam_files, "need the bam files for each sample"
    for sample in bam_files:
        if sample.endswith("yaml"):
            continue
        dt = {}
        dt['name'] = splitext_plus(op.basename(sample))[0].replace("-ready", "").replace("_data", "")
        dt['config'] = config
        dt['bam'] = op.abspath(sample)
        if vcf_files:
            dt['vcf'] = _find_bam(vcf_files, sample)
        if fastqc_files:
            dt['fastqc'] = _find_fastqc(fastqc_files, sample)
        data.append([dt])
    return data

def calculate_cg_depth_coverage(data, args):
    safe_makedir(args.out)
    resources = {'name': 'vcf_stats', 'mem': 2, 'cores': 1}
    data = _update_algorithm(data, resources)
    cluster.send_job(calc_variants_stats, data, args, resources)

def bias_exome_coverage(data, args):
    safe_makedir(args.out)
    resources = {'name': 'bias', 'mem': 1, 'cores': 1}
    data = _update_algorithm(data, resources)
    cluster.send_job(calculate_bias_over_multiple_regions, data, args, resources)

def average_exome_coverage(data, args):
    # dfs = [_calc_total_exome_coverage(bam, bed_file) for bam in in_bams]
    safe_makedir(args.out)
    resources = {'name': 'bedtools', 'mem': 8, 'cores': 1}
    data = _update_algorithm(data, resources)
    cluster.send_job(_calc_total_exome_coverage, data, args, resources)
    # df = rbind(dfs)
    # df.to_csv(tx_out_file, mode='a', index=False, header=["r10", "r25", "r50", "region", "size", "sample"])

def bcbio_metrics(args):
    """
    parse project.yaml file to get metrics for each bam
    """
    project = yaml.load(open(args.bams[0]))
    out_dir = safe_makedir(args.out)
    out_file = op.join(out_dir, "metrics.tsv")
    dt_together = []
    with file_transaction(out_file) as out_tx:
        for s in project['samples']:
            m = s['summary']['metrics']
            for me in m:
                if isinstance(m[me], list):
                    m[me] = ";".join(m[me])
            dt = pd.DataFrame.from_dict(m)
            dt.columns = [k.replace(" ", "_").replace("(", "").replace(")", "") for k in dt.columns]
            dt['sample'] = s['description']
            dt_together.append(dt)
        dt_together = rbind(dt_together)
        dt_together.to_csv(out_tx, index=False, sep="\t")

def report(out_dir):
    """
    create rmd template
    """
    out_dir = safe_makedir(out_dir)
    template = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../ecov/report.Rmd"))
    content = open(template).read()
    out_content = string.Template(content).safe_substitute({'path_results': op.abspath(".")})
    out_file = op.join(out_dir, "report-ready.Rmd")
    with open(out_file, 'w') as out_handle:
        print >>out_handle, out_content

def complete(args):
    """
    Run all modules
    """
    assert args.reference, "need the reference genome"
    assert args.bams, "no files detected. Add vcf and bam files"
    assert args.region, "need region bed file"

    data = _prepare_samples(args)
    vcf = [d[0]['vcf'] for d in data]
    bam = [d[0]['bam'] for d in data]
    fastqc = [d[0]['fastqc'] for d in data]
    yaml_file = [fn for fn in args.bams if fn.endswith("yaml")]

    assert len(vcf) == len(bam), "no paired bam/vcf files found. %s %s" % (vcf, bam)
    assert yaml_file, "No bcbio yaml file found."
    assert fastqc, "No fastqc files"

    cluster = []
    if args.scheduler:
        cluster = ['-n', args.numcores, '-s', args.scheduler, '-q', args.queue, '-p', args.tag, '-t', args.paralleltype]
        if args.resources:
            cluster += ['-r'] + args.resources
    cluster = map(str, cluster)

    print "doing basic-bam"
    new_args = ['--run', 'basic-bam', '--out', 'basic-bam'] + bam
    new_args = params().parse_args(new_args)
    calculate_bam(new_args)

    print "doing metrics"
    new_args = ['--run', 'metrics', '--out', 'metrics', yaml_file[0]]
    new_args = params().parse_args(new_args)
    bcbio_metrics(new_args)

    print "doing fastqc parsing"
    new_args = ['--run', 'fastqc', '--out', 'fastqc'] + fastqc
    new_args = params().parse_args(new_args)
    merge_fastq(data, new_args)

    print "doing stats-coverage"
    new_args = ['--run', 'stats-coverage', '--out', 'coverage', '--region', args.region] + bam + cluster
    new_args = params().parse_args(new_args)
    average_exome_coverage(data, new_args)

    print "doing bias-coverage"
    new_args = ['--run', 'bias-coverage', '--out', 'bias', '--region', args.region] + bam + cluster
    new_args = params().parse_args(new_args)
    bias_exome_coverage(data, new_args)

    print "doing cg-depth in vcf files"
    new_args = ['--run', 'cg-vcf', '--out', 'cg', '--region', args.region, '--reference', args.reference] + bam + vcf + cluster
    new_args = params().parse_args(new_args)
    calculate_cg_depth_coverage(data, new_args)

    print "doing report"
    report("report")

def _new_complete(args):
    data = _read_final(args.bams[0])

    assert args.reference, "need the reference genome"
    assert args.bams, "no files detected. Add vcf and bam files"
    assert args.region, "need region bed file"

    vcf_type = data.values()[0]['vcf'].keys()[0]
    vcf = [d['vcf'][vcf_type] for d in data.values() if 'vcf' in d]
    bam = [d['bam']['ready'] for d in data.values() if 'bam' in d]
    fastqc = [d['qc']['fastqc'] for d in data.values() if 'qc' in d]
    yaml_file = args.bams[0]

    assert len(vcf) == len(bam), "no paired bam/vcf files found. %s %s" % (vcf, bam)
    assert yaml_file, "No bcbio yaml file found."
    assert fastqc, "No fastqc files"

    cluster = []
    if args.scheduler:
        cluster = ['-n', args.numcores, '-s', args.scheduler, '-q', args.queue, '-p', args.tag, '-t', args.paralleltype]
        if args.resources:
            cluster += ['-r'] + args.resources
    cluster = map(str, cluster)
    galaxy = []
    if args.galaxy:
        galaxy = ['--galaxy', args.galaxy]

    print "copy qsignature"
    fn = glob.glob(op.join(_get_final_folder(yaml_file)['upload'], "*/mixup_check/qsignature.ma"))
    if file_exists(fn[0]) and not file_exists("qsignature.ma"):
        shutil.copy(fn[0], "qsignature.ma")

    print "doing basic-bam"
    new_args = ['--run', 'basic-bam', '--out', 'basic-bam'] + galaxy + bam
    new_args = params().parse_args(new_args)
    calculate_bam(new_args)

    print "doing metrics"
    new_args = ['--run', 'metrics', '--out', 'metrics', yaml_file] + galaxy
    new_args = params().parse_args(new_args)
    bcbio_metrics(new_args)

    print "doing fastqc parsing"
    new_args = ['--run', 'fastqc', '--out', 'fastqc'] + galaxy + fastqc + bam
    new_args = params().parse_args(new_args)
    data = _prepare_samples(new_args)
    merge_fastq(data, new_args)

    print "doing stats-coverage"
    new_args = ['--run', 'stats-coverage', '--out', 'coverage', '--region', args.region] + galaxy + bam + cluster
    new_args = params().parse_args(new_args)
    data = _prepare_samples(new_args)
    average_exome_coverage(data, new_args)

    # print "doing bias-coverage"
    # new_args = ['--run', 'bias-coverage', '--out', 'bias', '--region', args.region, '--n_sample', str(args.n_sample)] + galaxy + bam + cluster
    # new_args = params().parse_args(new_args)
    # data = _prepare_samples(new_args)
    # bias_exome_coverage(data, new_args)

    print "doing cg-depth in vcf files"
    new_args = ['--run', 'cg-vcf', '--out', 'cg', '--region', args.region, '--reference', args.reference] + galaxy + bam + vcf + cluster
    new_args = params().parse_args(new_args)
    data = _prepare_samples(new_args)
    calculate_cg_depth_coverage(data, new_args)

    print "doing report"
    report("report")

def _read_qc_files(qc_dir):
    """
    get the fastqc files from sample
    """
    qc = {}
    for fn in glob.glob(op.join(qc_dir, '*/*data.txt')):
        qc_fn = op.relpath(fn, qc_dir)
        qc_type = qc_fn.split(os.sep)[0]
        if qc_type == "fastqc":
            qc[qc_type] = fn
    return qc

def _get_final_folder(yaml_file):
    project = yaml.load(open(yaml_file))
    return project

def _read_final(yaml_file):
    """
    Get files for each sample from the bcbio upload folder
    """
    project = _get_final_folder(yaml_file)
    final = project['upload']
    samples =  [s['description'] for s in project['samples']]
    data = defaultdict(dict)
    print "bcbio results at %s" % final
    for fn in glob.glob(op.join(final, '*/*')):
        if fn.endswith('tbi') or fn.endswith('bai'):
            continue
        rel_path = op.relpath(fn, final)
        sample = rel_path.split(os.sep)[0]
        if sample not in samples:
            continue
        fn_type, ext = splitext_plus(rel_path.split(os.sep)[1].replace(sample + "-", ""))
        if fn_type == "qc":
            data[sample]["qc"] = _read_qc_files(fn)
            continue
        ext = ext.replace(".gz", "")[1:]
        if ext not in data[sample]:
            data[sample][ext] = {}
        data[sample][ext].update({fn_type: fn})
    return data

def params():
    parser = ArgumentParser(description="Create file with coverage of a region")
    parser = arguments.myargs(parser)
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("--reference", help="genome fasta file.")
    parser.add_argument("--out", help="output file.")
    parser.add_argument("bams", nargs="*", help="Bam files.")
    parser.add_argument("--run", required=True, help="type of analysis", choices=['complete', 'report', 'metrics', 'basic-bam', 'cg-vcf', 'stats-coverage', 'tstv', 'bias-coverage', 'plot', 'fastqc', 'final'])
    parser.add_argument("--n_sample", default=1000, help="sample bed files with this number of lines")
    parser.add_argument("--seed", help="replication of sampling")
    return parser


if __name__ == "__main__":
    parser = params()
    args = parser.parse_args()

    if args.run == "stats-coverage":
        data = _prepare_samples(args)
        average_exome_coverage(data, args)
    elif args.run == "bias-coverage":
        data = _prepare_samples(args)
        bias_exome_coverage(data, args)
    elif args.run == "tstv":
        calculate_tstv(args)
    elif args.run == "basic-bam":
        calculate_bam(args)
    elif args.run == "metrics":
        bcbio_metrics(args)
    elif args.run == "cg-vcf":
        data = _prepare_samples(args)
        calculate_cg_depth_coverage(data, args)
    elif args.run == "plot":
        save_multiple_regions_coverage(args.bams, args.out, args.region)
    elif args.run == "fastqc":
        data = _prepare_samples(args)
        merge_fastq(data, args)
    elif args.run == "report":
        report(args)
    elif args.run == "complete":
        complete(args)
    elif args.run == "final":
        _new_complete(args)
