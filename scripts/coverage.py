"""
calculate coverage across a list of regions
"""
import os
import glob
import shutil

from argparse import ArgumentParser
import yaml
import os.path as op

from bcbio import log
from bcbio.log import logger
from bcbio.install import _get_data_dir
from bcbio import utils
from bcbio.utils import file_exists, splitext_plus, safe_makedir, rbind, chdir
from bcbio.bam import is_bam
from bcbio.bam.fastq import is_fastq, combine_pairs
from bcbio.distributed import clargs, resources, prun
from bcbio.provenance import system
from bcbio.distributed.transaction import file_transaction
from bcbio.install import _get_data_dir

from bcbiocov.prepare import bcbio_complete, report

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
    return system_config

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

def bias_exome_coverage(data, args):
    safe_makedir(args.out)
    resources = {'name': 'bias', 'mem': 1, 'cores': 1}
    data = _update_algorithm(data, resources)
    cluster.send_job(calculate_bias_over_multiple_regions, data, args, resources)

def _complete(args):
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

def _read_final(yaml_file, args, config):
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
        data[sample]['name'] = sample
        fn_type, ext = splitext_plus(rel_path.split(os.sep)[1].replace(sample + "-", ""))
        if fn_type == "qc":
            is_there =  _read_qc_files(fn)
            if is_there:
                data[sample]["qc"] = is_there
            continue
        ext = ext.replace(".gz", "")[1:]
        if ext not in data[sample]:
            data[sample][ext] = {}
        data[sample][ext].update({fn_type: fn})
    for sample in data:
        data[sample]['config'] = config
        data[sample]["region"] = op.abspath(args.region) if args.region else ""
        data[sample]["reference"] = op.abspath(args.reference) if args.reference else ""
    return [[data[sample]] for sample in data]

def params():
    parser = ArgumentParser(description="Create file with coverage of a region")
    parser.add_argument("--region", help="bed file with regions.")
    parser.add_argument("--reference", help="genome fasta file.")
    parser.add_argument("--out", help="output file.")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("--run", required=True, help="type of analysis", choices=['all', 'report', 'bcbio'])
    parser.add_argument("--n_sample", default=1000, help="sample bed files with this number of lines")
    parser.add_argument("--seed", help="replication of sampling")
    parser.add_argument("-n", "--numcores", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("-c", "--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("-m", "--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("--retries", default=0, type=int,
                        help=("Number of retries of failed tasks during "
                              "distributed processing. Default 0 "
                              "(no retries)"))
    parser.add_argument("-s", "--scheduler", help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge", "pbspro"])
    parser.add_argument("-r", "--resources", help="Extra scheduler resource flags.", default=[], action="append")
    parser.add_argument("-q", "--queue", help="Queue to submit jobs to.")
    parser.add_argument("-p", "--tag", help="Tag name to label jobs on the cluster", default="bcb-cov")
    parser.add_argument("-t", "--paralleltype",
                        choices=["local", "ipython"],
                        default="local", help="Run with iptyhon")
    parser.add_argument("--galaxy", help="galaxy file.")


    return parser

if __name__ == "__main__":
    parser = params()
    args = parser.parse_args()
    assert args.files, "no files detected. Add vcf and bam files"

    system_config = _config(args)
    with open(system_config) as in_handle:
        config = yaml.load(in_handle)
        config["log_dir"] = os.path.join(os.path.abspath(os.getcwd()), "log")
        config["algorithm"] = {"num_cores": 1}

    parallel = clargs.to_parallel(args)
    samples = _read_final(args.files[0], args, config)
    project = _get_final_folder(args.files[0])
    final_folder = project['upload']
    summary = yaml.load(open(args.files[0]))
    qsignature = glob.glob(op.join(summary['upload'], "*/mixup_check/qsignature.ma"))
    qsignature = qsignature[0] if qsignature else None
    if args.run == "report":
        report(args)
    # elif args.run == "all":
    #    complete(args)
    elif args.run == "bcbio":

        parallel.update({'progs': ['samtools']})
        parallel = log.create_base_logger(config, parallel)
        log.setup_local_logging(config, parallel)
        dirs = {'work': os.path.abspath(os.getcwd())}
        system.write_info(dirs, parallel, config)
        sysinfo = system.machine_info()[0]
        parallel = resources.calculate(parallel, [samples], sysinfo, config)

        with prun.start(parallel, samples, config, dirs) as run_parallel:
            bcbio_complete(run_parallel, samples, summary, dirs, qsignature, args.region)
