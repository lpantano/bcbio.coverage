This package has a set of modules to get coverage statistics for targeted DNA sequencing.

Main modules:

* stats-coverage:  

`coverage.py --run stats-coverage --out out_dir --region region_file.bed *bam`

This will generate two files per sample in the out_dir folder:
`*_cov.csv` will have the summary data and `*dat` the raw coverage information 

* bias-coverage: 

`coverage.py --run bias-coverage --region ../regions.bed --out bias ../sample1.bam --n_sample 20 --seed 42`

this will generate the coverage bias information inside the `bias` folder. Each sample
will have `*_bias.tsv` file.

* cg-vcf:

`coverage.py --run cg-vcf --region ../regions.bed --out cg ../sample1.vcf.gz --reference  REFERENCE.FA`

this will generate `*cg-depth-parse.tsv` for each sample inside `cg` folder.


# Installation

if you already have a bcbio installation working:

`git clone https://github.com/lpantano/exomeCov.git`
`git clone https://github.com/lpantano/ich-wrapper.git`

export PYTHONPATH=$PATH_to_exomeCov:$PATH_to_ich-wrapper

# Usage
make sure to use the bcbio python:

`python /$PATH_to_exomeCov/scripts/coverage.py`

You can do an exec and add it to your path modifying the file `scripts/coverage`
