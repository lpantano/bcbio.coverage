set -o pipefail  # trace ERR through pipes
set -o errtrace  # trace ERR through 'time command' and other functions
set -o nounset   ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

REF=$1
# ~/orch/groups/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa 

rm -rf work

mkdir work
cd work

coverage.py --run stats-coverage --out coverage --region ../regions.bed ../sample1.bam
coverage.py --run bias-coverage --region ../regions.bed --out bias ../sample1.bam --n_sample 20 --seed 42
coverage.py --run cg-vcf --region ../regions.bed --out cg ../sample1.vcf.gz  ../sample1.bam --reference $REF
coverage.py --run basic-bam --out basic-bam ../sample1.bam
coverage.py --run metrics ../final/2015-01-01_test/project-summary.yaml --out metrics 
coverage.py --run report --out report

cd ..
rm -rf work_complete
mkdir work_complete
cd work_complete

coverage.py --run complete ../sample1.bam ../sample1.vcf.gz ../final/2015-01-01_test/project-summary.yaml --out test --reference $REF --region ../regions.bed 
coverage.py --run fastqc --out fastqc ../sample1_data.txt ../sample1.bam 
