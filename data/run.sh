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

coverage.py --run stats-coverage --out stats --region ../regions.bed ../sample1.bam
coverage.py --run bias-coverage --region ../regions.bed --out bias ../sample1.bam --n_sample 20 --seed 42
coverage.py --run cg-vcf --region ../regions.bed --out cg ../sample1.vcf.gz --reference $REF
