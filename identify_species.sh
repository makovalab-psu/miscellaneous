#!/bin/bash
#runBlasr.sh mitochondrion_reference reads_1
set -e
set -x
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -t 0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=biomonika@psu.edu

#source paths if on server nn15
if [[ `uname -a` == *nn15* ]]; then 
    source /galaxy/home/biomonika/.bash_profile
fi

reads_1=$2
reads_2=$3

reference_mitochondrion_genome=$1

if [ ! -f ${reference_mitochondrion_genome} ]; then
    echo "Mitochondrial not provided!"
    exit
fi

if [ ! -f ${reads_1} ]; then
    echo "Library not provided!"
    exit
fi

hash blastn 2>/dev/null || { echo "I require blastn in path but it's not installed.  Aborting."; exit 1; }
hash makeblastdb 2>/dev/null || { echo "I require makeblastdb in path but it's not installed.  Aborting."; exit 1; }
hash blasr 2>/dev/null || { echo "I require blasr in path but it's not installed.  Aborting."; exit 1; }
hash seqtk 2>/dev/null || { echo "I require seqtk in path but it's not installed.  Aborting."; exit 1; }
hash samtools 2>/dev/null || { echo "I require samtools in path but it's not installed.  Aborting."; exit 1; }
path_to_picard=/nfs/brubeck.bx.psu.edu/scratch1/monika/prog/picard-tools-1/picard-tools-1.106/
path_to_seqtk=/galaxy/home/biomonika/seqtk/

reference_name=`basename $reference_mitochondrion_genome`
reads_name=`basename $reads_1`

SAM=${reads_name}_${reference_name}.sam
BAM=${reads_name}_${reference_name}.bam
SORTED=${reads_name}_${reference_name}.sorted
OUTPUT=${reads_name}_${reference_name}_mito_hits

if [ ! -f ${reference_mitochondrion_genome}.bwt ]; then
    echo "File not indexed! Indexing."
    bwa index $reference_mitochondrion_genome
fi

if [ ! -f `echo $reference_name`.nsq ]; then
    echo "Blast database doesn't exist. Creating one."
    makeblastdb -in $reference_mitochondrion_genome -out $reference_mitochondrion_genome -dbtype nucl -parse_seqids
fi

#bwa mem -t 63 -k 15 $reference_mitochondrion_genome $reads_1 $reads_2 >$SAM
blasr -nproc 63 -bestn 1 -sam -minPctIdentity 85 $reads_1 $reference_mitochondrion_genome -sa `echo $reference_name | sed 's/.fasta//'`.sa >$SAM

echo "Generated SAM file."

#remove unmapped reads and based on alignment quality
samtools view -bhS -F 4 -q 30 $SAM > $BAM
echo "Generated BAM file."

samtools sort $BAM $SORTED
echo "BAM file sorted."

wait;

#remove supplementary alignments
samtools view -bh -F 2048 ${SORTED}.bam >$BAM 
rm $SAM
rm ${SORTED}.bam

java -jar ${path_to_picard}/SamToFastq.jar I=$BAM F=${OUTPUT}.fastq
echo "Single-end (forward if paired-end) reads converted to fastq."

${path_to_seqtk}/seqtk seq -A ${OUTPUT}.fastq >${OUTPUT}.fasta

echo "Total number of sequences with mitochondrial hits: " >RESULTS_$OUTPUT
grep ">" ${OUTPUT}.fasta | wc -l >>RESULTS_$OUTPUT
echo "---" >>RESULTS_$OUTPUT

echo "BLAST ANALYSIS" >>RESULTS_$OUTPUT

#length of blast alignment at least 1000 bp, percentage identity at least 85
blastn -db $reference_name -query ${OUTPUT}.fasta -outfmt 6 -perc_identity 85 | sed 's/gi.*|//' | awk '{if ($4>1000) print $0}' >blast_$OUTPUT
cat blast_$OUTPUT | cut -f2 | sort | uniq -c | sort -rgk1 | head | tr -s " " | sed 's/^ //g' >candidates

cat candidates >>RESULTS_$OUTPUT; echo -e "" >>RESULTS_$OUTPUT

wait;

echo "---" >>RESULTS_$OUTPUT
echo "MAPPING ANALYSIS" >>RESULTS_$OUTPUT

samtools index $BAM
samtools idxstats $BAM | sed 's/gi.*|//' | awk '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}' | awk '{print $2 " " $1}' | sort -rgk1 | head >candidates
cat candidates >>RESULTS_$OUTPUT; echo -e "" >>RESULTS_$OUTPUT

rm candidates
rm ${OUTPUT}.fastq

echo "Done."

