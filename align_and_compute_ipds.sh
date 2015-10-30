#align_and_compute_ipds.sh reference bax_file

reference=$1
bax_file=$2
alignment=${reference}_${bax_file}.cmp.h5

source /nfs/brubeck.bx.psu.edu/scratch4/software/smrtanalysis/current/etc/setup.sh #load path to the software needed

#align reads
pbalign -vv ${bax_file} ${reference} ${alignment} --forQuiver --metrics IPD,DeletionQV,DeletionTag,InsertionQV,MergeQV,SubstitutionQV

#call IPDs
ipdSummary.py --reference ${reference} ${alignment} --outfile ipds_${alignment}