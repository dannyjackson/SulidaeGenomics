# RAxML of Candidate Introgressed Genes
cd /xdisk/mcnew/dannyjackson/sulidae/analyses/genetrees

GFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/bMorBas.EGAPx.gff

# GPR148
mkdir GPR148
cd GPR148
REGION=$(grep 'GPR148' "$GFF" | grep 'ID=gene' | awk '{print $1 ":" $4 "-" $5}')


TSPAN12
CTDP1