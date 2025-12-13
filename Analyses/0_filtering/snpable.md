# snpable markdown

sbatch 4.1_snpable_ref.sh

PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta

awk '{print $1}' ${REF}.fai > ${PROJDIR}/referencelists/SCAFFOLDS.txt

grep 'CM' ${PROJDIR}/referencelists/SCAFFOLDS.txt > ${PROJDIR}/referencelists/CONTIGS.txt

# Z chromosome: 	CM062600.1
grep 'CM062600' ${PROJDIR}/referencelists/CONTIGS.txt > ${PROJDIR}/referencelists/Z_chrom.txt
# MT genome: 		CM062610.1
grep 'CM062610' ${PROJDIR}/referencelists/CONTIGS.txt > ${PROJDIR}/referencelists/MTgenome.txt

grep -vE 'CM062600|CM062610' ${PROJDIR}/referencelists/CONTIGS.txt > ${PROJDIR}/referencelists/AUTOSOMES.txt



# Next, you'll have to turn this into a "mappability mask" to be used with MSMC. To do this, edit the makeMappabilityMask.py script. You'll need to edit the paths indicated on lines 26 and 30. 
# On line 26, direct the script to the masked fasta file created by the run_snpable2 script (note: I recommend using the full path).
# On line 30, you need to specify the output location for the individual scaffold masks. In the path, use curly braces {} to indicate where the name of the scaffold should go. 

sbatch 4.2_snpable_mapabilitymask.sh

sbatch 4.3_phasevcf.sh


sbatch --array=1-29 4.2callvariants.sh 

