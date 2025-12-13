#!/bin/bash

#SBATCH --job-name=SLURM_fasterqdump
#SBATCH --ntasks=8
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.SLURM_fasterqdump.%j

# downloading fastas from SRA

module load sratoolkit/3.1.1

echo -e 'SRR19149590\nSRR19149589\nSRR19149578\nSRR19149567\nSRR19149562\nSRR19149561\nSRR19149560\nSRR19149559\nSRR19149558\nSRR19149557\nSRR19149588\nSRR19149587\nSRR19149586\nSRR19149585\nSRR19149584\nSRR19149583\nSRR19149582\nSRR19149581\nSRR19149580\nSRR19149579\nSRR19149577\nSRR19149576\nSRR19149575\nSRR19149574\nSRR19149573\nSRR19149572\nSRR19149571\nSRR19149570\nSRR19149569\nSRR19149568\nSRR19149566\nSRR19149565\nSRR19149564\nSRR19149563' > /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/filenames_SRA.txt

cd /xdisk/mcnew/dannyjackson/sulidae/raw_sequences/

while read -r bird; do
fasterq-dump $bird 
done < filenames_SRA.txt
