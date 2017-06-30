


module load Molgenis-Compute/v16.11.1-Java-1.8.0_74
#module load Molgenis-Compute/v16.05.1-Java-1.8.0_45

sh convert.sh parameters.csv parameters.converted.csv

sh $EBROOTMOLGENISMINCOMPUTE/molgenis_compute.sh \
    --backend slurm \
    --generate \
    -p parameters.converted.csv \
    -p samplesheetPrePhasing.csv \
    -p BIOS_phasing/chromosomesAutosomalOnly.csv \
    -p individual_bam_link.txt \
    -w BIOS_phasing/workflow.csv \
    -rundir ./jobs/phasing/ \
    --weave 



cd ./jobs/phasing/
mkdir readbackedPhasing
for i in {22..1}; 
do 
  mkdir readbackedPhasing/chr$i
  echo $i; 
  grep -l chr$i Readback*sh | xargs mv -t readbackedPhasing/chr$i;  
done
cd -
#
#perl -pi -e 's/partition=dev/partition=long/g' BeagleGenotypeCalling*.sh
#perl -pi -e 's/partition=dev/partition=medium/g' Convert*.sh
#
#cd -





#--header BIOS_phasing/templates/slurm/header.ftl \


