set -e
module load tabix/0.2.6-foss-2015b
module load VCFtools/0.1.14-foss-2015b-Perl-5.22.0-bare 

projectDir=/groups/umcg-bios/scr01/umcg-ndeklein/
rnaseq_rare_variants=${projectDir}RNA-seq_rare_variants/
comparison_files=${rnaseq_rare_variants}/comparison_files/
output_folder=${comparison_files}outputs_QC_filter_RNA_seq_pass_filtered/
wgs_folder=${rnaseq_rare_variants}GoNL_WGS_calls/
wgs_postfix=.release5.NoMAFSelection
GQ=10
callrate=25
input_dir=filtered_vcf_GQ${GQ}_callrate${callrate}

# INCASE IT IS RERUN, REMOVE PREVIOUS RESULT AS NEW DATA GETS APPENDED TO IT
#	merged calling
rm -f ${output_folder}RNA_GQ10_call25_count.txt
rm -f ${output_folder}skipped*
rm -f ${output_folder}WGS_filtered_count.txt


mkdir -p ${output_folder}

for chr in {1..22}
	
do
    echo "Starting CHR ${chr}"
	#	tabix
    if [ ! -f ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf.gz ]; then
        echo "bgzip ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf.gz"
	    bgzip -c ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf  > ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf.gz
    	tabix -p vcf ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf.gz
    fi
	#	WGS
    echo "unzip ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf.gz, grep on PASS, write to ${output_folder}WGS_filtered_count.txt"
	gunzip -c ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix}.vcf.gz | grep 'PASS' | wc -l >> ${output_folder}WGS_filtered_count.txt

	##	Step 2: Read how many variants in WGS and RNA-seq
    echo "counting how many variants in ${input_dir}/lld_plus_gonl_chr${chr}.gg.PASSED.vcf.gz"
	gunzip -c ${input_dir}/lld_plus_gonl_chr${chr}.gg.PASSED.vcf.gz | grep 'GT:AD' | wc -l >> ${output_folder}/RNA_GQ${GQ}_call${callrate}_count.txt
	
	
	##	Step 3: Compare
	#	without individual calling
	java -jar //groups/umcg-bios/scr01/umcg-ndeklein/CompareGenotypeCalls.jar \
		-d1 ${wgs_folder}gonl-abc_samples.chr${chr}${wgs_postfix} \
		-D1 VCF \
		-d2 ${input_dir}/lld_plus_gonl_chr10.gg.PASSED.vcf.gz \
		-D2 VCF \
		-o ${output_folder}intersected_filter_GQ${GQ}_call${callrate}_chr${chr} \
		-s ${rnaseq_rare_variants}linking_file.txt 2>&1 | tee ${output_folder}intersected_GQ${GQ}_call${callrate}_chr${chr}.log


    awk -v x=$chr '/Number of SNPs compared:/{print x, "compared", $NF} /Skipped vars due to incosistant alleles:/{print x, "skipped", $NF}' ${output_folder}intersected_GQ${GQ}_call${callrate}_chr${chr}.log  >>  ${output_folder}skipped_merged_no_filter.txt

done
##############################################################################################################################################


#	Step 4: append output files:

#Rscript script_concatenate_variant_test.R

#	Step 5: Use R script to.gen.gzerate report plots and files

#Rscript script_rMarkdown_render.R


