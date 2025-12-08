# CHM13 edge liftover GRCh38 pos
screen out independent CHM13 graph edges and liftover to GRCh38 

## run pipeline
nextflow run main.nf \
-ansi-log false \
-resume \
-profile mccleary \
--edge_file data/chr{chr}.edge_list.txt \
--pangenie_decom_bi_vcf data/chr{chr}.var.vcf \
--liftover_vcf data/chr{chr}.liftover.vcf


