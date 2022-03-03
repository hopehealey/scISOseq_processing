# scISOseq_processing
Scripts to process scISOseq data and run scRNAseq analysis.

## scISOr-Seq Workflow for final annotation 
![alt text](https://github.com/hopehealey/scISOseq_processing/blob/d861c92879b85fe2f6a1c182dab2536bfd21f989/scISOrSeq_paper_workflow.png)





## Folders 
Custom_Scripts: contains scripts written for the processing of scISOr-Seq data and assigning ensembl ids after tama merge.

scRNAseq_cellranger_and_analysis: contains scripts to run cell ranger and complete Seurat analysis

IsoSeq_processing: contains scripts written to process Iso-Seq and scISOr-Seq data using the Pacbio SMRT PacBio SMRTAnalysis software, minimap2,  cDNA cupcake, and sqanti3

Tama: scripts used to run tama merge and convert to gtf

## Important Notes

scISOr_Seq_processing.py requires that the input file be in fastq format. Additionally, if using the single cell flag, it requires a list of the cell barcodes. These barcodes can be taken from the barcovdes.tsv.gz file in the filtered feature bc matrix folder of a seurat output by using the below command.

```
zcat barcodes.tsv.gz | cut -f 1 -d "-" > barcodes_nodash.tsv
```
