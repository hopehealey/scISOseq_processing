#!/usr/bin/env python3

'''
The purpose of this script is to take in genome annotation files and a sam file with reads (designed for scRNAseq reads, but could work for any reads).
It examines how many reads align to each gene (this is useful for scRNAseq because their output generally is organized by cell).
For each gene, it finds the exons of each transcript and finds the read counts of each exon. It gives information in two categories, read counts in the final exon and in all of the exons.
Since some genes have multiple exons, it only keeps the information that has the highest value for each gene. 
'''

from types import new_class
import pysam

### Alignment file
samfile = pysam.AlignmentFile("/Users/hopehealey/Dropbox (University of Oregon)/Hope_Dissertation_Folder/SingleCell/Gac_single_cell/70hpf_pilot_trial/Cell_ranger_ensGenome/210609_cellranger/outs/possorted_genome_bam.bam", "rb")
#gene_annots = "../../Gac_single_cell/70hpf_pilot_trial/Cell_ranger_ensGenome/genome_files/Gasterosteus_aculeatus.BROADS1.104.gtf"
#output_file = "genome_annotation_results.tsv"

import argparse
import regex as re

### argpase statements, need an annotation file, name of an output, and a word that describes the type of annotation file
parser=argparse.ArgumentParser("A program to have fun assessing where reads are landing in genes")
parser.add_argument("-annot1", "--annotation_file1", help="gtf file", required=True)
parser.add_argument("-out", "--output_file", help="name of the output file", required=True)
parser.add_argument("-name", "--annot_name", help="name that describes what annotation the gtf comes from", required=True)

args = parser.parse_args()


### generates some initial information
print("Total mapped reads", samfile.mapped)
print("Total unmapped reads", samfile.unmapped)


### setting up dictionaries for future use
genes_reads_max = {}
genes_reads_final_exon = {}
genes = []

### reading the gtf file, pulling out exon information, and saving it
with open(args.annotation_file1, "r") as a1:
    for line in a1:
        line = line.strip("/n")
        original_line = line
        
        ### only calling in lines with sequence information
        if "#!" not in line:
            line = line.split()
        #now keeping relevant information   
            
           
            ### isolating only current lines with transcript identity
            if line[2]=='transcript':
                    ### getting the geneid of the current line
                    geneid = line[9]
                    
                    ### (if reading comments for first time, go to the end of this if statement to understand it better)
                    ### making sure that the exons of the previous transcript have been loaded in the dictionary
                    if len(genes)!=0:
                        ### gathering some relevant information for future use
                        print(transcript)
                        orientation = transcript[6]
                        chromosome = transcript[0]

                        geneid = transcript[9]
                        geneid = geneid.strip('"')
                        geneid = geneid.strip(";")
                        geneid = geneid.strip('"')
                    
                        ex_count=0
                        ex_num=len(exons[tran_id])
                        ex_total=0
                        #ex_total_length=0
                        #print(tran_id)
                        #print(ex_num)                    
                        ### loop through each exon in the transcript
                        for j in exons[tran_id]:
                            print(j)
                            ### gathering more information
                            ex_count+=1
                            ex_start = int(j[3]) 
                            ex_end = int(j[4])
                            #length = abs(ex_end-ex_start)
                            
                            if orientation=="+" or orientation=="-": 
                            ### if the current exon is the LAST exon 
                                print("start", tran_id, geneid,ex_count, ex_num, ex_total)
                                print(chromosome, ex_start, ex_end, samfile.count(chromosome, ex_start, ex_end))
                                if ex_count==ex_num:
                                    ### update total read count in each exon and the length
                                    ex_total+=samfile.count(chromosome, ex_start, ex_end)
                                    #ex_total_length+=length
                                    
                                    ### updating the two dictionaries
                                    ### adding genes in there for the first time with the current counts
                                    ### if it is already there, then only changing the dictionary entry if the current count is higher
                                    if geneid not in genes_reads_final_exon: 
                                        genes_reads_final_exon[geneid]=samfile.count(chromosome, ex_start, ex_end)
                                    else: 
                                        
                                        current = genes_reads_final_exon[geneid]
                                        new = samfile.count(chromosome, ex_start, ex_end)
                                        print("Another final exon: current", current, "and now new", new)
                                        #print(ex_start, ex_end)
                                        if new > current:
                                            genes_reads_final_exon[geneid]=samfile.count(chromosome, ex_start, ex_end)
                                            print(tran_id,"replacing", current, "with", new)
                                        
                                    print(tran_id, ex_total)
                                    if geneid not in genes_reads_max:
                                        genes_reads_max[geneid]=ex_total
                                    else: 
                                        print("Another transcript: current",  genes_reads_max[geneid], "and now new", ex_total)
                                        current_total = genes_reads_max[geneid]
                                        if ex_total > current_total:
                                            genes_reads_max[geneid] = ex_total
                                            print(tran_id,"replacing", current_total, "with", ex_total)
                                ### if the current exon is NOT the final one, then just update transcript total counts and lengths
                                
                                else:
                                    ex_total+=samfile.count(chromosome, ex_start, ex_end)
                                    #ex_total_length+=length

                    #         elif orientation == "-": 
                    #         ### if the current exon is the LAST exon 
                    #             print("start", tran_id, geneid,ex_count, ex_num, ex_total)
                    #             print(chromosome, ex_start, ex_end, samfile.count(chromosome, ex_start, ex_end))
                    #             if ex_count==1:
                    #                 ### update total read count in each exon and the length
                    #                 ex_total+=samfile.count(chromosome, ex_start, ex_end)
                    #                 #ex_total_length+=length
                                    
                    #                 ### updating the two dictionaries
                    #                 ### adding genes in there for the first time with the current counts
                    #                 ### if it is already there, then only changing the dictionary entry if the current count is higher
                    #                 if geneid not in genes_reads_final_exon: 
                    #                     genes_reads_final_exon[geneid]=samfile.count(chromosome, ex_start, ex_end)
                    #                 else: 
                                        
                    #                     current = genes_reads_final_exon[geneid]
                    #                     new = samfile.count(chromosome, ex_start, ex_end)
                    #                     #print("Another transcript: current", current, "and now new", new)
                    #                     #print(ex_start, ex_end)
                    #                     if new > current:
                    #                         genes_reads_final_exon[geneid]=samfile.count(chromosome, ex_start, ex_end)
                    #                         #print("replacing", current, "with", new)
                    #                 if ex_count==ex_num:
                    #                     print("putting in reads dic", tran_id, ex_total)
                    #                     if geneid not in genes_reads_max:
                    #                             genes_reads_max[geneid]=ex_total
                    #                     else: 
                    #                         current_total = genes_reads_max[geneid]
                    #                         if ex_total > current_total:
                    #                             genes_reads_max[geneid]=ex_total
                                       
                    #             ### if the current exon is NOT the final one, then just update transcript total counts and lengths
                                
                    #             elif ex_count<ex_num:
                    #                 ex_total+=samfile.count(chromosome, ex_start, ex_end)
                    #                 #ex_total_length+=length
                    #             elif ex_count==ex_num:
                    #                 ex_total+=samfile.count(chromosome, ex_start, ex_end)

                                   
                    #                 if geneid not in genes_reads_max:
                    #                         genes_reads_max[geneid]=ex_total
                    #                         print("putting in reads dic", tran_id, ex_total)
                    #                 else: 
                    #                     current_total = genes_reads_max[geneid]
                    #                     if ex_total > current_total:
                    #                         genes_reads_max[geneid]=ex_total
                    #                         print("putting in reads dic", tran_id, ex_total)
                                
                    #             print("end", tran_id, geneid,ex_count, ex_total, samfile.count(chromosome, ex_start, ex_end))
        
                    # ### after going through the previous transcript, the next transcript is initialized
                    transcript = line 
                    genes.append(geneid)
                    exons = {}

                    # #print(len(genes))
                    # if (len(genes)%1000)==0:
                    #     print("1000 genes in")
                    #     #print(genes_reads_max)
                    #     #print(genes_reads_final_exon)

            ### saving the information of each exon by the transcript id        
            elif line[2]=='exon':
                tran_id = re.search('transcript_id ["A-Za-z0-9_.]+', original_line)
                tran_id = tran_id.group()

                #geneid = re.search('gene_id ["A-Za-z0-9_.]+', original_line)
                #geneid = geneid.group()
                
                if tran_id not in exons:
                    exons[tran_id]=[line]
                else:
                    exons[tran_id].append(line)

            # elif line[2]=='three_prime_utr':
            #     print(line)
            #     num_ex = len(exons[tran_id])
            #     print(num_ex)
            #     exons[tran_id][num_ex-1]=line
            #     print(exons[tran_id][num_ex-1])



### outputting the information into tsvs
with open(args.output_file, "w") as out:
    out.write("geneid" + "\t" + "maximum_reads_per_transcript" + "\t" + "maximum_reads_per_final_exon" + "\t" + "annotation" + "\n")
    for record in genes_reads_max:
        out.write(str(record) + "\t" + str(genes_reads_max[record]) + "\t" + str(genes_reads_final_exon[record]) + "\t" + str(args.annot_name) +"\n")



print(genes_reads_max)
print(genes_reads_final_exon)

samfile.close()

