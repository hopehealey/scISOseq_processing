#!/usr/bin/env python
print("Hi hope")

import argparse

parser=argparse.ArgumentParser("A program to have fun with ISOseq and single cell isoseq by removing primers and polyA/T tails and dedupping")
parser.add_argument("-singlecell", "--single_cell_flag", help="says whether the reads are from a single cell experiment", required=False, type=bool)
parser.add_argument("-blast", "--blast_output", help="a file that has blast information of the reads against transcripts", required=False)
parser.add_argument("-seq", "--input_sequences", help="the actual read sequences", required=True)
parser.add_argument("-barcode", "--cell_barcodes", help="the cell barcodes from a previous cell ranger run", required=False)
parser.add_argument("-out", "--output_start", help="name to start the output files", required=False)
parser.add_argument("-3primer", "--p3_primer", help="the forward primer sequence", required=True)
parser.add_argument("-5primer", "--p5_primer", help="the reverse primer sequence", required=True)

args = parser.parse_args()


#reading in all the files
seqs=args.input_sequences
#seqs = "Gac_pacbio_test.fastq"
#seqs="4786.Gac_70dpf_PacBio--Gac_70dpf_PacBio.fastq"
if args.single_cell_flag:
    barcodes = args.cell_barcodes
if args.blast_output:
    blast_file = args.blast_output
    print("yes")

print(args.single_cell_flag)

#setting up output files
seqs_to_keep = args.output_start + ".fa"
#seqs_to_keep = "SCRIPTOUTPUT_4786.Gac_70dpf_PacBio--Gac_70dpf_PacBio.fasta"
if args.blast_output:
    if args.single_cell_flag:
        blast_cb_seqs_savers = args.output_start + "_cb+blast.fa"
if args.single_cell_flag:
    cell_barcode_savers = args.output_start + "_cb.fa"
    my_csv_table=  args.output_start + ".csv"
#my_csv_table = "SCRIPTOUTPUT_4786.Gac_70dpf_PacBio--Gac_70dpf_PacBio.csv"
    print(cell_barcode_savers)

#importing modules
from os import read
from Bio.Seq import Seq

#lengths of barcodes
cellB_len = 16
UMI_len = 12
Min_poly=20

#inputting primers
Seq_3p = args.p3_primer
Seq_3p_rev = str(Seq(Seq_3p).reverse_complement())
Seq_5p = args.p5_primer
Seq_5p_rev =  str(Seq(Seq_5p).reverse_complement())
Tstring = "TTTTTTTTTTTTTTTTTTTT"
Astring = "AAAAAAAAAAAAAAAAAAAA"

#saving stuff 
mysequences = {}
mysequence_info = {}

#counters 
num_reads_total = 0


#reading in and saving cell barcodes (note these are from a run with the ensembl version of the genome)
if args.single_cell_flag:
    print("here")
    molecule_counter=0
    cellb_save={}
    blast_save={}   
    cellbarcodes = {}
    with open(barcodes, "r") as cb:
        for line in cb:
            line=line.strip("\n")
    #lineR = Seq(line).reverse_complement()
            cellbarcodes[line]=1
    #cellbarcodes[str(lineR)]=1

###reading through blast results and saving the relevant components
if args.blast_output:
    my_blast_results = {}
    with open(blast_file, "r") as blastoff:
        for records in blastoff:
            records = records.strip()
            records = records.split()
            #print(records)
            #defining the key components (record, the transcript it's hitting to, the percent identity, the evalue, and where the hit starts and ends)
            record_id = records[0]
            transcript_id =records[1]
            percent_id = float(records[2])
            evalue = float(records[10])
            start=float(records[6])
            end=float(records[7])
            aligned_length=end-start #getting overall aligned length

            # #if it is present and the transcript is the same then there's a bit of work
            if record_id in my_blast_results and my_blast_results[record_id][0]==transcript_id:
                #defining the new vs old records
                old_record = my_blast_results[record_id]
                new_record = [transcript_id, percent_id, aligned_length, evalue, start, end]
                #if they don't have the same start and end site then proceed to update the start and end sites of the saved record
                if old_record[4]!= new_record[4] and old_record[5]!= new_record[5]:
                    #if the old records start value is larger than the new records start then 
                    #make mystart equal to the new record's start
                    #otherwise, keep the old records start
                    if old_record[4] > new_record[4]:
                        mystart=start
                    else:
                        mystart=old_record[4]
                    #if the old records end site is less than the new records end site
                    #then make myend equal to the new records end
                    #otherwise keep the old records end site
                    if old_record[5]< new_record[5]:
                        myend=end
                    else:
                        myend=old_record[5]
                    #averaging percentage identities and evalues
                    new_entry_per = (old_record[1]+percent_id)/2
                    #length is the smallest start and largest end site
                    new_entry_len = myend-mystart
                    new_entry_eval = (old_record[3]+evalue)/2
                    #defining new record
                    newest_record = [old_record[0], new_entry_per, new_entry_len, new_entry_eval, mystart, myend]
                    #adding newest record to the dictionary
                    my_blast_results[record_id]=newest_record
                
            #  #if it is present and the transcript is different then check for evalue and keep the one with the lower evalue
            elif record_id in my_blast_results and my_blast_results[record_id][0]!=transcript_id:
                current_record = my_blast_results[record_id]
                current_eval = current_record[3]
                #print("oldeval", current_eval)
                #if the new evalue is less than current then replace the current record
                if current_eval > evalue:
                    #print("old", record_id, my_blast_results[record_id])
                    my_blast_results[record_id]=[transcript_id, percent_id, aligned_length, evalue,start,end]
                    #print("new", record_id, my_blast_results[record_id])

            #if the record id isn't in the dictionary then put it in
            elif record_id not in my_blast_results:
                my_blast_results[record_id]=[transcript_id, percent_id, aligned_length, evalue, start,end]

print("whew, made it through loading in the blast information")

### reading through fastq (aka doing the fun stuff now)
count=0
with open(seqs, 'r') as fastq:
    for lines in fastq:
        count+=1
        lines=lines.strip()
        header = lines
        
        #checkpoint so you know the script is actually running
        if count%100000==0:     
            print("at", count, "lines out of 8,679,384")
            print("total reads", len(mysequences)+len(cellb_save)+len(blast_save))
            print("total 5p-3p", len(mysequences))
            if args.single_cell_flag:
                print("total cell barcode saviors", len(cellb_save))
            if args.blast_output:
                print("total reads to keep from blast", len(blast_save))

        #getting the header
        if count%4==1:
            lines=lines.strip("@")
            id = lines
            read_id=">" + lines
            #print(read_id)

        #isolating sequence line, sequence comes after header so we are ready to go once we have the sequence
        elif count%4==2:
            sequence = lines
            num_reads_total+=1

            #now checking for the 5p and 3p primers AND making sure that they only occur once
            if sequence.count(Seq_5p)==1 or sequence.count(Seq_5p_rev)==1:
                if sequence.count(Seq_3p)==1 or sequence.count(Seq_3p_rev)==1:
                    #print(len(lines))

                    ###finding the primers in the sequence
                    Seq_3p_start_site = sequence.find(Seq_3p)
                    #print(Seq_3p_start_site)
                    Seq_3p_rev_start_site = sequence.find(Seq_3p_rev)
                    #print(Seq_3p_rev_start_site)
                    Seq_5p_start_site = sequence.find(Seq_5p)
                    #print(Seq_5p_start_site)
                    Seq_5p_rev_start_site = sequence.find(Seq_5p_rev)
                    #print(Seq_5p_rev_start_site)

                    ###removing the primers from the sequence (doing it for each orientation separately)
                    ###captured POLYA TRANSCRIPTS
                    if Seq_3p_start_site==-1 and len(sequence)>(len(Seq_3p)+len(Seq_5p)):
                        #finding primers and removing them
                        remove_indices = Seq_3p_rev_start_site 
                        remove_indices2 = Seq_5p_rev_start_site + len(Seq_5p_rev)
                        primersremoved_sequence = sequence[remove_indices2:remove_indices]
                        #print(len(primersremoved_sequence))
                        #print("orignal", sequence)
                        #print("modified", primersremoved_sequence)
                        #print(primersremoved_sequence)

                        #figuring out barcodes based on the distance to the end (bc polyA)
                        if args.single_cell_flag:
                            cell_barcode_rev=primersremoved_sequence[-cellB_len:]
                            cell_barcode = Seq(cell_barcode_rev).reverse_complement()
                            #print(cell_barcode_rev, len(cell_barcode_rev))
                            UMI_rev=primersremoved_sequence[-(cellB_len+UMI_len):-cellB_len]
                            UMI = Seq(UMI_rev).reverse_complement()
                            #print(UMI_rev, len(UMI_rev))

                        #removing polyA tail (20bp minimum) and removing
                        
                        #print(len(Astring))
                        potential_end=primersremoved_sequence.find(Astring)
                        #print(primersremoved_sequence.rfind(Tstring))
                        final_sequence = primersremoved_sequence[:potential_end]
                        #print(final_sequence)

                    ###captured POLYT TRANSCRIPTS###
                    elif Seq_3p_rev_start_site==-1 and len(sequence)>(len(Seq_3p)+len(Seq_5p)): 
                        #removing primers
                        remove_indices= Seq_3p_start_site + len(Seq_3p)
                        remove_indices2= Seq_5p_start_site 
                        primersremoved_sequence = sequence[remove_indices:remove_indices2]
                        #print("length primers removed strategy", len(primersremoved_sequence))   
                        #print(primersremoved_sequence)
                        #print("cell", cell_barcode, len(cell_barcode))
                        #print("UMI", UMI, len(UMI))
                        #finding polyT tail and removing then reverse complementing the transcript
                        potential_end=primersremoved_sequence.rfind(Tstring)+20
                        #print(primersremoved_sequence.rfind(Tstring))
                        final_sequence = Seq(primersremoved_sequence[potential_end:]).reverse_complement()
                        
                        if args.single_cell_flag:
                            #finding cell barcode and UMI
                            cell_barcode=primersremoved_sequence[0:cellB_len]
                            cell_barcode_rev = Seq(cell_barcode).reverse_complement()
                            UMI=primersremoved_sequence[cellB_len:(UMI_len+cellB_len)]
                            UMI_rev = Seq(UMI).reverse_complement()

                    ### ALL DONE = putting stuff in the dictionaries (one that will make the fasta file and one that will make the csv file)
                    if len(final_sequence)>50:
            
                        if args.single_cell_flag:
                            dic_key=cell_barcode+UMI
                            if dic_key not in mysequence_info and len(final_sequence)>50:
                                mysequences[read_id]=str(final_sequence)
                                mysequence_info[dic_key]=["molecule/"+str(molecule_counter), str(UMI), str(UMI_rev), str(cell_barcode), str(cell_barcode_rev), len(final_sequence), 1]
                                molecule_counter+=1
                            elif dic_key in mysequence_info and len(final_sequence)>50:
                                print(read_id, "UMI:", str(UMI), "Cellbarcode:", str(cell_barcode), "umi and cb seen before")
                                mylist=mysequence_info[dic_key]
                                newcount=mylist[6]+1
                                mysequence_info[dic_key]=[mylist[0], str(UMI), str(UMI_rev), str(cell_barcode), str(cell_barcode_rev), len(final_sequence), newcount]
           
           ### Approach2  - retaining by cell barcodes
           #should only take reads that weren't captured by approach 1
            elif read_id not in mysequence_info and args.single_cell_flag==True:
                #looping through the cell barcodes and checking for their presence in the reads
                for i in cellbarcodes:
                    lineR = Seq(i).reverse_complement()
                    #if the cell barcode or its reverse comp is in the sequence then proceed and find it
                    if i in sequence or str(lineR) in sequence:
                        location=sequence.find(i)
                        locationrev=sequence.find(str(lineR))
                        #retaining reads with polyA (also making sure there are more polyA than T)
                        if locationrev > 1 and sequence.find(Astring)>1 and sequence.count(Astring)>sequence.count(Tstring):
                            cell_barcode=Seq(str(lineR)).reverse_complement()
                            cell_barcode_rev=str(lineR)
                            cell_bar_loc = sequence.find(cell_barcode_rev)
                            # finding UMI in relation to the cell barcode
                            UMI_rev=sequence[(cell_bar_loc-UMI_len):(cell_bar_loc)]
                            UMI=Seq(UMI_rev).reverse_complement()
                            # print(sequence)
                            # print(cell_barcode_rev)
                            # print(cell_bar_loc)
                            # print(UMI_rev, len(UMI_rev))
                            #finding polyA tail 
                            polyAloc = sequence.find(Astring)
                            potential_end=sequence.find(Astring)
                            #print("new end", sequence.find(Tstring))
                            final_sequence = Seq(sequence[len(Seq_5p):potential_end]) ###NOT BEST WAY BUT NEED TO ENSURE THAT 5p primer is gone
                            
                            #if the transcript isn't in my dictionary then add it as long as it's not 0bp long post clipping
                            dic_key=cell_barcode+UMI
                            if dic_key not in mysequence_info and len(final_sequence)>50 and ((sequence.rfind(Astring)+19)-sequence.find(UMI_rev))==0:
                                #print(read_id)
                                molecule_counter+=1
                                #mysequences[read_id]=str(final_sequence)
                                cellb_save[read_id]=str(final_sequence)
                                mysequence_info[dic_key]=["molecule/"+str(molecule_counter), str(UMI), str(UMI_rev), str(cell_barcode), str(cell_barcode_rev), len(final_sequence), 1]
                            #if it is there then update the molecule count and tell me more about the sequence so I can confirm the barcodes all line up
                            elif dic_key in mysequence_info and len(final_sequence)>50 and ((sequence.rfind(Astring)+19)-sequence.find(UMI_rev))==0:
                                print(read_id, "UMI:", str(UMI), "Cellbarcode:", str(cell_barcode), "umi and cb seen before")
                                mylist=mysequence_info[dic_key]
                                newcount=mylist[6]+1
                                mysequence_info[dic_key]=[mylist[0], mylist[1], mylist[2], mylist[3], mylist[4], mylist[5], newcount]
                                #print(mysequence_info[read_id])
                                #print(cell_barcode)
                                #print(len(final_sequence))
                                #print(UMI)
                                #print(mysequences[read_id])
                                #print(sequence)
                                #print(final_sequence)

                        #retaining reads with polyT
                        elif locationrev==-1 and sequence.rfind(Tstring)>1 and sequence.count(Tstring)>sequence.count(Astring): 
                            #finding barcodes and umis
                            cell_barcode=i
                            cell_barcode_rev=Seq(i).reverse_complement()
                            UMI=sequence[(location+cellB_len):(location+cellB_len+UMI_len)]
                            #print("loc", (location+cellB_len+UMI_len))
                            UMI_rev=Seq(UMI).reverse_complement()
                            #print(sequence)
                            #print("cellb", cell_barcode, len(cell_barcode))
                            #print("umi", UMI, len(UMI))
                            polyTloc = sequence.rfind(Tstring)
                            potential_end=sequence.rfind(Tstring)+20
                            #print("new end", sequence.find(Tstring))
                            #taking reverse complement of final sequence
                            final_sequence = Seq(sequence[potential_end:(-len(Seq_5p))]).reverse_complement()
                            #print(final_sequence)

                            dic_key=cell_barcode+UMI
                            if dic_key not in mysequence_info and len(final_sequence)>50 and sequence.find(Tstring)==(location+cellB_len+UMI_len):
                                #print(read_id)
                                #mysequences[read_id]=str(final_sequence)
                                molecule_counter+=1
                                cellb_save[read_id]=str(final_sequence)
                                mysequence_info[dic_key]=["molecule/"+str(molecule_counter), str(UMI), str(UMI_rev), str(cell_barcode), str(cell_barcode_rev), len(final_sequence), 1]
                            elif dic_key in mysequence_info and len(final_sequence)>50 and sequence.find(Tstring)==(location+cellB_len+UMI_len):
                                print(read_id, "UMI:", str(UMI), "Cellbarcode:", str(cell_barcode), "umi and cb seen before")
                                mylist=mysequence_info[dic_key]
                                newcount=mylist[6]+1
                                mysequence_info[dic_key]=[mylist[0], mylist[1], mylist[2], mylist[3], mylist[4], mylist[5], newcount]
                                #print(mysequence_info[read_id])
                                #print(cell_barcode)
                                #print(len(final_sequence))
                                #print(UMI)
                                #print(mysequences[read_id])
                                #print(sequence)
                                #print(final_sequence)

                #### Approach3 - blast comparison 
                ##only focusing on reads where they have a blast hit and the read id isn't found with approach 1 or 2
                if read_id not in cellb_save and read_id not in mysequences:
                    if args.blast_output:
                        if read_id.strip(">") in my_blast_results:
                            #grabbing parameters for quality filtering
                            #print(read_id.strip(">"))
                            current_record = my_blast_results[read_id.strip(">")]
                            #print(current_record)
                            current_percent = current_record[1]
                            transcript_id = current_record[0]
                            current_eval = current_record[3]
                            aligned_length = current_record[2]
                            read_length = len(sequence)
                            #taking ratio of the length that aligned to a transcript vs the read length to avoid chimeras 
                            ratio_lengths = (aligned_length/read_length)
                            #main quality filter: using percent identity, e-value, and the ratio of lengths
                            if current_percent > 95 and current_eval<(1e-20) and ratio_lengths > .8:
                                #identifying polyA transcripts
                                if sequence.count(Astring)>sequence.count(Tstring):
                                    #polyA
                                    polyAloc = sequence.find(Astring)
                                    potential_end=sequence.find(Astring)
                                    #print("new end", sequence.find(Tstring))
                                    final_sequence = Seq(sequence[len(Seq_5p):potential_end]) ###NOT BEST WAY BUT NEED TO ENSURE THAT 5p primer is gone
                                    #only keeping the fasta file (no contribution to the csv bc these shouldn't have cell barcodes recognized by cell ranger, could still have barcodes though)
                                    blast_save[read_id]=final_sequence
                                #identifying polyT transcripts
                                elif sequence.count(Tstring)>sequence.count(Astring):
                                    #polyT
                                    polyTloc = sequence.rfind(Tstring)
                                    potential_end=sequence.rfind(Tstring)+20
                                    #print("new end", sequence.find(Tstring))
                                    #saving the reverse complement
                                    final_sequence = Seq(sequence[potential_end:(-len(Seq_5p))]).reverse_complement()
                                    blast_save[read_id]=final_sequence
                                else:
                                    pass #don't want it without the polyA or T tail
                                
        

#printing out counters and stuff
print("Total number of reads:", num_reads_total)
print("Number of reads with 5p and 3p primers:", len(mysequences))
if args.blast_output:
    print("Reads to keep from blast", len(blast_save))
if args.single_cell_flag:
    print("Reads kept from cell barcode presence", len(cellb_save))
if args.single_cell_flag==True:
    if args.blast_output:
        print("Number of retained sequences", len(mysequences)+len(cellb_save)+len(blast_save))
#print(mysequences)
#print(mysequence_info)

#shuttling to output files
with open(seqs_to_keep, "w") as fasta:
    for key in mysequences:
        fasta.write(key + "\n" + mysequences[key] + "\n")

if args.single_cell_flag:
    with open(my_csv_table, "w") as csv:
        csv.write("sequence_id" + "\t" + "UMI" + "\t" + "UMIrev" + "\t" + "BC" + "\t" + "BCrev" + "\t" + "length" + "\t" + "count" + "\n")
        for key in mysequence_info:
            temp = mysequence_info[key]
            csv.write(temp[0].strip(">") + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\t" + temp[4] + "\t" + str(temp[5]) + "\t" + str(temp[6]) + "\n")

    with open(cell_barcode_savers, "w") as fasta:
        for key in cellb_save:
            fasta.write(key + "\n" + cellb_save[key] + "\n")
        for key in mysequences:
            fasta.write(key + "\n" + mysequences[key] + "\n")

if args.single_cell_flag==True:
    if args.blast_output:
        with open(blast_cb_seqs_savers, "w") as fasta:
            for key in cellb_save:
                fasta.write(key + "\n" + cellb_save[key] + "\n")
            for key in mysequences:
                fasta.write(key + "\n" + mysequences[key] + "\n")
            for key in blast_save:
                fasta.write(key + "\n" + str(blast_save[key]) + "\n")