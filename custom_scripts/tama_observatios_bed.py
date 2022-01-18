#!/usr/bin/env python
print("Hi!")

import pysam

gtf_name="Gac_white_70hpf_ens_as_1.bed"
#ens_g_names="ens_id_as1_w_G_names_locs.txt"
new="Gac_white_70hpf_ens_as_1_distributions.csv"
samfile = pysam.AlignmentFile("/Users/hopehealey/Dropbox (University of Oregon)/Hope_Dissertation_Folder/SingleCell/Gac_single_cell/70hpf_pilot_trial/Cell_ranger_ensGenome/210609_cellranger/outs/possorted_genome_bam.bam", "rb")

geneids = {}
transids = {}
#getting ids

import regex as re
problem_children={}
problem_children_locations = {}
problem_children_found = {}


#sizes
size_lib = []
end_lib = []
start_lib = []

with open (gtf_name, "r") as names:
    for line in names:
        line=line.strip()
        line_in = line.split("\t")
        moreline=line_in[3]
        moreline=(moreline.split(";")) 

        chro=line_in[0]
        start=line_in[1]
        end = line_in[2]

        #ex start sites
        print("line", line)
        ex_st_sites = line_in[11].split(',')
        
        #print(ex_st_sites)

        #ex lengths
        ex_length = line_in[10].split(',')
        #print(ex_length)
        #print(ex_st_sites[0], line_in[1])
        if line_in[5]=="+":
            last_exon_start = int(ex_st_sites[-1]) + int(line_in[1])
            last_exon_end = last_exon_start + int(ex_length[-1])
            first_exon_start =  int(ex_st_sites[0]) + int(line_in[1])
            first_exon_end = first_exon_start + int(ex_length[0])
            print(first_exon_start, first_exon_end)

        elif line_in[5]=="-":
            last_exon_start = int(ex_st_sites[0]) + int(line_in[1])
            last_exon_end = last_exon_start + int(ex_length[0])

            print("last start", last_exon_start, "last end", last_exon_end)
            
            first_exon_start =  int(ex_st_sites[-1]) + int(line_in[1])
            first_exon_end = first_exon_start + int(ex_length[-1])
            print("first start", first_exon_start, "first end",first_exon_end)

        #print(last_exon_start, last_exon_end)

        gname = moreline[0]
        tname = moreline[1]

        if len(moreline)>2:
            ensgname=moreline[2]
            enstname=moreline[3]
            if gname not in geneids:
                geneids[gname]=[chro, start, end, gname, tname,ensgname, enstname, last_exon_start, last_exon_end, first_exon_start, first_exon_end]
            elif gname in geneids:
                if geneids[gname][5]!= ensgname:
                    print("conflict", geneids[gname][5], geneids[gname][0], geneids[gname][1], geneids[gname][2], ensgname, chro, start, end)
                    if gname not in problem_children:
                        problem_children[gname]=tname
                        problem_children_locations[gname]=[[chro, start, end, gname, tname, ensgname, enstname]]
                        problem_children_locations[gname].append([geneids[gname][0],geneids[gname][1], geneids[gname][2],geneids[gname][3], geneids[gname][4], geneids[gname][5], geneids[gname][6]])
                        print(problem_children_locations[gname])
                    if gname in problem_children:
                        problem_children_locations[gname].append([chro, start, end, gname, tname, ensgname, enstname])
                    #problem_children_locations[(geneids[gname][0],geneids[gname][1], geneids[gname][2])]=[geneids[gname][3], geneids[gname][4], geneids[gname][5], geneids[gname][6]]




print(len(problem_children))
#print(problem_children_locations)
print(len(problem_children_locations))
print("done file 1")

c=0
a=0
transcripts=0

with open(gtf_name, "r") as gtf, open(new, "w") as stats:
    stats.write("gene_name" + ","+ "size_dif" + "," + "start_diff" + ","+ "end_diff"+ ","+"iso_final_minus_ens_final" + "," + "iso_first_minus_ens_first" + "\n")
    for line in gtf:
        line=line.strip()
        safeline=line
        line=line.split()
        orientation = line[5]
        #print(line[0])
        #print(line)
        names = line[3].split(';')
        #print(gname)
        gname = names[0]
        #print(gname)

        #ex start sites
        ex_st_sites = line[11].split(',')
        #print(ex_st_sites)

        #ex lengths
        ex_length = line[10].split(',')
        #print(ex_length)

        chromosome = line[0]

        if gname not in problem_children and gname in geneids and len(names)<3:
            a+=1
            # print("ADDING")
            # print(names)
            #pulling info [chro, start, end, gname, tname,ensgname, enstname, last_exon_start, last_exon_end]


            #info is ensmbl
            myinfo=geneids[gname]
            #isoseq driven annotation - ensembl
            #
            if orientation=="+":
                size_diff =  abs(int(line[2])-int(line[1])) - abs(int(myinfo[2])-int(myinfo[1]))
                #size_lib.append(size_diff)

                start_diff = int(line[1]) - int(myinfo[1])
                #start_lib.append(start_diff)
                
                end_diff = int(line[2]) - int(myinfo[2])
                #end_lib.append(end_diff)
                ####getting READ COUNT difference for LAST exons in ISO-ENS
                start = int(line[1])
                last_exon_start = int(ex_st_sites[-1]) + int(line[1])
                last_exon_end = last_exon_start + int(ex_length[-1])
                # print("iso last exon", chromosome, last_exon_start, last_exon_end)
                # print("ens last exon", myinfo[0], myinfo[7], myinfo[8])
                counts = samfile.count(chromosome, last_exon_start, last_exon_end)

                ens_counts = samfile.count(myinfo[0], myinfo[7], myinfo[8])

                ####getting READ COUNT difference for FIRST exons in ISO-ENS
                first_exon_start =  int(ex_st_sites[0]) + int(line[1])
                first_exon_end = first_exon_start + int(ex_length[0])
                first_counts = samfile.count(chromosome, first_exon_start, first_exon_end)
                first_ens_counts = samfile.count(myinfo[0], myinfo[9], myinfo[10])


                # print("last counts", counts)
                # print("last ens counts", ens_counts)
                # print("iso first exon", chromosome, first_exon_start, first_exon_end)
                # print("ens first exon", myinfo[0], myinfo[9], myinfo[10])
                # print("first counts", first_counts)
                # print("first ens counts", first_ens_counts)                

            if orientation=="-":
                #print(line)
                size_diff =  -abs(int(line[2])-int(line[1])) + abs(int(myinfo[2])-int(myinfo[1]))
                end_diff = -int(line[1]) + int(myinfo[1])
                start_diff = -int(line[2]) + int(myinfo[2])

                last_exon_start = int(ex_st_sites[0]) + int(line[1])
                #print("last start", last_exon_start)
                last_exon_end = last_exon_start + int(ex_length[0])

                counts = samfile.count(chromosome, last_exon_start, last_exon_end)

                ens_counts = samfile.count(myinfo[0], myinfo[7], myinfo[8])

                # print("iso last exon", chromosome, last_exon_start, last_exon_end)
                # print("ens last exon", myinfo[0], myinfo[7], myinfo[8])
                first_exon_start =  int(ex_st_sites[-1]) + int(line[1])
                first_exon_end = first_exon_start + int(ex_length[-1])
                # print("iso first exon", chromosome, first_exon_start, first_exon_end)
                # print("ens first exon", myinfo[0], myinfo[9], myinfo[10])
                first_counts = samfile.count(chromosome, first_exon_start, first_exon_end)
                first_ens_counts = samfile.count(myinfo[0], myinfo[9], myinfo[10])

                
                # print("counts", counts)
                # print("ens counts", ens_counts)


                # print("first counts", first_counts)
                # print(" first ens counts", first_ens_counts)                
                # print(end_diff)
                # print(start_diff)
                # print(size_diff)


            #print(line)
            #print(myinfo)

            # print(end_diff)
            # print(start_diff)
            # print(size_diff)

            end_lib.append(end_diff)
            start_lib.append(start_diff)
            size_lib.append(size_diff)
            count_dif = counts - ens_counts
            first_count_dif = first_counts - first_ens_counts
            stats.write(line[3] + "," + str(size_diff) + "," + str(start_diff) + "," + str(end_diff) + "," + str(count_dif) + "," + str(first_count_dif) + "\n")



print(sum(size_lib)/len(size_lib))
print(min(size_lib))
print(max(size_lib))

print(sum(start_lib)/len(start_lib))
print(min(start_lib))
print(max(start_lib))

print(sum(end_lib)/len(end_lib))
print(min(end_lib))
print(max(end_lib))


