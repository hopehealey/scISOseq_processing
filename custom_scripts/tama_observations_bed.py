#!/usr/bin/env python
print("Hi!")

gtf_name="Gac_white_70hpf_ens_as_1.bed"
ens_g_names="ens_id_as1_w_G_names_locs.txt"
new="Gac_white_70hpf_ens_as_1_distributions.csv"

geneids = {}
transids = {}
#getting ids

import regex as re
genes_w_multiple_ens_ids_for_one_g_id={}
genes_w_multiple_ens_ids_for_one_g_id_locations = {}
genes_w_multiple_ens_ids_for_one_g_id_found = {}


#sizes
size_lib = []
end_lib = []
start_lib = []

with open (ens_g_names, "r") as names:
    for line in names:
        line=line.strip()
        line_in = line.split("\t")
        moreline=line_in[3]
        moreline=(moreline.split(";")) 

        chro=line_in[0]
        start=line_in[1]
        end = line_in[2]
        gname = moreline[0]
        tname = moreline[1]

        if len(moreline)>2:
            ensgname=moreline[2]
            enstname=moreline[3]
            if gname not in geneids:
                geneids[gname]=[chro, start, end, gname, tname,ensgname, enstname]
            elif gname in geneids:
                if geneids[gname][5]!= ensgname:
                    print("conflict", geneids[gname][5], geneids[gname][0], geneids[gname][1], geneids[gname][2], ensgname, chro, start, end)
                    if gname not in genes_w_multiple_ens_ids_for_one_g_id:
                        genes_w_multiple_ens_ids_for_one_g_id[gname]=tname
                        genes_w_multiple_ens_ids_for_one_g_id_locations[gname]=[[chro, start, end, gname, tname, ensgname, enstname]]
                        genes_w_multiple_ens_ids_for_one_g_id_locations[gname].append([geneids[gname][0],geneids[gname][1], geneids[gname][2],geneids[gname][3], geneids[gname][4], geneids[gname][5], geneids[gname][6]])
                        print(genes_w_multiple_ens_ids_for_one_g_id_locations[gname])
                    if gname in genes_w_multiple_ens_ids_for_one_g_id:
                        genes_w_multiple_ens_ids_for_one_g_id_locations[gname].append([chro, start, end, gname, tname, ensgname, enstname])
                    #genes_w_multiple_ens_ids_for_one_g_id_locations[(geneids[gname][0],geneids[gname][1], geneids[gname][2])]=[geneids[gname][3], geneids[gname][4], geneids[gname][5], geneids[gname][6]]




print(len(genes_w_multiple_ens_ids_for_one_g_id))
#print(genes_w_multiple_ens_ids_for_one_g_id_locations)
print(len(genes_w_multiple_ens_ids_for_one_g_id_locations))
print("done file 1")

c=0
a=0
transcripts=0

with open(gtf_name, "r") as gtf, open(new, "w") as stats:
    stats.write("gene_name" + ","+ "size_dif" + "," + "start_diff" + ","+ "end_diff"+ "\n")
    for line in gtf:
        line=line.strip()
        safeline=line
        line=line.split()
        #print(line[0])
        #print(line)
        names = line[3].split(';')
        #print(gname)
        gname = names[0]
        #print(gname)
        if gname not in genes_w_multiple_ens_ids_for_one_g_id and gname in geneids and len(names)<3:
            a+=1
            #print(names)
            #pulling info [chro, start, end, gname, tname,ensgname, enstname]


            #info is ensmbl
            myinfo=geneids[gname]
            #isoseq driven annotation - ensembl
            #
            size_diff =  abs(int(line[2])-int(line[1])) - abs(int(myinfo[2])-int(myinfo[1]))
            #size_lib.append(size_diff)

            start_diff = int(line[1]) - int(myinfo[1])
            #start_lib.append(start_diff)
            
            end_diff = int(line[2]) - int(myinfo[2])
            #end_lib.append(end_diff)

            print(line)
            print(myinfo)

            print(end_diff)
            print(start_diff)
            print(size_diff)

            end_lib.append(end_diff)
            start_lib.append(start_diff)
            size_lib.append(size_diff)
                
            stats.write(line[3] + "," + str(size_diff) + "," + str(start_diff) + "," + str(end_diff) + "\n")



print(sum(size_lib)/len(size_lib))
print(min(size_lib))
print(max(size_lib))

print(sum(start_lib)/len(start_lib))
print(min(start_lib))
print(max(start_lib))

print(sum(end_lib)/len(end_lib))
print(min(end_lib))
print(max(end_lib))


