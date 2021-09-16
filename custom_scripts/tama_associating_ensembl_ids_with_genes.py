#!/usr/bin/env python
print("Hi!")

gtf_name="Gac_white_70hpf_ens_as_1.bed"
ens_g_names="ens_id_as1_w_G_names_locs.txt"
new="Gac_white_70hpf_ens_as_1_modnames.bed"

geneids = {}
transids = {}
#getting ids

import regex as re
genes_w_multiple_ens_ids_for_one_g_id={}
genes_w_multiple_ens_ids_for_one_g_id_locations = {}
genes_w_multiple_ens_ids_for_one_g_id_found = {}

#reading in gene names to get the corresponding ensembl id for each gene name as well as it's location
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

        #if the gene id has an associated ensembl id then save it in a dictionary for future use
        if len(moreline)>2:
            ensgname=moreline[2]
            enstname=moreline[3]
            #if the gene id is not in the dictionary already then add it in
            if gname not in geneids:
                geneids[gname]=[chro, start, end, gname, tname,ensgname, enstname]
            #if the gene id is already in the dictionary (i.e. when there are multiple splice variants for one gene id) 
            #then check if the ensembl id associated with that gene name is the same as the ensembl id associated with the current line
            #if it is the same, then we don't need to change anything (this dictionary will just be used to associated ensembl names with gene ids)
            #if it is different, then print out the potential conflict (there will be some cases of conflict)
            #add the gene information to a dictionary which will be examined later and include the start and end sites, also adding in the same information from gene id that it originally matched to
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

#now opening up the real bed file and the bed file we are creating
#reading through each line of the bed file 
with open(gtf_name, "r") as gtf, open(new, "w") as newby:
    for line in gtf:
        line=line.strip()
        safeline=line
        line=line.split()
        #print(line[0])
        #print(line)
        gname = line[3].split(';')
        gname = gname[0]
        #print(gname)

        #If the gene id is not in the odd cases (multiple ensembl ids for one gene id) and the gene id does have an ensembl id associated with it
        #then proceed
        #replace the gene id with the ensembl id 
        #and print the line out to the new file
        if gname not in genes_w_multiple_ens_ids_for_one_g_id and gname in geneids:
            a+=1
            #pulling info [chro, start, end, gname, tname,ensgname, enstname]
            myinfo=geneids[gname]
            #print(line)
            #print(re.sub('(gname)', 'myinfo[5]', safeline))
            name=gname
            #print(name)
            line[3]=line[3].replace(name, myinfo[5])
            #print(line)
            s = "\t".join(line)
            newby.write(s + "\n")
            #print("done")

        #if there are multiple ensembl ids for the one gene id...
        #try to find the ensembl id that is the closest distance to the particular transcript
        #looking at both start and end sites and picking the ensembl id that is the closest for both
        #however, we added in a second filter to prevent spurious naming
        #we allowed only best matches within a certain distance of the original ensembl annotation to take on that id
        #Because 3' ends of transcripts seem to be more variable with Iso-Seq data, we allowed a distance of 600bp from the end of the ensembl annotation
        #Since there is less variability at the 5' ends, we only allowed 50bp of difference between the ensembl model and this model
        #we then repeated writing out to the new file as done in the previous if statement
        elif gname in genes_w_multiple_ens_ids_for_one_g_id:
             c+=1
             #print(gname, line)
             chro=line[0]
             start=float(line[1])
             end=float(line[2])
             num_trans = len(genes_w_multiple_ens_ids_for_one_g_id_locations[gname])
             #print(gname, chro, start, end)
            
            #finding the gene that this record is the best match to probably
             bestmatch=''
             bestmatchdist=0
             bestmatchdist2=0
             for i in range(num_trans):
                #print(i, "it starts", genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][1], "gene nane", genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][5])
                dist=abs(float(genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][1])-start)
                dist2=abs(float(genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][2])-end)
                #print(genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i])
                #print(i, dist)
                if i==0:
                    bestmatchdist=dist
                    bestmatchdist2=dist2
                    bestmatchg=genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][5]
                    bestmatcht=genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][6]
                    #print("starting", bestmatchdist)
                elif dist<bestmatchdist and dist2<bestmatchdist2:
                    bestmatchdist=dist
                    bestmatchdist2=dist2
                    bestmatchg=genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][5]
                    bestmatcht=genes_w_multiple_ens_ids_for_one_g_id_locations[gname][i][6]
                    #print("best fit", start, bestmatchdist)
                else:
                    #print("not best match", dist)
                    pass
             #print(line[3], genes_w_multiple_ens_ids_for_one_g_id_locations[gname])
            # print("FINAL PICK", bestmatcht, bestmatchdist)

             if line[0]=="MT":
                 name="MT"+gname
             else:
                 name=gname
            #print(name)
             line[3]=line[3].replace(name, bestmatchg)
             if bestmatchdist<50 and bestmatchdist2<600:
                s = "\t".join(line)
                newby.write(s + "\n")
             #print("NEWLINE", line)

        #if the gene id does not have an ensembl id that matches to it (i.e. novel genes)
        #then the line will stay the same 
        elif gname not in geneids or gname in genes_w_multiple_ens_ids_for_one_g_id:
            #print(safeline)
            c+=1
            newby.write(safeline + "\n")
            


print("lines that have an ensembl id associated with the gene id",a)
print("genes_w_multiple_ens_ids_for_one_g_id_locations+normal lines without gene ids", c)
print("genes_w_multiple_ens_ids_for_one_g_id_locations", len(genes_w_multiple_ens_ids_for_one_g_id))
#print("transcripts saved", transcripts)

#print(genes_w_multiple_ens_ids_for_one_g_id)