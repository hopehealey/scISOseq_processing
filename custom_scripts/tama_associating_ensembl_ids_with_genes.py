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
        if gname not in genes_w_multiple_ens_ids_for_one_g_id and gname in geneids:
            a+=1
            #pulling info [chro, start, end, gname, tname,ensgname, enstname]
            myinfo=geneids[gname]
            #print("yes ma'am", gname, myinfo[5])
            #print(line)
            #print(re.sub('(gname)', 'myinfo[5]', safeline))
            name=gname
            #print(name)
            line[3]=line[3].replace(name, myinfo[5])
            #print(line)
            s = "\t".join(line)
            newby.write(s + "\n")
            #print("done")

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


        elif gname not in geneids or gname in genes_w_multiple_ens_ids_for_one_g_id:
            #print(safeline)
            c+=1
            newby.write(safeline + "\n")
            


print("safe lines",a)
print("problems lines+normal lines without gene ids", c)
print("problem children", len(genes_w_multiple_ens_ids_for_one_g_id))
#print("transcripts saved", transcripts)

#print(genes_w_multiple_ens_ids_for_one_g_id)