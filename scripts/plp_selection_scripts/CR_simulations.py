from __future__ import division
import csv
import random
from datetime import datetime
import ast

genes_length=".../450k/regions/1929_genes_average_max_length.txt"
panels_list = ".../450k/regions/gene-panel-gencode-v34.txt"
panels_length=".../450k/regions/panel_length_average_based.txt"

consanguinity = ".../450k/CR/new_gene_names/new_freq/ARC_per_gene_consanguinity_no_sharedED.txt"
new_simu=".../450k/CR/new_gene_names/new_freq/CR_simulations_no_sharedED"
consan_ratio=".../450k/CR/new_gene_names/new_freq/CR_scores_no_sharedED.txt"

#minimum simulations analysis:
consan_ratio_min="../450k/CR_new/CR_scores_min_20per.txt"
new_simu_min="../450k/CR_new/CR_20_per_pop_1000_simulations_all_vars"



def simulations():
    with open (panels_list, "rt") as original, open (consanguinity, "rt") as consan,open (genes_length, "rt") as length,open (panels_length, "rt") as plength, open (new_simu+".txt", "w") as final:
        panels=csv.reader(original, delimiter=",")
        con=csv.reader(consan,delimiter=",")
        leng=csv.reader(length,delimiter="\t")
        pleng=csv.reader(plength,delimiter="\t")
        writer=csv.writer(final,delimiter="\t")

        next(con)
        next(leng)

        d_len={}
        d_genes_num={}
        d_con = {}
        genes_list = {}
                
        for g in pleng:
            d_genes_num[g[0].strip()]=float(g[1])

        for r in con:
            d_con[r[0]] = [r[1],r[2]]

        for r in leng:
            d_len[r[0]] = float(r[2])

        n=0        
        for r in panels:
            n+=1
            genes_list[n]=r[0]

        for p in d_genes_num.keys():
            row_p=[p]
            for i in range(10000):
                print(i)
                tmp_len_p = [0]
                chosen_genes=[]
                while len(chosen_genes)<1929:
                    random.seed(datetime.now())
                    rand=random.randint(1,1929)                    
                    if rand not in chosen_genes:
                        chosen_genes.append(rand)
                        if tmp_len_p[0] <d_genes_num[p]:
                            tmp_len_p[0]+=d_len[genes_list[rand]]
                            tmp_len_p.append(genes_list[rand])                                                
                            
                c_non=0
                c_con=0
                for gene in tmp_len_p[1:]:                                    
                    c_non+= float(d_con[gene][0])
                    c_con+= float(d_con[gene][1])
                if c_non!=0:
                    tmp_len_p.append(c_con/c_non)
                else:
                    tmp_len_p.append(0)
               
                row_p.extend([tmp_len_p[-1]])
                              
            writer.writerow(row_p)
                
##simulations()     

def analysis(new_simu):
    with open (new_simu+".txt", "rt") as infile, open (consan_ratio,"rt") as CFs, open (new_simu+"_analysis.txt","w") as final:
        reader=csv.reader(infile,delimiter="\t")
        CF=csv.reader(CFs, delimiter="\t")
        writer=csv.writer(final,delimiter="\t")

        d_CF={}
        dC={}
##        next(CF)

        for r in CF:
            d_CF[r[0].strip()]=float(r[1])
            dC[r[0].strip()]=0
                  
        for r in reader:
            dC[r[0]]=[0,0]
            for f in r[1:]:
                if float(f)>=d_CF[r[0]]:
                    dC[r[0]][0]+=1
                if float(f)<d_CF[r[0]]:
                    dC[r[0]][1]+=1                               
        
        for k,v in dC.items():
            if type(v)!=int:
                writer.writerow([k,v[0], v[1]]) 
                                     
analysis(new_simu)

def analysis_min(new_simu_min):
    with open (new_simu+".txt", "rt") as infile, open (consan_ratio_min,"rt") as CFs, open (new_simu_min+"_analysis.txt","w") as final:
        reader=csv.reader(infile,delimiter="\t")
        CF=csv.reader(CFs, delimiter="\t")
        writer=csv.writer(final,delimiter="\t")

        d_CF={}
        dC={}
##        next(CF)

        for r in CF:
            d_CF[r[0].strip()]=float(r[1])
            dC[r[0].strip()]=0
                  
        for r in reader:
            dC[r[0]]=[0,0]
            for f in r[1:]:
                if float(f)>=d_CF[r[0]]:
                    dC[r[0]][0]+=1
                if float(f)<d_CF[r[0]]:
                    dC[r[0]][1]+=1                               
        
        for k,v in dC.items():
            if type(v)!=int:
                writer.writerow([k,v[0], v[1]]) 
                                     

#analysis_min(new_simu_min)      




new_simu_ID_vs_all=".../450k/CR/new_gene_names/CR_simulations_ID_vs_all"
consan_ratio_ID_vs_all=".../450k/CR/new_gene_names/CR_scores_ID_vs_all.txt"

def simulations_ID_vs_all():
    with open (panels_list, "rt") as original, open (consanguinity, "rt") as consan,open (genes_length, "rt") as length,open (panels_length, "rt") as plength, open (new_simu_ID_vs_all+".txt", "w") as final:
        panels=csv.reader(original, delimiter=",")
        con=csv.reader(consan,delimiter=",")
        leng=csv.reader(length,delimiter="\t")
        pleng=csv.reader(plength,delimiter="\t")
        writer=csv.writer(final,delimiter="\t")

        next(con)
        next(leng)

        d_len={}
        d_genes_num={}
        d_con = {}
        genes_list = {}

        d_genes_num["all_others"]=0
                
        for g in pleng:
            if g[0]=="ID-total":
                d_genes_num[g[0].strip()]=float(g[1])
            else:
                d_genes_num["all_others"]+=float(g[1])

        for r in con:
            d_con[r[0]] = [r[1],r[2]]

        for r in leng:
            d_len[r[0]] = float(r[2])

        n=0        
        for r in panels:
            n+=1
            genes_list[n]=r[0]

        for p in d_genes_num.keys():
            row_p=[p]
            for i in range(10000):
                print(i)
                tmp_len_p = [0]
                chosen_genes=[]
                while len(chosen_genes)<1929:
                    random.seed(datetime.now())
                    rand=random.randint(1,1929)                    
                    if rand not in chosen_genes:
                        chosen_genes.append(rand)
                        if tmp_len_p[0] <d_genes_num[p]:
                            tmp_len_p[0]+=d_len[genes_list[rand]]
                            tmp_len_p.append(genes_list[rand])                                                
                            
                c_non=0
                c_con=0
                for gene in tmp_len_p[1:]:                                    
                    c_non+= float(d_con[gene][0])
                    c_con+= float(d_con[gene][1])
                if c_non!=0:
                    tmp_len_p.append(c_con/c_non)
                else:
                    tmp_len_p.append(0)
               
                row_p.extend([tmp_len_p[-1]])
                              
            writer.writerow(row_p)
                
#simulations_ID_vs_all()     

def analysis_ID_vs_all(new_simu_ID_vs_all):
    with open (new_simu_ID_vs_all+".txt", "rt") as infile, open (consan_ratio_ID_vs_all,"rt") as CFs, open (new_simu_ID_vs_all+"_analysis.txt","w") as final:
        reader=csv.reader(infile,delimiter="\t")
        CF=csv.reader(CFs, delimiter="\t")
        writer=csv.writer(final,delimiter="\t")

        d_CF={}
        dC={}
##        next(CF)

        for r in CF:
            d_CF[r[0].strip()]=float(r[1])
            dC[r[0].strip()]=0
                  
        for r in reader:
            dC[r[0]]=[0,0]
            for f in r[1:]:
                if float(f)>=d_CF[r[0]]:
                    dC[r[0]][0]+=1
                if float(f)<d_CF[r[0]]:
                    dC[r[0]][1]+=1                               
        
        for k,v in dC.items():
            if type(v)!=int:
                writer.writerow([k,v[0], v[1]]) 
                                     
#analysis_ID_vs_all(new_simu_ID_vs_all)

