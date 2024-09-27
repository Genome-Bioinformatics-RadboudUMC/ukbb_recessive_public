from __future__ import division
import csv
import random
from datetime import datetime
from collections import Counter
import ast
import math
import statistics

plps_file=".../450k/plp_selection/basic/new_gene_names/new_freq/all_chr_total_presumable_plps_HFE_final_sorted.txt"
plps_file_geno=".../450k/RAP_output_per_chr/filtered_plps/basic/new_gene_names/new_freq/all_chr.all_parts_final.csv"
genes_file=".../450k/regions/gene-panel-gencode-v34.txt"
mild_vars=".../450k/plp_selection/basic/new_gene_names/new_freq/all_chr_total_presumable_plp_mild_var.txt"
CR_scores=".../450k/CR/new_gene_names/new_freq/CR_scores.txt"
output_folder=".../450k/CR/new_gene_names/new_freq/"

#examples:
#per=20/10/75
#analysis_type= genes_per_panel/pop
#variants_source= no_decipher/no_decipher_clinvar_rare_LOF
#change function in the activation function=make_CR_for_given_samples/make_CR_for_given_genes

per=75
analysis_type= "genes_per_panel"
variants_source= "all_vars"

file_name= output_folder+"CR_"+str(per)+"_per_"+analysis_type+"_1000_simulations_"+variants_source

#samples_num=40621
#samples_num=81242
samples_num=406213

def make_CR_for_given_samples(infile):
    '''
    This functions randomly choose XX samples and calculate CRs for this group of samples
    '''
    with open (infile, 'rt') as plpfile,open (genes_file,'rt') as genesfile,open (mild_vars,'rt') as mildvars:
        reader= csv.reader(plpfile, delimiter="\t")
        genes= csv.reader(genesfile, delimiter=",")
        mild= csv.reader(mildvars, delimiter="\t")

        samples_list=[]
        d_gene_freq={}
        d_arc={}
        mvars=[]
        d_CR={}

        for g in genes:
            d_gene_freq[g[0]]=[0,0,0,0]             #het-severe,hom-severe,het-mild,hom-mild
            d_arc[g[0]]=[0,0,g[1]]
            d_CR[g[1]]=[0,0,0]

        for v in mild:
            mvars.append(v[0]+"_"+v[1]+"_"+v[2]+"_"+v[3])

        for i in range(1,samples_num+1):
            random.seed(datetime.now())
            rand=random.randint(1,406213)
            samples_list.append(rand)

        for row in reader:
            tmpl=[]
            gene=row[4]
            for i in samples_list:
                if row[i+19]=="0/1" or row[i+19]=="1/0":
                    tmpl.append(1)
                elif row[i+19]=="1/1":
                    tmpl.append(2)

            if 1 in tmpl or 2 in tmpl:
                c=Counter(tmpl)
                if row[0]+"_"+row[1]+"_"+row[2]+"_"+row[3] in mvars:                    
                    d_gene_freq[gene][2]+=c[1]
                    d_gene_freq[gene][3]+=c[2]
                else:
                    d_gene_freq[gene][0]+=c[1]
                    d_gene_freq[gene][1]+=c[2]

        for k,v in d_gene_freq.items():
            sAC=v[0]+(2*v[1])
            mAC=v[2]+(2*v[3])
            sQ= sAC/(samples_num*2)
            mQ= mAC/(samples_num*2)
            sP=1-sQ
            mP=1-mQ
            sPQ2=2*sQ*sP
            mPQ2=2*mQ*mP
            d_arc[k][0]= (sPQ2*sPQ2 + 2*sPQ2*mPQ2)* samples_num
            d_arc[k][1]= sPQ2*(1/8)* samples_num
        
        for k,v in d_arc.items():
            d_CR[v[2]][0]+=v[0]
            d_CR[v[2]][1]+=v[1]

        for k,v in d_CR.items():
            if v[0]==0:
                d_CR[k][2]=0
            else:            
                d_CR[k][2]=v[1]/v[0]


        return d_CR

def make_CR_for_given_genes(infile):
    '''
    This functions randomly choose XX of genes per panel and and calculate CRs for this group of genes
    '''
    with open (infile, 'rU') as plpfile,open (genes_file,'rt') as genesfile,open (mild_vars,'rt') as mildvars:
        reader= csv.reader(plpfile, delimiter="\t")
        genes= csv.reader(genesfile, delimiter=",")
        mild= csv.reader(mildvars, delimiter="\t")

        samples_list=[]
        d_gene_freq={}
        d_arc={}
        panel_genes={}
        mvars=[]
        d_CR={}

        for g in genes:
            d_gene_freq[g[0]]=[0,0,0,0]             #het-severe,hom-severe,het-mild,hom-mild
            d_arc[g[0]]=[0,0,g[1]]
            d_CR[g[1]]=[0,0,0]
            if g[1] not in panel_genes.keys():
                panel_genes[g[1]]=[g[0]]
            else:
                panel_genes[g[1]].append(g[0])

        for v in mild:
            mvars.append(v[0]+"_"+v[1]+"_"+v[2]+"_"+v[3])

        for k,v in panel_genes.items():
            for i in range(int(math.ceil(0.75*len(v)))):
                random.seed(datetime.now())
                rand=random.randint(1,len(v))
                samples_list.append(v[rand-1])

        for row in reader:
            tmpl=[]
            gene=row[4]
            if gene in samples_list:
                if row[0]+"_"+row[1]+"_"+row[2]+"_"+row[3] in mvars:                    
                    d_gene_freq[gene][2]+=int(row[22])
                    d_gene_freq[gene][3]+=int(row[23])
                else:
                    d_gene_freq[gene][0]+=int(row[22])
                    d_gene_freq[gene][1]+=int(row[23])

        for k,v in d_gene_freq.items():
            sAC=v[0]+(2*v[1])
            mAC=v[2]+(2*v[3])
            sQ= sAC/(samples_num*2)
            mQ= mAC/(samples_num*2)
            sP=1-sQ
            mP=1-mQ
            sPQ2=2*sQ*sP
            mPQ2=2*mQ*mP
            d_arc[k][0]= (sPQ2*sPQ2 + 2*sPQ2*mPQ2)* samples_num
            d_arc[k][1]= sPQ2*(1/8)* samples_num
        
        for k,v in d_arc.items():
            d_CR[v[2]][0]+=v[0]
            d_CR[v[2]][1]+=v[1]

        for k,v in d_CR.items():
            if v[0]==0:
                d_CR[k][2]=0
            else:            
                d_CR[k][2]=v[1]/v[0]


        return d_CR



def activate_permutations():
    '''
    This function activate 1000 permutations for the given function (edit function before use)
    '''
    with open (genes_file,'rt') as genesfile, open (file_name+".txt", 'w')as final:
        genes= csv.reader(genesfile, delimiter=",")
        writer=csv.writer(final, delimiter="\t")

        d_CR_final={}

        for g in genes:
            d_CR_final[g[1]]=[]   

        #*******CHANGE_FUNCTION********
        for i in range(1000):
            d=make_CR_for_given_samples(plps_file_geno)
            for k,v in d.items():
                d_CR_final[k].append(v[2])
            print (i)

        for k,v in d_CR_final.items():
            writer.writerow([k,v])
                

activate_permutations()


def analyze():
    with open (file_name+".txt", 'rt') as analysis, open (CR_scores,'rt') as CRscores, open (file_name+"_analysis.txt",'w') as final:
        reader= csv.reader(analysis,delimiter="\t")
        scores= csv.reader(CRscores,delimiter="\t")
        writer= csv.writer(final,delimiter="\t") 

        CR={}

        for s in scores:
            CR[s[0]]=s[1]
        
        writer.writerow(["Panel","True CR", "Simulations- Min","Simulations- Max","Simulations- Mean","SD","Coefficient of variation"])    
        for r in reader:
            v= ast.literal_eval(r[1])
            if r[0]!="No_panel" and r[0]!="Tumor" and r[0]!="Cardiovascular":
                writer.writerow([r[0],CR[r[0]],min(v),max(v),sum(v)/len(v),statistics.stdev(v),statistics.stdev(v)/statistics.mean(v)])
            

analyze()

















            
