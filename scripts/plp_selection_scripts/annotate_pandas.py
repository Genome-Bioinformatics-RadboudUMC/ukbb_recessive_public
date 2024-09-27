from __future__ import division
import csv
import ast
import os
import io
import pandas as pd
import numpy as np
import gzip
import gc 

#needed files
input_file  = ".../450k/RAP_output_per_chr/files_for_annotation/chrnumber.all_parts_final_for_annotation_sorted.txt"
hcdiff= ".../450k/RAP_output_per_chr/files_for_annotation/chrnumber_no_genotypes.hcdiff_format_annotated.txt"
clinvar = "../450k/clinvar/clinvar_20230430.vcf.gz"
hgmd ="../450k/hgmd/hgmd_pro_2018.3_hg38.vcf"
vkgl ="../450k/vkgl/VKGL_public_consensus_apr2021_arranged.csv"
Intervar = "../450k/annovar/example/450k/myannonumber.hg38_multianno_for_annotation.txt.intervar"
Decipher = "../450k/decipher_data_230321/decipher-snvs-grch38-2021-03-21.txt"
mpc = "../450k/mpc_scores/mpc.hg38.vcf.vcf"
mpc_1929 = "../450k/mpc_scores/mpc_1929.hg38.vcf"
MOI = "../450k/MOI-Pred/MOI-Pred_hg38_allchr_recessive_only.txt"
ADAR_genes = ".../450k/regions/ADAR_genes_gencode-v34.txt"
AR_genes = ".../450k/regions/AR_genes_gencode-v34.txt"
HGMD_transcripts = "../450k/HGMD_transcripts_with_matching_ENS.txt"
regions_1929 = "../450k/regions/1929_genes/transcripts_exons_hg38_merged_10bp_with_chr.bed"
output_file = ".../450k/plp_selection/chrnumber_1929_genes_annotated.txt"


#load files into dataframes
annotated_file_header= ["chr","position","ref","alt","gene","region","synonymous","gnomad-E_AF","loftee","cDNA","protein","variant_type","cadd_score","mpc_score","MPI-pred_score","decipher","vkgl","hgmd-DM","clinvar","clinvar_stars","intervar","hets","homs"]

def read_vcf(path):
    try:
        with gzip.open(path, 'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t',
            usecols =['#CHROM','POS','REF','ALT','INFO']
        ).rename(columns={'#CHROM': 'CHROM'})
    except gzip.BadGzipFile:
        with open(path, 'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t',
            usecols =['#CHROM','POS','REF','ALT','INFO']
        ).rename(columns={'#CHROM': 'CHROM'})

df_input = pd.read_csv(input_file, sep="\t", names=["chr","position","ref","alt","hets","homs"])
df_inhouse = pd.read_csv(hcdiff, sep="\t", usecols=["Chromosome","Start position", "Reference", "Variant", "Gene name", "Gene component", "Synonymous", "Hgvsc", "Hgvsp", "Loftee", "CADD_PHRED", "gnomAD-E AF","Variant type","Protein Effect","All transcripts"])
df_clinvar = read_vcf(clinvar)
df_hgmd = read_vcf(hgmd)
df_vkgl = pd.read_csv(vkgl, sep=",",usecols=["chromosome","start","ref","alt","classification"])
df_intervar = pd.read_csv(Intervar, sep="\t", usecols=["#Chr","Start", "Ref", "Alt", "Ref.Gene", " InterVar: InterVar and Evidence "])
df_decipher = pd.read_csv(Decipher, sep="\t", skiprows=1, usecols =["chr","start", "ref_allele", "alt_allele", "pathogenicity"])
##os.system("module load bioinf/bedtools/2.25.0")
##os.system("bedtools intersect -a " + mpc + " -b " + regions_1929 + " > "+ mpc_1929)
df_mpc = pd.read_csv(mpc_1929, sep="\t", usecols=[0,1,3,4,7],names=["CHROM","POS","REF","ALT","INFO"])
df_moi = pd.read_csv(MOI, sep="\t", usecols=[0,1,2,3,4,7],names=["CHROM","POS","REF","ALT","GENE","MOI-Pred_score"])
df_AR = pd.read_csv(AR_genes, sep="\t", header=None)
df_ADAR = pd.read_csv(ADAR_genes, sep="\t", header=None)

##edit files_functions
##general functions
def add_chr(value):
    return "chr"+str(value)

def add_key(chromosome,position,ref,alt):
    if chromosome !=None:
        return chromosome+"-"+str(position)+"-"+ref+"-"+alt

def add_key_gene(chromosome,position,ref,alt,gene):
    if chromosome !=None:
        return chromosome+"-"+str(position)+"-"+ref+"-"+alt+"-"+gene
    
#positions editing
def edit_pos_input(row):
    position, ref, alt = row['position'], row['ref'], row['alt']
    if len(ref)==1 and len(alt)==1:
        return row
    elif len(ref) > len(alt):
        row['position'] = str(int(position)+1)
        row['ref'] = ref[1:]
        row['alt'] = "."        
        return row
    elif len(ref) < len(alt):
        row['ref'] = "."
        row['alt'] = alt[1:]
        return row

def edit_pos_vkgl(row):
    position, ref, alt = row['start'], row['ref'], row['alt']
    if len(ref)==1 and len(alt)==1:
        return row
    elif len(ref) > len(alt):
        row['position'] = str(int(position)+1)
        row['ref'] = ref[1:]
        row['alt'] = "."        
        return row
    elif len(ref) < len(alt):
        row['ref'] = "."
        row['alt'] = alt[1:]
        return row

def edit_pos_clinvar_hgmd(row):
    position, ref, alt = row['POS'], row['REF'], row['ALT']
    if len(ref)==1 and len(alt)==1:
        return row
    elif len(ref) > len(alt):
        row['POS'] = str(int(position)+1)
        row['REF'] = ref[1:]
        row['ALT'] = "."        
        return row
    elif len(ref) < len(alt):
        row['REF'] = "."
        row['ALT'] = alt[1:]
        return row

def edit_pos_inhouse(row):
    position, ref, alt = row['Start position'], row['Reference'], row['Variant']
    if str(alt)=="nan":
        row['Variant'] = "."        
        return row
    elif str(ref)=="nan":
        row['Reference'] = "."
        return row
    elif len(ref)==1 and len(alt)==1:
        return row

def edit_pos_intervar(row):
    ref, alt = row['Ref'], row['Alt']
    if len(ref)==1 and len(alt)==1 and ref!="-" and alt!="-":
        return row
    elif len(ref) > len(alt) or (len(ref)==1 and len(alt)==1 and alt=="-"):
        row['Alt'] = "."        
        return row
    elif len(ref) < len(alt) or (len(ref)==1 and len(alt)==1 and ref=="-"):
        row['Ref'] = "."
        return row
    else:
        return row

def edit_pos_decipher(row):
    position, ref, alt = row['start'], row['ref_allele'], row['alt_allele']
    if len(ref)==1 and len(alt)==1:
        return row
    elif len(ref) > len(alt):
        row['start'] = str(int(position)+1)
        row['ref_allele'] = ref[1:]
        row['alt_allele'] = "."           
        return row
    elif len(ref) < len(alt):
        row['ref_allele'] = "."
        row['alt_allele'] = alt[1:]
        return row

##inpout_df
#edit hets and homs as int
def edit_het_hom(value):
##    try:
        return int(value)
##    except ValueError:
##        return value

    
##inhouse_df
#edit apperance of cDNA and protein level
def edit_c_p(value):
    try:
        return value.split(":")[1]
    except IndexError:
        return value
  
#create a dictionary of transcripts per gene
def gene_transcript():
    with open (HGMD_transcripts,'rt') as transcripts:
        hgmd_trans= csv.reader(transcripts, delimiter='\t')

        d_names={"HJV":"HFE2","PJVK":"DFNB59"}
        d_ens={}

        for w in hgmd_trans:
            if w[0] not in d_ens.keys():
                if w[0] in d_names.keys():
                    d_ens[d_names[w[0]]]=ast.literal_eval(w[2])
                else:
                    d_ens[w[0]]=ast.literal_eval(w[2])
            else:
                if w[0] in d_names.keys():
                    d_ens[d_names[w[0]]].extend(ast.literal_eval(w[2]))
                else:
                    d_ens[w[0]].extend(ast.literal_eval(w[2]))
        return d_ens

#used as a part of edit_consequence() end counts the #num of lof/missense/other for each variant
def count_vars_type(x):
    lof=0
    missense=0
    other=0
    LOF=["stop_gained" ,"start_lost" , "frameshift_variant" ,"splice_donor_variant" , "splice_acceptor_variant"]
    Missense= ["missense_variant"]
    for j in x:
        if j in LOF:
            lof+=1
        elif j in Missense:
            missense+=1
        else:
            other+=1
    return(lof,missense,other)

#used as a part of edit_consequence() end determine the final variant outcome based on variant types counts
def determine(lof,missense,other,proteinEffect):
    if (lof>0 and missense==0 and other==0) or (lof==0 and missense>0 and other==0) or (lof==0 and missense==0 and other>0):
        return proteinEffect.split(";")[0]
    else:
        perlof= lof/(lof+missense+other)*100
        permis= missense/(lof+missense+other)*100
        perother= other/(lof+missense+other)*100
        if perlof>50:
            return proteinEffect.split(";")[0]
        elif permis>50:
            return "missense_variant"
        elif perother>50:
            return "other"
        else:
            return "no"

#choose variant outcome based on transcripts      
def edit_consequence(d_ens,geneName,proteinEffect,allTranscripts):
    tmpL=allTranscripts.split(" ")
    varList=[]
    rightV=[]

    for u in tmpL:
        if len(u)>1:
            ens=u.split("(")[0].split(".")[0]
            vType=u.split("(")[1][:-1]
            varList.append([ens,vType])

    if geneName in d_ens.keys():   
        for pair in varList:
            if pair[0] in d_ens[geneName]:
                rightV.append(pair[1])
                                
        if len(set(rightV))==1 and rightV!=[]:        
            var = rightV[0]
            return var
        elif rightV!=[]:
            
            lof,missense,other = count_vars_type(rightV)
            var = determine(lof,missense,other,proteinEffect)
            if var=="no":
                return rightV
            else:
                return var
        else:
            tmp=[]
            for pair in varList:
                tmp.append(pair[1])
            lof,missense,other = count_vars_type(tmp)
            var = determine(lof,missense,other,proteinEffect)
            if var=="no":
                for pair in varList:
                    rightV.append(pair[1])
                return rightV
            else:
                return var
    else:
        tmp=[]
        for pair in varList:
            tmp.append(pair[1])
        lof,missense,other = count_vars_type(tmp)
        var = determine(lof,missense,other,proteinEffect)
        if var=="no":
            for pair in varList:
                rightV.append(pair[1])
            return rightV
        else:
            return var

    
##clinvar
def edit_clinvar(row):
    value = row['INFO']
    if value != None:
        inf= value.split(';')
        if 'CLNSIGINCL' in inf:
            clin_clas="-"
            clin_rev="-"
        else:
            clin_clas="-"
            clin_rev="-"
            for i in inf:
                if "CLNSIG=" in i:
                    clin_clas= str(i.split('=')[1])
                elif "CLNREVSTAT=" in i:
                    tmp_rev= i.split('=')[1]
                    if "no_assertion" in tmp_rev:
                        clin_rev= 0
                    elif "single" in tmp_rev or "conflicting" in tmp_rev:
                        clin_rev= 1
                    elif "multiple" in tmp_rev:
                        clin_rev= 2
                    elif "expert" in tmp_rev:
                        clin_rev= 3
                    elif "practice" in tmp_rev:
                        clin_rev= 4
                    else:
                        clin_rev= "-"
        return pd.Series([clin_clas, clin_rev])


##hgmd
def edit_hgmd(value):
    for inf in value.split(";"):
        if "CLASS" in inf and inf.split("=")[1]=="DM":
            return "Y"
        else:
            return None


##intervar
def edit_intervar(value):
    tmp = value.split(": ")[1].split(" PVS1")[0]
    if "significance" in tmp:
        return "Uncertain_significance"
    elif "benign" in tmp:
        return "Likely_benign"
    elif "Benign" in tmp:
        return "Benign"
    elif "pathogenic" in tmp:
        return "Likely_pathogenic"
    elif "Pathogenic" in tmp:
        return "Pathogenic"

##mpc
def edit_mpc(row):
    value = row['INFO']
    if value != None:
        inf= value.split('|')
        return pd.Series([inf[4],inf[-1]]) 


#activate editing
df_input['true_key'] = df_input.apply(lambda x: add_key(x['chr'],x['position'],x['ref'],x['alt']), axis=1)
df_input = df_input.apply(edit_pos_input, axis=1)
df_input['hets'] = df_input['hets'].astype('Int64')
df_input['homs'] = df_input['homs'].astype('Int64')
df_input['key'] = df_input.apply(lambda x: add_key(x['chr'],x['position'],x['ref'],x['alt']), axis=1)
                                               
df_inhouse['Hgvsc'] = df_inhouse['Hgvsc'].astype(str)
df_inhouse['Hgvsp'] = df_inhouse['Hgvsp'].astype(str)
df_inhouse['Hgvsc'] = df_inhouse['Hgvsc'].apply(edit_c_p)
df_inhouse['Hgvsp'] = df_inhouse['Hgvsp'].apply(edit_c_p)
d_ens = gene_transcript()
df_inhouse['Protein Effect'] = df_inhouse.apply(lambda x: edit_consequence(d_ens, x['Gene name'],x['Protein Effect'],x['All transcripts']), axis=1)
df_inhouse = df_inhouse.apply(edit_pos_inhouse, axis=1)
df_inhouse.columns = df_inhouse.columns.str.replace('gnomAD-E AF', 'gnomAD-E_AF')
df_inhouse.columns = df_inhouse.columns.str.replace('Gene name', 'gene')
df_inhouse['gene'] = df_inhouse['gene'].astype(str)
df_inhouse.columns = df_inhouse.columns.str.replace('Gene component', 'region')
df_inhouse.columns = df_inhouse.columns.str.replace('Synonymous', 'synonymous')
df_inhouse.columns = df_inhouse.columns.str.replace('Variant type', 'variant_type')
df_inhouse.columns = df_inhouse.columns.str.replace('Loftee', 'loftee')
df_inhouse.columns = df_inhouse.columns.str.replace('CADD_PHRED', 'CADD_score')
df_inhouse['key'] = df_inhouse.apply(lambda x: add_key(x['Chromosome'],x['Start position'],x['Reference'],x['Variant']), axis=1)
df_inhouse['key_gene'] = df_inhouse.apply(lambda x: add_key_gene(x['Chromosome'],x['Start position'],x['Reference'],x['Variant'],x['gene']), axis=1)

df_clinvar[['clinvar','clinvar_stars']] = df_clinvar.apply(edit_clinvar,axis=1)
df_clinvar['CHROM'] = df_clinvar['CHROM'].apply(add_chr)
df_clinvar = df_clinvar.apply(edit_pos_clinvar_hgmd, axis=1)
df_clinvar.drop('INFO', axis=1, inplace=True)
df_clinvar['clinvar_stars'] = df_clinvar['clinvar_stars'].astype('Int64')
df_clinvar['key'] = df_clinvar.apply(lambda x: add_key(x['CHROM'],x['POS'],x['REF'],x['ALT']), axis=1)

df_hgmd['hgmd-DM'] = df_hgmd['INFO'].apply(edit_hgmd)
df_hgmd['CHROM'] = df_hgmd['CHROM'].apply(add_chr)
df_hgmd = df_hgmd.apply(edit_pos_clinvar_hgmd, axis=1)
df_hgmd['key'] = df_hgmd.apply(lambda x: add_key(x['CHROM'],x['POS'],x['REF'],x['ALT']), axis=1)

df_vkgl['chromosome'] = df_vkgl['chromosome'].apply(add_chr)
df_vkgl = df_vkgl.apply(edit_pos_vkgl, axis=1)
df_vkgl.columns = df_vkgl.columns.str.replace('classification', 'vkgl')
df_vkgl['start'] = df_vkgl['start'].astype('Int64')
df_vkgl['key'] = df_vkgl.apply(lambda x: add_key(x['chromosome'],x['start'],x['ref'],x['alt']), axis=1)

df_intervar.rename(columns={'#Chr': 'Chr',' InterVar: InterVar and Evidence ':'intervar','Ref.Gene':'gene'},inplace=True)
df_intervar['gene'] = df_intervar['gene'].astype(str)
df_intervar['Chr'] = df_intervar['Chr'].apply(add_chr)
df_intervar['intervar'] = df_intervar['intervar'].apply(edit_intervar)
df_intervar = df_intervar.apply(edit_pos_intervar, axis=1)
df_intervar['Start'] = df_intervar['Start'].astype('Int64')
df_intervar['key_gene'] = df_intervar.apply(lambda x: add_key_gene(x['Chr'],x['Start'],x['Ref'],x['Alt'],x['gene']), axis=1)

df_decipher['chr'] = df_decipher['chr'].apply(add_chr)
df_decipher = df_decipher.apply(edit_pos_decipher, axis=1)
df_decipher.columns = df_decipher.columns.str.replace('pathogenicity', 'decipher')
df_decipher['key'] = df_decipher.apply(lambda x: add_key(x['chr'],x['start'],x['ref_allele'],x['alt_allele']), axis=1)

df_mpc[['gene','mpc_score']] = df_mpc.apply(edit_mpc,axis=1)
df_mpc.drop('INFO', axis=1, inplace=True)
df_mpc['gene'] = df_mpc['gene'].astype(str)
df_mpc['key_gene'] = df_mpc.apply(lambda x: add_key_gene(x['CHROM'],x['POS'],x['REF'],x['ALT'],x['gene']), axis=1)

df_moi['CHROM'] = df_moi['CHROM'].apply(add_chr)
df_moi['key'] = df_moi.apply(lambda x: add_key(x['CHROM'],x['POS'],x['REF'],x['ALT']), axis=1)


#merging all database into one df
##print("input")
df_combined = df_input[['key','true_key','hets','homs']].drop_duplicates().merge(
    df_inhouse[['key','key_gene','variant_type','gene','region','Protein Effect', 
                'synonymous','Hgvsc','Hgvsp','loftee','gnomAD-E_AF', 'CADD_score']], how='left', on='key')

##print("clinvar")
df_combined = df_combined.merge(df_clinvar[['key','clinvar','clinvar_stars']].drop_duplicates(), how='left', on='key')

##print("vkgl")
df_combined = df_combined.merge(df_vkgl[['key','vkgl']].drop_duplicates(), how='left', on='key')

##print("hgmd")
df_combined = df_combined.merge(df_hgmd[['key','hgmd-DM']].drop_duplicates(), how='left', on='key')

##print("intervar")
df_combined = df_combined.merge(df_intervar[['key_gene','intervar']].drop_duplicates(), how='left', on='key_gene')

##print("decipher")
df_combined = df_combined.merge(df_decipher[['key','decipher']].drop_duplicates(), how='left', on='key')

##print("mpc")
df_combined = df_combined.merge(df_mpc[['key_gene','mpc_score']].drop_duplicates(), how='left', on='key_gene')

df_combined = df_combined.merge(df_moi[['key','MOI-Pred_score']].drop_duplicates(), how='left', on='key')




#edit apperance of combined df
df_combined[['chr','position','ref','alt']]=df_combined.true_key.str.split('-',expand=True)
df_combined.drop('key', axis=1, inplace=True)
df_combined.drop('true_key', axis=1, inplace=True)
df_combined.drop('key_gene', axis=1, inplace=True)
cols=['chr','position','ref','alt','gene','region','synonymous','Hgvsc','Hgvsp','variant_type','Protein Effect','gnomAD-E_AF','loftee','mpc_score','CADD_score', 'MOI-Pred_score','decipher','vkgl','hgmd-DM','clinvar','clinvar_stars','intervar','hets','homs']
df_combined= df_combined[cols]
df_combined.to_csv(output_file, index=False, sep ='\t')


