from __future__ import division
import csv
import os
import pandas as pd
import numpy as np

#enter annotated file name
annotated_file = ".../450k/plp_selection/chrnumber_1929_genes_annotated.txt"
output_file = ".../450k/plp_selection/basic/new_gene_names/new_freq/chrnumber_total_presumable_plps.txt"
output_file_leftovers = ".../450k/plp_selection/basic/new_gene_names/new_freq/chrnumber_total_presumable_plps_frequency_leftovers.txt"

#enter num of samples
samples_num= 406213

#needed files
ADAR_genes = ".../450k/regions/ADAR_genes_gencode-v34.txt"
AR_genes = ".../450k/regions/AR_genes_gencode-v34.txt"

#header
annotated_file_header=['chr','position','ref','alt','gene','region','synonymous','Hgvsc','Hgvsp','variant_type','Protein Effect','gnomAD-E_AF','loftee','mpc_score','CADD_score','MOI-Pred_score','decipher','vkgl','hgmd-DM','clinvar','clinvar_stars','intervar','hets','homs']


#load files
df_AR = pd.read_csv(AR_genes, sep="\t", header=None)
df_ADAR = pd.read_csv(ADAR_genes, sep="\t", header=None)
df_annotated = pd.read_csv(annotated_file, sep="\t", engine='python')

#edit annotated file- convert strings into int/float in specific columns
def change_dtype(value):
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value

df_annotated.loc[:,'clinvar_stars']=df_annotated['clinvar_stars'].apply(change_dtype)
df_annotated['gnomAD-E_AF'] = df_annotated['gnomAD-E_AF'].astype(float)
df_annotated['position'] = df_annotated['position'].astype('Int64')
df_annotated['ref'] = df_annotated['ref'].astype(str)
df_annotated['alt'] = df_annotated['alt'].astype(str)
df_annotated['hets'] = df_annotated['hets'].astype('Int64')
df_annotated['homs'] = df_annotated['homs'].astype('Int64')
df_annotated['synonymous'] = df_annotated['synonymous'] = df_annotated['synonymous']=='TRUE'
df_annotated['CADD_score'] = df_annotated['CADD_score'].astype(float)
df_annotated['mpc_score'] = df_annotated['mpc_score'].astype(float)
df_annotated['MOI-Pred_score'] = df_annotated['MOI-Pred_score'].astype(float)
df_annotated['vkgl'] = df_annotated['vkgl'].astype(str)
df_annotated['hgmd-DM'] = df_annotated['hgmd-DM'].astype(str)
df_annotated['intervar'] = df_annotated['intervar'].astype(str)
df_annotated['variant_type'] = df_annotated['variant_type'].astype(str)
df_annotated['region'] = df_annotated['region'].astype(str)
df_annotated['loftee'] = df_annotated['loftee'].astype(str)


#conditions AR/ADAR subs/indels
cond_AR = (df_annotated['gene'].isin(df_AR[0]))
cond_ADAR = (df_annotated['gene'].isin(df_ADAR[0]))
cond_subs = (((df_annotated['ref'].str.len()==1) & (df_annotated['ref']!=".")) & ((df_annotated['alt'].str.len()==1) & (df_annotated['alt']!=".")))
cond_indels = (~(cond_subs))  

#frequency condition (hets<5%, homs<1%)
cond_freq_cohort = ((df_annotated['hets']<samples_num*0.05) & (df_annotated['homs']<samples_num*0.01))
cond_not_freq_cohort = (~(cond_freq_cohort))

#conditions tier 1- "safe tier"
cond_clinvar = ((df_annotated['clinvar'].str.contains("athogenic")) & (df_annotated[~(df_annotated['clinvar_stars']==" ")]['clinvar_stars']>=2))
cond_vkgl = (df_annotated['vkgl']=="LP")          

tier_1_cond = (cond_clinvar | cond_vkgl)

#conditions tier 2- "rare LOF"
cond_region = ((df_annotated['region'].str.contains("CANONICAL")) | (df_annotated['region']=="EXON_REGION"))
cond_synonymous = (~df_annotated['synonymous'])
variant_types = ["stop_gained","frameshift","splice_donor","splice_acceptor"]
cond_variant_type = (df_annotated['Protein Effect'].str.contains('|'.join(variant_types),na=False))

cond_lof = (cond_region & cond_synonymous & cond_variant_type)
cond_freq_var = (df_annotated['gnomAD-E_AF']<1)

cond_clinvar_not = ~(df_annotated['clinvar'].str.contains("enign") & (df_annotated[~(df_annotated['clinvar_stars']==" ")]['clinvar_stars']>=2))
cond_vkgl_not = ~(df_annotated['vkgl']=="LB")

tier_2_cond = (cond_lof & cond_freq_var & cond_clinvar_not & cond_vkgl_not)

#addition for indels
cond_indel_length = ~((df_annotated['ref'].str.len()>10) | (df_annotated['alt'].str.len()>10))
cond_indel_loftee = (df_annotated['loftee']=="HC")

cond_indels_extras = (cond_indel_length & cond_indel_loftee)

def remove_adjacent_indels(in_df):
    in_df.to_csv('../450k/tmp.vcf',sep="\t", index=False)
    os.system("echo end >> ../450k/tmp.vcf")
    with open ('../450k/tmp.vcf', "rt") as inPut, open ('../450k/tmp.txt', "w") as final:
        reader= csv.reader(inPut, delimiter="\t")
        writer= csv.writer(final, delimiter="\t")

        next(reader)
        nl1=next(reader)
        nl2=next(reader)

        while nl1[0]!="end" and nl2[0]!="end":
            if int(float(nl2[1]))-int(float(nl1[1]))>10:
                writer.writerow(nl1)
                nl1=nl2
                nl2=next(reader)
            else:
                nl1=next(reader)
                if nl1[0]!="end":
                    nl2=next(reader)
        if nl2[0]=="end":
            writer.writerow(nl1)
    return pd.read_csv('../450k/tmp.txt', sep="\t", names=annotated_file_header)

#conditions tier 3- "missed missense"
cond_missense = (df_annotated['Protein Effect'].str.contains("missense_variant"))
cond_clinvar_weak = ((df_annotated['clinvar'].str.contains("Pathogenic")) | (df_annotated['clinvar']=="Likely_pathogenic"))
cond_hgmd = (df_annotated['hgmd-DM']=="Y")
cond_intervar = (df_annotated['intervar'].str.contains("athogenic"))

tier_3_cond = ((cond_clinvar_not & cond_vkgl_not & cond_missense) & ((cond_clinvar_weak & cond_hgmd) | (cond_hgmd & cond_intervar) | (cond_intervar & cond_clinvar_weak)))
tier_3_cond_CADD_MPC = ((cond_clinvar_not & cond_vkgl_not & cond_missense) & (df_annotated['CADD_score']>25) & (df_annotated['mpc_score']>2))
tier_3_cond_CADD_MOI = ((cond_clinvar_not & cond_vkgl_not & cond_missense) & (df_annotated['CADD_score']>25) & (df_annotated['MOI-Pred_score']>0.9))


#combined conditions- choose which tier 3 is relevant
cond_AR_subs_plps  = (cond_AR & cond_subs & (tier_1_cond | tier_2_cond | tier_3_cond ))
#cond_AR_subs_plps  = (cond_AR & cond_subs & (tier_1_cond | tier_2_cond | (tier_3_cond | tier_3_cond_CADD_MPC)))
#cond_AR_subs_plps  = (cond_AR & cond_subs & (tier_1_cond | tier_2_cond | (tier_3_cond | tier_3_cond_CADD_MOI)))

cond_ADAR_subs_plps = (cond_ADAR & cond_subs & (tier_1_cond | tier_2_cond | tier_3_cond) & cond_lof)     #only LOF variants for ADAR genes
#cond_ADAR_subs_plps = (cond_ADAR & cond_subs & (tier_1_cond | tier_2_cond | (tier_3_cond | tier_3_cond_CADD_MPC)) & cond_lof)     #only LOF variants for ADAR genes
#cond_ADAR_subs_plps = (cond_ADAR & cond_subs & (tier_1_cond | tier_2_cond | (tier_3_cond | tier_3_cond_CADD_MOI)) & cond_lof)     #only LOF variants for ADAR genes

cond_AR_indels  = (cond_AR & cond_indels)
cond_ADAR_indels = (cond_ADAR & cond_indels & cond_lof)
cond_AR_indels_tier_1  = (cond_AR_indels & tier_1_cond )
cond_ADAR_indels_tier_1 = (cond_ADAR_indels & tier_1_cond)

#extra steps needed for indels tier2
df_AR_indels_tier2_prep = df_annotated[(cond_AR_indels & tier_2_cond & cond_indels_extras & cond_freq_cohort) & (~(cond_AR_indels_tier_1))]
df_ADAR_indels_tier2_prep = df_annotated[(cond_ADAR_indels & tier_2_cond & cond_indels_extras & cond_freq_cohort) & (~(cond_ADAR_indels_tier_1))]


#select plps

#subs:
df_AR_subs_plps = df_annotated[cond_AR_subs_plps & cond_freq_cohort]
df_ADAR_subs_plps = df_annotated[cond_ADAR_subs_plps & cond_freq_cohort]

#indels:
df_AR_indels_tier_1 = df_annotated[cond_AR_indels_tier_1 & cond_freq_cohort]
df_ADAR_indels_tier_1 = df_annotated[cond_ADAR_indels_tier_1 & cond_freq_cohort]
df_AR_indels_tier_2 = df_AR_indels_tier2_prep.pipe(remove_adjacent_indels)
df_ADAR_indels_tier_2 = df_ADAR_indels_tier2_prep.pipe(remove_adjacent_indels)

df_AR_indels_plps = pd.concat([df_AR_indels_tier_1, df_AR_indels_tier_2])
df_ADAR_indels_plps = pd.concat([df_ADAR_indels_tier_1, df_ADAR_indels_tier_2])

#all plps:
df_all_plps = pd.concat([df_AR_subs_plps, df_ADAR_subs_plps, df_AR_indels_plps, df_ADAR_indels_plps])
df_all_plps_sorted = df_all_plps.sort_values('position')
df_all_plps_sorted.to_csv(output_file ,sep="\t", index=False)


#dealing with frequency drop-outs
#subs:
df_AR_subs_plps_leftovers = df_annotated[cond_AR_subs_plps & cond_not_freq_cohort]
df_ADAR_subs_plps_leftovers = df_annotated[cond_ADAR_subs_plps & cond_not_freq_cohort]

#indels:
df_AR_indels_tier_1_leftovers = df_annotated[cond_AR_indels_tier_1 & cond_not_freq_cohort]
df_ADAR_indels_tier_1_leftovers = df_annotated[cond_ADAR_indels_tier_1 & cond_not_freq_cohort]
df_AR_indels_tier_2_leftovers = df_annotated[(cond_AR_indels & tier_2_cond & cond_indels_extras & cond_not_freq_cohort) & (~(cond_AR_indels_tier_1))]
df_ADAR_indels_tier_2_leftovers = df_annotated[(cond_ADAR_indels & tier_2_cond & cond_indels_extras & cond_not_freq_cohort) & (~(cond_ADAR_indels_tier_1))]

df_AR_indels_plps_leftovers = pd.concat([df_AR_indels_tier_1_leftovers, df_AR_indels_tier_2_leftovers])
df_ADAR_indels_plps_leftovers = pd.concat([df_ADAR_indels_tier_1_leftovers, df_ADAR_indels_tier_2_leftovers])

#all plps:
df_all_plps_leftovers = pd.concat([df_AR_subs_plps_leftovers, df_ADAR_subs_plps_leftovers, df_AR_indels_plps_leftovers, df_ADAR_indels_plps_leftovers])
df_all_plps_sorted_leftovers = df_all_plps_leftovers.sort_values('position')
df_all_plps_sorted_leftovers.to_csv(output_file_leftovers ,sep="\t", index=False)




