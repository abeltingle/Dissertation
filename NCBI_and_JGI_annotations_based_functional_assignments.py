#!/usr/bin/env python
# coding: utf-8

# In[43]:


###################
# EC_based_functional_assignments.py
# Copyright 2023, Abel Ingle and Daniel Noguera
# Revised: January 26, 2023
###################

# Import packages

import pandas as pd
from pathlib import Path
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
import re
# Import packages

import subprocess
from subprocess import run, Popen, PIPE

from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
#%matplotlib inline

ids=pd.read_table("codefiles/MAG_names_accessions_ids.txt", encoding='latin-1').fillna('')
ids=ids.set_index("Strain",drop=False)

# The paths of these directories should be consistent with your filenaming scheme

reactions_table=pd.read_table("codefiles/metabolic_reactions.txt", encoding='latin-1').fillna('')

# Create a dictionary of EC numbers associated with reach metabolic reaction

ecnumbers={}
for row, index in reactions_table.iterrows():
    if index["Enzyme comission number"]=='':
        pass
    else:
        ecnumbers.update({index["Reaction ID"]:index["Enzyme comission number"]})

# Create a dataframe to be populated with 0's and 1's associated with each reaction-genome combination

ec=pd.DataFrame(ecnumbers.keys(), index=ecnumbers.keys())
ec=ec.rename(columns=({0:"BiGG Models Name"}))

# Introduce file paths of each JGI annotation


JGI_annotation_directory = 'codefiles/JGI_annotations/'
JGI_annotation_files=Path(JGI_annotation_directory).glob('*.gff')

for jgi_anno in JGI_annotation_files:

    genome = (ids[ids['JGI_bin_ID'] == jgi_anno.stem].index.tolist())[0]

    # If the genome does not yet have a dataframe column, add one
    
    if genome not in ec.columns:
        ec.insert(len(ec.columns),column=genome,value=0)
    
    # Create dataframe of annotation to create list of all E.C. numbers
    
    JGI_gff_annotation = gffpd.read_gff3(jgi_anno)
    gff_annos=[]
    for line in JGI_gff_annotation.df["attributes"]:
        if "ec_number" in line:
            ec_nos=line.split("ec_number=EC:")[1].split(";")[0].split(",EC:")
            for ec_no in ec_nos:
                gff_annos.append(ec_no)
    
    final_list_ec_nos_per_genome=set(gff_annos)
    
    # Check if EC number for each reaction is in the annotation
    
    for key, value in ecnumbers.items():
        if "or" in value:
            ecn=value.split(" or ")
            enzyme_comission_threshold = any(e in final_list_ec_nos_per_genome for e in ecn)
            
        elif "and" in value:
            ecn=value.split(" and ")
            enzyme_comission_threshold = all(e in final_list_ec_nos_per_genome for e in ecn)
        
        elif value in final_list_ec_nos_per_genome:
            enzyme_comission_threshold = True
            
        else:
            enzyme_comission_threshold = False
                
                
        if enzyme_comission_threshold is True:
            ec.at[key,genome]=1
            
            
        else:
            pass
        print(enzyme_comission_threshold, genome,key,value)
        
id_to_reaction_dic={}
for row, index in reactions_table.iterrows():
    id_to_reaction_dic.update({index["Reaction ID"]:index["Reaction Name"]})
    
for key in ecnumbers.keys():
    ec.at[key,"BiGG Models Name"]=id_to_reaction_dic[key]

ec.set_index("BiGG Models Name",drop=False)
    
ec.to_csv("output/JGI_annotation_based_functional_assignment.txt",index=None,sep="\t",mode="w")

ec=pd.DataFrame(ecnumbers.keys(), index=ecnumbers.keys())
ec=ec.rename(columns=({0:"BiGG Models Name"}))

# Introduce file paths of each NCBI annotation

NCBI_annotation_directory='codefiles/NCBI_annotations/'
NCBI_annotation_files=Path(NCBI_annotation_directory).glob('*.gbff')

for ncbi_anno in NCBI_annotation_files:
    
    NCBI_record_IDs=list()
    for record in SeqIO.parse(ncbi_anno, 'genbank'):
        NCBI_record_IDs.append(record.id)
    ncbi_accession=NCBI_record_IDs[0][0:6]
    
    genome = (ids[ids['NCBI_genome_accession'] == ncbi_accession].index.tolist())[0]
    
    # If the genome does not yet have a dataframe column, add one
    
    if genome not in ec.columns:
        ec.insert(len(ec.columns),column=genome,value=0)
    
    # Create dataframe of annotation to create list of all E.C. numbers
    
    gbff_annos = []
    for record in SeqIO.parse(ncbi_anno, 'genbank'):
        record_ids=list()
        for f in record.features:
            if f.type == "CDS":
                if "EC_number" in f.qualifiers:
                    for number in f.qualifiers["EC_number"]:
                        gbff_annos.append(number)
                        
    final_list_ec_nos_per_genome=set(gbff_annos)

    
    # Check if EC number for each reaction is in the annotation
    
    for key, value in ecnumbers.items():
        if "or" in value:
            ecn=value.split(" or ")
            enzyme_comission_threshold = any(e in final_list_ec_nos_per_genome for e in ecn)
            
        elif "and" in value:
            ecn=value.split(" and ")
            enzyme_comission_threshold = all(e in final_list_ec_nos_per_genome for e in ecn)
        
        elif value in final_list_ec_nos_per_genome:
            enzyme_comission_threshold = True
            
        else:
            enzyme_comission_threshold = False
                
                
        if enzyme_comission_threshold is True:
            ec.at[key,genome]=1
            
            
        else:
            pass
        print(enzyme_comission_threshold, genome,key,value)
        
id_to_reaction_dic={}
for row, index in reactions_table.iterrows():
    id_to_reaction_dic.update({index["Reaction ID"]:index["Reaction Name"]})
    
for key in ecnumbers.keys():
    ec.at[key,"BiGG Models Name"]=id_to_reaction_dic[key]

ec.set_index("BiGG Models Name",drop=False)
    
ec.to_csv("output/NCBI_annotation_based_functional_assignment.txt",index=None,sep="\t",mode="w")

###################
# HOmology_Vs_Annotation.py
# Copyright 2023, Abel Ingle and Daniel Noguera
# Revised: January 26, 2023
###################



cMAG_directory= 'refgenomes'

cMAG_files=Path(cMAG_directory).glob('*.fasta')

mf_wo_ind=pd.read_csv("output/BLAST_based_functional_assignment.txt",sep="\t")
mf=pd.DataFrame(mf_wo_ind.set_index("BiGG Models Name",drop=False))

ncbi_ec_wo_ind=pd.read_csv("output/NCBI_annotation_based_functional_assignment.txt",sep="\t")
ncbi_ec=pd.DataFrame(ncbi_ec_wo_ind.set_index("BiGG Models Name",drop=False))

jgi_ec_wo_ind=pd.read_csv("output/JGI_annotation_based_functional_assignment.txt",sep="\t")
jgi_ec=pd.DataFrame(jgi_ec_wo_ind.set_index("BiGG Models Name",drop=False))

reactions_table=pd.read_table("codefiles/metabolic_reactions.txt", encoding='latin-1').fillna('')

ecnumbers={}
for row, ind in reactions_table.iterrows():
    if ind["Enzyme comission number"]=='':
        pass
    else:
        ecnumbers.update({ind["Reaction ID"]:ind["Enzyme comission number"]})


mk_Venndiagram_directory='mkdir output/Venn_diagrams'

subprocess.run([mk_Venndiagram_directory],shell=True)

cMAG_directory= 'refgenomes'

cMAG_files=Path(cMAG_directory).glob('*.fasta')

final=mf

for cmag in cMAG_files:
    
    #Ab = Contained in group A, but not B
    
    
    Ab=0

    #aB = Contained in group B, but not A
    aB=0

    #AB = Contained in both group A and B
    AB=0
    
    Bnj=0
    BNj=0
    BnJ=0
    bnJ=0
    bNj=0
    bNJ=0
    BNJ=0
    bnj=0



    
    for index in ncbi_ec.index:
        if ncbi_ec.at[index,cmag.stem] != 1 and jgi_ec.at[index,cmag.stem] != 1 and mf.at[index,cmag.stem] == 1:
            Bnj=Bnj+1
            ncbi_ec.at[index,cmag.stem]="B"
            final.at[index,cmag.stem]=1
        elif ncbi_ec.at[index,cmag.stem] == 1 and jgi_ec.at[index,cmag.stem] != 1 and mf.at[index,cmag.stem] == 1:
            BNj=BNj+1
            ncbi_ec.at[index,cmag.stem]="BN"
            final.at[index,cmag.stem]=1
        elif ncbi_ec.at[index,cmag.stem] != 1 and jgi_ec.at[index,cmag.stem] == 1 and mf.at[index,cmag.stem] == 1:
            BnJ=BnJ+1
            ncbi_ec.at[index,cmag.stem]="BJ"
            final.at[index,cmag.stem]=1
        elif ncbi_ec.at[index,cmag.stem] != 1 and jgi_ec.at[index,cmag.stem] == 1 and mf.at[index,cmag.stem] != 1:
            bnJ=bnJ+1
            ncbi_ec.at[index,cmag.stem]="J"
            final.at[index,cmag.stem]=1
        elif ncbi_ec.at[index,cmag.stem] == 1 and jgi_ec.at[index,cmag.stem] != 1 and mf.at[index,cmag.stem] != 1:
            bNj=bNj+1
            ncbi_ec.at[index,cmag.stem]="N"
            final.at[index,cmag.stem]=1
        elif ncbi_ec.at[index,cmag.stem] == 1 and jgi_ec.at[index,cmag.stem] == 1 and mf.at[index,cmag.stem] != 1:
            bNJ=bNJ+1
            ncbi_ec.at[index,cmag.stem]="NJ"
            final.at[index,cmag.stem]=1
        elif ncbi_ec.at[index,cmag.stem] == 1 and jgi_ec.at[index,cmag.stem] == 1 and mf.at[index,cmag.stem] == 1:
            BNJ=BNJ+1
            ncbi_ec.at[index,cmag.stem]="BNJ"
            final.at[index,cmag.stem]=1
        else:
            bnj=bnj+1
            ncbi_ec.at[index,cmag.stem]=""
            final.at[index,cmag.stem]=0
    # Label figures
    number_of_absent_reactions=str(bnj)
    
    venn3(subsets = (Bnj,bNj,BNj,bnJ,BnJ,bNJ,BNJ), set_labels = ('BLAST-based','NCBI annotation informed','JGI annotation informed'),
          set_colors=('red','blue','yellow'),
          alpha = 0.5);
    venn3_circles(subsets = (Bnj,bNj,BNj,bnJ,BnJ,bNJ,BNJ))
    plt.title(cmag.stem+" ("+number_of_absent_reactions+" reactions absent)")
    plt.savefig('output/Venn_diagrams/{0}.png'.format(cmag.stem), facecolor="white", format="png")
    plt.close()

ncbi_ec.to_csv("output/BLAST_EC_set_operations.txt",index=None,sep="\t",mode="w")
header = 'echo "B: BLAST based \nN: NCBI annotation informed \nJ: JGI annotation informed \nBlank: neither\n" | cat - output/BLAST_EC_set_operations.txt > temp && mv temp output/BLAST_EC_set_operations.txt'
header_process = subprocess.run(header,
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
final.to_csv("output/genomic_potential_reactions.txt", index=None,sep="\t",mode="w")

# Rename DataFrame 

gpr=final

# Make directory to store fermentative pathway assignments

subprocess.call(["mkdir", "output/BLAST_annotation_fermentations"])

# Excel file that summarizes various fermentation pathways of various substrates is read as DataFrame

fermentative_lifestyle_reqs=pd.read_excel("codefiles/fermentation_requirements.xlsx",sheet_name=None)

# For every substrate of interest, make a DataFrame of substrate-specific pathways, 
# and their respective biochemical reactions involved

for key in fermentative_lifestyle_reqs.keys():
    sole_e_donor_pathways=pd.DataFrame(fermentative_lifestyle_reqs[key].fillna(''))
    pathways=list(sole_e_donor_pathways.columns)
    fermentations={key+"_fermentation_pathways": pathways}
    
    # A DataFrame indexed with substrate-specific pathways
    # that summarizes the presence / absence of each genome's pathway-specific protein encoding regions
    
    fermentative_lifestyles=pd.DataFrame(fermentations, index=pathways)
    fermentative_lifestyles2=pd.DataFrame(fermentations, index=pathways)
    
    # For every pathway queried per substrate
    
    for pathway in sole_e_donor_pathways.columns:
        
        # For every step of every pathway, make reactions possible into a list
        
        for row, index in sole_e_donor_pathways.iterrows():
            """
            For some metabolic reactions within a pathway, 
            sufficient biochemistry can proceed with varying enzymes or complexes;
            hence, an "or" condition is used.
            """
            
            if sole_e_donor_pathways.at[row,pathway] == "":
                pass
            
            else:
                split_cell=sole_e_donor_pathways.at[row,pathway].split(" or ")
                sole_e_donor_pathways.at[row,pathway]=split_cell
                
    # Now, add columns with genomes as headers
    
    cMAG_directory= 'refgenomes'
    cMAG_files=Path(cMAG_directory).glob('*.fna')
    
    for cmag in cMAG_files:
        
        for pw in pathways:
            
            if cmag.stem not in fermentative_lifestyles.columns:
                
                fermentative_lifestyles.insert(len(fermentative_lifestyles.columns),column=cmag.stem,value=0)
                fermentative_lifestyles2.insert(len(fermentative_lifestyles2.columns),column=cmag.stem,value=0)
                fermentative_lifestyles.insert(len(fermentative_lifestyles.columns),column=cmag.stem+"_missing_reactions",value="")
                
            truth=list()
            
            missing_reactions=list()
            
            """
            Parse DataFrame, and check if the BLAST based functional assignment indicated a 1 or 0
            for every biochemical reaction for every genome. If 0's are involved, add biochemical 
            reaction to list of missing reactions.
            """
            
            for row,index in sole_e_donor_pathways.iterrows():
                if sole_e_donor_pathways.at[row,pw]=="":
                    pass
                else:
                    question=any(gpr.at[rxn,cmag.stem] == 1 for rxn in sole_e_donor_pathways.at[row,pw])
                    truth.append(question)
                    
                    # Append list of reactions missing for every substrate-specific pathway
                    
                    if question is False:
                        missing_reactions.append(sole_e_donor_pathways.at[row,pw])
                        
            fermentative_lifestyles.at[pw,cmag.stem+"_missing_reactions"]=missing_reactions
            
            # If all biochemical reactions required for the pathway have at least one enzyme assigned 
            # a '1' in the functional assignments, then the genome receives a '1' in the pathways DataFrame.
            
            if all(truths is True for truths in truth):
                fermentative_lifestyles.at[pw,cmag.stem]=1
                fermentative_lifestyles2.at[pw,cmag.stem]=1
                print(cmag.stem,"can perform",pw,"fermentation using",key,"as an electron donor.")
               
    fermentative_lifestyles.to_csv("output/BLAST_annotation_fermentations/{}_fermentations_missing_reactions.txt".format(key), index=None,sep="\t",mode="w")
    fermentative_lifestyles3=fermentative_lifestyles2.transpose()
    fermentative_lifestyles3.to_csv("output/BLAST_annotation_fermentations/{}_fermentations.txt".format(key),sep="\t",header=None,mode="w")
print("Done predicting substrate-specific fermentative pathways.")


# In[ ]:




