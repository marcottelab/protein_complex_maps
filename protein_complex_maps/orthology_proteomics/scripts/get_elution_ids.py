from __future__ import print_function
import pandas as pd
import sys
import argparse
import logging
def create_tables(species, level, experiment, orthology_file, elution_file, peptide_file, contam_file):
    '''
    Break this  into multiple functions
    Break into 
    '''

    LOG_FILENAME = species + "_" +  experiment + "_" + level + '.log'
    print("LOGFILENAME", LOG_FILENAME)
    #print("SPECIES", species)
    #print("Experiment", experiment)
    #print("LEVE:", level)

    logger = logging.getLogger()
    filehandler = logging.FileHandler(LOG_FILENAME, mode="w")
    formatter = logging.Formatter('%(message)s')
    filehandler.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(filehandler)


    logger.setLevel(logging.INFO)


    #experiment = elution_file.replace(".csv", "")
    #:experiment = experiment.split("/")[-1]

       
    peps = pd.DataFrame(pd.read_csv(peptide_file)) 
    print("Peptide file loaded")
    #Remove undistinguishable isoleucine/leucine with J
    peps['Peptide'] = peps['Peptide'].str.replace('I', 'J')
    peps['Peptide'] = peps['Peptide'].str.replace('L', 'J')


    peps = peps.set_index(['ProteinID'])
 
    print(peps)

    contam = pd.DataFrame(pd.read_csv(contam_file))
    contam['Peptide'] = contam['Peptide'].str.replace('I', 'J')
    contam['Peptide'] = contam['Peptide'].str.replace('L', 'J')
    contam_list = contam['Peptide'].tolist()
  
    #print(contam_list)

    
    frac = pd.DataFrame(pd.read_csv(elution_file))
    #Remove undistinguishable isoleucine/leucine
    frac['Peptide'] = frac['Peptide'].str.replace('I', 'J')
    frac['Peptide'] = frac['Peptide'].str.replace('L', 'J')
    frac['Peptide'] = frac['Peptide'].str.replace('*', '')
    contam=contam.set_index(['Peptide'])



    #print(frac.shape)
    #print("pre")
    #print(frac[frac['Peptide']=="TNAENEFVTJK"])
    #print(contam[contam['Peptide']=="TNAENEFVTJK"])
    frac = frac[~frac.isin(contam_list)]
    #print("post")
    #print(frac[frac['Peptide']=="TNAENEFVTJK"])
    #print(frac.shape)

    print("find pep")
    
    frac = frac.set_index(['Peptide'])

    #print(contam)
    #print(frac)

    #join_contam = contam.join(frac, how="inner")
    #print(join_contam)
    #print(join_contam.shape)
    #print("DONE")

   
    



    prot = pd.DataFrame(pd.read_csv(orthology_file, sep="\t"))
    prot = prot.set_index(['ProteinID'])
    
    
    print(peps)
   
    
    group_pep = prot.join(peps, how = "left")

    msg = species + "\t" + experiment + "\t" + "total" +"\t" +"proteome_peptides" + "\t"   + str(group_pep.shape[0])
    logger.info(msg)

    group_pep = group_pep.reset_index()
    
    group_pep_cols = ['GroupID', 'Peptide']
    
    group_pep = group_pep[group_pep_cols]
    
    
    #Get unique Rows, as one peptide can occur multiples in one group
    uniq_group_pep = group_pep.drop_duplicates()
    msg = species + "\t" + experiment + "\t" + level +"\t" + "nonredundant_peptides" + "\t" + str(uniq_group_pep.shape[0])
    logger.info(msg)


    #formerly nonredundant_orthogroup
    nonred_group = "peptide_assignments/" + species + "/nonredundant_orthogroup_peptides_" + species + "_" + level + ".csv"
    uniq_group_pep.to_csv(nonred_group)
    
    #uniq_group_pep = final_group_pep.set_index(['Peptide'])
    #Get Peptides which only occur in one group
    final_group_pep = uniq_group_pep.drop_duplicates(subset=['Peptide'], keep=False) #current Docs/version have subset

    msg = species + "\t" + experiment + "\t" + level + "\t" + "unique_peptides" + "\t" + str(final_group_pep.shape[0])
    logger.info(msg)
    #This is currently same between group and protein identifications.
 
    #formerly identifying_orthogroup
    ident_group = "peptide_assignments/" + species + "/orthogroup_unique_peptides" + species + "_" + level + ".csv"
    final_group_pep.to_csv(ident_group)
    
    final_group_pep = final_group_pep.set_index(['Peptide'])
    
    #Join peptides in experiment to unique Group Identifying peptides
    frac_group = frac.join(final_group_pep, how='left')

    msg = species + "\t" + experiment + "\t" + level + "\t" + "identified_peptides" + "\t" + str(frac_group.shape[0])
    logger.info(msg)

   
    frac_group = frac_group.reset_index()
    group_identified_peptides = frac_group[['Peptide']].drop_duplicates()

    #formerly identifiedpeps
    identified_group = "peptide_assignments/" + species + "/exp_orthogroup_identified_peptides" + "_" + experiment + "_" + species + "_" + level + ".csv"
    group_identified_peptides.to_csv(identified_group)
    
        
    frac_group = frac_group[['ExperimentID', 'FractionID', 'GroupID', 'PeptideCount']]
    
    grouped_frac_group = frac_group.groupby(by=['ExperimentID', 'FractionID', 'GroupID'])['PeptideCount'].sum()
    
    final_frac_group =  grouped_frac_group.reset_index(name='Total_SpecCounts')
   
    msg = species + "\t" + experiment + "\t" + level + "\t" +"spectral_counts" + "\t" + str(final_frac_group['Total_SpecCounts'].sum())
    logger.info(msg)



    elut_group = "identified_elutions/" + species + "/" + experiment +"_elution_" + species + "_" + level + ".csv"
    final_frac_group.to_csv(elut_group)


    #ungrouped analysis    
    protein_pep = prot.join(peps, how = "left")
    #msg = species + "\t" + experiment + "\t" + "protein_pep" + "\t" + str(protein_pep.shape[0])
    #logger.info(msg)

    protein_pep = protein_pep.reset_index()
    
    
    protein_pep_cols = ['ProteinID', 'Peptide']
    
    protein_pep = protein_pep[protein_pep_cols]
    
    
    #Get unique Rows, as one peptide can occur multiples in one protein
    uniq_protein_pep = protein_pep.drop_duplicates()
    msg = species + "\t" + experiment + "\t" + "protein" + "\t" + "nonredundant_peptides" + "\t" + str(uniq_protein_pep.shape[0])
    logger.info(msg)

    nonred_prot = "peptide_assignments/" + species +"/nonredundant_protein_peptides" + species + ".csv"
    uniq_protein_pep.to_csv(nonred_prot)

    ###This is just dropping subsequent appearances. It's not removing anything with a duplicate. 

  
    #Get Peptides which only occur in one protein
    final_protein_pep = uniq_protein_pep.drop_duplicates(subset=['Peptide'], keep=False) #current Docs/version have subset
    msg = species + "\t" + experiment + "\t" + "protein" +"\t" +"unique_peptides" + "\t" + str(final_protein_pep.shape[0])
    logger.info(msg)


    ident_prot = "peptide_assignments/" + species + "/protein_unique_peptides_" + species + ".csv"
  
    final_protein_pep.to_csv(ident_prot)
    
    final_protein_pep = final_protein_pep.set_index(['Peptide'])
    
    #Join peptides in experiment to unique protein Identifying peptides
    frac_protein = frac.join(final_protein_pep, how='left')
    msg = species + "\t" + experiment + "\t" + "protein" +"\t" +"identified_peptides" + "\t" + str(frac_protein.shape[0])
    logger.info(msg)
  
    frac_protein = frac_protein.reset_index()

    temp_filename = "peptide_assignments/" + species + "/elution_peptides" + "_" + experiment + "_" + species+ ".csv"
    frac_protein.to_csv(temp_filename, index=False)

    protein_identified_peptides = frac_protein[['Peptide']].drop_duplicates()
    identified_prot = "peptide_assignments/" + species + "/exp_protein_identified_peptides" + "_" + experiment + "_" + species + ".csv"

    protein_identified_peptides.to_csv(identified_prot)
    
        
    frac_protein = frac_protein[['ExperimentID', 'FractionID', 'ProteinID', 'PeptideCount']]

    

    
    grouped_frac_protein = frac_protein.groupby(by=['ExperimentID', 'FractionID', 'ProteinID'])['PeptideCount'].sum()
    
    final_frac_protein =  grouped_frac_protein.reset_index(name='Total_SpecCounts')
    msg = species + "\t" + experiment + "\t" + "protein" + "\t" +"spectral_counts" + "\t" + str(final_frac_protein['Total_SpecCounts'].sum())
    logger.info(msg)


    elut_prot = "identified_elutions/" + species + "/" + experiment +"_elution_" + species + "_proteins.csv"
    final_frac_protein.to_csv(elut_prot)
    
 
    
    #Test identifiable peptides
    
    #No difference in arabidopsis? 
    #Check every step
    #Grouping condences 28000 identifications to 22000
    
    #claire@Oppenheimer:~/for_mySQL$ awk -F',' '{print $4}' singleIdentified.csv | sort -u | wc -l
    #10381
    #claire@Oppenheimer:~/for_mySQL$ awk -F',' '{print $4}' GroupIdentified.csv | sort -u | wc -l
    #6015
    
    #So, grouping is getting fewer groups, but average large size groups. Grouping get 20% more identifying peptides
    
    
    
    
    
    
parser = argparse.ArgumentParser(description='Interpret mass spec experiments using orthologous groups of proteins to make identifications')

parser.add_argument('species_code', action="store", type=str)
parser.add_argument('phylogenetic_level', action="store", type=str)
parser.add_argument('experiment', action="store", type=str)
parser.add_argument('orthology_file', action="store", type=str)
parser.add_argument('elution_file', action="store", type=str)
parser.add_argument('peptides_file', action="store", type=str)
parser.add_argument('contam_file', action="store", type=str)


inputs = parser.parse_args()

print(inputs.species_code)
print(inputs.phylogenetic_level)
print("experiment", inputs.experiment)
print(inputs.orthology_file)

create_tables(inputs.species_code, inputs.phylogenetic_level, inputs.experiment, inputs.orthology_file, inputs.elution_file, inputs.peptides_file, inputs.contam_file)






    
