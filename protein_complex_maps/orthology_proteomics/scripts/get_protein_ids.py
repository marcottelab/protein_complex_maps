from __future__ import print_function
import pandas as pd
import sys
import argparse
import logging

#This needs some work

def create_tables(species, experiment, elution_file, peptide_file):
    '''
    Break this  into multiple functions
    Break into 
    '''

    LOG_FILENAME = species + "_" +  experiment + "_" + '.log'
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
    uniq_protein_pep = peps.drop_duplicates()
    final_protein_pep = uniq_protein_pep.drop_duplicates(subset=['Peptide'], keep=False) #current Docs/version have subset
    print(uniq_protein_pep)
    print(final_protein_pep) 
     
    frac = pd.DataFrame(pd.read_csv(elution_file))
    #Remove undistinguishable isoleucine/leucine
    frac['Peptide'] = frac['Peptide'].str.replace('I', 'J')
    frac['Peptide'] = frac['Peptide'].str.replace('L', 'J')
    frac['Peptide'] = frac['Peptide'].str.replace('*', '')

    frac = frac.set_index(['Peptide'])
    
    #prot = pd.DataFrame(pd.read_csv(orthology_file, sep="\t"))
    #prot = prot.set_index(['ProteinID'])
    
    final_protein_pep = final_protein_pep.reset_index()
    final_protein_pep = final_protein_pep.set_index(['Peptide'])
    
    #Join peptides in experiment to unique protein Identifying peptides
    frac_protein = frac.join(final_protein_pep, how='left')
    print(final_protein_pep)
    
    print(frac_protein)

    msg = species + "\t" + experiment + "\t" + "protein" +"\t" +"identified_peptides" + "\t" + str(frac_protein.shape[0])
    logger.info(msg)
  
    frac_protein = frac_protein.reset_index()

    temp_filename = "peptide_assignments/" + species + "/elution_peptides" + "_" + experiment + "_" + species+ ".csv"
    frac_protein.to_csv(temp_filename, index=False)

    protein_identified_peptides = frac_protein[['Peptide']].drop_duplicates()
    #identified_prot = "peptide_assignments/" + species + "/exp_protein_identified_peptides" + "_" + experiment + "_" + species + ".csv"

    #protein_identified_peptides.to_csv(identified_prot)
    
    
    
    
    #ungrouped analysis    
    #protein_pep = prot.join(peps, how = "left")
    #msg = species + "\t" + experiment + "\t" + "protein_pep" + "\t" + str(protein_pep.shape[0])
    #logger.info(msg)

    #protein_pep = protein_pep.reset_index()
    
    
    #protein_pep_cols = ['ProteinID', 'Peptide']
    
    #protein_pep = protein_pep[protein_pep_cols]
    
    
    #Get unique Rows, as one peptide can occur multiples in one protein
    #msg = species + "\t" + experiment + "\t" + "protein" + "\t" + "nonredundant_peptides" + "\t" + str(uniq_protein_pep.shape[0])
    #logger.info(msg)

    #nonred_prot = "peptide_assignments/" + species +"/nonredundant_protein_peptides" + species + ".csv"
    #uniq_protein_pep.to_csv(nonred_prot)

    ###This is just dropping subsequent appearances. It's not removing anything with a duplicate. 

  
    #Get Peptides which only occur in one protein
    #msg = species + "\t" + experiment + "\t" + "protein" +"\t" +"unique_peptides" + "\t" + str(final_protein_pep.shape[0])
    #logger.info(msg)


    #ident_prot = "peptide_assignments/" + species + "/protein_unique_peptides_" + species + ".csv"
  
    #final_protein_pep.to_csv(ident_prot)
    
       
    #frac_protein = frac_protein[['ExperimentID', 'FractionID', 'ProteinID', 'PeptideCount']]

    

    
    #grouped_frac_protein = frac_protein.groupby(by=['ExperimentID', 'FractionID', 'ProteinID'])['PeptideCount'].sum()
    
    #final_frac_protein =  grouped_frac_protein.reset_index(name='Total_SpecCounts')
    #msg = species + "\t" + experiment + "\t" + "protein" + "\t" +"spectral_counts" + "\t" + str(final_frac_protein['Total_SpecCounts'].sum())
    #logger.info(msg)


    #elut_prot = "identified_elutions/" + species + "/" + experiment +"_elution_" + species + "_proteins.csv"
    #final_frac_protein.to_csv(elut_prot)
    
 
    
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
parser.add_argument('experiment', action="store", type=str)
parser.add_argument('elution_file', action="store", type=str)
parser.add_argument('peptides_file', action="store", type=str)
 

inputs = parser.parse_args()

print(inputs)

create_tables(inputs.species_code, inputs.experiment, inputs.elution_file, inputs.peptides_file)






    
