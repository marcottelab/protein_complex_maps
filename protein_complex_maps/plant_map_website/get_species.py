import pandas as pd
import argparse


def identify_species(input_list, conversion_tbl):

    Species_dict={'allplants':'All plants', 'arath':'Arabidopsis thaliana', 'braol':'Brassica oleraceae', 'chlre':'Chlamydomonas reinhardtii', 'orysj':'Oryza sativa', 'traes':'Tricicum aestivum', 'selml':'Selaginella moellendorfii', 'cerri':'Ceratopteris ricardii'}
   
     
    #print conversion_tbl.head
    list_species = []
    #print(input_list)    
    for ID in input_list:
        output = conversion_tbl[conversion_tbl['ProteinID'].str.contains(ID)]
        print(output)
        spec = output['Species'].tolist()[0]
        list_species.append(spec)

    if len(list(set(list_species))) > 1:
         print "All query proteins must be from the same species"
         print list_species         

         uniprot_query_str = ("+OR+").join(input_list)
         uniprot_query_str = "http://www.uniprot.org/uniprot/?query=" + uniprot_query_str + "&sort=score&columns=id%2Centry%20name%2Corganism" 

       #ihttp://www.uniprot.org/uniprot/?query=yourlist:M20170315F725F458AC8690F874DD868E4ED79B88A95DDD1&sort=yourlist:M20170315F725F458AC8690F874DD868E4ED79B88A95DDD1&columns=yourlist(M20170315F725F458AC8690F874DD868E4ED79B88A95DDD1),id%2Centry%20name%2Cgenes%28PREFERRED%29%2Cgenes%2Ccomment%28DISRUPTION%20PHENOTYPE%29
         print uniprot_query_str
      

     

    elif len(list(set(list_species))) == 0:
         print "Could not determine species of input queries"
         print list_species
  

    else:
         print Species_dict[spec] + " query"

         return spec, Species_dict[spec]

 
if __name__ == '__main__':

    input_list=["WRK58_ARATH", "Q9SR92", "W5B5G3"]
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")

    identify_species(input_list, conversion_tbl)
