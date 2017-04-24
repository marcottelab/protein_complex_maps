import pandas as pd
import logging
import protein_complex_maps.complex_map_website.plot_corum_dists as pcd
#import matplotlib.pyplot as plt
import sys
#import seaborn as sns
import protein_complex_maps.complex_map_website.complex_db as cdb
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_
import argparse


db = cdb.get_db()



#app = cdb.get_app()

#from flask.ext.wtf import Form
#from flask_wtf import Form
#from wtforms.fields import StringField, SubmitField

def make_data_frame(query, columns):
    """
    Takes a sqlalchemy query and a list of columns, returns a dataframe.
    Move this to other file cdb.make_data_frame
    """
    def make_row(x):
        return dict([(c, getattr(x, c)) for c in columns])
    return pd.DataFrame([make_row(x) for x in query])


def get_protein_interactions(bait):
    '''
    This query takes a list of bait IDs in a 'known' complex and finds protins which fit in with them
    '''
    #print(bait)
    try:
          
        bg_interactions = db.session.query(cdb.Edge).filter((cdb.Edge.protein_key2.in_(bait) | cdb.Edge.protein_key.in_(bait))).all()
        bg_query  =  make_data_frame(bg_interactions, ['protein_key', 'protein_key2', 'in_complex', 'score'])
        #print("pull down")
        #print(bg_query)

        #print("median score from query", bg_query['score'].median())

        #print("each protein median correlation")
        all_prots = pd.concat([bg_query[['protein_key', 'score']], bg_query[['protein_key2', 'score']]])
        all_prots['protein_key'] = all_prots['protein_key'].fillna(all_prots['protein_key2'])
        #print(all_prots)
        df_all_prots = pd.DataFrame(bg_query.groupby(['protein_key'])['score'].median()).sort_values(by='score', ascending=False)
        #print(df_all_prots)
        return(df_all_prots)

    except NoResultFound:
        print "NoResultFound"
     
def get_baitbait_interactions(bait):
    '''
    Find interactions between input bait complex proteins
    '''
    try:
        print(bait) 
        bait_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.protein_key2.in_(bait), cdb.Edge.protein_key.in_(bait)).all()

        bait_query  =  make_data_frame(bait_interactions, ['protein_key', 'protein_key2', 'in_complex', 'score'])
        #print("interactions among baits")
        #print(bait_query)

        bait_median  = bait_query['score'].median()
        bait_std  = bait_query['score'].std()
        bait_mean = bait_query['score'].mean() 
        #print("bait median and standard deviation")
        #print(bait_median, bait_std) 
        lower_one_std_bound = bait_median - bait_std
        #print("lower bound")
        #print(lower_one_std_bound)

        return bait_query, lower_one_std_bound, bait_median, bait_std, bait_mean
    except NoResultFound:
        print "NoResultFound"
    except Exception as e: 
        print(e)

def annotate_baitbait(bait_query, df_predicted):
        '''
        Annotate interactions that are bait-bait
        '''
        print("annotating bait_query")
        bait_query['bait_bait'] = True
        bait_query_annot = bait_query[['bait_bait', 'protein_key','protein_key2', 'score']]
        print(bait_query_annot)
        bait_query_annot = bait_query_annot.set_index(['protein_key', 'protein_key2', 'score'])

        df_predicted = df_predicted.set_index(['protein_key', 'protein_key2', 'score'])

        baitbait_annotated = df_predicted.join(bait_query_annot, how="outer")
        baitbait_annotated = baitbait_annotated.reset_index()
        print(baitbait_annotated)
        return baitbait_annotated

def predict_complex_members(lower_one_std_bound, df_all_prots):
        '''
        Get interactions with score above median- 1 std of bait-bait interactions
        '''

        predicted_members = df_all_prots[df_all_prots['score'] >= lower_one_std_bound]
        print(predicted_members) 

        annots = list(predicted_members.index.astype(str))
        print(annots)
       
        pooled_interactions = bait + annots 
        print(pooled_interactions)

        final_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.protein_key2.in_(pooled_interactions), cdb.Edge.protein_key.in_(pooled_interactions)).all()
      
        pd_total= make_data_frame(final_interactions,  ['protein_key', 'protein_key2', 'in_complex', 'score'])
        print(pd_total)
        #df_predicted = pd_total[pd_total.values in bait]
        #print(df_predicted)


        #Get high scoring interactions
        final_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.protein_key2.in_(bait) | cdb.Edge.protein_key.in_(bait)). filter(cdb.Edge.score >= lower_one_std_bound).all()

        df_predicted= make_data_frame(final_interactions,  ['protein_key', 'protein_key2', 'in_complex', 'score'])
        print(df_predicted)
        print(bait)
        return df_predicted
    
def make_annots():


     annotation_query = db.session.query(cdb.Conversion).all()
     annots = make_data_frame(annotation_query, ['genename', 'proteinname', 'gene_id', 'uniprot_acc'])
     annots = annots.convert_objects(convert_numeric=True)

     annots = annots[~annots['gene_id'].isnull()]
     return annots

def annotate_table(baitbait_annotated, annots):

        '''
        So, database held in SQL. Pandas for working with it
        SQLalchemy output needs to be parsed into table after join
        Too much hassle
        '''
       
        print(annots['gene_id'])

        #cdm Can also choose proteinname, uniprot_acc
        annots = annots[['gene_id', 'genename']]
        print(annots['gene_id'].dtype) 
        #cdm gene_id = Entrez gene id
        annots= annots.set_index(['gene_id'])
        #cdm Index column needs to match target
        annots.index.names =["protein_key"]
        print(annots)
        print(baitbait_annotated['protein_key'].dtype) 
        baitbait_annotated = baitbait_annotated.set_index(['protein_key'])
        print(baitbait_annotated)
    
        protein_key1_annot = baitbait_annotated.join(annots, how="left")
        protein_key1_annot = protein_key1_annot.reset_index()
        print(protein_key1_annot) 
        #cdm Same thing for the second interaction partner
        protein_key1_annot = protein_key1_annot.set_index(['protein_key2'])
        print("start pandas")

        annots.index.names = ['protein_key2']
        protein_key2_annot = protein_key1_annot.join(annots, how = "left", rsuffix="2")
        protein_key2_annot = protein_key2_annot.reset_index()
 
        #cdm Order columns
        final_annotated = protein_key2_annot[['score','genename','genename2', 'protein_key', 'protein_key2', 'bait_bait']]
       
        final_annotated = final_annotated.sort_values(by='score', ascending = False)
        return final_annotated

def load_args():

    #cdm This is for the website input
    #enrichment = request.args.get('enrichment')
    #form = SearchForm()
    #print enrichmenti
   parser = argparse.ArgumentParser(description="Adds positive and negative labels to feature matrix")
   parser.add_argument("--bait_complex", action="store", dest="bait", required=True, 
                  help="A space separating string of geneIDs to use as a bait complex")
   parser.add_argument("--format", action="store", dest="id_format", required=True, 
                  help="gene_id, genename, or uniprot_acc")
   parser.add_argument("--logname", action="store", dest="logname", required=False, default='shared_bait_feature.log',
                                    help="filename for logging, default=shared_bait_feature.log")

   parser.add_argument("--save_stats_file", action="store", dest="stats_file", required=False)
   parser.add_argument("--corum_scores", action="store", dest="corum_scores", required=True)


   args = parser.parse_args()
   print(args)
   print(args.bait)
   return args

def run_process(bait, annots):

 
    #Get all interactions containing at least one bait protein
    df_all_prots = get_protein_interactions(bait)

    #Get all interactions between two bait proteins
    baitbait_analysis = get_baitbait_interactions(bait)

    bait_query = baitbait_analysis[0]



    #Get lower limit for predictions (median of baitbait minus 1 stdev
 
    lower_one_std_bound = baitbait_analysis[1]
    median = baitbait_analysis[2]
    stdev = baitbait_analysis[3]
    mean = baitbait_analysis[4]
     



    #Predict interactions above lower score limit for predictions
    df_predicted = predict_complex_members(lower_one_std_bound, df_all_prots)

    #Annotate all interactions by whether they are bait-bait
    baitbait_annotated = annotate_baitbait(bait_query, df_predicted)
    #Annotate final table with genenames
    final_annotated = annotate_table(baitbait_annotated, annots)


    return final_annotated
  
def sampling_process(bait, stats_file=None):
    print(bait)
   
    #max_interactions = len(bait) * len(bait)
    #print(max_interactions)

    
    baitbait_analysis = get_baitbait_interactions(bait)

    bait_query = baitbait_analysis[0]

    lower_one_std_bound = baitbait_analysis[1]
    median = baitbait_analysis[2]
    stdev = baitbait_analysis[3]
    mean = baitbait_analysis[4]
  

  
    if stats_file:
        stats = list(get_baitbait_interactions(bait)[1:])
        log_data = str(stats).replace("[", "").replace("]", "")
        print("log_data", log_data)
        stat_log = open(stats_file, "a")
        
        stat_log.write(log_data +"\n")

    #HERE will be plot of distribution of gold standard lower limits/median 
    #Make with corum
    #Fixed image? A key for quality of input complex
    #print("lower bound", baitbait_analysis[1])
    #print("median", baitbait_analysis[2])
    #print("stdev", baitbait_analysis[3])
    #if baitbait_analysis[1] < 0.5:
    #      print("Low quality input complex", baitbait_analysis[1])
    #      print("Distribution of gold standard lows")
    #      print("and an arrow point to where input dist is")

    lower_list = []
    median_list = []
    stdev_list = []
    mean_list = []
  
    baitbait_dict={}
    zeropoint2_dict={}

    total_int = get_protein_interactions(bait) 
    print(len(total_int.index))
    trimlist = []
    #full_query = len(bait_query.index)
    #print(len(bait_query.index))
    for index in range(len(bait)):
        #Not using pop because not deleting item permanently
        print("gene ID", bait[index])
        bait_subset = bait[:index] + bait[index+1 :]
        #print(bait_subset)
        samp_total_int = get_protein_interactions(bait_subset) 
        lost =len(samp_total_int.index)
        totsbait =  len(total_int.index)
        prop_02_interactions = 1 - float(lost)/len(total_int.index)
        num_02_interactions = totsbait - lost 
        print(num_02_interactions, "what", totsbait, lost)
        #Bar graph of this
        #print("involvement in +0.2 interactions", prop_02_interactions)
        
        zeropoint2_dict[bait[index]] = num_02_interactions
        print(zeropoint2_dict)

        if prop_02_interactions == 0.0:
           trimlist.append(bait[index])

        bait_sample_analysis = get_baitbait_interactions(bait_subset)
        try:
            bait_sampling = bait_sample_analysis[0]
        except:
            continue
        lower_list.append(bait_sample_analysis[1])
        median_list.append(bait_sample_analysis[2])
        stdev_list.append(bait_sample_analysis[3])
        mean_list.append(bait_sample_analysis[4])

        diff_median = median -bait_sample_analysis[2]
        if diff_median < 0:
           print("removing ", bait[index], " would improve median by ", diff_median)
  
        lower_list.append(bait_sample_analysis[1])
        median_list.append(bait_sample_analysis[2])
        stdev_list.append(bait_sample_analysis[3])
        mean_list.append(bait_sample_analysis[4])




        lost = (len(bait_sampling.index))
        #print(lost, len(bait_query))
        prop_bait_interactions = 1 - float(lost)/len(bait_query)
        num_bait_interactions = len(bait_query.index) - lost
        #Bar graph of this
        #print("involvement in bait bait interactions", prop_bait_interactions)
        baitbait_dict[bait[index]] = num_bait_interactions

    for prot in trimlist:
        bait.remove(prot) 
    print(trimlist, "had no interactions above score 0.2")
    print(zeropoint2_dict)
    print(baitbait_dict)
    pcd.corum_plots(args.corum_scores, lower_one_std_bound, median, stdev, mean, lower_list, median_list, stdev_list, mean_list)
   
    pcd.interaction_share(baitbait_dict, "Number of bait-bait interactions", zeropoint2_dict, "Number of interactions above score 0.2")
    #pcd.interaction_share(zeropoint2_dict, "Involvement in interactions above score 0.2")



    return(bait)
     

def convert_id(bait, annots, id_format, target_id):
    
    bait_ids = annots[annots[id_format].isin(bait)]
    print(bait_ids)
    bait_vector = bait_ids[target_id]
    print(bait_vector)
    bait= bait_vector.tolist()
    print(bait)
    return(bait)



def setup_log(logname):

    LOG_FILENAME = logname
    global logger
    logger = logging.getLogger()
    filehandler = logging.FileHandler(LOG_FILENAME, mode="w")
    formatter = logging.Formatter('%(message)s')
    filehandler.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(filehandler)


    logger.setLevel(logging.INFO)


if __name__ == "__main__":

    db.create_all()  # make our sqlalchemy tables
    args = load_args()

    setup_log(args.logname)
    bait  = args.bait.split(",")
    annots = make_annots()
  
    print(bait) 
    if args.id_format != 'gene_id':
       try:
          bait = convert_id(bait, annots, args.id_format, "gene_id")
       except Exception as e:
           print(e)
    print("marker")
    bait_genenames = convert_id(bait, annots, "gene_id", "genename")

    print(bait)
    if args.stats_file:  
        sampling_process(bait, args.stats_file)

   
    bait = sampling_process(bait)
    final_annotated = run_process(bait, annots)

    full_network = final_annotated[['genename','genename2','score']]
    pcd.draw_network(full_network, bait_genenames)
    

    final_annotated.to_csv("query_output.txt", sep = "\t", index= False)







   #app.run()

    #OFD subcomplex
    #bait = ['8481', '5108', '79848', '197335', '64770'] 
 # '8481 5108 79848 197335 64770'
    #JBW full list
    #This really is a process for user determined gold standard proteins
    #bait = ['123872', '352909', '161582', '54919', '25804', '1981', '10146', '10963', '9987', '3189', '283237', '150275', '63892', '10985', '8382', '26150', '120935', '10963']
#'123872,352909,161582,54919,25804,1981,10146,10963,9987,3189,283237,150275,63892,10985,8382,26150,120935,10963'

    #EIF3 full complex members A-M
    #bait = ['10480', '27335', '3646', '51386', '72868', '8661', '8662', '8663', '8664', '8665', '8666', '8667', '8668', '8669']
#'10480,27335,3646,51386,72868,8661,8662,8663,8664,8665,8666,8667,8668,8669'

