import pandas as pd
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import complex_db as cdb
#import GroupID_complex_maps.complex_map_website.complex_db as cdb
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_
import argparse
from get_groups_from_prots import get_groups
import plot_corum_dists as pcd


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


def annotate_nodes(nodes, annots, bait, df_all_prots):
     '''bait in terms of euNOG
     '''
     #print "starting annotate_nodes"
     #nodes = nodes.set_index(['GroupID'])
     annots = annots.set_index(['GroupID'])
     #bait_series = pd.Series(bait)

     bait_df = pd.DataFrame(True, index=bait, columns = ['Bait'])
     #/print(bait_df)

     nodes_tmp = nodes.join(annots, how='left')
     #print "first nodes_tmp"
     #print nodes_tmp
     nodes_tmp = nodes_tmp.join(bait_df, how='left')
     #print "second one"
     #print nodes_tmp
     nodes_scores = nodes_tmp.join(df_all_prots, how='outer', rsuffix="bla")
     #print "df_all_prots"
     #print(df_all_prots)
     nodes_scores.index.name = 'GroupID'
     nodes_scores = nodes_scores.reset_index()
     #print(nodes_scores)
     nodes_mean = nodes_scores.groupby(['GroupID'])['score'].mean()
     nodes_annotated =nodes_tmp.join(nodes_mean, how="left")
     nodes_annotated = nodes_annotated.sort_values(['Degree', 'score'], ascending=[False, False])
     nodes_annotated = nodes_annotated.drop_duplicates()
     nodes_annotated['score'][nodes_annotated.Bait == True] = "-"
     nodes_annotated['Bait'][nodes_annotated.Bait != True] = "-"

     #print(nodes_annotated)
     return nodes_annotated



def get_group_interactions(bait):
    '''
    This query takes a list of bait IDs in a 'known' complex and finds protins which fit in with them
    '''
    #print(bait)
    try:
          
        bg_interactions = db.session.query(cdb.Edge).filter((cdb.Edge.GroupID_key2.in_(bait) | cdb.Edge.GroupID_key.in_(bait))).all()
        bg_query  =  make_data_frame(bg_interactions, ['GroupID_key', 'GroupID_key2', 'in_complex', 'score'])
        #print("pull down")
       #print(bg_query)

       #print("median score from query", bg_query['score'].median())

       #print("each GroupID median correlation")
        all_prots = pd.concat([bg_query[['GroupID_key', 'score']], bg_query[['GroupID_key2', 'score']]])
        all_prots['GroupID_key'] = all_prots['GroupID_key'].fillna(all_prots['GroupID_key2'])
        #print(all_prots)
        df_all_prots = pd.DataFrame(bg_query.groupby(['GroupID_key'])['score'].median()).sort_values(by='score', ascending=False)
        #print(df_all_prots)
        return(df_all_prots, bg_query)

    except Exception as E:
       print "NoResultFound"
     
def get_baitbait_interactions(bait):
    '''
    Find interactions between input bait complex GroupIDs
    '''
    try:
        
        bait_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.GroupID_key2.in_(bait), cdb.Edge.GroupID_key.in_(bait)).all()
        #print bait_interactions
        bait_query  =  make_data_frame(bait_interactions, ['GroupID_key', 'GroupID_key2', 'in_complex', 'score'])
       #print("interactions among baits")
       #print(bait_query)

        return bait_query
    except Exception as E:
        return NoResultFound


def stats_interactions(inp_query):
        try:
            inp_median  = inp_query['score'].median()
            inp_std  = inp_query['score'].std()
            inp_mean = inp_query['score'].mean() 
        #print("inp median and standard deviation")
        #print(inp_median, inp_std) 
            lower_one_std_bound = inp_median - inp_std
         #print("lower bound")
        #print(lower_one_std_bound)
        except Exception as E:
              lower_one_std_bound = 0
              inp_median=0
              inp_std= 0
              inp_mean=0   
        return lower_one_std_bound, inp_median, inp_std, inp_mean


def annotate_baitbait(bait_query, df_predicted):
        '''
        Annotate interactions that are bait-bait
        '''
        #print("annotating bait_query")
        bait_query['bait_bait'] = True
        bait_query_annot = bait_query[['bait_bait', 'GroupID_key','GroupID_key2', 'score']]
        #print(bait_query_annot)
        bait_query_annot = bait_query_annot.set_index(['GroupID_key', 'GroupID_key2', 'score'])

        df_predicted = df_predicted.set_index(['GroupID_key', 'GroupID_key2', 'score'])

        baitbait_annotated = df_predicted.join(bait_query_annot, how="outer")
        baitbait_annotated = baitbait_annotated.reset_index()
        #print(baitbait_annotated)
        return baitbait_annotated



def predict_complex_members(bait, lower_one_std_bound, df_all_prots):
        '''
        Get interactions with score above median- 1 std of bait-bait interactions
        '''

        predicted_members = df_all_prots[df_all_prots['score'] >= lower_one_std_bound]
        #print(predicted_members) 

        annots = list(predicted_members.index.astype(str))
        #print(annots)
       
        pooled_interactions = bait + annots 
        #print(pooled_interactions)

        final_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.GroupID_key2.in_(pooled_interactions), cdb.Edge.GroupID_key.in_(pooled_interactions)).all()
      
        pd_total= make_data_frame(final_interactions,  ['GroupID_key', 'GroupID_key2', 'in_complex', 'score'])
        #print(pd_total)
        #df_predicted = pd_total[pd_total.values in bait]
        #print(df_predicted)


        #Get high scoring interactions
        final_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.GroupID_key2.in_(bait) | cdb.Edge.GroupID_key.in_(bait)). filter(cdb.Edge.score >= lower_one_std_bound).all()

        df_predicted= make_data_frame(final_interactions,  ['GroupID_key', 'GroupID_key2', 'in_complex', 'score'])
        #print(df_predicted)
        #print(bait)
        return df_predicted
    
def make_annots():


     annotation_query = db.session.query(cdb.Conversion).all()
     annots = make_data_frame(annotation_query, ['GroupID', 'ProteinID', 'Species'])
     annots = annots.convert_objects(convert_numeric=True)

     #Maybe not needed
     #annots = annots[~annots['gene_id'].isnull()]
     return annots

def annotate_table(baitbait_annotated):

        '''
        So, database held in SQL. Pandas for working with it
        SQLalchemy output needs to be parsed into table after join
        Too much hassle
        '''
        eggnog_annots = pd.read_csv("all_annotations.csv")      

        eggnog_annots = eggnog_annots[['GroupID', 'Annotation']]
        eggnog_annots= eggnog_annots.set_index(['GroupID'])
        #cdm Index column needs to match target
        eggnog_annots.index.names =["GroupID_key"]
        #print(eggnog_annots)
        #print(baitbait_annotated['GroupID_key'].dtype) 
        baitbait_annotated = baitbait_annotated.set_index(['GroupID_key'])
        #print(baitbait_annotated)
    
        GroupID_key1_annot = baitbait_annotated.join(eggnog_annots, how="left")
        GroupID_key1_annot = GroupID_key1_annot.reset_index()
        #print(GroupID_key1_annot) 
        #cdm Same thing for the second interaction partner
        GroupID_key1_annot = GroupID_key1_annot.set_index(['GroupID_key2'])
        #print("start pandas")

        eggnog_annots.index.names = ['GroupID_key2']
        GroupID_key2_annot = GroupID_key1_annot.join(eggnog_annots, how = "left", rsuffix="2")
        GroupID_key2_annot = GroupID_key2_annot.reset_index()
 
        #cdm Order columns
        final_annotated = GroupID_key2_annot[['score','Annotation','Annotation2', 'GroupID_key', 'GroupID_key2', 'bait_bait']]
       
        final_annotated = final_annotated.sort_values(by='score', ascending = False)
        final_annotated['bait_bait'][final_annotated.bait_bait != True] = "-"

        return final_annotated

def load_args():

    #cdm This is for the website input
    #enrichment = request.args.get('enrichment')
    #form = SearchForm()
    #print enrichmenti
   parser = argparse.ArgumentParser(description="Adds positive and negative labels to feature matrix")
   parser.add_argument("--bait_complex", action="store", dest="bait", required=True, 
                  help="A space separating string of ids to use as a bait complex")
   parser.add_argument("--format", action="store", dest="id_format", required=True, 
                  help="groupid")

   parser.add_argument("--logname", action="store", dest="logname", required=False, default='shared_bait_feature.log',
                                    help="filename for logging, default=shared_bait_feature.log")

   parser.add_argument("--save_stats_file", action="store", dest="stats_file", required=False)

   args = parser.parse_args()
   return args

def run_process(bait, annots):

    #print("Run process")
    #Get all interactions containing at least one bait GroupID
    df_all_prots, df_all_interactions = get_group_interactions(bait)
    #print(df_all_prots)
    #Get all interactions between two bait GroupIDs
    try:
        bait_query = get_baitbait_interactions(bait)
    
        #Get lower limit for predictions (median of baitbait minus 1 stdev
        lower_one_std_bound = stats_interactions(bait_query)[0]
        lower_one_std_bound = get_baitbait_interactions(bait)
    
    
        #Predict interactions above lower score limit for predictions
        df_predicted = predict_complex_members(bait, lower_one_std_bound, df_all_prots)
    
        #Annotate all interactions by whether they are bait-bait
        baitbait_annotated = annotate_baitbait(bait_query, df_predicted)
    except Exception as e:
        #print e
        baitbait_annotated = annotate_baitbait(bait_query, df_all_interactions)
 

       
    #Annotated with eggNOG annotations
    final_annotated =  annotate_table(baitbait_annotated)
    return final_annotated, df_all_prots
  
def sampling_process(bait, stats_file=None):
   #print "START SAMPLING"
   #print(bait)

    #max_interactions = len(bait) * len(bait)
    #print(max_interactions)


    bait_query = get_baitbait_interactions(bait)
 
#CDM
    lower_one_std_bound, median, stdev, mean = stats_interactions(bait_query)
    #lower_one_std_bound = baitbait_analysis[0]
    #median = baitbait_analysis[1]
    #stdev = baitbait_analysis[2]
    #mean = baitbait_analysis[3]

    if stats_file:
        stats = list(stats_interactions(bait_query))
        log_data = str(stats).replace("[", "").replace("]", "")
        #print("log_data", log_data)
        stat_log = open(stats_file, "a")

        stat_log.write(log_data +"\n")

    lower_list = []
    median_list = []
    stdev_list = []
    mean_list = []

    baitbait_dict={}
    zeropoint2_dict={}

    total_int = get_group_interactions(bait)
    #print(len(total_int.index))
    trimlist = []
    #full_query = len(bait_query.index)
    #print(len(bait_query.index))
    suggestions = []
    for index in range(len(bait)):
        #Not using pop because not deleting item permanently
        #print("gene ID", bait[index])
        bait_subset = bait[:index] + bait[index+1 :]
        #print(bait_subset)

       # samp_total_int = get_group_interactions(bait_subset)[0]
       # lost =len(samp_total_int.index)
       # totsbait =  len(total_int.index)
       # prop_02_interactions = 1 - float(lost)/len(total_int.index)
       # num_02_interactions = totsbait - lost
        #print(num_02_interactions, "what", totsbait, lost)
        #Bar graph of this
        #print("involvement in +0.2 interactions", prop_02_interactions)

        #zeropoint2_dict[bait[index]] = num_02_interactions
        #print(zeropoint2_dict)

        #if prop_02_interactions == 0.0:
        #   trimlist.append(bait[index])

        bait_subset_query = get_baitbait_interactions(bait_subset)
        sample_analysis = stats_interactions(bait_subset_query)

        lower_list.append(sample_analysis[0])
        median_list.append(sample_analysis[1])
        stdev_list.append(sample_analysis[2])
        mean_list.append(sample_analysis[3])

        diff_median = median - sample_analysis[1]
        if diff_median < 0:
           suggestion = "Remove "+ str(bait[index]) + " to raise median by "+ str(abs(diff_median))
           suggestions.append(suggestion)
        diff_mean = mean - sample_analysis[3]
        if diff_mean < 0:
           suggestion = "Remove "+ str(bait[index]) + " to raise mean by "+ str(abs(diff_mean))
           suggestions.append(suggestion)


        #lower_list.append(bait_sample_analysis[1])
        #median_list.append(bait_sample_analysis[2])
        #stdev_list.append(bait_sample_analysis[3])
        #mean_list.append(bait_sample_analysis[4])




        #lost = (len(bait_sampling.index))
        #print(lost, len(bait_query))
        #prop_bait_interactions = 1 - float(lost)/len(bait_query)
        #num_bait_interactions = len(bait_query.index) - lost
        #Bar graph of this
        #print("involvement in bait bait interactions", prop_bait_interactions)
        #baitbait_dict[bait[index]] = num_bait_interactions

    for prot in trimlist:
        bait.remove(prot)
    #print(trimlist, "had no interactions above score 0.2")
    #print(zeropoint2_dict)
    #print(baitbait_dict)
    corum_scores = "mammal_corum_scores.txt"
    pcd.corum_plots_bokeh(corum_scores, lower_one_std_bound, median, stdev, mean, lower_list, median_list, stdev_list, mean_list)
    return median, mean, suggestions
    #pcd.interaction_share(baitbait_dict, "Number of bait-bait interactions", zeropoint2_dict, "Number of interactions above score 0.2")

    #pcd.interaction_share(zeropoint2_dict, "Involvement in interactions above score 0.2")
     

def convert_id(bait, annots, id_format):
    
    bait_ids = annots[annots[id_format].isin(bait)]
    #print(bait_ids)
    bait_vector = bait_ids['gene_id']
    #print(bait_vector)
    bait= bait_vector.tolist()
    #print(bait)
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

    bait  = args.bait.split(" ")
   #print bait
    #GroupID ProteinID Species
    annots = make_annots()

    if args.id_format != 'GroupID':
       try:
          #print "converting..."
          #print bait
          bait = get_groups(bait, annots) 
          #print bait
          #bait = convert_id(bait, annots, args.id_format)
       except Exception as e:
          print(e)
    #bait should now be in GroupID formate
    #print bait  

    postsamp_bait = sampling_process(bait, args.stats_file)
    #final_annotated = run_process(postsamp_bait, annots)
    

    #final_annotated.to_csv("query_output.txt", sep = "\t", index= False)



