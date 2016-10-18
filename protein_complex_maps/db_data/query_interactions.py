import pandas as pd
import matplotlib as plt
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
        
        bait_interactions = db.session.query(cdb.Edge).filter(cdb.Edge.protein_key2.in_(bait), cdb.Edge.protein_key.in_(bait)).all()

        bait_query  =  make_data_frame(bait_interactions, ['protein_key', 'protein_key2', 'in_complex', 'score'])
        #print("interactions among baits")
        #print(bait_query)

        bait_median  = bait_query['score'].median()
        bait_std  = bait_query['score'].std()
 
        #print("bait median and standard deviation")
        #print(bait_median, bait_std) 
        lower_one_std_bound = bait_median - bait_std
        #print("lower bound")
        #print(lower_one_std_bound)

        return bait_query, lower_one_std_bound, bait_median, bait_std
    except NoResultFound:
        print "NoResultFound"


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
 

   args = parser.parse_args()
   return args

def run_process(bait, annots):

 
    #Get all interactions containing at least one bait protein
    df_all_prots = get_protein_interactions(bait)

    #Get all interactions between two bait proteins
    bait_query = get_baitbait_interactions(bait)[0]

    #Get lower limit for predictions (median of baitbait minus 1 stdev
    lower_one_std_bound = get_baitbait_interactions(bait)[1]

    #Predict interactions above lower score limit for predictions
    df_predicted = predict_complex_members(lower_one_std_bound, df_all_prots)

    #Annotate all interactions by whether they are bait-bait
    baitbait_annotated = annotate_baitbait(bait_query, df_predicted)
    #Annotate final table with genenames
    final_annotated = annotate_table(baitbait_annotated, annots)


    return final_annotated
  
def sampling_process(bait):
    print(bait)
   
    #max_interactions = len(bait) * len(bait)
    #print(max_interactions)
    baitbait_analysis = get_baitbait_interactions(bait)

    bait_query = baitbait_analysis[0]
  

    #HERE will be plot of distribution of gold standard lower limits/median 
    #Make with corum
    #Fixed image? A key for quality of input complex
    print("lower bound", baitbait_analysis[1])
    print("median", baitbait_analysis[2])
    print("stdev", baitbait_analysis[3])
    if baitbait_analysis[1] < 0.5:
          print("Low quality input complex", baitbait_analysis[1])
          print("Distribution of gold standard lows")
          print("and an arrow point to where input dist is")
  

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
        #print(lost, len(total_int.index))
        prop_02_interactions = 1 - float(lost)/len(total_int.index)
 
        #Bar graph of this
        print("involvement in +0.2 interactions", prop_02_interactions)

        if prop_02_interactions == 0.0:
           trimlist.append(bait[index])

        bait_sampling = get_baitbait_interactions(bait_subset)[0]
        lost = (len(bait_sampling.index))
        #print(lost, len(bait_query))
        prop_bait_interactions = 1 - float(lost)/len(bait_query)
        #Bar graph of this
        print("involvement in bait bait interactions", prop_bait_interactions)
    for prot in trimlist:
        bait.remove(prot) 
    print(trimlist, "had no interactions above score 0.2")
    return(bait)
     

def convert_id(bait, annots, id_format):
    
    bait_ids = annots[annots[id_format].isin(bait)]
    print(bait_ids)
    bait_vector = bait_ids['gene_id']
    print(bait_vector)
    bait= bait_vector.tolist()
    print(bait)
    return(bait)





if __name__ == "__main__":

    db.create_all()  # make our sqlalchemy tables
    args = load_args()

    bait  = args.bait.split(" ")

    annots = make_annots()
    if args.id_format != 'gene_id':
       try:
          bait = convert_id(bait, annots, args.id_format)
       except Exception as e:
           print(e)
  
    bait = sampling_process(bait)
    final_annotated = run_process(bait, annots)
    

    final_annotated.to_csv("query_output.txt", sep = "\t", index= False)











   #app.run()

    #OFD subcomplex
    #bait = ['8481', '5108', '79848', '197335', '64770'] 

    #JBW full list
    #This really is a process for user determined gold standard proteins
    #bait = ['123872', '352909', '161582', '54919', '25804', '1981', '10146', '10963', '9987', '3189', '283237', '150275', '63892', '10985', '8382', '26150', '120935', '10963']
#'123872 352909 161582 54919 25804 1981 10146 10963 9987 3189 283237 150275 63892 10985 8382 26150 120935 10963'

    #EIF3 full complex members A-M
    #bait = ['10480', '27335', '3646', '51386', '72868', '8661', '8662', '8663', '8664', '8665', '8666', '8667', '8668', '8669']
# 10480 27335 3646 51386 72868 8661 8662 8663 8664 8665 8666 8667 8668 8669 

def annotate_from_sql():
        #CDM trying to get protein name annotations... stuck
        tester2 = db.session.query(cdb.Edge, cdb.Protein).filter(cdb.Edge.protein_key2.in_(bait), cdb.Edge.protein_key.in_(bait)).all()
        #pd_tester2 = make_data_frame(tester2, [ 'protein_key', 'protein_key2', 'in_complex', 'score'])
        for row in tester2:
            print(row.Protein.proteinname)
            print(row.Edge.protein_key)
            print(row.Edge.protein_key2)
            print(row.Edge.score)
        #print(tester2.__table__.columns)
        #print(pd_tester2)


        #pd_tester2 = make_data_frame(tester2, [ 'protein_key', 'protein_key2', 'in_complex', 'score'])

        print(pd_tester2)

        #genenames_etc = db.session.filter(cdb.Protein).filter(cdb.Protein.gene_id=="3024")
        #print(genenames_etc)

 
        #print(bg_query.sort(['score'], ascending = 0).head)

       #complexes = []
    #return render_template('index.html', form=form, complexes=complexes)



def annotate_table2(baitbait_annotated):
        '''
        Working away from pandas joins to sql query
        Nope, this is way too slow.
        Pandas Joins.
        '''   
        print("start sql")
        joined_key1 = db.session.query(cdb.Edge, cdb.Conversion).filter(cdb.Conversion.gene_id==cdb.Edge.protein_key).all()#.add_columns(cdb.Conversion.proteinname)
        print("end sql")
        #print(dir(joined_key1))
        #for g in joined_key1:
        #    print "", g
        #    for m in g.Edge:
        #        print 2 * " ", m
        #        for i in m.Conversion:
        ##            print 4 * " ", i
        #for x in joined_key1:
            #print(dir(x.Edge))
          #  print(x.Conversion.proteinname, x.Edge.protein_key)
         #   #print(dir(x.
            #print(x.Edge.genename)      

        df_joined= make_data_frame(joined_key1, ['protein_key', 'protein_key2','proteinname', 'in_complex', 'score'])
        print(df_joined)
        # .group_by(cdb.Edge.gene_id.first()) 
        joined_key2 = db.session.query(cdb.Edge, cdb.Conversion).filter(cdb.Conversion.gene_id==cdb.Edge.protein_key2).all()#.add_columns(cdb.Conversion.proteinname)
        print(dir(joined))


