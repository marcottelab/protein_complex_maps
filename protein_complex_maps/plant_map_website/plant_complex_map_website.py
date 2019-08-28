#!/usr/bin/env python
import random
from bs4 import BeautifulSoup
from collections import OrderedDict
#import protein_complex_maps.plant_map_website.complex_db as cdb
import plant_complex_db as cdb
import pandas as pd
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_
#import networkx as nx

db = cdb.get_db()
app = cdb.get_app()

from flask_wtf import FlaskForm  #Changed from just Form to avoid deprecation warning
from wtforms.fields import StringField, SubmitField, SelectField
from flask import make_response


from flask import render_template
from flask import url_for, redirect, request, jsonify

app.jinja_env.add_extension('jinja2.ext.do') #Allow appending to list in html 


def getitem(obj, item, default):
    if item not in obj:
        return default
    else:
        return obj[item]

def getitems(obj, item, default):
    if item not in obj:
        return default
    else:
        print obj.values()
        print obj
        return obj.getlist(item) 
       


class SearchForm(FlaskForm):
    species_list = [("arath", "Arabidopsis"), 
                    ("braol", "Broccoli"), 
                    ("chlre", "Chlamydomonas"),
                    ("cocnu", "Coconut"), 
                    ("cansa", "Hemp"), 
                    ("sollc", "Tomato"), 
                    ("maize", "Maize"), 
                    ("chqui", "Quinoa"), 
                    ("orysj","Rice"), 
                    ("soybn","Soy"), 
                    ("selml", "Spikemoss"), 
                    ("wheat", "Wheat"), 
                    ("cerri", "C-Fern")]

    OrthogroupID = StringField(u'virNOG Orthogroup ID ex: ENOG411DWGM :')
    ProteinID = StringField(u'Protein ID ex. F4KCR6, SC13A_ARATH, AT1G02090')
    Species = SelectField(u'Species', choices = species_list, default = 'arath')
    #enrichment = StringField(u'Enrichment (ex. cilium):')
    #submit = SubmitField(u'Search')
    submit = SubmitField(u'Search complexes')
   
    submitinteractions = SubmitField(u'Get top interactions')



@app.route("/")
def root(complexes=[]):
    print(complexes)
    #complexes = cdb.Complex.query.all()
    form = SearchForm()
    return render_template('index.html', form = form, complexes = complexes)

def OrthogroupQuery(Input_OrthogroupID, Species, error, cdb, template):
    OrthogroupID = db.session.query(cdb.Orthogroup).filter((func.upper(cdb.Orthogroup.OrthogroupID) == func.upper(Input_OrthogroupID))).first()
    print(OrthogroupID)
    #if len(OrthogroupIDs) == 0:
    if not OrthogroupID:
        #kdrew: input genename is not valid, flash message
        error = "Could not find given virNOG orthogroup ID: %s" % Input_OrthogroupID

        return render_template(template, form = form, complexes = [], error = error, Species = Species)
    return(OrthogroupID)

def ProteinQuery(Input_ProteinID, error, cdb, template):
    ProteinID = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.ProteinID) == func.upper(Input_ProteinID))).first()
    print(ProteinID)
    print(ProteinID.orthogroups)
    if not ProteinID.orthogroups:
        #kdrew: input ProteinID is not valid, flash message
        error = "Could not find orthogroup for given Protein ID: %s.\n Try using a Uniprot.org Accession" % Input_ProteinID

        return render_template(template, form = form, complexes = [], error = error)


    OrthogroupID_string = ProteinID.orthogroups.OrthogroupID
    Species = ProteinID.Spec

    OrthogroupID = db.session.query(cdb.Orthogroup).filter((func.upper(cdb.Orthogroup.OrthogroupID) == func.upper(OrthogroupID_string))).first()

    return ProteinID, OrthogroupID, Species

def troubleshoot_clusters(orthogroup_clusters):
                #Keep for trouble shooting syntax
                #for prot in orthogroup_clusters.hiercomplexes:
                #      print(prot)
                #      print(dir(prot))                            
                #      for bla in prot.orthogroups:
                #             print(dir(bla))
                #             print(bla.Orthoannots)
                             #print(dir(bla.Orthoannots))
                             #print(dir(bla))
                             #print(bla.OrthogroupID)
                #             print("Protein object", bla.Proteins)
                #             for p in bla.Proteins:
                #                  print(p.ProteinID, p.Spec)
                             #print(stop)             
            return 0


 
@app.route("/displayComplexesForOrthogroupID")
def displayComplexesForOrthogroupID():
    Species = request.args.get('Species')
    Input_OrthogroupID = request.args.get('OrthogroupID').strip().upper()
    form = SearchForm()
    error=None

    #CDM: See if orthogroup is a valid orthogroup ID
    OrthogroupID = OrthogroupQuery(Input_OrthogroupID, Species, error, cdb,'index.html')
    print(OrthogroupID.OrthogroupID)
    complexes = []
    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
    # Get hier complexes associated with an orthogroup ID, Then get member proteins of each cluster level
    complexes = orthogroup_clusters.hiercomplexes
    
    #troubleshoot_clusters(orthogroup_clusters)               
    if len(orthogroup_clusters.hiercomplexes) == 0:
        error = "No complexes found for given virNOG orthogroup ID: %s" % Input_OrthogroupID

    return render_template('getcomplexes.html', form = form, complexes = complexes, Species = Species, Input_OrthogroupID = Input_OrthogroupID, error = error)

@app.route("/displayComplexesForProteinID")
def displayComplexesForProteinID():
    Input_ProteinID = request.args.get('ProteinID').strip().upper()
    form = SearchForm()
    error=None 

    ProteinID, OrthogroupID, Species = ProteinQuery(Input_ProteinID, error, cdb, "index.html")
    complexes = []
    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
    complexes = orthogroup_clusters.hiercomplexes

    if len(orthogroup_clusters.hiercomplexes) == 0:
        error = "No complexes found for given Protein ID %s and its virNOG orthogroup ID: %s" % (Input_ProteinID, OrthogroupID.OrthogroupID)
    return render_template('getcomplexes.html', form = form, complexes = complexes, Species = Species, Input_ProteinID = Input_ProteinID, error = error)


@app.route("/getInteractionsForOrthogroupID")
def getInteractionsForOrthogroupID():
    Species = request.args.get('Species')
    Input_OrthogroupID = request.args.get('OrthogroupID').strip().upper()
    form = SearchForm()
    error=None

    OrthogroupID = OrthogroupQuery(Input_OrthogroupID, Species, error, cdb,'index.html')
    interactions = []
    for score in OrthogroupID.Scores:
        print(score.ScoreVal, score.InteractionID)
        interaction = db.session.query(cdb.Score).filter(cdb.Score.InteractionID == score.InteractionID).all()
        interactions.append(interaction)
    return render_template('getinteractions.html', form = form, interactions = interactions,  Species = Species, Input_OrthogroupID = Input_OrthogroupID, error = error)

@app.route("/getInteractionsForProteinID")
def getInteractionsForProteinID():
    Input_ProteinID = request.args.get('ProteinID').strip().upper()
    form = SearchForm()
    error=None

    ProteinID, OrthogroupID, Species = ProteinQuery(Input_ProteinID, error, cdb, "index.html")

    interactions = []
    for score in OrthogroupID.Scores:
        print(score.ScoreVal, score.InteractionID)
        interaction = db.session.query(cdb.Score).filter(cdb.Score.InteractionID == score.InteractionID).all()
        interactions.append(interaction)
    return render_template('getinteractions.html', form = form, interactions = interactions,  Species = Species, Input_ProteinID = Input_ProteinID, error = error)

def OrthogroupQuery(Input_OrthogroupID, Species, error, cdb, template):
    OrthogroupID = db.session.query(cdb.Orthogroup).filter((func.upper(cdb.Orthogroup.OrthogroupID) == func.upper(Input_OrthogroupID))).first()
    print(OrthogroupID)
    #if len(OrthogroupIDs) == 0:
    if not OrthogroupID:
        #kdrew: input genename is not valid, flash message
        error = "Could not find given virNOG orthogroup ID: %s" % Input_OrthogroupID

        return render_template(template, form = form, complexes = [], error = error, Species = Species)
    return(OrthogroupID)

def ProteinQuery(Input_ProteinID, error, cdb, template):
    ProteinID = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.ProteinID) == func.upper(Input_ProteinID))).first()
    print(ProteinID)
    print(ProteinID.orthogroups)
    if not ProteinID.orthogroups:
        #kdrew: input ProteinID is not valid, flash message
        error = "Could not find orthogroup for given Protein ID: %s" % Input_ProteinID

        return render_template(tmplate, form = form, complexes = [], error = error)


    OrthogroupID_string = ProteinID.orthogroups.OrthogroupID
    Species = ProteinID.Spec

    OrthogroupID = db.session.query(cdb.Orthogroup).filter((func.upper(cdb.Orthogroup.OrthogroupID) == func.upper(OrthogroupID_string))).first()

    return ProteinID, OrthogroupID, Species

def troubleshoot_clusters(orthogroup_clusters):
                #Keep for trouble shooting syntax
                #for prot in orthogroup_clusters.hiercomplexes:
                #      print(prot)
                #      print(dir(prot))                            
                #      for bla in prot.orthogroups:
                #             print(dir(bla))
                #             print(bla.Orthoannots)
                             #print(dir(bla.Orthoannots))
                             #print(dir(bla))
                             #print(bla.OrthogroupID)
                #             print("Protein object", bla.Proteins)
                #             for p in bla.Proteins:
                #                  print(p.ProteinID, p.Spec)
                             #print(stop)             
            return 0


 
#@app.route("/displayComplexesForOrthogroupID")
#def displayComplexesForOrthogroupID():
#    Species = request.args.get('Species')
#    Input_OrthogroupID = request.args.get('OrthogroupID')
#    form = SearchForm()
#    error=None
#
#    #CDM: See if orthogroup is a valid orthogroup ID
#    OrthogroupID = OrthogroupQuery(Input_OrthogroupID, Species, error, cdb,'index.html')
#    print(OrthogroupID.OrthogroupID)
#    complexes = []
#    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
#    # Get hier complexes associated with an orthogroup ID, Then get member proteins of each cluster level
#    complexes = orthogroup_clusters.hiercomplexes
#    
#    #troubleshoot_clusters(orthogroup_clusters)               
#    if len(orthogroup_clusters.hiercomplexes) == 0:
#        error = "No complexes found for given virNOG orthogroup ID: %s" % Input_OrthogroupID
#
#    return render_template('index.html', form = form, complexes = complexes, Species = Species, Input_OrthogroupID = Input_OrthogroupID, error = error)
#
#@app.route("/displayComplexesForProteinID")
#def displayComplexesForProteinID():
#    Input_ProteinID = request.args.get('ProteinID')
#    form = SearchForm()
#    error=None 
#
#    ProteinID, OrthogroupID, Species = ProteinQuery(Input_ProteinID, error, cdb, "index.html")
#    complexes = []
#    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
#    complexes = orthogroup_clusters.hiercomplexes
#
#    if len(orthogroup_clusters.hiercomplexes) == 0:
#        error = "No complexes found for given Protein ID %s and its virNOG orthogroup ID: %s" % (Input_ProteinID, OrthogroupID.OrthogroupID)
#    return render_template('index.html', form = form, complexes = complexes, Species = Species, Input_ProteinID = Input_ProteinID, error = error)


#@app.route("/getInteractionsForOrthogroupID")
#def getInteractionsForOrthogroupID():
#    print("GET INTERACTIONS")
#    Species = request.args.get('Species')
#    Input_OrthogroupID = request.args.get('OrthogroupID')
#    form = SearchForm()
#    error=None
#
#    # Need to get OrthogroupID propagating to Score
#    OrthogroupID = OrthogroupQuery(Input_OrthogroupID, Species, error, cdb,'index.html')
#    print(OrthogroupID.OrthogroupID)
#    print(dir(OrthogroupID))
#    print(OrthogroupID.Scores)
#    print(dir(OrthogroupID.Scores))
#    edgecodes = []
#    interactions = []
#    for score in OrthogroupID.Scores:
#        print(score.ScoreVal, score.InteractionID)
#        interaction = db.session.query(cdb.Score).filter(cdb.Score.InteractionID == score.InteractionID).all()
#        interactions.append(interaction)
#        #for p in pair:
#        #   print(p.scores.OrthogroupID)        
#
#          
#
#    #orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
# 
#
#    
#    #CDM: See if orthogroup is a valid orthogroup ID
#    #OrthogroupID = OrthogroupQuery(Input_OrthogroupID, Species, error, cdb,'getinteractions.html')
#    #print(dir(OrthogroupID))
#    return render_template('getinteractions.html', form = form, interactions = interactions,  Species = Species, Input_OrthogroupID = Input_OrthogroupID, error = error)
#




@app.route(u'/search', methods=[u'POST'])
def searchComplexes():
    form = SearchForm()
    complexes = []
    print("formvalues")
    print(form.submit, form.submit.data)
    print(form.submitinteractions, form.submitinteractions.data)

    if form.validate_on_submit():

        if form.submit.data == True:
            if len(form.OrthogroupID.data) > 0:
                if len(form.Species.data) > 0:
                   return redirect(url_for('displayComplexesForOrthogroupID', OrthogroupID = form.OrthogroupID.data, Species = form.Species.data))
    
                else:
                   return redirect(url_for('displayComplexesForOrthogroupID', OrthogroupID = form.OrthogroupID.data))
            elif len(form.ProteinID.data) > 0:
                return redirect(url_for('displayComplexesForProteinID', ProteinID = form.ProteinID.data))
        
        if form.submitinteractions.data == True:
            if len(form.OrthogroupID.data) > 0:
                if len(form.Species.data) > 0:
                   return redirect(url_for('getInteractionsForOrthogroupID', OrthogroupID = form.OrthogroupID.data, Species = form.Species.data))
    
                else:
                   return redirect(url_for('getInteractionsForOrthogroupID', OrthogroupID = form.OrthogroupID.data))
            elif len(form.ProteinID.data) > 0:
                return redirect(url_for('getInteractionsForProteinID', ProteinID = form.ProteinID.data))
 
    #kdrew: added hoping it would fix redirect problem on stale connections
    return render_template('index.html', form = form, complexes = complexes)









# Unused below here
@app.route("/displayComplexesForProtein")
def displayComplexesForProtein():
    protein_search = request.args.get('protein')
    form = SearchForm()
    error = None
    #print protein
    #kdrew: do error checking
    complexes = []
    try:
        proteins = db.session.query(cdb.Protein).filter(cdb.Protein.proteinname.like('%'+protein_search+'%')).all()
        for protein in proteins:
            print protein
            complexes = complexes + protein.complexes

        #kdrew: remove redudant complexes
        complexes = list(set(complexes))

    except NoResultFound:
        complexes = []

    if len(complexes) == 0:
        error = "No complexes found for given search term: %s" % protein_search

    return render_template('index.html', form=form, complexes=complexes, error=error)





#testing viz


#Breaking up functions


@app.route('/finder', methods=['GET'])
def finding():
    """ Query the network
    """
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")

    # Grab the inputs arguments from the URL
    args = request.args
    # Get all the form arguments in the url with defaults
    bait =  getitem(args, 'bait', 'F4JY76_ARATH EIF3K_ARATH B3H7J6_ARATH')
    degree = getitem(args, 'deg', 2)

    #reselect = getitems(args, 'reselect', 1)
 
    reselect =   request.args.getlist('reselect')


    print("test reselect")
    print(reselect)

    bait_list = bait.split(" ")

    #check that proteinID inputs are valid
    bait_list = valid_query(bait_list, conversion_tbl)

    #Identify species of input proteins
    input_protein_species, longform_species = identify_species(bait_list, conversion_tbl)


    #Map from protein to group    
    protein_table, group_table = make_conversion_tables(bait_list, conversion_tbl, input_protein_species)


    #Get groups from a protein list
    #Somewhat duplicated in make_conversion_tables...
 
    #This line only for when there are ortholog groups in the mix  
    bait = get_groups(bait_list, conversion_tbl)

    bait_str = " ".join(bait)

    #Check for check marks ticked
    print("reselect", reselect)
    #checks = getitems('reselect')
    if reselect :
        bait_list = reselect
        bait_str = ' '.join(reselect)
        bait=bait_list

    print(bait) 
    try:
          final_annotated, df_all_prots =run_process(bait, conversion_tbl)
          med_score, mean_score, suggestions = sampling_process(bait)

          suggestion_str = "\n".join(suggestions)
          suggestion_df  = pd.DataFrame({'Suggestions': suggestions})
          print(suggestion_df)
          suggestion_html = suggestion_df.to_html(classes='SuggestionTbl', index=False) 
         

    except Exception as E:
        print(Exception)
        print("Failed to get results")
        med_score = "0"
        mean_score = "0"
        suggestion_html = "No results found"
        final_annotated = pd.DataFrame(columns = ['score','Annotation','Annotation2', 'GroupID_key', 'GroupID_key2', 'bait_bait'])
 
    print("Draw network")

    try:
        print(bait, degree)
        nodes_table, network_script, network_div = pcd.networking(final_annotated, bait, degree, df_all_prots)


        formData = request.values if request.method == "POST" else request.values
        response = "Form Contents <pre>%s</pre>" % "<br/>\n".join(["%s:%s" % item for item in formData.items()] )



    except Exception as e:
        print e
        nodes_table = ""
        results_table = ""
        #results_table = df_all_interactions.to_html(classes='ResultsTbl', index=False) 
 
        network_script = ""
        network_div = ""
    # Create a polynomial line graph with those arguments
    #x = list(range(_from, to + 1))
    #fig = figure(title="Polynomial")
    #fig.line(x, [i ** 2 for i in x], color=color, line_width=2)
    results_table = final_annotated.to_html(classes='tablesorter" id = "my_id', index=False) 
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    bait_plus = "+".join(bait)

    spec_plus = longform_species.replace(" ", "+")

    elution_link = "http://127.0.0.1:5000/complexdisplay?proteins={}&species_longform={}".format(bait_plus, spec_plus)


    #This is temporary until I can get a better network plotter
    network_div = network_div.replace('<div class="bk-plotdiv"', '')
    network_div = network_div.replace('</div>', '', 1)
    network_div = network_div.replace('>', '', 1)
    html = render_template(
        'finder.html',
        #plot_script=script,
        #plot_div=div,
        #exp_link=expression_link,
        elut_link=elution_link,
        net_script=network_script,
        net_div=network_div,
        prot_tbl=protein_table,
        #proteins=proteins,
        nodes_tbl=nodes_table,
        results_tbl=results_table,
        deg=degree,
        js_resources=js_resources,
        median=med_score,
        mean=mean_score,
        suggest=suggestion_html,
        bt=bait_str,
        css_resources=css_resources,
    )
    return encode_utf8(html)

@app.route("/about")
def displayAbout():
    return render_template('about.html')

@app.route("/download")
def displayDownload():
    return render_template('download.html')


if __name__ == "__main__":
    db.create_all()  # make our sqlalchemy tables
    app.run(threaded=True, debug=True)
    

#@app.route('/proteinquery')
#def protein_query():
#    """ Query a list of proteins against out data
#    """
#    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")
#
#    #Save this in a file somewhere
#    spec_conv_dict = {"Arabidopsis thaliana":"arath", "Brassica oleraceae":"braol", "Chlamydomonas reinhardtii":"chlre", "Oryza sativa":"orysj",
#"Selaginella moellendorfii":"selml", "Triticum aestivum":"traes", "Ceratopteris richardii":"cerri", "All plants":"allplants"}  
#
#    # Grab the inputs arguments from the URL
#    args = request.args
#
#    proteins =  getitem(args, 'proteins', "WRK58_ARATH Q9SR92")
#    species_longform =  getitem(args, 'species_longform', "All plants")
#
#
#    print species_longform
#
#    species = spec_conv_dict[species_longform]
#    print species
#    protein_list = proteins.split(" ")
#    protein_list = valid_query(protein_list, conversion_tbl)
#    protein_table, group_table = make_conversion_tables(protein_list, conversion_tbl, species)
#
#    input_protein_species = identify_species(protein_list, conversion_tbl)[0]
#    input_protein_species_longform = [x for x in spec_conv_dict.keys() if spec_conv_dict[x] == input_protein_species][0]
#
#
#    proteins_plus = proteins.replace(" ", "+")
#    species_longform_plus =species_longform.replace(" ", "+")
#    display_link = "http://127.0.0.1:5000/proteinquery?proteins={}&species_longform={}".format(proteins_plus, species_longform_plus)
#     
# 
#    #js and cs needed for the plot
#    js_resources = INLINE.render_js()
#    css_resources = INLINE.render_css()
#
#    #Get plot over to the html webpage
#    html = render_template(
#        'proteinquery.html',
#        js_resources=js_resources,
#        css_resources=css_resources,
#        proteins=proteins,
#        disp_link= display_link 
#    )
#    return encode_utf8(html)
#

#@app.route('/complexdisplay')
#def display_complex():
#    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")
#
#    #Save this in a file somewhere
#    spec_conv_dict = {"Arabidopsis thaliana":"arath", "Brassica oleraceae":"braol", "Chlamydomonas reinhardtii":"chlre", "Oryza sativa":"orysj",
#"Selaginella moellendorfii":"selml", "Triticum aestivum":"traes", "Ceratopteris richardii":"cerri", "All plants":"allplants"}  
#
#    # Grab the inputs arguments from the URL
#    args = request.args
#
#    proteins =  getitem(args, 'proteins', "WRK58_ARATH Q9SR92")
#    species_longform =  getitem(args, 'species_longform', "All plants")
#
#
#    print species_longform
#
#    species = spec_conv_dict[species_longform]
#    print species
#    protein_list = proteins.split(" ")
#    protein_table, group_table = make_conversion_tables(protein_list, conversion_tbl, species)
#
#    input_protein_species = identify_species(protein_list, conversion_tbl)[0]
#    input_protein_species_longform = [x for x in spec_conv_dict.keys() if spec_conv_dict[x] == input_protein_species][0]
#
#    print(species)
#    #There's some redundancy with make_conversion_table pulling groupIDs
#    #Going to be sql'd anyway later
#    complexes = get_groups(protein_list, conversion_tbl)
#    print complexes
#     
# 
#    #js and cs needed for the plot
#    js_resources = INLINE.render_js()
#    css_resources = INLINE.render_css()
#
#    try:
#       results1 = make_protein_sparklines(complexes, species, species_longform, conversion_tbl)
#       script1, div1 = results1
#    except Exception as e:
#        print e
#        print "Something went wrong"
#        script1="No input proteins observed in proteomics data for " + species_longform
#        div1="none"
#
#    #Get plot over to the html webpage
#    html = render_template(
#        'complexdisplay.html',
#        js_resources=js_resources,
#        css_resources=css_resources,
#        complexes=complexes,
#        proteins=proteins,
#        prot_tbl=protein_table,
#        group_tbl=group_table,
#        species_longform=species_longform,
#        inp_species_longform=input_protein_species_longform,
#        plot_script1=script1,
#        plot_div1=div1,
# 
#    )
# 
#    return encode_utf8(html)
#
#@app.route('/sparkline_simple')
#def spark():
#    """ Simple sparkline
#    """
#
#    # Grab the inputs arguments from the URL
#    args = request.args
#
#    # Get all the form arguments in the url with defaults
#    #color = colors[getitem(args, 'color', 'Black')]
#    #_from = int(getitem(args, '_from', 0))
#    #to = int(getitem(args, 'to', 10))
#    complexes = getitem(args, 'complexes', "ENOG410IDXB ENOG410IDXU")
# 
#
#    # Create a polynomial line graph with those arguments
#    #x = list(range(_from, to + 1))
#    #fig = figure(title="Polynomial")
#    #fig.line(x, [i ** 2 for i in x], color=color, line_width=2)
#
#    results = make_sparklines(complexes)
#
#    js_resources = INLINE.render_js()
#    css_resources = INLINE.render_css()
#
#    script, div = results
#    html = render_template(
#        'viewer.html',
#        plot_script=script,
#        plot_div=div,
#        js_resources=js_resources,
#        css_resources=css_resources,
#        complexes=complexes
#    )
#    return encode_utf8(html)
#
#


