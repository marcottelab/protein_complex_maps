import random
from bs4 import BeautifulSoup
from collections import OrderedDict
#import protein_complex_maps.plant_map_website.complex_db as cdb
import complex_db as cdb
import pandas as pd
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_
import networkx as nx

db = cdb.get_db()
app = cdb.get_app()

#from flask.ext.wtf import Form
from flask_wtf import Form
from wtforms.fields import StringField, SubmitField
#testing viz
from model import InputForm
from compute import compute
from flask import make_response
from lineplot import make_protein_sparklines
from validate_query import valid_query
from get_species import identify_species
from get_groups_from_prots import get_groups
from make_conv_tables import make_conversion_tables
from get_distributions import sampling_process, run_process, annotate_nodes
import plot_corum_dists as pcd

from bokeh.resources import INLINE
from bokeh.util.string import encode_utf8

#PLOT_OPTIONS = dict(plot_width=800, plot_height=300)
#SCATTER_OPTIONS = dict(size=12, alpha=0.5)

colors = {
    'Black': '#000000',
    'Red':   '#FF0000',
    'Green': '#00FF00',
    'Blue':  '#0000FF',
}

def getitem(obj, item, default):
    if item not in obj:
        return default
    else:
        return obj[item]


class SearchForm(Form):
    complex_id = StringField(u'Complex ID:')
    genename = StringField(u'Gene Name (ex. OFD1):')
    enrichment = StringField(u'Enrichment (ex. cilium):')
    protein = StringField(u'Protein (ex. Centrosomal protein):')
    submit = SubmitField(u'Search')

from flask import render_template
from flask import url_for, redirect, request, jsonify


@app.route("/")
def root(complexes=[]):
    print complexes
    #complexes = cdb.Complex.query.all()
    form = SearchForm()
    return render_template('index.html', form=form, complexes=complexes)

@app.route("/displayComplexesForGeneName")
def displayComplexesForGeneName():
    genename = request.args.get('genename')
    form = SearchForm()
    #kdrew: do error checking
    error=None

    #kdrew: tests to see if genename is a valid genename
    #protein = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.genename) == func.upper(genename))).one()
    genes = db.session.query(cdb.Gene).filter((func.upper(cdb.Gene.genename) == func.upper(genename))).all()

    if len(genes) == 0:
        #kdrew: input genename is not valid, flash message
        error = "Could not find given genename: %s" % genename

        return render_template('index.html', form=form, complexes=[], error=error)


    complexes = []
    for gene in genes:
        try:
            proteins = db.session.query(cdb.Protein).filter((cdb.Protein.gene_id == gene.gene_id)).all()

        except NoResultFound:
            #kdrew: input genename is not valid, flash message
            error = "Could not find given genename: %s" % genename

            return render_template('index.html', form=form, complexes=[], error=error)

        for protein in proteins:
            try:
                complexes = complexes + protein.complexes
            except NoResultFound:
                continue

    if len(complexes) == 0:
        error = "No complexes found for given genename: %s" % genename

    complexes = list(set(complexes))

    return render_template('index.html', form=form, complexes=complexes, error=error)

@app.route("/displayComplexesForEnrichment")
def displayComplexesForEnrichment():
    enrichment = request.args.get('enrichment')
    form = SearchForm()
    error=None
    #print enrichment
    #kdrew: do error checking
    try:
        enrichment_complex_keys_query = db.session.query(cdb.ComplexEnrichment.complex_key).filter(
                    ( cdb.ComplexEnrichment.t_name.like('%'+enrichment+'%')) | (cdb.ComplexEnrichment.term_id.like('%'+enrichment+'%' ) ) )
        enrichment_complex_keys = enrichment_complex_keys_query.all()
        enrichment_complex_keys_set = set([x[0] for x in enrichment_complex_keys])
        if len(enrichment_complex_keys_set) == 0:
            complexes = []
        else:
            complexes = db.session.query(cdb.Complex).filter(cdb.Complex.id.in_(enrichment_complex_keys_set)).all()
    except NoResultFound:
        complexes = []

    if len(complexes) == 0:
        error = "No complexes found for given enrichment term: %s" % enrichment

    return render_template('index.html', form=form, complexes=complexes, error=error)

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

@app.route("/displayComplexes")
def displayComplexes():
    complex_key = request.args.get('complex_key')
    form = SearchForm()
    error=None
    #kdrew: do error checking
    try:
        comp = db.session.query(cdb.Complex).filter_by(complex_id=complex_key).one()
    except NoResultFound:
        comp = None

    if comp == None:
        error = "No complexes found: %s" % complex_key

    return render_template('complex.html', form=form, comp=comp, error=error)


@app.route(u'/search', methods=[u'POST'])
def searchComplexes():
    form = SearchForm()
    complexes = []
    if form.validate_on_submit():
        if len(form.genename.data) > 0:
            return redirect(url_for('displayComplexesForGeneName', genename=form.genename.data))
        elif len(form.enrichment.data) > 0:
            return redirect(url_for('displayComplexesForEnrichment', enrichment=form.enrichment.data))
        elif len(form.protein.data) > 0:
            return redirect(url_for('displayComplexesForProtein', protein=form.protein.data))


    #kdrew: added hoping it would fix redirect problem on stale connections
    return render_template('index.html', form=form, complexes=complexes)


#testing viz


@app.route('/finder', methods=['GET', 'POST'])
def polynomial():
    """ Very simple embedding of a polynomial chart
    """
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")
    # Grab the inputs arguments from the URL
    args = request.args
    # Get all the form arguments in the url with defaults
    proteins =  getitem(args, 'proteins', 'F4JY76_ARATH EIF3K_ARATH B3H7J6_ARATH')
    degree = getitem(args, 'deg', 2)
    name1 = getitem(args, 'name1', 1)
    print name1
    protein_list = proteins.split(" ")
    protein_list = valid_query(protein_list, conversion_tbl)
    input_protein_species = identify_species(protein_list, conversion_tbl)
    protein_table, group_table = make_conversion_tables(protein_list, conversion_tbl, input_protein_species)
    complexes = get_groups(protein_list, conversion_tbl)
    print "COMPLEXES"   
    print complexes
    complex_str = " ".join(complexes)
    checks = request.form.getlist('name1')
    print checks
    if checks:
        complexes = checks  
        print complexes

      
    
    final_annotated, df_all_prots =run_process(complexes, conversion_tbl)
    med_score, mean_score, suggestions = sampling_process(complexes)

    suggestion_str = "\n".join(suggestions)
   
    print("Draw network")
    #print(bait_group_annotations)
    #full_network = final_annotated[['Annotation','Annotation2','score']]
    #G=pcd.get_network(full_network, 'Annotation', 'Annotation2')
    #pcd.draw_network(G, bait_group_annotations, 2)

    try:

        full_network = final_annotated[['GroupID_key','GroupID_key2','score']]
        G = pcd.get_network(full_network, 'GroupID_key', 'GroupID_key2')


        print(G.degree())
        nodes = pd.DataFrame.from_dict(G.degree(), orient='index')
        nodes.columns = ['Degree']
        print nodes
        annots = pd.read_csv("all_annotations.csv")
        nodes_annotated = annotate_nodes(nodes, annots, complexes, df_all_prots)

        #print(nodes_annotated)
       
        nodes_annotated = nodes_annotated.sort_values(['Degree', 'score'], ascending=False)
        nodes_annotated = nodes_annotated.reset_index()

        #nodes_annotated_dict = nodes_annotated[['index', 'Annotation']].to_dict()

        nodes_annotated_dict =  dict(zip(nodes_annotated['index'], nodes_annotated['Annotation']))
        #print(nodes_annotated_dict)

        #print(nodes_annotated)
        nodes_table = nodes_annotated.to_html(classes='ResultsTbl', index=False)      
        nodes_soup = BeautifulSoup(nodes_table, 'html.parser')
        #print nodes_soup
        #print "SOUP"
        #print nodes_soup
        result = ""
        result = result + '<table border="1" class="dataframe ResultsTbl"><thead><tr style="text-align: right;"><th>index</th><th>Degree</th><th>Code</th><th>Annotation</th><th>Bait</th><th>score</th><th>select</th></tr></thead><tbody>'

        rows = nodes_soup.findAll('tr')
        header=True
        for row in rows:
            cols = row.findAll('td')
            #print cols
            n = 0
            #if header==True:
               #to_add = "<th>Select</th>"
               #result = result + to_add
               #header =  False
 
            result = result  + "<tr>"     
            for col in cols:
                result = result + str(col)
                if n ==0:
                     vallabel = col.text
                if n == 5:
                    #if header== True:
                    #    to_add = "<th>Select</th>"
                    #    result = result + to_add
                    #    header =  False
                    #else:
                     #print dir(col) 
                     #vallabel = col.text
                     if vallabel in complexes:
                         add_check = ' checked="Checked"' 
                     else:
                         add_check = ''
                     print vallabel
                     to_add = '<td> <input type="checkbox" name="name1" value="' + vallabel + '"' + add_check + ' />&nbsp; </td>'
                     result = result + to_add +"</tr>"
                n = n + 1           
        result = result + '</tbody></table>'
        #print result
        nodes_table = result

        G = pcd.filter_nodes(G, degree)
        weights,pos, colorvalues =  pcd.draw_network(G, complexes, degree)   
        #print weights, pos, colorvalues

        network_script, network_div = pcd.bokeh_network(G, weights, pos, colorvalues, nodes_annotated_dict)
    
       # print request.form['reselect_query']   
        ##
        #if request.method == 'POST':

        formData = request.values if request.method == "POST" else request.values
        response = "Form Contents <pre>%s</pre>" % "<br/>\n".join(["%s:%s" % item for item in formData.items()] )

        print response
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
    results_table = final_annotated.to_html(classes='ResultsTbl', index=False) 
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    html = render_template(
        'finder.html',
        #plot_script=script,
        #plot_div=div,
        net_script=network_script,
        net_div=network_div,
        prot_tbl=protein_table,
        nodes_tbl=nodes_table,
        results_tbl=results_table,
        deg=degree,
        js_resources=js_resources,
        median=med_score,
        mean=mean_score,
        suggest=suggestion_str,
        bait=complex_str,
        css_resources=css_resources,
    )
    return encode_utf8(html)

@app.route('/viewer')
def spark():
    """ Simple sparkline
    """

    # Grab the inputs arguments from the URL
    args = request.args

    # Get all the form arguments in the url with defaults
    #color = colors[getitem(args, 'color', 'Black')]
    #_from = int(getitem(args, '_from', 0))
    #to = int(getitem(args, 'to', 10))
    complexes = getitem(args, 'complexes', "ENOG410IDXB ENOG410IDXU")
 

    # Create a polynomial line graph with those arguments
    #x = list(range(_from, to + 1))
    #fig = figure(title="Polynomial")
    #fig.line(x, [i ** 2 for i in x], color=color, line_width=2)

    results = make_sparklines(complexes)

    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    script, div = results
    html = render_template(
        'viewer.html',
        plot_script=script,
        plot_div=div,
        js_resources=js_resources,
        css_resources=css_resources,
        complexes=complexes
    )
    return encode_utf8(html)

@app.route('/proteinquery')
def protein_query():
    """ Query a list of proteins against out data
    """
    conversion_tbl = pd.read_csv("all_tophits_protlength.txt", sep="\t")

    #Save this in a file somewhere
    spec_conv_dict = {"Arabidopsis thaliana":"arath", "Brassica oleraceae":"braol", "Chlamydomonas reinhardtii":"chlre", "Oryza sativa":"orysj",
"Selaginella moellendorfii":"selml", "Triticum aestivum":"traes", "Ceratopteris richardii":"cerri", "All plants":"All plants"}  

    # Grab the inputs arguments from the URL
    args = request.args

    proteins =  getitem(args, 'proteins', "WRK58_ARATH Q9SR92")
    species_longform =  getitem(args, 'species_longform', "All plants")


    print species_longform

    species = spec_conv_dict[species_longform]
    print species
    protein_list = proteins.split(" ")
    protein_list = valid_query(protein_list, conversion_tbl)
    protein_table, group_table = make_conversion_tables(protein_list, conversion_tbl, species)

    input_protein_species = identify_species(protein_list, conversion_tbl)
    input_protein_species_longform = [x for x in spec_conv_dict.keys() if spec_conv_dict[x] == input_protein_species][0]

    print(species)
    #There's some redundancy with make_conversion_table pulling groupIDs
    #Going to be sql'd anyway later
    complexes = get_groups(protein_list, conversion_tbl)
    print complexes
     
 
    #js and cs needed for the plot
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    try:
       results1 = make_protein_sparklines(complexes, species, species_longform, conversion_tbl)
       script1, div1 = results1
    except Exception as e:
        print e
        print "Something went wrong"
        script1="No input proteins observed in proteomics data for " + species_longform
        div1="none"

    #Get plot over to the html webpage
    html = render_template(
        'proteinquery.html',
        js_resources=js_resources,
        css_resources=css_resources,
        complexes=complexes,
        proteins=proteins,
        prot_tbl=protein_table,
        group_tbl=group_table,
        species_longform=species_longform,
        inp_species_longform=input_protein_species_longform,
        plot_script1=script1,
        plot_div1=div1,
 
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
    

