#!/usr/bin/env python
import random
#from bs4 import BeautifulSoup
#from collections import OrderedDict
import plant_complex_db as cdb
#import pandas as pd
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_

db = cdb.get_db()
app = cdb.get_app()

from flask_wtf import FlaskForm  #Changed from just Form to avoid deprecation warning
from wtforms.fields import StringField, SubmitField, SelectField
from flask import make_response
from flask import render_template
from flask import url_for, redirect, request, jsonify
import re
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

    OrthogroupID = StringField(u'virNOG Orthogroup:')
    ProteinID = StringField(u'Protein ID')
    Species = SelectField(u'Species', choices = species_list, default = 'arath')
    submit = SubmitField(u'Search complexes')
    submitinteractions = SubmitField(u'Search top interactions')
 
   
@app.route("/")
def root(complexes=[]):
    print(complexes)
    form = SearchForm()
    return render_template('index.html', form = form, complexes = complexes)


def OrthogroupQuery(Input_OrthogroupID):
    OrthogroupID = db.session.query(cdb.Orthogroup).filter((func.upper(cdb.Orthogroup.OrthogroupID) == func.upper(Input_OrthogroupID))).first()
    print(OrthogroupID)
    return(OrthogroupID)

def FullTextQuery(Input_FullText):

    OrthogroupIDs = []
    OrthogroupIDs.append(db.session.query(cdb.Orthogroup).filter((func.upper(cdb.Orthoannot.ArathGO.like('%'+func.upper(Input_GO.rstrip("*"))+'%'))).all()))
    #ProteinIDs = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.ProteinID).like('%'+func.upper(Input_ProteinID.rstrip("*"))+'%'))).all()
    if len(ProteinIDs) > 1:
        for prot in ProteinIDs:
           OID = OrthogroupQuery(prot.orthogroups.OrthogroupID)
           if OID.Scores:
                with_scores.append(OID)
                ProteinID = prot
                OrthogroupID = OID
           OrthogroupIDs.append(OID)
    if len(Orthogroups) > 1: 
        return render_template('resolveambiguityortho.html', ogs = OrthogroupIDs, form = form, error = error)

    return render_template('resolveambiguityortho.html', ogs = OrthogroupIDs, form = form, error = error)

def ProteinQuery(Input_ProteinID, Species):
    print(Input_ProteinID)
    if Species == "orysj":
       Input_ProteinID = re.sub('^LOC_', '', Input_ProteinID)
    print(Input_ProteinID)

    ProteinIDs = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.ProteinID) == func.upper(Input_ProteinID))).all()
    OrthogroupIDs = []
    for prot in ProteinIDs:
        OrthogroupID = OrthogroupQuery(prot.orthogroups.OrthogroupID)
        OrthogroupIDs.append(OrthogroupID)
    print(Species)
    OrthogroupIDs = list(set(OrthogroupIDs))
    if len(OrthogroupIDs) > 0:
       return(OrthogroupIDs)
    #If not getting exact matches, try wildcard search
    else:
       ProteinIDs = db.session.query(cdb.Protein).filter(func.upper(cdb.Protein.ProteinID).like('%'+func.upper(Input_ProteinID.rstrip("*"))+'%')).filter(func.upper(cdb.Protein.Spec) == func.upper(Species)).all()
       OrthogroupIDs = []
       for prot in ProteinIDs:
           OrthogroupID = OrthogroupQuery(prot.orthogroups.OrthogroupID)
           OrthogroupIDs.append(OrthogroupID)

    OrthogroupIDs = list(set(OrthogroupIDs))
    print(OrthogroupIDs)
    return(OrthogroupIDs)


def ProteinQuery2(Input_ProteinID, Species):
    #cdb.ComplexEnrichment.t_name.like('%'+enrichment+'%')
    ProteinIDs = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.ProteinID) == func.upper(Input_ProteinID))).all()
    if len(ProteinIDs) == 1:
       return(ProteinIDs)
    else:
       #ProteinIDs = db.session.query(cdb.Protein).filter(func.upper(cdb.Protein.ProteinID).like('%'+func.upper(Input_ProteinID.rstrip("*"))+'%')).filter(cdb.Protein.Spec == Species).all()
       ProteinIDs = db.session.query(cdb.Protein).filter(func.upper(cdb.Protein.ProteinID).like('%'+func.upper(Input_ProteinID.rstrip("*"))+'%')).all()
    return(ProteinIDs)

@app.route("/displayComplexesForProteinID")

def displayComplexesForProteinID():
    Species = request.args.get('Species')
    Input_ProteinID = request.args.get('ProteinID').strip().upper()
    form = SearchForm()
    error=None 
    print(Species)
   
    OrthogroupIDs = ProteinQuery(Input_ProteinID, Species)
    #ProteinIDs = ProteinQuery(Input_ProteinID)
    #ProteinIDs.sort(key=lambda x: x.orthogroups.OrthogroupID)

    if len(OrthogroupIDs) > 1:
        return render_template('resolveambiguityortho.html', ogs = OrthogroupIDs, Species = Species, form = form, error = error)


    # Check quality of Input_ProteinID query


    if not OrthogroupIDs:##[0].orthogroups:
        #kdrew: input ProteinID is not valid, flash message
        error = "Could not find orthogroup for given Protein ID: %s.\n Try using a Uniprot.org Accession" % Input_ProteinID
        return render_template('index.html', form = form, complexes = [], error = error)
    else:
        #ProteinIDs = ProteinQuery(Input_ProteinID, Species)
        #ProteinID = ProteinIDs[0]
        OrthogroupID  = OrthogroupIDs[0]
        #OrthogroupID_string = ProteinID.orthogroups.OrthogroupID
        #OrthogroupID = OrthogroupQuery(OrthogroupID_string)

    #OrthogroupID_string = ProteinID.orthogroups.OrthogroupID

    #Species = ProteinID.Spec
 

    complexes = []
    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
    complexes = orthogroup_clusters.hiercomplexes

    if len(orthogroup_clusters.hiercomplexes) == 0:
        error = "No complexes found for given Protein ID %s and its virNOG orthogroup ID: %s. Try for Protein Interactions" % (Input_ProteinID, OrthogroupID.OrthogroupID)
        #return render_template('index.html', form = form, complexes = [], Species = Species, Query = OrthogroupID, Input_ProteinID = Input_ProteinID, error = error)
       
        
    return render_template('getcomplexes.html', form = form, complexes = complexes, Species = Species, Query = OrthogroupID, Input_ProteinID = Input_ProteinID, error = error)




def displayComplexesForProteinID2():
    Input_ProteinID = request.args.get('ProteinID').strip().upper()
    form = SearchForm()
    error=None 
   
    ProteinIDs = ProteinQuery(Input_ProteinID)
    #ProteinIDs = ProteinQuery(Input_ProteinID)
    ProteinIDs.sort(key=lambda x: x.orthogroups.OrthogroupID)


    # Check quality of Input_ProteinID query
    if len(ProteinIDs) > 1:
        with_scores = []
        for prot in ProteinIDs:
           OID = OrthogroupQuery(prot.orthogroups.OrthogroupID)
           if OID.Scores:
                with_scores.append(OID)
                ProteinID = prot
                OrthogroupID = OID
        if len(with_scores) > 1:
   
            return render_template('resolveambiguity.html', prots = ProteinIDs, form = form, error = error)


    if not ProteinIDs:##[0].orthogroups:
        #kdrew: input ProteinID is not valid, flash message
        error = "Could not find orthogroup for given Protein ID: %s.\n Try using a Uniprot.org Accession" % Input_ProteinID
        return render_template('index.html', form = form, complexes = [], error = error)

    else:
        ProteinID = ProteinIDs[0]
        OrthogroupID_string = ProteinID.orthogroups.OrthogroupID
        OrthogroupID = OrthogroupQuery(OrthogroupID_string)

    OrthogroupID_string = ProteinID.orthogroups.OrthogroupID

    Species = ProteinID.Spec
 

    complexes = []
    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()
    complexes = orthogroup_clusters.hiercomplexes

    if len(orthogroup_clusters.hiercomplexes) == 0:
        error = "No complexes found for given Protein ID %s and its virNOG orthogroup ID: %s. Try for Protein Interactions" % (Input_ProteinID, OrthogroupID.OrthogroupID)
        #return render_template('index.html', form = form, complexes = [], Species = Species, Query = OrthogroupID, Input_ProteinID = Input_ProteinID, error = error)
       
        
    return render_template('getcomplexes.html', form = form, complexes = complexes, Species = Species, Query = OrthogroupID, Input_ProteinID = Input_ProteinID, error = error)


@app.route("/displayComplexesForOrthogroupID")
def displayComplexesForOrthogroupID():
    Species = request.args.get('Species')
    Input_OrthogroupID = request.args.get('OrthogroupID').strip().upper()
    form = SearchForm()
    error=None
    print(Species)
    #CDM: See if orthogroup is a valid orthogroup ID
    OrthogroupID = OrthogroupQuery(Input_OrthogroupID)

    if not OrthogroupID:
        #kdrew: input genename is not valid, flash message
        error = "Could not find given virNOG orthogroup ID: %s" % Input_OrthogroupID

        return render_template('index.html', form = form, complexes = [], error = error, Species = Species)
    complexes = []
    orthogroup_clusters = (db.session.query(cdb.Orthogroup).filter(cdb.Orthogroup.id == OrthogroupID.id)).first()

    # Get hier complexes associated with an orthogroup ID, Then get member proteins of each cluster level
    complexes = orthogroup_clusters.hiercomplexes
    
    if len(orthogroup_clusters.hiercomplexes) == 0:
        error = "No complexes found for given virNOG orthogroup ID: %s" % Input_OrthogroupID

    return render_template('getcomplexes.html', form = form, complexes = complexes, Species = Species, Query = OrthogroupID, current_OrthogroupID = Input_OrthogroupID, Input_OrthogroupID = Input_OrthogroupID, error = error)

@app.route("/getInteractionsForOrthogroupID")
def getInteractionsForOrthogroupID():
    Species = request.args.get('Species')
    Input_OrthogroupID = request.args.get('OrthogroupID').strip().upper()
    form = SearchForm()
    error=None

    OrthogroupID = OrthogroupQuery(Input_OrthogroupID)
    if not OrthogroupID:
        #kdrew: input genename is not valid, flash message
        error = "Could not find given virNOG orthogroup ID: %s" % Input_OrthogroupID
        return render_template('index.html', form = form, complexes = [], error = error, Species = Species)

    interactions = []
    for score in OrthogroupID.Scores:
        interaction = db.session.query(cdb.Score).filter(cdb.Score.InteractionID == score.InteractionID).all()
        interactions.append(interaction)
    return render_template('getinteractions.html', form = form, interactions = interactions,  Species = Species, Query = OrthogroupID, Input_OrthogroupID = Input_OrthogroupID, current_OrthogroupID = Input_OrthogroupID, error = error)


@app.route("/getInteractionsForProteinID")
def getInteractionsForProteinID():
    Input_ProteinID = request.args.get('ProteinID').strip().upper()
    form = SearchForm()
    error=None

    ProteinIDs = ProteinQuery(Input_ProteinID)

    # Check quality of Input_ProteinID query
    if len(ProteinIDs) > 1:
        with_scores = []
        for prot in ProteinIDs:
           OID = OrthogroupQuery(prot.orthogroups.OrthogroupID)
           if OID.Scores:
                with_scores.append(OID)
                ProteinID = prot
                OrthogroupID = OID
        if len(with_scores) > 1:
            
            return render_template('resolveambiguity.html', prots = ProteinIDs, form = form, error = error)

    if not ProteinIDs:##[0].orthogroups:
        #kdrew: input ProteinID is not valid, flash message
        error = "Could not find orthogroup for given Protein ID: %s.\n Try using a Uniprot.org Accession" % Input_ProteinID
        return render_template('index.html', form = form, complexes = [], error = error)

    else:
        ProteinID = ProteinIDs[0]
        OrthogroupID_string = ProteinID.orthogroups.OrthogroupID
        OrthogroupID = OrthogroupQuery(OrthogroupID_string)

    Species = ProteinID.Spec

    interactions = []
    if not OrthogroupID.Scores:
        error = "No strong interactions found for given protein ID: %s" % Input_ProteinID
        return render_template('index.html', form = form, complexes = [], error = error)
      
    else:   
        for score in OrthogroupID.Scores:
            interaction = db.session.query(cdb.Score).filter(cdb.Score.InteractionID == score.InteractionID).all()
            interactions.append(interaction)
        return render_template('getinteractions.html', form = form, interactions = interactions,  Species = Species, Query = OrthogroupID, Input_ProteinID = Input_ProteinID, current_OrthogroupID = OrthogroupID_string, error = error)



@app.route(u'/search', methods = ['POST'])
def searchComplexes():
    form = SearchForm()
    complexes = []
    
    if form.validate_on_submit():
 
        if form.submit.data == True:
            if len(form.OrthogroupID.data) > 0:
                if len(form.Species.data) > 0:
                   return redirect(url_for('displayComplexesForOrthogroupID', OrthogroupID = form.OrthogroupID.data, Species = form.Species.data))
    
                else:
                   return redirect(url_for('displayComplexesForOrthogroupID', OrthogroupID = form.OrthogroupID.data))
            elif len(form.ProteinID.data) > 0:
                return redirect(url_for('displayComplexesForProteinID', ProteinID = form.ProteinID.data, Species = form.Species.data))
        
        if form.submitinteractions.data == True:
            if len(form.OrthogroupID.data) > 0:
                if len(form.Species.data) > 0:
                   return redirect(url_for('getInteractionsForOrthogroupID', OrthogroupID = form.OrthogroupID.data, Species = form.Species.data))
    
                else:
                   return redirect(url_for('getInteractionsForOrthogroupID', OrthogroupID = form.OrthogroupID.data))
            elif len(form.ProteinID.data) > 0:
                return redirect(url_for('getInteractionsForProteinID', ProteinID = form.ProteinID.data, Species = form.Species.data))
 
    #kdrew: added hoping it would fix redirect problem on stale connections
    return render_template('index.html', form = form, complexes = complexes)



@app.route("/about")
def displayAbout():
    return render_template('about.html')

@app.route("/download")
def displayDownload():
    return render_template('download.html')


if __name__ == "__main__":
    db.create_all()  # make our sqlalchemy tables
    app.run(threaded=True, debug=True)
 

