
import protein_complex_maps.complex_map_website.complex_db as cdb
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_

import mpmath as mpm
#import scipy.misc as misc

db = cdb.get_db()
app = cdb.get_app()

#from flask.ext.wtf import Form
from flask_wtf import Form
from wtforms.fields import StringField, SubmitField, TextAreaField

#kdrew: from http://stackoverflow.com/questions/33468821/is-scipy-misc-comb-faster-than-an-ad-hoc-binomial-computation
def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

#kdrew: hypergeometric test, k = overlap, n = total number of genes input, m = total number of genes in complex, N = total number of genes in all complexes
def pval(k,n,m,N):
    pv = 0.0                                                 
    #N_choose_m = 1.0*misc.comb(N,m)
    for i in range(k,int(min(m,n)+1)):
        pi = ( mpm.binomial(n,i) * mpm.binomial((N-n), (m-i)) ) / mpm.binomial(N,m)
        #pi = ( misc.comb(n,i) * misc.comb((N-n), (m-i)) ) /  N_choose_m
        #kdrew: somethings wrong with this implementation, returns 0 always (or the cases I've tried)
        #pi = ( choose(n,i) * choose((N-n), (m-i)) ) / choose(N,m)
        pv += pi
    return pv


class SearchForm(Form):
    complex_id = StringField(u'Complex ID:')
    #genename = StringField(u'Gene Name (ex. OFD1):')
    listOfGenenames = TextAreaField(u'List of Gene Names (ex. OFD1 PCM1 CSPP1):')
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

#@app.route("/displayComplexesForGeneName")
#def displayComplexesForGeneName():
#    genename = request.args.get('genename')
#    form = SearchForm()
#
#    complexes, error = getComplexesForGeneName(genename)
#    return render_template('index.html', form=form, complexes=complexes, error=error)
#
def getComplexesForGeneName(genename):
    #kdrew: do error checking
    error = ""
    complexes = []

    #kdrew: tests to see if genename is a valid genename
    #protein = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.genename) == func.upper(genename))).one()
    genes = db.session.query(cdb.Gene).filter((func.upper(cdb.Gene.genename) == func.upper(genename))).all()

    if len(genes) == 0:
        #kdrew: input genename is not valid, flash message
        error = "Could not find given genename: %s" % genename
        return complexes, error

    for gene in genes:
        try:
            proteins = db.session.query(cdb.Protein).filter((cdb.Protein.gene_id == gene.gene_id)).all()
        except NoResultFound:
            #kdrew: input genename is not valid, flash message
            error = error + "Could not find given genename: %s" % genename

        for protein in proteins:
            try:
                complexes = complexes + protein.complexes
            except NoResultFound:
                continue

    if len(complexes) == 0:
        error = "No complexes found for given genename: %s" % genename

    complexes = list(set(complexes))

    return complexes, error

@app.route("/displayComplexesForListOfGeneNames")
def displayComplexesForListOfGeneNames():
    listOfGenenames = request.args.get('listOfGenenames').split()
    form = SearchForm()
    #kdrew: do error checking
    error=None

    #print listOfGenenames

    all_genes = []
    for genename in listOfGenenames:
        #print genename
        #kdrew: tests to see if genename is a valid genename
        genes = db.session.query(cdb.Gene).filter((func.upper(cdb.Gene.genename) == func.upper(genename))).all()

        if len(genes) == 0:
            #kdrew: input genename is not valid, flash message
            if error == None:
                error = ""
            error = error + "Could not find genename: %s<br>" % genename

            #return render_template('index.html', form=form, complexes=[], error=error)

        all_genes = all_genes + genes

    #print [g.genename for g in all_genes]

    complexes = []
    all_proteins = []
    error_proteins = []
    for gene in all_genes:
        try:
            proteins = db.session.query(cdb.Protein).filter((cdb.Protein.gene_id == gene.gene_id)).all()

        except NoResultFound:
            if error == None:
                error = ""
            #kdrew: input genename is not valid, flash message
            error = error + "Could not find genename: %s<br>" % gene.genename

            #return render_template('index.html', form=form, complexes=[], error=error)

        all_proteins = all_proteins + proteins

        error_proteins_current = []
        for protein in proteins:
            if len(protein.complexes) == 0:
                error_proteins_current.append(protein)
            try:
                complexes = complexes + protein.complexes
                #kdrew: if found complexes for any of the proteins attached to gene then do not report error
                error_proteins_current = []
            except NoResultFound:
                continue

        #kdrew: add proteins from this iteration that did not have complexes to the list of error proteins
        error_proteins = error_proteins + error_proteins_current

    if len(complexes) == 0:
        if error == None:
            error = ""
        error = error + "No complexes found for genenames: %s<br>" % ', '.join([g.genename for g in all_genes])
    if len(error_proteins) > 0:
        if error == None:
            error = ""
        error = error + "No complexes found for genenames: %s<br>" % ', '.join([p.genename() for p in error_proteins])


    n = len(all_proteins)
    N = db.session.query(cdb.Protein).distinct(cdb.Protein.gene_id).count()
    pvalue_dict = dict()
    for c in set(complexes):
        k = complexes.count(c)
        m = len(c.proteins) 
        print "complex: %s k: %s n: %s m: %s N: %s" % (c.complex_id,k,n,m,N)
        pvalue = pval(k=k,n=n,m=m,N=N)
        pvalue_dict[c] = pvalue
        

    #complexes = list(set(complexes))
    #complexes = [x[1] for x in sorted(((complexes.count(e), e) for e in set(complexes)), reverse=True)]
    complexes = [x[1] for x in sorted((((complexes.count(e), -1*e.top_rank), e) for e in set(complexes)), reverse=True)]



    #print [p.id for p in all_proteins]
    return render_template('index.html', form=form, complexes=complexes, prot_ids=[p.id for p in all_proteins], pvalue_dict=pvalue_dict, error=error)

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
        #if len(form.genename.data) > 0:
        #    return redirect(url_for('displayComplexesForGeneName', genename=form.genename.data))
        #elif len(form.listOfGenenames.data) > 0:
        if len(form.listOfGenenames.data) > 0:
            return redirect(url_for('displayComplexesForListOfGeneNames', listOfGenenames=form.listOfGenenames.data))
        elif len(form.enrichment.data) > 0:
            return redirect(url_for('displayComplexesForEnrichment', enrichment=form.enrichment.data))
        elif len(form.protein.data) > 0:
            return redirect(url_for('displayComplexesForProtein', protein=form.protein.data))


    #kdrew: added hoping it would fix redirect problem on stale connections
    return render_template('index.html', form=form, complexes=complexes)

@app.route("/about")
def displayAbout():
    return render_template('about.html')

@app.route("/download")
def displayDownload():
    return render_template('download.html')



if __name__ == "__main__":
    db.create_all()  # make our sqlalchemy tables
    app.run(threaded=True)


