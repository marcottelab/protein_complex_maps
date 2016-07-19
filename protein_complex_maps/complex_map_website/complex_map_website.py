
import protein_complex_maps.complex_map_website.complex_db as cdb
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_

db = cdb.get_db()
app = cdb.get_app()

#from flask.ext.wtf import Form
from flask_wtf import Form
from wtforms.fields import StringField, SubmitField


class SearchForm(Form):
    complex_id = StringField(u'Complex ID:')
    genename = StringField(u'Gene Name (ex. OFD1):')
    enrichment = StringField(u'Enrichment (ex. cilium):')
    protein = StringField(u'Protein (ex. Centrosomal protein):')
    submit = SubmitField(u'Search')

from flask import render_template
from flask import url_for, redirect, request

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
    try:
        protein = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.genename)==func.upper(genename))).one()
        complexes = protein.complexes.all()
    except NoResultFound:
        complexes = []

    return render_template('index.html', form=form, complexes=complexes)

@app.route("/displayComplexesForEnrichment")
def displayComplexesForEnrichment():
    enrichment = request.args.get('enrichment')
    form = SearchForm()
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
        print "NoResultFound"
        complexes = []

    return render_template('index.html', form=form, complexes=complexes)

@app.route("/displayComplexesForProtein")
def displayComplexesForProtein():
    protein_search = request.args.get('protein')
    form = SearchForm()
    #print protein
    #kdrew: do error checking
    complexes = []
    try:
        proteins = db.session.query(cdb.Protein).filter(cdb.Protein.proteinname.like('%'+protein_search+'%')).all()
        for protein in proteins:
            print protein
            complexes = complexes + protein.complexes.all()

        #kdrew: remove redudant complexes
        complexes = list(set(complexes))

    except NoResultFound:
        print "NoResultFound"
        complexes = []

    return render_template('index.html', form=form, complexes=complexes)

@app.route("/displayComplexes")
def displayComplexes():
    complex_key = request.args.get('complex_key')
    form = SearchForm()
    #kdrew: do error checking
    try:
        complexes = db.session.query(cdb.Complex).filter_by(complex_id=complex_key).all()
    except NoResultFound:
        complexes = []

    return render_template('complex.html', form=form, complexes=complexes)


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



if __name__ == "__main__":
    db.create_all()  # make our sqlalchemy tables
    app.run()


