
import protein_complex_maps.complex_map_website.complex_db as cdb

db = cdb.get_db()
app = cdb.get_app()

from flask.ext.wtf import Form
from wtforms.fields import StringField, SubmitField


class SearchForm(Form):
    complex_id = StringField(u'Complex ID:')
    genename = StringField(u'Gene Name:')
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
    print "in displayComplexes"
    genename = request.args.get('genename')
    form = SearchForm()

    protein = db.session.query(cdb.Protein).filter_by(genename=genename).one()
    complexes = protein.complexes.all()

    return render_template('index.html', form=form, complexes=complexes)


@app.route(u'/search', methods=[u'POST'])
def searchComplexes():
    form = SearchForm()
    complexes = []
    if form.validate_on_submit():
        return redirect(url_for('displayComplexesForGeneName', genename=form.genename.data))


if __name__ == "__main__":
    db.create_all()  # make our sqlalchemy tables
    app.run()


