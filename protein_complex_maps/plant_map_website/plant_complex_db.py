from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func, or_, and_, UniqueConstraint

app = Flask(__name__)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///plant.db'

app.config['SECRET_KEY'] = 'please, tell nobody'

db = SQLAlchemy(app)

def drop_table(tablename):
    tablename.drop()

def get_db():
    return db

def get_app():
    return app

def get_or_create(db, model, **kwargs):
    session = db.session
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance
    else:
        instance = model(**kwargs)
        #session.add(instance) # unnecessary?
        #session.commit() # unnecessary?
    return instance


class Hiercomplex(db.Model):
    """A single complex"""
    # Each complex contains multiple orthogroups
    # One orthogroup can belong to many complexes
    # Many complexes to Many orthogroups
    id = db.Column(db.Integer, primary_key=True)
    clustid = db.Column(db.Integer, unique = False)
    clustid_set = db.Column(db.String(63))
    orthogroups = db.relationship('Orthogroup', secondary = 'orthogroup_complex_mapping', backref = db.backref('hiercomplexes'), lazy = 'select') # or 'lazy = 'dynamic'
  

class Orthogroup(db.Model):
    """A single group"""
    id = db.Column(db.Integer, primary_key=True) #, autoincrement = True)
    OrthogroupID = db.Column(db.String(63), unique = True)
    __table_args__ = (UniqueConstraint('OrthogroupID'),)


class Protein(db.Model):
    """A single group"""
    # Each protein is only in one orthogroup, but each orthogroup has many proteins
    # One orthogroup to Many proteins
    id = db.Column(db.Integer, primary_key=True) #, autoincrement=True)
    ProteinID = db.Column(db.String(63))
    Spec = db.Column(db.String(63))
    IDtype = db.Column(db.String(63)) 
    #orthogroups = db.relationship('Orthogroup', secondary = 'orthogroup_protein_mapping', backref = db.backref('Proteins'), lazy = 'select') # or 'lazy = 'dynamic'
    #__table_args__ = (UniqueConstraint('ProteinID'),)
    OrthogroupID_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
    orthogroups = db.relationship("Orthogroup", backref = "Proteins", lazy = "select")

class Orthoannot(db.Model):
    """A single group"""
    # One orthoannotation to One orthogroup
    id = db.Column(db.Integer, primary_key=True)
    EggnogAnnot = db.Column(db.String(63))
    Tair = db.Column(db.String(63))
    ArathGenenames = db.Column(db.String(63))
    ArathProtnames = db.Column(db.String(63))
    ArathGO = db.Column(db.String(63))
    Counts = db.Column(db.Integer)
    # ADD ArathProteinNames
    OrthogroupID_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
    orthogroups = db.relationship("Orthogroup", backref = "Orthoannots", lazy = "select")
    #orthogroups = db.relationship('Orthogroup', secondary = 'orthogroup_annot_mapping', backref = db.backref('Orthoannots'), lazy = 'select') # or 'lazy = 'dynamic'
 
class Score(db.Model):
    """A orthogroup orthogroup edge with score"""
    # One orthogroup has Many scores
   # Each score has two orthogroups
    id = db.Column(db.Integer, primary_key=True)
    InteractionID = db.Column(db.Integer)
    OrthogroupID_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
    ScoreVal = db.Column(db.Float) # redundant
    orthogroups = db.relationship("Orthogroup", backref = "Scores", lazy = "select") # varname should be orthogroups

class OrthogroupComplexMapping(db.Model):
    """A mapping between groups and complexes"""
    orthogroup_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id'), primary_key=True)
    hiercomplex_key = db.Column(db.Integer, db.ForeignKey('hiercomplex.id'), primary_key=True)

# To Do delete this class
#class OrthogroupProteinMapping(db.Model):
#    """A mapping between groups and complexes"""
#    orthogroup_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id'), primary_key=True)
#    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'), primary_key=True)


#class OrthogroupAnnotMapping(db.Model):
#    """A mapping between groups and complexes"""
#    orthogroup_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id'), primary_key=True)
#    orthoannot_key = db.Column(db.Integer, db.ForeignKey('orthoannot.id'), primary_key=True)


#    OrthogroupID1_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
#    OrthogroupID2_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
#    ScoreVal = db.Column(db.Float)
     

#class Score(db.Model):
#    """A orthogroup orthogroup edge with score"""
    # One orthogroup has Many scores
    # Each score has two orthogroups
#    id = db.Column(db.Integer, primary_key=True)
#    OrthogroupID1_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
#    OrthogroupID2_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
#    ScoreVal = db.Column(db.Float)
    # Not sure if this will work
#    orthogroup1s = db.relationship("Orthogroup", backref = "Scores", lazy = "select", foreign_keys = 'Score.OrthogroupID1_key')
#    orthogroup2s = db.relationship("Orthogroup", backref = "Scores", lazy = "select", foreign_keys = 'Score.OrthogroupID2_key')

    #Payments = db.relationship('Payment', backref = 'payer', lazy = 'dynamic', foreign_keys = 'Payment.uidPayer')
    #Received = db.realtionship('Payment', backref = 'Receiver', lazy = 'dynamic, foreign_keys = 'Payment.uidReceiver') 

 
#class Edge(db.Model):
#    """A orthogroup orthogroup edge"""
#    id = db.Column(db.Integer, primary_key=True)
#    OrthogroupID_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
#    OrthogroupID_key2 = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
#    in_complex = db.Column(db.Integer)
#    score =db.Column(db.Float)
#
#    evidences = db.relationship('Evidence', lazy='dynamic')
#    #orthogroup_key_genename = db.relationship('Edge', primaryjoin='Edge.orthogroup_key'=='Conversion.genename', foreign_keys='Conversion.genename')
#    #orthogroup_key2_genename = db.relationship('Edge', primaryjoin='Edge.orthogroup_key2'=='Conversion.genename', foreign_keys='Conversion.genename')
#
#    def get_orthogroups(self,):
#        orthogroup1 = db.session.query(Orthogroup).filter(Orthogroup.id==self.orthogroup_key).first()
#        orthogroup2 = db.session.query(Orthogroup).filter(Orthogroup.id==self.orthogroup_key2).first()
#        return (orthogroup1, orthogroup2)
#
#

    


