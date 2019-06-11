
from flask import Flask
#from flask.ext.sqlalchemy import SQLAlchemy
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func, or_, and_

app = Flask(__name__)
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
        session.add(instance)
        session.commit()
    return instance


class Hiercomplex(db.Model):
    """A single complex"""
    id = db.Column(db.Integer, primary_key=True)
    clustid_1 = db.Column(db.Integer, unique = False)
    clustid_2 = db.Column(db.Integer, unique = False)
    clustid_3 = db.Column(db.Integer, unique = False)
    clustid_4 = db.Column(db.Integer, unique = False)
   

class Conversion(db.Model):
    """A mapping between different OrthogroupID types"""
    id = db.Column(db.Integer, primary_key=True)
    OrthogroupID = db.Column(db.String(63))
    Spec = db.Column(db.String(63))
    ProteinID = db.Column(db.String(63)) 

class Orthogroup(db.Model):
    """A single group"""
    id = db.Column(db.Integer, primary_key=True)
    OrthogroupID = db.Column(db.String(63))
    #complexes = db.relationship('Complex', secondary='protein_complex_mapping',  back_populates='groups', lazy='dynamic')

#class Edge(db.Model):
#    """A group group edge"""
#    id = db.Column(db.Integer, primary_key=True)
#    OrthogroupID_key = db.Column(db.Integer, db.ForeignKey('group.id') )
#    OrthogroupID_key2 = db.Column(db.Integer, db.ForeignKey('group.id') )
#    in_complex = db.Column(db.Integer)
#    score =db.Column(db.Float)
#
#    evidences = db.relationship('Evidence', lazy='dynamic')
#    #group_key_genename = db.relationship('Edge', primaryjoin='Edge.group_key'=='Conversion.genename', foreign_keys='Conversion.genename')
#    #group_key2_genename = db.relationship('Edge', primaryjoin='Edge.group_key2'=='Conversion.genename', foreign_keys='Conversion.genename')
#
#    def get_groups(self,):
#        group1 = db.session.query(Orthogroup).filter(Orthogroup.id==self.group_key).first()
#        group2 = db.session.query(Orthogroup).filter(Orthogroup.id==self.group_key2).first()
#        return (group1, group2)
#
#

    


class OrthogroupComplexMapping(db.Model):
    """A mapping between groups and complexes"""
    orthogroup_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id'), primary_key=True)
    hiercomplex_key = db.Column(db.Integer, db.ForeignKey('hiercomplex.id'), primary_key=True)


