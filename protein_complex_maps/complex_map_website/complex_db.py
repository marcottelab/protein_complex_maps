
from flask import Flask
#from flask.ext.sqlalchemy import SQLAlchemy
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import func, or_, and_

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'

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


class Conversion(db.Model):
    """A mapping between different ID types"""
    id = db.Column(db.Integer, primary_key=True)
    groupname = db.Column(db.String(63))
    genename = db.Column(db.String(63)) 
    uniprot_acc = db.Column(db.String(63))
    group_id = db.Column(db.String(63)) 
    group_annotation = db.Column(db.String(63)) 
    entryname = db.Column(db.String(63)) 

       

class Group(db.Model):
    """A single group"""
    id = db.Column(db.Integer, primary_key=True)
    group_id = db.Column(db.String(63))
    uniprot_acc = db.Column(db.String(63))
    genename = db.Column(db.String(255))
    proteinname = db.Column(db.String(255))
    uniprot_url = db.Column(db.String(255))
    #complexes = db.relationship('Complex', secondary='group_complex_mapping',  back_populates='groups', lazy='dynamic')

    def uniprot_link(self,):
        if self.genename == "":
            retstr = self.group_id
        else:
            retstr = "<a href=%s target=\"_blank\">%s</a>" % (self.uniprot_url, self.genename)
        return retstr

class Edge(db.Model):
    """A group group edge"""
    id = db.Column(db.Integer, primary_key=True)
    group_key = db.Column(db.Integer, db.ForeignKey('group.id') )
    group_key2 = db.Column(db.Integer, db.ForeignKey('group.id') )
    in_complex = db.Column(db.Integer)
    score =db.Column(db.Float)

    #evidences = db.relationship('Evidence', lazy='dynamic')
    #group_key_genename = db.relationship('Edge', primaryjoin='Edge.group_key'=='Conversion.genename', foreign_keys='Conversion.genename')
    #group_key2_genename = db.relationship('Edge', primaryjoin='Edge.group_key2'=='Conversion.genename', foreign_keys='Conversion.genename')

    def get_groups(self,):
        grp1 = db.session.query(Group).filter(Group.id==self.group_key).first()
        grp2 = db.session.query(Group).filter(Group.id==self.group_key2).first()
        return (prot1, prot2)


