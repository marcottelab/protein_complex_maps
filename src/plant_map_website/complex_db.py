
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


class Complex(db.Model):
    """A single complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_id = db.Column(db.Integer, unique=True)
    #groups = db.relationship('Group', secondary='protein_complex_mapping', back_populates='complexes', lazy='dynamic')
    enrichments = db.relationship('ComplexEnrichment', lazy='dynamic')

    def complex_link(self,):
        retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.complex_id, self.complex_id)
        return retstr

    def edges(self,):
        es = []
        pairs = set()
        for group1 in self.groups:
            for group2 in self.groups:
                #kdrew: check to see if the pair has already been entered
                if frozenset((group1,group2)) not in pairs:
                    edge = db.session.query(Edge).filter((and_(Edge.group_key == group1.id, Edge.group_key2 == group2.id) | and_(Edge.group_key == group2.id,Edge.group_key2 == group1.id))).first()
                    #es = es + edge
                    if edge != None:
                        es.append(edge)
                        #kdrew: add pair of prot ids to keep track of edges already added
                        pairs.add(frozenset((group1,group2)))

        return sorted(list(set(es)), key=lambda es: es.score, reverse=True)

class Conversion(db.Model):
    """A mapping between different GroupID types"""
    id = db.Column(db.Integer, primary_key=True)
    GroupID = db.Column(db.String(63))
    Species = db.Column(db.String(63))
    ProteinID = db.Column(db.String(63)) 
  

class Group(db.Model):
    """A single group"""
    id = db.Column(db.Integer, primary_key=True)
    GroupID = db.Column(db.String(63))
    #complexes = db.relationship('Complex', secondary='protein_complex_mapping',  back_populates='groups', lazy='dynamic')

class Edge(db.Model):
    """A group group edge"""
    id = db.Column(db.Integer, primary_key=True)
    GroupID_key = db.Column(db.Integer, db.ForeignKey('group.id') )
    GroupID_key2 = db.Column(db.Integer, db.ForeignKey('group.id') )
    in_complex = db.Column(db.Integer)
    score =db.Column(db.Float)

    evidences = db.relationship('Evidence', lazy='dynamic')
    #group_key_genename = db.relationship('Edge', primaryjoin='Edge.group_key'=='Conversion.genename', foreign_keys='Conversion.genename')
    #group_key2_genename = db.relationship('Edge', primaryjoin='Edge.group_key2'=='Conversion.genename', foreign_keys='Conversion.genename')

    def get_groups(self,):
        group1 = db.session.query(Group).filter(Group.id==self.group_key).first()
        group2 = db.session.query(Group).filter(Group.id==self.group_key2).first()
        return (group1, group2)

'''class Edge2(db.Model):
    """Not necessarily complexes group group edge threshold svm score 0.2"""
    id = db.Column(db.Integer, primary_key=True)
    group_key = db.Column(db.Integer, db.ForeignKey('group.id') )
    group_key2 = db.Column(db.Integer, db.ForeignKey('group.id') )
    #score = db.Column(db.Float)
    svmscore =db.Column(db.Float)

    #evidences = db.relationship('Evidence', lazy='dynamic')

    def get_groups(self,):
        group1 = db.session.query(Group).filter(Group.id==self.group_key).first()
        group2 = db.session.query(Group).filter(Group.id==self.group_key2).first()
        return (group1, group2)
'''



class Evidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    evidence_type = db.Column(db.String(255))
    


class GroupComplexMapping(db.Model):
    """A mapping between groups and complexes"""
    group_key = db.Column(db.Integer, db.ForeignKey('group.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)

class ComplexEnrichment(db.Model):
    """Annotation Enrichment for a Complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'))
    #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term GroupID     t type  t group    t name and depth in group        Q&T list
    corr_pval = db.Column(db.Float)
    t_count = db.Column(db.Integer)
    q_count = db.Column(db.Integer)
    qandt_count = db.Column(db.Integer)
    qandt_by_q = db.Column(db.Float)
    qandt_by_t = db.Column(db.Float)
    term_id = db.Column(db.String(255))
    t_type = db.Column(db.String(63))
    t_group = db.Column(db.Integer)
    t_name = db.Column(db.String(255))
    depth_in_group = db.Column(db.Integer)
    qandt_list = db.Column(db.String(255))

    def get_groups(self,):
        groups = []
        for acc in self.qandt_list.split(','):
            prot = db.session.query(Group).filter_by(uniprot_acc=acc).first()
            groups.append(prot)
        return groups





