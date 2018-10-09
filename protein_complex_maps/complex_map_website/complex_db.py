
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
    orthogroups = db.relationship('Orthogroup', secondary='protein_complex_mapping', back_populates='complexes', lazy='dynamic')
    enrichments = db.relationship('ComplexEnrichment', lazy='dynamic')

    def complex_link(self,):
        retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.complex_id, self.complex_id)
        return retstr

    def edges(self,):
        es = []
        pairs = set()
        for orthogroup1 in self.orthogroups:
            for orthogroup2 in self.orthogroups:
                #kdrew: check to see if the pair has already been entered
                if frozenset((orthogroup1,orthogroup2)) not in pairs:
                    edge = db.session.query(Edge).filter((and_(Edge.orthogroup_key == orthogroup1.id, Edge.orthogroup_key2 == orthogroup2.id) | and_(Edge.orthogroup_key == orthogroup2.id,Edge.orthogroup_key2 == orthogroup1.id))).first()
                    #es = es + edge
                    if edge != None:
                        es.append(edge)
                        #kdrew: add pair of prot ids to keep track of edges already added
                        pairs.add(frozenset((orthogroup1,orthogroup2)))

        return sorted(list(set(es)), key=lambda es: es.score, reverse=True)

class Conversion(db.Model):
    """A mapping between different OrthogroupID types"""
    id = db.Column(db.Integer, primary_key=True)
    OrthogroupID = db.Column(db.String(63))
    Species = db.Column(db.String(63))
    ProteinID = db.Column(db.String(63)) 
  

class Orthogroup(db.Model):
    """A single orthogroup"""
    id = db.Column(db.Integer, primary_key=True)
    OrthogroupID = db.Column(db.String(63))
    #complexes = db.relationship('Complex', secondary='protein_complex_mapping',  back_populates='orthogroups', lazy='dynamic')

class Edge(db.Model):
    """A orthogroup orthogroup edge"""
    id = db.Column(db.Integer, primary_key=True)
    OrthogroupID_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
    OrthogroupID_key2 = db.Column(db.Integer, db.ForeignKey('orthogroup.id') )
    in_complex = db.Column(db.Integer)
    score =db.Column(db.Float)

    evidences = db.relationship('Evidence', lazy='dynamic')
    #orthogroup_key_genename = db.relationship('Edge', primaryjoin='Edge.orthogroup_key'=='Conversion.genename', foreign_keys='Conversion.genename')
    #orthogroup_key2_genename = db.relationship('Edge', primaryjoin='Edge.orthogroup_key2'=='Conversion.genename', foreign_keys='Conversion.genename')

    def get_orthogroups(self,):
        orthogroup1 = db.session.query(Orthogroup).filter(Orthogroup.id==self.orthogroup_key).first()
        orthogroup2 = db.session.query(Orthogroup).filter(Orthogroup.id==self.orthogroup_key2).first()
        return (orthogroup1, orthogroup2)


class Evidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    evidence_type = db.Column(db.String(255))
    


class OrthogroupComplexMapping(db.Model):
    """A mapping between orthogroups and complexes"""
    orthogroup_key = db.Column(db.Integer, db.ForeignKey('orthogroup.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)

class ComplexEnrichment(db.Model):
    """Annotation Enrichment for a Complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'))
    #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term OrthogroupID     t type  t orthogroup    t name and depth in orthogroup        Q&T list
    corr_pval = db.Column(db.Float)
    t_count = db.Column(db.Integer)
    q_count = db.Column(db.Integer)
    qandt_count = db.Column(db.Integer)
    qandt_by_q = db.Column(db.Float)
    qandt_by_t = db.Column(db.Float)
    term_id = db.Column(db.String(255))
    t_type = db.Column(db.String(63))
    t_orthogroup = db.Column(db.Integer)
    t_name = db.Column(db.String(255))
    depth_in_orthogroup = db.Column(db.Integer)
    qandt_list = db.Column(db.String(255))

    def get_orthogroups(self,):
        orthogroups = []
        for acc in self.qandt_list.split(','):
            prot = db.session.query(Orthogroup).filter_by(uniprot_acc=acc).first()
            orthogroups.append(prot)
        return orthogroups





