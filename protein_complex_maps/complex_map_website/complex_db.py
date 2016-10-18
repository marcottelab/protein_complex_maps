
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
    proteins = db.relationship('Protein', secondary='protein_complex_mapping', back_populates='complexes', lazy='dynamic')
    enrichments = db.relationship('ComplexEnrichment', lazy='dynamic')

    def complex_link(self,):
        retstr = "<a href=displayComplexes?complex_key=%s>%s</a>" % (self.complex_id, self.complex_id)
        return retstr

    def edges(self,):
        es = []
        pairs = set()
        for prot1 in self.proteins:
            for prot2 in self.proteins:
                #kdrew: check to see if the pair has already been entered
                if frozenset((prot1,prot2)) not in pairs:
                    edge = db.session.query(Edge).filter((and_(Edge.protein_key == prot1.id, Edge.protein_key2 == prot2.id) | and_(Edge.protein_key == prot2.id,Edge.protein_key2 == prot1.id))).first()
                    #es = es + edge
                    if edge != None:
                        es.append(edge)
                        #kdrew: add pair of prot ids to keep track of edges already added
                        pairs.add(frozenset((prot1,prot2)))

        return sorted(list(set(es)), key=lambda es: es.score, reverse=True)

class Conversion(db.Model):
    """A mapping between different ID types"""
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(63))
    proteinname = db.Column(db.String(63))
    genename = db.Column(db.String(63)) 
    uniprot_acc = db.Column(db.String(63))
       

class Protein(db.Model):
    """A single protein"""
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(63))
    uniprot_acc = db.Column(db.String(63))
    genename = db.Column(db.String(255))
    proteinname = db.Column(db.String(255))
    uniprot_url = db.Column(db.String(255))
    complexes = db.relationship('Complex', secondary='protein_complex_mapping',  back_populates='proteins', lazy='dynamic')

    def uniprot_link(self,):
        if self.genename == "":
            retstr = self.gene_id
        else:
            retstr = "<a href=%s target=\"_blank\">%s</a>" % (self.uniprot_url, self.genename)
        return retstr

class Edge(db.Model):
    """A protein protein edge"""
    id = db.Column(db.Integer, primary_key=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id') )
    protein_key2 = db.Column(db.Integer, db.ForeignKey('protein.id') )
    in_complex = db.Column(db.Integer)
    score =db.Column(db.Float)

    evidences = db.relationship('Evidence', lazy='dynamic')
    #protein_key_genename = db.relationship('Edge', primaryjoin='Edge.protein_key'=='Conversion.genename', foreign_keys='Conversion.genename')
    #protein_key2_genename = db.relationship('Edge', primaryjoin='Edge.protein_key2'=='Conversion.genename', foreign_keys='Conversion.genename')

    def get_proteins(self,):
        prot1 = db.session.query(Protein).filter(Protein.id==self.protein_key).first()
        prot2 = db.session.query(Protein).filter(Protein.id==self.protein_key2).first()
        return (prot1, prot2)

'''class Edge2(db.Model):
    """Not necessarily complexes protein protein edge threshold svm score 0.2"""
    id = db.Column(db.Integer, primary_key=True)
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id') )
    protein_key2 = db.Column(db.Integer, db.ForeignKey('protein.id') )
    #score = db.Column(db.Float)
    svmscore =db.Column(db.Float)

    #evidences = db.relationship('Evidence', lazy='dynamic')

    def get_proteins(self,):
        prot1 = db.session.query(Protein).filter(Protein.id==self.protein_key).first()
        prot2 = db.session.query(Protein).filter(Protein.id==self.protein_key2).first()
        return (prot1, prot2)
'''



class Evidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    edge_key = db.Column(db.Integer, db.ForeignKey('edge.id'))
    evidence_type = db.Column(db.String(255))
    


class ProteinComplexMapping(db.Model):
    """A mapping between proteins and complexes"""
    protein_key = db.Column(db.Integer, db.ForeignKey('protein.id'), primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'), primary_key=True)

class ComplexEnrichment(db.Model):
    """Annotation Enrichment for a Complex"""
    id = db.Column(db.Integer, primary_key=True)
    complex_key = db.Column(db.Integer, db.ForeignKey('complex.id'))
    #  signf   corr. p-value   T   Q   Q&T Q&T/Q   Q&T/T   term ID     t type  t group    t name and depth in group        Q&T list
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

    def get_proteins(self,):
        proteins = []
        for acc in self.qandt_list.split(','):
            prot = db.session.query(Protein).filter_by(uniprot_acc=acc).first()
            proteins.append(prot)
        return proteins





