#!/usr/bin/env python
'''
Define data objects of pixsim.
'''
import json
import numpy
import io

from datetime import datetime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey, Table
import sqlalchemy.types as types

from sqlalchemy import Column, Integer, String, Enum
Base = declarative_base()

class NumpyArray(types.TypeDecorator):

    impl = types.BLOB

    def process_bind_param(self, value, dialect):
        out = io.BytesIO()
        numpy.save(out, value)
        out.seek(0)
        return out.read()

    def process_result_value(self, value, dialect):
        out = io.BytesIO(value)
        out.seek(0)
        return numpy.load(out, allow_pickle=True)

    def copy(self, **kw):
        return NumpyArray(self.impl.length)    

class JSONBLOB(types.TypeDecorator):

    impl = String

    def process_bind_param(self, value, dialect):
        return json.dumps(value)

    def process_result_value(self, value, dialect):
        return json.loads(value)

    def copy(self, **kw):
        return JSONBLOB(self.impl.length)    

#todo: Simplify these 
result_types = [
    'geo',
    'boundary',
    'raster',
    'drift',
    'volume',
    'evaluate',
    'points',
    'step'
]

array_types = [
    'tuples',       # N_tuples X N_tuple_length.  List of tuples
    'points',       # N_points X 3.  List of (x,y,z) points
    'rays',         # N_rays X 6.  List of (x1,y1,z1, x2,y2,z2) endpoints
    'indices',      # N_groups X N_groupsize.  Indices into some other array, in fixed sized groups
    'scalar',       # N_scalars.  One scalar value per something (eg, per point)
    'linspace',
    'mgrid',
    'gscalar',
    'gvector'
]

subs = Table('subs',
    Base.metadata,
    Column("result_id", Integer, ForeignKey("result.id")),
    Column("array_id", Integer, ForeignKey("array.id"))
)

class Result(Base):
    __tablename__ = 'result'
    
    id = Column(Integer, primary_key=True)
    parent_id = Column(Integer, ForeignKey('result.id'))
    name = Column(String, default='')
    typename = Column(Enum(*result_types))
    #params = Column(JSONBLOB, default='[]')
    created = Column(types.TIMESTAMP, default=datetime.now)

    data   = relationship("Array", secondary=subs, backref="arrays")
    parent = relationship("Result", remote_side=[id])

    def array_data_by_type(self):
        '''
        Return dictionary key'ed by array type.  No check for uniqueness.
        '''
        return {a.type:a.data for a in self.data}
    def array_data_by_name(self):
        '''
        Return dictionary key'ed by array name.  No check for uniqueness.
        '''
        return {a.name:a.data for a in self.data}

    def triplets(self):
        '''
        Return (typename, name data) triplet for each array, in order.
        '''
        return [ (a.typename,a.name,a.data) for a in self.data ]
        

class Array(Base):
    __tablename__ = 'array'

    id = Column(Integer, primary_key=True)
    name = Column(String, default='')
    typename = Column(Enum(*array_types))
    data = Column(NumpyArray)

        
