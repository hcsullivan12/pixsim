#!/usr/bin/env python
'''
Interface for persistance.
'''

import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker
from sqlalchemy import desc
import pixsim.models
import numpy as np

import sqlalchemy.exc
IntegrityError = sqlalchemy.exc.IntegrityError

def session(database = None):
    assert(database is not None)

    if ":///" not in database:      
        database = "sqlite:///" + database

    engine = sa.create_engine(database)
    pixsim.models.Base.metadata.create_all(engine)
    Session = sa.orm.sessionmaker(engine)
    return Session()
    
def results(ses):
    return ses.query(pixsim.models.Result)

def arrays(ses):
    return ses.query(pixsim.models.Array)

def get_array(ses, typename=None, id=None):
    """Return array matching typename or id."""
    if typename is not None:
        return arrays(ses).filter_by(typename=typename).order_by(desc(pixsim.models.Array.created)).first()
    if id is None:
        return None

    if id < 0:
        return arrays(ses).order_by(desc(pixsim.models.Array.id)).first()
    else:
        return arrays(ses).get(id)

def dump_array(ses, arr_id):
    print 'Contents of',arrays(ses).get(arr_id).name,'...'
    arr = arrays(ses).get(arr_id).data
    dim = list(arr.shape)
    print 'Shape =',arr.shape
    if len(dim) == 4:
        print 'x =',arr[0,:,:,:].reshape(arr[0].size)
        print 'y =',arr[1,:,:,:].reshape(arr[1].size)
        print 'z =',arr[2,:,:,:].reshape(arr[2].size)
    else:
        print arr

def dump_table(ses):
    """Dump table."""
    print 'Results...'
    for res in results(ses):
        print 'id: %-2s  name: %-10s  typename: %-10s  data: %-2s' % (res.id, res.name, res.typename, len(res.data))
    print 'Arrays...'
    for arr in arrays(ses):
        print 'id: %-2s  name: %-10s  typename: %-10s  shape: %-10s' % (arr.id, arr.name, arr.typename, arr.data.shape)

def dump(ses, arr_id):
    """Dump contents."""
    if arr_id is not None:
        dump_array(ses, arr_id)
    else:
        dump_table(ses)
    
def get_result(ses, typename=None, id=None):
    """Return result matching type or id."""
    if typename is not None:
        return results(ses).filter_by(typename=typename).order_by(desc(pixsim.models.Result.created)).first()
    if id is None:
        return None

    if id < 0:
        return results(ses).order_by(desc(pixsim.models.Result.id)).first()
    else:
        return results(ses).get(id)