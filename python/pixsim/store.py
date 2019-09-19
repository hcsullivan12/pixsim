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

def get_array(ses, name=None, id=None):
    """Return array matching name or id."""
    if name is not None:
        return arrays(ses).filter_by(name=name).order_by(desc(pixsim.models.Array.created)).first()
    if id is None:
        return None

    if id < 0:
        return arrays(ses).order_by(desc(pixsim.models.Array.id)).first()
    else:
        return arrays(ses).get(id)

def dump_array(ses, arr):
    print 'Contents of',arr.name,'...'
    arr = arr.data
    dim = list(arr.shape)
    print 'Shape =',arr.shape
    if len(dim) == 4:
        print 'x =',arr[0,:,:,:].reshape(arr[0].size)
        print 'y =',arr[1,:,:,:].reshape(arr[1].size)
        print 'z =',arr[2,:,:,:].reshape(arr[2].size)
        print 'Average x component =',sum(arr[0,:,:,:].reshape(arr[0].size))/arr[0].size
    else:
        print arr

def dump_arrays(ses):
    print 'Arrays...'
    for arr in arrays(ses):
        print 'id: %-2s  name: %-10s  typename: %-10s  shape: %-10s' % (arr.id, arr.name, arr.typename, arr.data.shape)

def dump_result(ses, res):
    print 'id: %-2s  name: %-10s  typename: %-10s  data: %-2s  parent: %-2s' % (res.id, res.name, res.typename, len(res.data), res.parent_id)

def dump_results(ses):
    print 'Results...'
    for res in results(ses):
        dump_result(ses,res)

def dump_table(ses):
    """Dump table."""
    dump_results(ses)
    dump_arrays(ses)

def dump(ses, res_id=None, arr_id=None, res=None, arrs=None):
    """Dump contents."""
    if res_id:
        dump_result(ses, results(ses).get(res_id))
    if arr_id:
        dump_array(ses, arrays(ses).get(arr_id))
    if res:
        dump_results(ses)
    if arrs:
        dump_arrays(ses)
    dumpall = [res_id,arr_id,res,arrs]
    if not any(dumpall):
        dump_table(ses)

def get_last_ids(ses):
    return {'array':arrays(ses).order_by(desc(pixsim.models.Array.id)).first().id,
            'result':results(ses).order_by(desc(pixsim.models.Result.id)).first().id }

def interpret(s):
    """Interpret string as int or string"""
    id,name = None, None
    try:
        int(s)
        id=s
    except ValueError:
        name=s
    return id,name

def get_result(ses, source=None, name=None, id=None):
    """Return result matching type or id."""
    checkid = id
    checkname = name
    if source is not None:
        checkid,checkname = interpret(source)
    
    if checkid is not None:
        if checkid < 0:
            return results(ses).order_by(desc(pixsim.models.Result.checkid)).first()
        else:
            return results(ses).get(checkid)
    elif checkname is not None:
        return results(ses).filter_by(name=checkname).order_by(desc(pixsim.models.Result.created)).first()
    else:
        return None
