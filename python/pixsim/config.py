import os
import numpy as np

def convert(sec):
    '''
    Convert section values to string, float, etc.
    '''
    newsec = dict()
    for k,v in sec.iteritems():
        if '(' in v:
            newsec[k] = np.asarray(eval(v))
        elif '\'' in v or '\"' in v:
            newsec[k]=str(v[1:-1])
        elif '/' in v:
            newsec[k]=v
        else:
            newsec[k]=eval(v)
    return newsec

def parse(filepath):
    '''
    Read in configuration file.
    '''
    from ConfigParser import SafeConfigParser
    cfg = SafeConfigParser()
    
    if hasattr(filepath, "readline"):
        cfg.readfp(filepath)
    else:
        filepath = os.path.expanduser(os.path.expandvars(filepath))
        if not os.path.exists(filepath):
            raise ValueError('No such filepath: %s' % filepath)
        cfg.read(filepath)
    
    ret = dict()
    for secname in cfg.sections():
        sec = {kv[0]:kv[1] for kv in cfg.items(secname)}
        ret[secname] = convert(sec)
    return ret
