from rpy2.robjects.packages import importr
import rpy2.robjects as robj

import numpy as np

def convert(a):
    x = np.array(a)
    nr, nc = x.shape
    xvec = robj.FloatVector(list(x.transpose().reshape((x.size))))
    xr = robj.r.matrix(xvec, nrow=nr, ncol=nc)
    return xr

exact = importr('exact2x2')

def get_p_value_null(table):
    return exact.exact2x2(convert(table), robj.r("NULL"), 1/1.25, alternative='greater')[0][0]

def get_p_value_positive(table):
    return exact.exact2x2(convert(table), robj.r("NULL"), 1/1.25, alternative='less')[0][0]
