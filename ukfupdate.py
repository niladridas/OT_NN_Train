import numpy as np
from numpy.linalg import inv


def ukfupdate(xsigmapts, ysigmapts, yobs, sigw):
    # Calculating the mean of xSigmapts
    xm = xsigmapts.mean(1)
    # Calculating the mean of ySigmapts
    ym = ysigmapts.mean(1)
    # Calculating Pxx
    l1 = np.shape(xsigmapts)[1]
    pxx = np.zeros(l1, l1)
    pyy = np.zeros(l1, l1)
    pxy = np.zeros(l1, l1)
    for i in range(0, 2*l1):
        pxx = pxx + sigw.wc[i]*np.matmul(xsigmapts[:,i]-xm,(xsigmapts[:,i]-xm).transpose())
        pyy = pyy + sigw.wc[i]*np.matmul(ysigmapts[:,i]-ym,(ysigmapts[:,i]-ym).transpose())
        pxy = pxy + sigw.wc[i]*np.matmul(xsigmapts[:,i]-xm,(ysigmapts[:,i]-ym).transpose())
    K = np.matmul(pxy,inv(pyy))
    xmpost = xm + np.matmul(K, (yobs-ym))
    xcpost = pxx - np.matmul(np.matmul(K, pyy), K.transpose())
    updatedata = {"xmpost": xmpost, "xcpost": xcpost}
    return updatedata
