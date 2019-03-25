import numpy as np
from numpy.linalg import inv

def ukfupdate(xSigmapts,ySigmapts,yObs):
    # Calculating the mean of xSigmapts
    xM = xSigmapts.mean(1)
    # Calculating the mean of ySigmapts
    yM = ySigmapts.mean(1)
    # Calculating Pxx
    L = np.shape(xSigmapts)[1]
    Ly = np.shape(ySigmapts)[1]
    Pxx = np.zeros(L,L)
    Pyy = np.zeros(L,L)
    Pxy = np.zeros(L,L)
    for i in range (0,2*L):
        Pxx = Pxx + sigW.wc[i]*np.matmul(xSigmapts[:,i]-xM,(xSigmapts[:,i]-xM).transpose())
        Pyy = Pyy + sigW.wc[i]*np.matmul(ySigmapts[:,i]-yM,(ySigmapts[:,i]-yM).transpose())
        Pxy = Pxy + sigW.wc[i]*np.matmul(xSigmapts[:,i]-xM,(ySigmapts[:,i]-yM).transpose())
    K = np.matmul(Pxy,inv(Pyy))
    xMpost =  xM + np.matmul(K,(yObs-yM))
    xCpost = Pxx - np.matmul(np.matmul(K,Pyy),K.transpose())
    updatedata = {"xMpost": xMpost, "xCpost": xCpost}
    return updatedata
