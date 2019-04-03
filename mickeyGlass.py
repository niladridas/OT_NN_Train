# Vedang Deshpande 02 April 2019
import os
os.chdir("C:/_WORK_SPACE/test_py/neural-networks")
from neuralnetwork import *
# from neuralnetworkTests import *

'''Define functions for Mickey-Glass Sequence generation'''
def xdot(x_t, x_tau, a,b,n):
    return (-b*x_t + a*x_tau/(1.0 + x_tau**n))

def xNext(x_t, x_tau, dt, a,b,n):
    k1 = dt*xdot(x_t, x_tau, a,b,n)
    k2 = dt*xdot(x_t+0.5*k1, x_tau, a,b,n)
    k3 = dt*xdot(x_t+0.5*k2, x_tau, a,b,n)
    k4 = dt*xdot(x_t+k3, x_tau, a,b,n)
    return (x_t + k1/6 + k2/3 + k3/3 + k4/6)

def generateSequence(a,b,n,tau,dt,num,x0):
    xSequence = (x0,)
    for i in range(1,num):
        if i>=tau:
            x_tau = xSequence[i-tau]
        else:
            x_tau = 0
        x_i =  xNext(xSequence[i-1], x_tau, dt, a,b,n)
        xSequence = xSequence + (x_i,)
    return xSequence

def createTrainingSet(xSequence, tau):
    num = len(xSequence)
    trainSet = []
    for i in range(1,tau):
        example = ( ((0,)*(tau-i) + xSequence[0:i]), xSequence[i] )
        trainSet.append(example)
    for i in range(tau,num):
        example = ((xSequence[i-tau:i]), xSequence[i] )
        trainSet.append(example)
    return trainSet

''' Generate Mickey-Glass Sequence & trainig data set'''
a=0.3
b=0.1
n=10
dt=1
num=300
x0=0.5
tau=30
xSeq = generateSequence(a,b,n,tau,dt,num,x0)
print(xSeq)

trainSet = createTrainingSet(xSeq, tau)
print(trainSet)

'''Create and train neural network'''
numInputs = tau
numLayers = 2
numNodes = 10
network = makeNetwork(numInputs, numLayers, numNodes)
network.train(trainSet, learningRate=0.25, maxIterations=10000)

err = 0
for xHist, x_next in trainSet:
    # print('Error for {} is {:.4f}. Output was:{:.4f}'.format(xHist, x_next - network.evaluate(xHist), network.evaluate(xHist)))
    err +=  (x_next - network.evaluate(xHist))**2
rms = (err/len(trainSet))**(0.5)
print(rms)
