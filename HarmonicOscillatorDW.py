import numpy as np
import math
import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d import Axes3D

#mass H = 1.6727346e-27
#masss e- = 9.109e-31 kg
numAtoms = 1
initialWalkers = 2000
wvnmbr = 4.55634e-6
omega = 2000.0000*wvnmbr #in atomic units
dimensions = 1
mass = 1836.35/2
deltaT = 5.0 #simulation parameter - adjustable
alpha = 1.0/(2.0*deltaT) #simulation parameter - adjustable
D = 0.5
sigma = math.sqrt((2*D*deltaT)/mass)

class Walker:
    def __init__(self):
        self.coords = np.zeros((numAtoms,dimensions)) #1d surface
        self.WalkerV = 0.0

myWalkers = [Walker() for r in range(initialWalkers)]
whoFrom = np.arange(initialWalkers)

def getPotentialForWalkers(): #use coordinates of walkers to get V
    omsqd = math.pow(omega,2)
    prefactor = 0.50000*mass*omsqd
    for d in range(len(myWalkers)):
        crds = float(myWalkers[d].coords)
        crdssq = math.pow(crds,2)
        myWalkers[d].WalkerV = (prefactor*crdssq)

def giveBirth(singleWalker):
    myWalkers.append(singleWalker)

def deleteWalker(walkerIndex,whoFr,Desc):
    del myWalkers[walkerIndex]
    if Desc:
        whoFr = np.delete(whoFr,walkerIndex)
    return whoFr

def birthOrDeath(vref,whoFr,Desc):
    deathAr = []
    birthAr = []
    birthArInd = []

    for y in range(len(myWalkers)):
        Rng = np.random.random()
        curV = myWalkers[y].WalkerV

        if curV > Vref:
            P = -1*(curV - Vref)*deltaT
            exP = math.exp(P)
            if exP < Rng:
                deathAr.append(y)
        else:
            P = -1*(curV - Vref) * deltaT
            exP = math.exp(P) - 1
            if exP > Rng:
                singleWalker = copy.deepcopy(myWalkers[y])
                birthAr.append(singleWalker)
                if Desc:
                    #get the index from where the pregnant walker is from
                    birthArInd.append(whoFr[y])

    if deathAr: #If it's not empty
        for k in reversed(deathAr):
            whoFr = deleteWalker(k,whoFr,Desc)
    if birthAr:
        for j in birthAr:
            giveBirth(j)
        if Desc:
            whoFr = np.concatenate((whoFr, birthArInd))


    del deathAr[:]
    del birthAr[:]

    return whoFr


def moveRandomly():
    #choose a random number from gaussian distribution (1/2pisig)(e^(-dx^2/1sig^2))
    for p in range(len(myWalkers)):
        gaussStep = np.random.normal(loc = 0.00000000000000,scale=sigma)
        myWalkers[p].coords = myWalkers[p].coords + gaussStep

def getVref(): #Use potential of all walkers to calculate vref
    varray = np.array([k.WalkerV for k in myWalkers])
    Vbar = np.average(varray)
    vref = Vbar - alpha*((float(len(myWalkers))-float(initialWalkers))/float(initialWalkers))
    return vref

#Start!
vrefAr = np.zeros(1000)
xc = []
popAr = np.zeros(1000)
Vref = 10000
DW = False
for i in range(1000):
    moveRandomly()
    getPotentialForWalkers()
    if i==0:
        Vref = getVref()
    if i>=950:
        DW = True
        if i==950:
            whoFrom = np.arange(len(myWalkers))
            oneHundo = copy.deepcopy(myWalkers)
    whoFrom = birthOrDeath(Vref,whoFrom,DW)
    getPotentialForWalkers()
    Vref = getVref()
    #Plotting business
    vrefAr[i] = Vref
    popAr[i] = len(myWalkers)

dWeights = np.zeros(len(oneHundo))
print "Final NumWalkers = ",len(myWalkers)

unique,counts = np.unique(whoFrom,return_counts=True)
z=0
for u in unique: #teehee
    dWeights[u]=counts[z]
    z+=1

for wlker in range(len(oneHundo)):
    xc.append(float(oneHundo[wlker].coords))

bns = 30

#Plotting against analytical solution. Can comment out.
xx = np.linspace(-3,3,bns)
psisqHO = math.sqrt(mass*omega/math.pi)*np.exp(-1*2*mass*omega*np.square(xx) / (2))

fff,ax1 = plt.subplots()
r, binz = np.histogram(xc, bins=bns, normed=True,range=[-3, 3], weights=dWeights)
r2, binz = np.histogram(xc, bins=bns, normed=True,range=[-3, 3])
r2 = r2**2
r2 = r2*1.2
s, bz = np.histogram(dWeights, bins=100, normed=True,range=[0, 100])
print 'descendent weights:',dWeights
print 'sum',np.sum(dWeights)
plt.plot(np.linspace(-3,3,bns),r)
plt.plot(xx,psisqHO)
plt.plot(np.linspace(-3,3,bns),r2)
plt.show()
plt.plot(np.linspace(0,100,100),s)
plt.show()
print np.average(vrefAr[200:])
print omega/2
var = [h for h in vrefAr]
x = np.linspace(0,1000,1000)
plt.plot(x,vrefAr)
om = [omega/2]*len(x)
plt.plot(x,om)
plt.show()
plt.plot(x,popAr)
plt.show()