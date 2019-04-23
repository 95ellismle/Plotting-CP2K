import numpy as np
import matplotlib.pyplot as plt

nAtom = 4 
nRep = 3 
colors = ['r', 'k', (1, 0.5, 0.5)] 
 
def gaussianCTMQC(RIv, RJv, sigma):  
    """  
    This creates the data of a gaussian function. RIv is the  
    position the gaussian is evalutated at, RJv is the mean  
    and sigma the width.  
    """  
    normConst = (2*np.pi)**(-0.5)/sigma  
    exponent = -((RIv - RJv)**2/(2*(sigma**2)))  
    data = normConst * np.exp(exponent) 
    return data 


allPos = [np.random.rand(nAtom) for J in range(nRep)]
allPos = np.array(allPos)

allGauss = np.ones((nRep, nRep, nAtom))
allProdGauss = np.ones((nRep, nRep))

for I in range(nRep):
   for J in range(nRep):
      for v in range(nAtom):
         g = gaussianCTMQC(allPos[I][v], allPos[J][v], 0.2)
         allGauss[I, J, v] = g
         allProdGauss[I, J] *= g

nuclDens = np.ones(nRep)
for I in range(nRep):
   for J in range(nRep):
        nuclDens[I] += allProdGauss[I, J]
nuclDens /= nRep

for J in range(nRep):
   plt.plot(allPos[J],
            np.zeros(len(allPos[J])), 'o',
            color=colors[J], ms=9)

J = 0
v = 0
for v in range(nAtom):
   for J in range(nRep):
      for I in range(nRep):
          plt.plot(allPos[I, v], allGauss[I, J, v], '.',
                   color=colors[I], ms=9)
      
          gaussX = np.arange(-0.1,1.1, 0.01)
          gaussY = gaussianCTMQC(gaussX, allPos[J, v], 0.2)
          plt.plot(gaussX, gaussY, color=colors[I], lw=0.5,
                   ls='-')
   
plt.show()
 
