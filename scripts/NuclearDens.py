import numpy as np
import matplotlib.pyplot as plt

nPos = 4 
nRep = 3 
colors = ['r', 'k', (1, 0.5, 0.5)] 
 
def gaussianCTMQC(RIv, RJv, sigma):  
    """  
    This creates the data of a gaussian function. RIv is the  
    position the gaussian is evalutated at, RJv is the mean  
    and sigma the width.  
    """  
    normConst = 1
    exponent = -((RIv - RJv)**2/(2*(sigma**2)))  
    data = normConst * np.exp(exponent) 
    return data
 
 
gaussX = np.arange(-nPos, nPos*11, 0.01*nPos) 
avgGauss = np.ones(len(gaussX)) 
allXPos = [] 

f, ax = plt.subplots(1, 2) 
for y in range(0, 3*nRep, 3): 
    xpos = np.random.rand(nPos)*nPos*10 
    ypos = np.ones(nPos)*y 

    prodGauss = np.ones(len(gaussX)) 
    for x in xpos: 
        allXPos.append(x)
        gaussY = gaussianCTMQC(gaussX, x, nPos) 
        prodGauss *= gaussY 
        # plot atomic gaussians 
        ax[0].plot(gaussX, gaussY+y, color=colors[y//3], lw=0.5, ls='--') 

    # plot nuclear density 
    avgGauss += prodGauss
    ax[0].plot(gaussX, (prodGauss)+y+1.5, color=colors[y//3]) 

    # plot atoms 
    ax[0].plot(xpos, ypos, 'o', color=colors[y//3], ms=9) 
     
    # Nuclear Density annoatation     
    ax[0].annotate(r"|$\chi^{%i}$|^2" % (y//3), (max(gaussX)*0.9, ypos[0]+1.6), 
                fontsize=23, color=colors[y//3]) 

    #ax[0].axhline(ypos[0] + 2.7, ls='-', lw=0.5, color='k')


nuclDens = (avgGauss/np.max(prodGauss))
QM = -np.gradient(nuclDens)/nuclDens

qmAx = ax[1].twinx()
ax[1].plot(gaussX, nuclDens-np.min(nuclDens), color='k') 
ax[1].plot(allXPos, np.zeros(len(allXPos)), 'ko', ms=9)
qmAx.plot(gaussX, QM, '--', lw=0.7)

for axI in ax:
   axI.set_xticks([]) 
   axI.set_yticks([]) 
 
   axI.set_xlabel("Nuclear Geom") 

qmAx.set_xticks([]) 
qmAx.set_yticks([]) 

# Fill in the pos and neg cols
qmAx.fill_between(gaussX[QM > 0], QM[QM > 0], color='g', alpha=0.1) 
qmAx.fill_between(gaussX[QM < 0], QM[QM < 0], color='r', alpha=0.1) 

qmAx.set_ylabel(r"$\mathcal{Q}$")
ax[0].set_ylabel("Replica")
ax[1].set_ylabel(r"$\frac{1}{N_{tr}} \sum_{J} |\chi^{J}|^2$")
#ax[2].set_ylabel(r"$\mathcal{Q}$")

plt.show()
