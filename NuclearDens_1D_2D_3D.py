import matplotlib.pyplot as plt
import numpy as np
import mayavi.mlab as mlab
from load import load_pos


plotDimension = 3

# 1D Params
multiAtom = False

# 2D Params

# 3D Params
loadReal = False


def gaussian1D(x, mu, sigma):
    x -= mu
    exponent = -(x**2)/(2*(sigma**2))
    prefact = 1/(2*np.pi*sigma**2)**(0.5)
    return prefact * np.exp(exponent)


def gaussian2D(x, y, mu, sigma):
    x -= mu[0]
    y -= mu[1]
    exponent = -(x*x + y*y)/(2*(sigma**2))
    prefact = 1/(2*np.pi*sigma**2)**(0.5)
    return prefact * np.exp(exponent)


def gaussian3D(x, y, z, mu, sigma):
    x -= mu[0]
    y -= mu[1]
    z -= mu[2]
    exponent = -(x*x + y*y + z*z)/(2*(sigma**2))
    prefact = 1/(2*np.pi*sigma**2)**(0.5)
    return prefact * np.exp(exponent)


if plotDimension == 1:
    if multiAtom is False:
        f, a = plt.subplots(ncols=3)

        x = np.arange(0, 10, 0.001)
        for i, ax in enumerate(a):
            width = 1  # *(i+1)
            mu = 2 - 0.5*(i+1)
            y1 = gaussian1D(x.copy(), 5+mu, width)
            y2 = gaussian1D(x.copy(), 5-mu, width)

            nuclDens = y1 + y2
            QM = -np.gradient(nuclDens)/(nuclDens)

            ax.plot(x, y1, lw=0.7, alpha=0.5)
            ax.plot(x, y2, lw=0.7, alpha=0.5)
            ax.plot(x, nuclDens, 'k--', label=r"$|\chi|^2$")

            axQM = ax.twinx()
            axQM.plot(x, QM, 'r-', lw=2, label=r"$\mathcal{Q}$")

            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            axQM.set_yticks([])
            if i != 0:
                ax.spines['left'].set_visible(False)
                axQM.spines['left'].set_visible(False)
            ax.set_xlabel("Nuclear Positions")

        #    ax.set_title("Width = %.2g" % width, fontsize=18)
        #    ax.legend(loc='upper right')

        plt.suptitle("Nuclear Density and Quantum Momentum -1D, 1 atom 2 Traj",
                     fontsize=20)

        plt.show()

    else:
        # pos shape = (nRep , nAtom)
        pos = np.array([[5], [5]])
        
            


elif plotDimension == 2:

    x, y = np.mgrid[0:10:100j, 0:10:100j]

    mlab.figure(size=(1000, 1000))

    mu = [7, 7]
    sigma = 0.9
    gauss1 = gaussian2D(x.copy(), y.copy(), mu, sigma)
#    gauss1 = (gauss1/np.max(gauss1)) * 3

    mu = [3, 3]
    sigma = 2
    gauss2 = gaussian2D(x.copy(), y.copy(), mu, sigma)
#    gauss2 = (gauss2/np.max(gauss2)) * 3

    nuclDens = gauss1 + gauss2
    nuclDens /= np.max(nuclDens)
    QM = -np.array(np.gradient(nuclDens))/(nuclDens)

    sampleDown = 5
    x_qm = x[::sampleDown, ::sampleDown]
    y_qm = y[::sampleDown, ::sampleDown]
    z_qm = np.ones(x_qm.shape) - 3
    QMx = QM[0, ::sampleDown, ::sampleDown]
    QMy = QM[1, ::sampleDown, ::sampleDown]
    QMz = np.zeros(QMx.shape)

    nuclPts = mlab.surf(x, y, nuclDens, opacity=0.2)
    QMPts = mlab.quiver3d(x_qm, y_qm, z_qm, QMx, QMy, QMz,
                          scale_factor=3)

    nFrames = 800
    numZeros = len(str(nFrames))

    @mlab.animate(delay=10)
    def plotPic():
        global nuclPts
        for i in range(nFrames):
            muX = 3  # 5 + np.sin(0.01*i)*4
            muY = 3  # 5 + np.cos(0.01*i)*3
            sigma = 1.5  # - np.sin(0.01 * i)
            gaussA = gaussian2D(x.copy(), y.copy(),
                                [muX, muY], sigma)

            muX = 7  # 5 - 3*np.cos(0.01*i)
            muY = 7  # 5 - 4*np.sin(0.01*i)
            sigma = 1 + np.sin(0.01 * i)/2
            gaussB = gaussian2D(x.copy(), y.copy(),
                                [muX, muY], sigma)

            nuclDens = gaussA + gaussB
            QM = -np.array(np.gradient(nuclDens))/(nuclDens)

            QMx = QM[0, ::sampleDown, ::sampleDown]
            QMy = QM[1, ::sampleDown, ::sampleDown]

            nuclPts.mlab_source.set(scalars=nuclDens)
            QMPts.mlab_source.set(u=QMx,
                                  v=QMy)
            #saveFolder = "/homes/mellis/Documents/Graphs/CTMQC/New_QM/tmpImg"
            #fname = "%s/%s%i.png" % (saveFolder, "0"*(numZeros-len(str(i))), i)
            #print(fname)
            #mlab.savefig(fname)
            mlab.view(0, 0)
            yield

    plotPic()
    mlab.show()

elif plotDimension == 3:
    if loadReal is True:
        folder = "/scratch/mellis/flavoured-cptk/200Rep_2mol3"
        allPos = load_pos.load_all_pos_in_folder(folder,
                                                 reps=[1, 51, 101, 201])
        pKeys = list(allPos.keys())
        nrep = len(allPos)
        nstep = len(allPos[pKeys[0]][0])
        allPos = [[allPos[i][0][allPos[i][1] != 'Ne']] for i in allPos]
        natom = np.shape(allPos)[2]//nstep
        allPos = np.reshape(allPos, (nstep, nrep, natom, 3))
        minX, maxX = np.min(allPos[:, :, :, 0]), np.max(allPos[:, :, :, 0]),
        minY, maxY = np.min(allPos[:, :, :, 1]), np.max(allPos[:, :, :, 1]),
        minZ, maxZ = np.min(allPos[:, :, :, 2]), np.max(allPos[:, :, :, 2]),
        np.max(allPos[:, :, :, 2])

        x, y, z = np.mgrid[minX:maxX:20j, minY:maxY:20j, minZ:maxZ:20j]

        pos = allPos[0]
    else:
        pos = [
               [  # Replica 1
                [0, 0, 0],  # atom 1
               ],
               [  # Replica 2
                [10, 10, 10],  # atom 1
               ],
               [  # Replica 3
                [10, 0, 10],  # atom 1
               ],
               [  # Replica 4
                [0, 10, 10],  # atom 1
               ],
              ]
        x, y, z = np.mgrid[0:10:20j, 0:10:20j, 0:10:20j]

    sigma = 1.5

    pos = np.array(pos).astype(float)
    pos0 = pos.copy()  # Initial positions
    allGauss = [[gaussian3D(x.copy(), y.copy(), z.copy(),
                            pos[J, v], sigma)
                 for v in range(len(pos[0]))]
                for J in range(len(pos))]

    nuclDens = np.product(allGauss, axis=1)  # Take product over atoms
    nuclDens = np.sum(nuclDens, axis=0)  # Sum over products
    QM = -np.array(np.gradient(nuclDens))/(nuclDens)

    colors = [(0.5, 0.5, 0.5), (1, 1, 1), (0, 0, 0), (1, 0, 0), (0, 1, 0),
              (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (0.5, 1, 0.5),
              (0.5, 0.5, 1), (0.5, 0.5, 0), (0.5, 0, 0.5), (0.5, 0.5, 1),
              (0, 0.5, 0.5),
              ]

    nuclPts = []
    for v in range(len(pos[0])):
        nuclPts.append(mlab.points3d(pos[:, v, 0], pos[:, v, 1], pos[:, v, 2],
                                     resolution=12, color=colors[v],
                                     scale_factor=sigma/2))

    qmPts = mlab.quiver3d(x, y, z, QM[0, :, :],
                          QM[1, :, :], QM[2, :, :])

    @mlab.animate(delay=10)
    def anim():
        """
        A generator function for animating the 3D visualisation
        """
        global pos, allPos, loadReal, x, y, z, sigma
        nFrames = 629
        if loadReal is True:
            nFrames = len(allPos)
        i = 0
        stepSize = 1
        while i < nFrames:
            if loadReal is True:
                pos = allPos[i]
            else:
                pos[0, 0, :] = pos0[0, 0, :] + (-np.cos(0.01 * i) + 1) * 2.5
                pos[1, 0, :] = pos0[1, 0, :] - (-np.cos(0.01 * i) + 1) * 2.5

            allGauss = [[gaussian3D(x.copy(), y.copy(), z.copy(),
                                    pos[I, v], sigma)
                         for v in range(len(pos[0]))]
                        for I in range(len(pos))]
            nuclDens = np.product(allGauss, axis=1)  # Take product over atoms
            nuclDens = np.sum(nuclDens, axis=0)  # Sum over products
            QM = -np.array(np.gradient(nuclDens))/(nuclDens)

            for v in range(len(pos[0])):
                nuclPts[v].mlab_source.set(x=pos[:, v, 0],
                                           y=pos[:, v, 1],
                                           z=pos[:, v, 2])
            qmPts.mlab_source.set(u=QM[0],
                                  v=QM[1],
                                  w=QM[2])
            i += stepSize
            yield

    anim()
    mlab.show()
