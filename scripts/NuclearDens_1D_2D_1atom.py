import matplotlib.pyplot as plt
import numpy as np
import mayavi.mlab as mlab

plotDimension = 2


def gaussian1D(x, mu, sigma):
    exponent = -(x-mu)**2/(2*(sigma**2))
    prefact = 1/(2*np.pi*sigma**2)**(0.5)
    return prefact * np.exp(exponent)


def gaussian2D(x, y, mu, sigma):
    x -= mu[0]
    y -= mu[1]
    exponent = -(x*x + y*y)/(2*(sigma**2))
    prefact = 1/(2*np.pi*sigma**2)**(0.5)
    return prefact * np.exp(exponent)


if plotDimension == 1:

    f, a = plt.subplots(ncols=3)

    x = np.arange(0, 10, 0.01)

    for i, ax in enumerate(a):
        width = 0.5*(i+1)
        mu = 3  # 0.5*(i+1)
        y1 = gaussian1D(x, 5+mu, width)
        y2 = gaussian1D(x, 5-mu, width)

        nuclDens = y1 + y2
        QM = -np.gradient(nuclDens)/(nuclDens)

        ax.plot(x, y1, lw=0.7, alpha=0.5)
        ax.plot(x, y2, lw=0.7, alpha=0.5)
        ax.plot(x, nuclDens, 'k--', label=r"$|\chi|^2$")

        ax.plot(x, QM, 'r-', lw=2, label=r"$\mathcal{Q}$")

        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        if i != 0:
            ax.spines['left'].set_visible(False)
        ax.set_xlabel("Nuclear Positions")

    #    ax.set_title("Width = %.2g" % width, fontsize=18)
    #    ax.legend(loc='upper right')

    plt.suptitle("Nuclear Density and Quantum Momentum -1D, 1 atom 2 Traj",
                 fontsize=20)

    plt.show()

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

    @mlab.animate(delay=50)
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
            saveFolder = "/homes/mellis/Documents/Graphs/CTMQC/New_QM/tmpImg"
            fname = "%s/%s%i.png" % (saveFolder, "0"*(numZeros-len(str(i))), i)
            print(fname)
            mlab.savefig(fname)
            mlab.view(0, 0)
            yield

    plotPic()
    mlab.show()
