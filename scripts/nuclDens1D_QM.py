import numpy as np
import matplotlib.pyplot as plt

PI = np.pi


# Calc the prod of gaussians
def gaussian(sigma_vJ, R_vI, R_vJ):
    """
    Will calculate the gaussian assigned to an atom with
    width sigma_vJ

    Inputs:
        * sigma_vJ => [float] the width parameter
        * R_vI => [float] the pos on traj I
        * R_vJ => [float] the pos on traj J

    Outputs:
        * single point in Gaussian distribution
    """
    RI_RJ = (R_vI - R_vJ)

    preFact = (1/(2*PI*(sigma_vJ**2)))**(1.5)
    exponent = -(RI_RJ**2) / (2 * (sigma_vJ**2))
    return preFact * np.exp(exponent)

sigma, mean = 2, 1
x = np.arange(-10, 10, 0.01)
nuclDens = gaussian(sigma, x, mean)

QM = -np.gradient(nuclDens)/nuclDens


fig, ax = plt.subplots(2)
ax[0].plot(x, nuclDens, 'b-')
ax[1].plot(x, QM, 'b-')

# Make thing pretty
ax[0].axvline(mean, lw=0.7, color='k')
ax[0].axhline(0, lw=0.7, color='k')

ax[1].axvline(mean, lw=0.7, color='k')
ax[1].axhline(0, lw=0.7, color='k')

ax[0].grid("off")
ax[1].grid("off")

# Put in the color
posX, posQM = x[QM > 0], QM[QM > 0]
negX, negQM = x[QM < 0], QM[QM < 0]
ax[1].fill_between(posX, posQM, color='r', alpha=0.3)
ax[1].fill_between(negX, negQM, color='g', alpha=0.3)

# Pop on the arrows
ax[1].arrow( (min(negX)+mean)/2-2,  0, 1, 0,
            head_length=0.5, color='k')
ax[1].annotate(r"$\mathbf{F}_{qm}$", ((min(negX)+mean)/2-1.8, min(negQM)/5), fontsize=20)

ax[1].arrow( (max(posX)+mean)/2+2,  0, -1, 0,
            head_length=0.5, color='k')
ax[1].annotate(r"$\mathbf{F}_{qm}$", ((max(posX)+mean)/2 +1, max(posQM)/5), fontsize=20)

ax[1].annotate(r"$\mathbf{F}_{qm} \propto -\mathcal{Q}$", (min(negX) , max(posQM)/1.2), fontsize=22)

ax[0].set_ylabel(r"$|\chi|^2$")
ax[1].set_ylabel(r"$\mathcal{Q}$")
ax[1].set_xlabel("Nucl. Geom")

ax[0].set_title("1D case of Quantum Momentum for Gaussian Nuclear Density")
plt.show()
