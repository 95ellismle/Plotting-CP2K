import numpy as np


def calc_Rabi_Osc(t, H):
   Hab = H[0, 1]
   delE = H[1, 1] - H[0, 0]
   Hab4 = 4*Hab**2
   
   omegaR = np.sqrt(delE**2 + Hab4)

   pops = Hab4
   pops /= (delE**2 + Hab4)
   pops *= np.sin(0.5 * omegaR * t)**2
   pops = 1 - pops
   
   return pops


class Rabi(object):
    """
    Will plot the Quantum Momentum (in the QM_0 format) against time.

    Inputs:
        axes  => a list of the axes to plot on. The first item should be the
                 widget axis the second will be the axis to plot the data.
    """
    def __init__(self, axes):
        if self.plot:
            Rabi.widget_ax = axes[0]
            Rabi.plot_ax = axes[1]

            Rabi.plot_all(self)

        Rabi.plot_ax.set_ylabel(r"$|u_{RABI}|^2$")

    @staticmethod
    def plot_all(self):
        """
        Will plot the quantum momentum for all the replicas.
        """
        for fname in self.all_ham_data:
            allH, cols, t = self.all_ham_data[fname]
            H = allH[0]
            E, U = np.linalg.eigh(H)
            delE = E[1] - E[0]
            pops1 = calc_Rabi_Osc(t, H)
            pops2 = 1-pops1

            Rabi.plot_ax.plot(t, pops1, lw=0.7, alpha=self.alpha,
                              color=self.colors[0])
            Rabi.plot_ax.plot(t, pops2, lw=0.7, alpha=self.alpha,
                              color=self.colors[1])
	


