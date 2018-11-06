import numpy as np
from scipy import stats
from copy import deepcopy as dc


class GKDE(object):
    """Returns an object of `scipy.stats.kde.gaussian_kde`.

    Parameters
    ----------
    data : array
        Data.

    Returns
    -------
    pdf : probability density function
        A `scipy.stats.kde.gaussian_kde` instance, containing:

    """


    def __init__(self, data):
        self.data = data
        """data"""
        self.gkde = stats.kde.gaussian_kde(data.reshape(-1)) #Make sure 1D array
        """`stats.kde.gaussian_kde(data)`"""
        # Cavariance factor is Scotts factor
        self.b_w = self.gkde.covariance_factor()
        """bandwith - `self.gkde.covariance_factor()`"""
        self.resample = self.gkde.resample
        """`self.gkde.resample`"""

    def adband(self, fac=1):
        """Adjust bandwith.

        Parameters
        ----------
        fac : scalar, optional; <1>
            Factor by which to adjust the bandwith.

        """

        def covariance_factor(self, b_w=self.b_w, fac=fac):
            """Define new bandwith."""
            return fac*b_w

        setattr(self.gkde, 'covariance_factor',
                covariance_factor.__get__(self.gkde,    type(self.gkde)))
        self.gkde._compute_covariance()

        return self.gkde

    def bins(self, bins=100):
        """Define bins from min(data) to max(data).

        Parameters
        ----------
        bin : integer, optional; <100>
            Number of bins.

        """

        return np.linspace(self.data.min(), self.data.max(), int(bins))

    def __call__(self, values):
        return self.gkde(values)


def rho_poup(data):
    """Resistivity with "Indonesia formula".

    Poupon and Leveaux (1971), "Indonesia formula", [TAI.54.Poupon]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por, vsh : scalar or vector
        - m_e : scalar or vector, optional; <2.>
        - a_f : scalar or vector, optional; <1.>
        - s_w : scalar or vector, optional; <1.>
        - n_e : scalar or vector, optional; <2.>

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    """
    tdata = dc(data)

    try:
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
        por   = tdata['por']
        vsh   = tdata['vsh']
    except NameError:
        raise
    m_e = tdata.get('m_e', np.array(2.))
    a_f = tdata.get('a_f', np.array(1.))
    s_w = tdata.get('s_w', np.array(1.))
    n_e = tdata.get('n_e', np.array(2.))

    alpha = 1. - vsh/2.

    tpor = np.array(por, dtype=float, copy=True, ndmin=1)
    tpor[tpor == 0.] = np.NaN

    rho_b = (np.sqrt(s_w**n_e)*(np.sqrt(tpor**m_e/(a_f*rho_f)) +
             vsh**alpha/np.sqrt(rho_s)))**(-2)

    return np.nan_to_num(rho_b)
