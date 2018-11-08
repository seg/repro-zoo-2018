import numpy as np


def eperm_ema(por, eperm1, eperm2, aratio, eguess, eperm3=None, sw=None):
    """Effective electric permittivity using EMA.

    Markov et al., 2012, Journal of Applied Geophysics, Eq. 1.

    Recursive formula.

    Parameters
    ----------
    por : float or array
        Concentration of constituent 1 (host, wetting phase).

    eperm1, eperm2, eperm3 : float or array
        Electric permittivity of constituent 1, 2, 3.
        eperm3 is optional and corresponds to oil; however, if provided one
        has also to provide sw.
        Generally, eperm1=matrix, eperm2=water, eperm3=hydrocarbon.

    sw : float or array
        Water saturation (-), only required if eperm3 is provided.

    aratio : array
        Aspect ratio.

    eguess: float or array
        Initial guess of dielectric permittivity for recursion.

    Returns
    -------
    eperm: float or array
        Effective dielectric permittivity.

    """
    # Check and cast input
    c = _conc(por, sw)
    e = _stack(eperm1, eperm2, eperm3)
    eg = np.atleast_2d(eguess)
    if e.shape[0] == 1:
        eg = eg.transpose()

    def recursive(eperm, guess, dpl, conc, tol=1e-7):
        """Recursion formula to solve for effective electric permittivity.

        Alternatively we could use a root-finding algorithm.
        """

        # Calculate effective electric permittivity for this guess
        R = np.sum(1/(dpl*eperm[:, None] + (1 - dpl)*guess), axis=1)
        effe = np.sum(conc*eperm*R)/np.sum(conc*R)

        # If error above tolerance, call it again with new guess
        if np.abs(guess - effe) > tol:
            effe = recursive(eperm, effe, dpl, conc, tol)

        return effe

    # Get depolarization factors
    dpl = fdepol(aratio)

    # Loop over porosities
    ema = np.zeros((c.shape[0], e.shape[0]), dtype=e.dtype)
    for ci in range(c.shape[0]):
        for ei in range(e.shape[0]):

            # Calculate effective electric permittivity
            ema[ci, ei] = recursive(e[ei, :], eg[ci, ei], dpl[ci, :], c[ci, :])

    return np.squeeze(ema)


def eperm_crim(por, eperm1, eperm2, eperm3=None, sw=None):
    """Effective electric permittivity after CRIM.

    Markov et al., 2012, Journal of Applied Geophysics, Eq. 7.

    Parameters
    ----------
    por: float or array
        Concentration of constituent 1 (host, wetting phase).

    eperm1, eperm2, eperm3 : float or array
        Electric permittivity of constituent 1, 2, 3.
        eperm3 is optional and corresponds to oil; however, if provided one
        has also to provide sw.
        Generally, eperm1=matrix, eperm2=water, eperm3=hydrocarbon.

    sw : float or array
        Water saturation (-), only required if eperm3 is provided.


    Returns
    -------
    eperm: float or array
        Effective dielectric permittivity.

    """
    c = _conc(por, sw)
    e = _stack(eperm1, eperm2, eperm3)

    return np.squeeze(np.dot(c, np.sqrt(e.transpose()))**2)


def aspratio(por, eq='RHG'):
    """Return aspect ratios from Markov et al., 2012.

    Return aspect ratios as given in Markov et al., 2012, Journal of Applied
    Geophysics, Appendix 4.

    Parameters
    ----------
    por : float or array
        Porosity.

    eq : string
        One of:
        - 'RHG': Raymer, Hunt, and Gardner.
        - 'Willie': Willie.
        - 'Markov': Markov et al., 2012.
        - 'spherical': Spherical inclusions.

    Returns
    -------
    aratio : array
        Aspect ratios.

    """
    conc = _conc(por)
    aratio = np.zeros((conc.shape[0], 4))
    p = conc[:, 1]

    if eq == 'RHG':  # Raymer, Hunt, Gardner
        aratio[:, 0] = 0.99
        aratio[:, 1] = -2.1122*p*p + 1.4149*p + 0.0126
        aratio[:, 2] = -1.9848*p*p + 1.1307*p + 0.0093
        aratio[:, 3] = -0.2288*p*p + 0.1971*p + 0.0049

    elif eq == 'Willie':  # Willie
        aratio[:, 0] = 0.99
        aratio[:, 1] = -2.6618*p*p + 1.824*p - 0.008
        aratio[:, 2] = -1.5033*p*p + 0.9166*p + 0.0032
        aratio[:, 3] = -0.3691*p*p + 0.2963*p + 0.0034

    elif eq == 'Markov':  # Markov et al., 2012
        aratio[:, 0] = 0.99
        aratio[:, 1] = -1.98*p*p + 1.38*p + 0.01
        aratio[:, 2] = -1.98*p*p + 1.13*p + 0.01
        aratio[:, 3] = -0.23*p*p + 0.20*p + 0.005

    else:  # Spherical Inclusions
        aratio[:, 0] = 0.99
        aratio[:, 1] = 0.985
        aratio[:, 2] = 0.99
        aratio[:, 3] = 0.985

    return aratio


def fdepol(aratio):
    """Depolarization factor of ellipsoid without sigma corr.

    Using Carlson's Integral (`carlson_rd`).

    Parameters
    ----------
    aration : array
        Aspect ratios

    Returns
    -------
    depol : array
        Depolarization factors

    """
    aratio = np.array(aratio, ndmin=2)
    dpl = np.zeros((aratio.shape[0], int(aratio.shape[1]/2), 3))

    for ci in range(dpl.shape[0]):
        for i in range(dpl.shape[1]):
            ar1 = aratio[ci, 2*i]
            ar2 = aratio[ci, 1+2*i]

            L1 = carlson_rd(ar1*ar1, ar2*ar2, 1.0)
            L2 = carlson_rd(1.0, ar2*ar2, ar1*ar1)
            L3 = carlson_rd(1.0, ar1*ar1, ar2*ar2)

            dpl[ci, i, :] = ar1*ar2/3*np.array([L1, L2, L3])

    return dpl


def carlson_rd(x, y, z, errtol=1e-4):
    r"""Computes Carlson's elliptic integral of the second kind, Rd(x,y,z).

    Adjusted and modified from https://github.com/nipy/dipy, file dki.py
    (reconst/dki.py); accessed 30/06/2017; BSD licensed.

    .. math::
        R_D = \frac{3}{2} \int_{0}^{\infty} (t+x)^{-\frac{1}{2}}
        (t+y)^{-\frac{1}{2}}(t+z)  ^{-\frac{3}{2}}

    Parameters
    ----------
    x, y, z : float
        First, second, and third independent variable of the integral.

    errtol : float
        Error tolerance. Integral is computed with relative error less in
        magnitude than the defined value.

    Returns
    -------
    rd : float
        Value of the incomplete second order elliptic integral.

    Note
    -----
    x, y, and z have to be non-negative and at most x or y is zero.
    """
    # Cast and copy arrays
    xt = np.array(x).copy()
    yt = np.array(y).copy()
    zt = np.array(z).copy()

    # Initialize parameters
    tsum = 0.
    fac = 1.
    ave = (xt + yt + 3*zt)/5
    Q = (errtol/4)**(-1/6)*np.max(np.abs([ave - xt, ave - yt, ave - zt]))

    # Loop until converges
    while fac*Q > abs(ave):

        # Roots
        xtroot = np.sqrt(xt)
        ytroot = np.sqrt(yt)
        ztroot = np.sqrt(zt)

        lambd = xtroot*(ytroot + ztroot) + ytroot*ztroot

        tsum += fac/(ztroot*(zt + lambd))

        fac /= 4
        xt = (xt + lambd)/4
        yt = (yt + lambd)/4
        zt = (zt + lambd)/4
        ave = (ave + lambd)/4

    # Post loop calculation
    dxy = (ave - xt)/ave*(ave - yt)/ave
    dz = (ave - zt)/ave
    A = dxy - dz*dz
    B = dxy - 6*dz*dz
    C = 2*A + B
    D = B*(-3/14 + 9/88*B - 9/52*dz*C)
    E = dz*(1/6*C + dz*(-9/22*A + 3/26*dz*dxy))

    return 3*tsum + fac*(1 + D + E)/(ave*np.sqrt(ave))


def _stack(matrix, fluid1, fluid2=None):
    """Stacks electric permeability of matrix and fluids.

    Creates matrix where each one is repeated so the output has all possible
    combinations.

    Parameters
    ----------
    matrix : float or array
        Electric permeability of matrix.

    fluid1, fluid2 : float or array; fluid2 is optional
        Electric permeability of fluids.

    Returns
    -------
    eperm : array
        Electric permeability
        Dimension: (x, y):
          - x = #matrix*#fluid1*#fluid2
          - y = 3 if fluid2, else 2

    """
    # Cast input
    m = np.array(matrix, ndmin=1)
    f1 = np.array(fluid1, ndmin=1)

    # If array is identical, reduce to one value
    if np.all(m[:-1]-m[1:] == 0):
        m = np.array([m[0]])
    if np.all(f1[:-1]-f1[1:] == 0):
        f1 = np.array([f1[0]])

    # Return depending if `fluid2` or not
    if fluid2:
        f2 = np.array(fluid2, ndmin=1)
        if np.all(f2[:-1]-f2[1:] == 0):
            f2 = np.array([f2[0]])
        eperm = np.column_stack((np.tile(m, f1.size*f2.size),
                                f1.repeat(m.size*f2.size),
                                np.tile(f2.repeat(m.size), f1.size)))
        return eperm
    else:
        return np.column_stack((np.tile(m, f1.size), f1.repeat(m.size)))


def _conc(por, sw=None):
    """Checks concentrations.

    Parameters
    ----------
    por : float or array
        Porosity.

    sw : float
        Water saturation.

    Returns
    -------
    conc : array (2D matrix)
        Concentrations:
        - out[:, 0]: matrix percentage
        - out[:, 1]: water percentage
        - out[:, 2]: hydrocarbon percentage, only if `sw` given

    """
    # Cast por
    p = np.array(por, ndmin=1)

    # Return depending if `sw` or not
    if sw is None:
        return np.column_stack((1-p, p))
    else:
        return np.column_stack((1-p, p*sw, p*(1-sw)))
