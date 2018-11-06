import numpy as np
from scipy.constants import epsilon_0


def eps_w_klsw(frequency, temperature, salinity, resistivity=None, e_8=4.9,
               alpha=0., dw=None):
    """Complex dielectric permittivity of water after Klein and Swift, 1977.

    Note
    ----
    Salinity OR resistivity can be a function:
        - salinity(resistivity, temperature) or
        - resistivity(salinity, temperature).
    If resistivity is not provided, salinity has to be provided as a value, and
    `res_sal_klsw` is used.

    Parameters
    ----------
    frequency : float or array
        Frequency (Hz).

    temperature : float or array
        Temperature (°C).

    salinity : float, array; or function
        Salinity (ppm). Can be one of the `sal_res_***`-functions.

    resistivity : float, array; or function
        Resistivity (Ohm.m). Can be one of the `res_sal_***`-functions.

    e_8 : float
        Relative permittivity at infinite frequency (-).
        Defaults to 4.9 as given in paper.

    alpha : float
        Relaxation exponent (i omega tau)^(1-alpha).
        Defaults to 0.

    dw : bool
        Flag if distilled water or not, for the calculation of the static
        electric permittivity. If True, neither salinity nor resistivity is
        required.
        Defaults to False.

    Returns
    -------
    e_w: float or array
        Dielectric permittivity of water.

    """
    # Cast frequency
    frequency = np.array(frequency, ndmin=1)

    # If salinity is a function, calculate salinity as a function of
    # resistivity and temperature.
    if callable(salinity):
        salinity = salinity(resistivity, temperature)

    # Resistivity
    if dw:
        # Distilled water
        resistivity = 10000
    elif callable(resistivity):
        # If resistivity is a function, calculate resistivity as a function of
        # salinity and temperature.
        resistivity = resistivity(salinity, temperature)
    elif not np.any(resistivity):
        # If resistivity is not provided, calculate resistivity as in the
        # original Klein and Swift model.
        resistivity = res_sal_klsw(salinity, temperature)

    # Convert salinity from ppm to ppk.
    tS = salinity/1000

    if dw:
        # Static electric permittivity e_s for distilled water; Equation 8
        e_s = 88.045 - 0.4147*temperature + 6.295e-4*temperature**2
        e_s += 1.075e-5*temperature**3
    else:
        # Static electric permittivity e_s; Equations 13-15
        e_sT = 87.134 - 1.949e-1*temperature - 1.276e-2*temperature**2
        e_sT += 2.491e-4*temperature**3
        a_ST = 1 + 1.613e-5*tS*temperature - 3.656e-3*tS
        a_ST += 3.21e-5*tS**2 - 4.232e-7*tS**3
        e_s = e_sT*a_ST

    # Tau; Equations 16-18
    tau_T = 1.768e-11 - 6.086e-13*temperature + 1.104e-14*temperature**2
    tau_T -= 8.111e-17*temperature**3

    b_ST = 1 + 2.282e-5*tS*temperature - 7.638e-4*tS
    b_ST += - 7.760e-6*tS**2 + 1.105e-8*tS**3
    tau = tau_T*b_ST

    # Complex dielectric constant of water; Equation 5
    omega = 2*np.pi*frequency
    e_w = e_8 + (e_s - e_8)/(1 + (1j*omega*tau)**(1-alpha))
    e_w -= 1j/(omega*epsilon_0*resistivity)

    return e_w


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


def eperm_maga(por, eperm1, eperm2, aratio, eperm3=None, sw=None):
    """Effective electric permittivity after Maxwell-Garnett.

    Seleznev et al. (2006), SPE, page 4.

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

    Returns
    -------
    eperm: float or array
        Effective dielectric permittivity.

    """

    # Get inclusion fractions and permittivities
    c = _conc(por, sw)
    e = _stack(eperm1, eperm2, eperm3)

    # Calculate CRIM background model
    ecri = np.atleast_1d(eperm_crim(por, eperm1, eperm2, eperm3, sw))

    # Get depolarization factors
    dpl = fdepol(aratio)

    # Calculate inner sums
    ne = (e-ecri[:, None])[:, :, None]*dpl[0, :, :][None, :, :]
    seen = np.sum(ecri[:, None, None]/(ecri[:, None, None] + ne), 2)
    snen = np.sum(dpl[0, :, :][None, :, :]/(ecri[:, None, None] + ne), 2)

    # Calculate outer sums
    num = np.sum(c*(e-ecri[:, None]) * seen, 1)
    den = 3-np.sum(c*(e-ecri[:, None]) * snen, 1)

    # Return effective permittivity
    return np.squeeze(ecri + num/den)


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


def sal_res_crain(resistivity, temperature):
    """Return salinity following Crain's petrophysical handbook.

    Salinity is calculated as a function of resistivity and temperature.

    Note
    ----
    Resistivity OR temperature can be arrays, but not both at the same time.

    Parameters
    ----------
    resistivity : float or array
        Resistivity in Ohm.m.

    temperature : float or array
        Temperature in °C.

    Returns
    -------
    salinity : float or array
        Salinity in parts per million (ppm).

    """
    return 400000/(temperature*9/5 + 32)/resistivity**1.14


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
