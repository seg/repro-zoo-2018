r"""
`vel2res` -- Velocity to resistivity
====================================

The approach for my project is generally from seismic information to
resistivity prediction, using porosity as link,

.. math:: \rho = f(\phi),\quad  \phi = g(v_p) \quad\Rightarrow\quad
    \rho(\phi[v_p]) \ .

However, many models are outlined the other way round, and this file hence
contains also cross-correlations in the manner of:

.. math:: v_p = h(\phi),\quad \phi = k(\rho) \ .

Furthermore, relations are commonly defined in terms of conductivities.
As we are interested in resistivities in hydrocarbon exploration, I
try to express relations in terms of resistivities.

:Note: This module could be much better annotated, with formulas for the
    different models. Also, most models I coded up, but never extensively
    used and hence tested it, except the ones I use in the thesis.


:1.A Resistivity from porosity:

::

    rho_arit        -- Resistivities with arithmetic mean.
    rho_geom        -- Resistivities with geometric mean.
    rho_harm        -- Resistivities with harmonic mean.
    rho_arch        -- Resistivity with Archie.
    rho_herm        -- Resistivity with Hermance.
    rho_glov        -- Resistivity with Glover.
    rho_self        -- Resistivity with the self-similar model.
    rho_crim        -- Resistivity with CRIM.
    rho_hsbs        -- Resistivity with Hashin-Shtrikman bounds.
    rho_hslb        -- Resistivity with Hashin-Shtrikman bounds - lower bound.
    rho_hsub        -- Resistivity with Hashin-Shtrikman bounds - upper bound.
    rho_poup        -- Resistivity with "Indonesia formula".
    rho_hsub2       -- Resistivity with HS lower bound, Berryman.

:1.B Velocity from porosity:

::

    vp_arit         -- P-wave velocity with arithmetic mean.
    vp_geom         -- P-wave velocity with geometric mean.
    vp_harm         -- P-wave velocity with harmonic mean (time-average eq.)
    vp_raym         -- P-wave velocity with Raymer.
    vp_hsbs         -- P-wave velocity with Hashin-Shtrikman bounds.
    vp_hslb         -- P-wave velocity with Hashin-Shtrikman lower bound.
    vp_hsub         -- P-wave velocity with Hashin-Shtrikman upper bound.
    vp_gass         -- P-wave velocity using Gassmann relations.
    vp_aff          -- P-wave velocity with AFF.

:2.A Porosity from resistivity:

::

    por_r_arit      -- Porosity with arithmetic mean from resistivities.
    por_r_geom      -- Porosity with geometric mean from resistivities.
    por_r_harm      -- Porosity with harmonic mean from resistivities.
    por_r_arch      -- Porosity with Archie from resistivities.
    por_r_herm      -- Porosity with Hermance from resistivities.
    por_r_hsbs      -- Porosity with Hashin-Shtrikman bounds from resistivity.
    por_r_hslb      -- Porosity with HS lower bound from resistivity.
    por_r_hsub      -- Porosity with HS upper bound from resistivity.
    por_r_hsub2     -- Porosity with HS lower bound from resistivity, Berryman.
    por_r_self      -- Porosity with the self-similar model from resistivities.
    por_r_crim      -- Porosity with CRIM from resistivities.
    por_r_dems      -- Porosity with DEM from resistivities.

:2.B Porosity from velocity:

::

    por_v_arit      -- Porosity with arithmetic mean from P-wave velocities.
    por_v_geom      -- Porosity with geometric mean from P-wave velocities.
    por_v_harm      -- Porosity with harmonic mean from P-wave velocities.
    por_v_raym      -- Porosity with Raymer from P-wave velocities.
    por_v_hsbs      -- Porosity with HS bounds from P-wave velocities.
    por_v_hslb      -- Porosity with HS lower bounds from P-wave velocities.
    por_v_hsub      -- Porosity with HS upper bounds from P-wave velocities.
    por_v_gass      -- Porosity with Gassmann from P-wave velocities.
    por_v_aff       -- Porosity with AFF from P-wave velocities.

:3. Cross-relations via porosity:

::

    in2por2out      -- Cross-property vel->por->res or res->por->vel.

:4. Cross-relations directly:

::

    rho_faus        -- Resistivity with Faust.
    vp_faus         -- P-wave velocity with Faust.

:5. Relations for other parameters:

::

    rhof_cec        -- Temperature dependent rho_f values.
    rhof_mol        -- Temperature and molatility dependent rho_f values.
    m_e_folke       -- Porosity corrected cementation-exponent.
    param_depth     -- Depth dependent parameters.

:6. Other necessary relations:

::

    den_bulk        -- Bulk density from fluid and grain densities.
    por_dens        -- Porosity from fluid and grain densities.
    vp_modu         -- P-wave velocity from bulk and shear moduli.
    k_mu_hsbs       -- Bulk and shear moduli with HS bounds.
    k_gass          -- Bulk modulus with Gassmann.
    k_mu_krie       -- Bulk and shear moduli with Krief.

:Variable units (unless otherwise noted):

::

  rho_b  (Omega.m) Resistivity
  rho_s  (Omega.m)                 of the matrix
  rho_f  (Omega.m)                 of the fluid

  vp_b   (km/s)    P-wave velocity
  vp_s   (km/s)                    of the matrix
  vp_f   (km/s)                    of the fluid

  k_b    (GPa)     Bulk modulus
  k_s    (GPa)                     of the solid
  k_f    (GPa)                     of the fluid

  mu_b   (GPa)     Shear modulus
  mu_s   (GPa)                     of the solid
  mu_f   (GPa)                     of the fluid

  den_b  (g/cm3)   Density
  den_s  (g/cm3)                   of the matrix
  den_f  (g/cm3)                   of the fluid

  m_e    (-)       m-exponent (Archie, Hermance, Self-similar, AFF)
  a_f    (-)       a-factor   (Archie)
  a_k    (-)       a-factor   (Krief)
  p_e    (-)       p-exponent (Glover)
  y_e    (-)       y-exponent (CRIM)

  por    (-)       Porosity

  depth  (km)      Depth


:References:

.. [B.AGU.95.Berryman] Berryman, J. G.,  1995, Mixture theory for rock
    properties, in Rock Physics \& Phase Relations: A Handbook of Physical
    Constants: AGU, 3,  205--228, http://dx.doi.org/10.1029/RF003.
.. [B.CUP.09.Mavko] Mavko, G., T. Mukerji, and J. Dvorkin,  2009, The Rock
    Physics Handbook: Cambridge University Press Cambridge,
    http://www.cambridge.org/9780521861366.
.. [B.S.07.Ellis] Ellis, D. V., and J. M. Singer,  2007, The Logging for
    Earth Scientists, 2 ed.: Springer, ISBN: 978-1-4020-3738-2.
.. [TAI.42.Archie] Archie, G. E.,  1942, The electrical resistivity log as an
    aid in determining some reservoir characteristics: Trans. AIME,  54--62,
    doi: 10.2118/942054-G.
.. [JPT.84.Clavier] Clavier, C., G. Coates, and J. Dumanoir,  1984, Theoretical
    and experimental bases for the dual-water model for interpretation of shaly
    sands: Journal of Petroleum Technology, 24, 153--168,
    http://dx.doi.org/10.2118/6859-PA.
.. [SEG.10.Engelmark] Engelmark, F.,  2010, Velocity to resistivity transform
    via porosity: SEG Technical Program Expanded Abstracts, 29, 2501--2505,
    http://library.seg.org/doi/abs/10.1190/1.3513358.
.. [GEO.53.Faust] Faust, L. Y.,  1953, A velocity function including lithologic
    variation: Geophysics, 18, 271--288, http://dx.doi.org/10.1190/1.1437869.
.. [AAPGB.92.Issler] Issler, D. R.,  1992, A new approach to shale compaction
    and stratigraphic restoration, Beaufort-Mackenzie Basin and Mackenzie
    Corridor, Northern Canada: American Association of Petroleum Geologists
    Bulletin, 76, 1170--1189.
.. [GEO.92.Sen] Sen, P. N., and P. A. Goode,  1992, Influence of temperature on
    electrical conductivity on shaly sands: Geophysics, 57, 89--96,
    http://dx.doi.org/10.1190/1.1443191.
.. [GEO.81.Sen] Sen, P. N., C. Scala, and M. H. Cohen,  1981, A self-similar
    model for sedimentary-rocks with application to the dielectric constant of
    fused glass-beads: Geophysics, 46, 781--795,
    http://dx.doi.org/10.1190/1.1441215.
.. [GRL.79.Hermance] Hermance, J. F.,  1979, The electrical conductivity of
    materials containing partial melt: A simple model from Archie's law:
    Geophysical Research Letters, 6, 613--616,
    http://dx.doi.org/10.1029/GL006i007p00613.
.. [EPSL.00.Glover] Glover, P. W. J., M. J. Hole, and J. Pous,  2000, A
    modified Archie’s law for two conducting phases: Earth and Planetary
    Science Letters, 180, 369--383,
    http://dx.doi.org/10.1016/S0012-821X(00)00168-0.
.. [GEO.1983.Bussian] Bussian, A. E.,  1983, Electrical conductance in a porous
    medium: Geophysics, 48, 1258--1268, http://dx.doi.org/10.1190/1.1441549.
.. [B.PER.96.Schon] Schön, J.H.,  1996, Physical properties of rocks:
    Fundamentals and principles of petrophysics: Pergamon Press, ISBN:
    978-0080410081.
.. [SPWLA.80.Raymer] Raymer, L. L., E. R. Hunt, and J. S. Gardner,  1980, An
    improved sonic transit time-to-porosity transform: Presented at the SPWLA
    21st Annual Logging Symposium, SPWLA.
.. [MPS.63.Hashin] Hashin, Z., and S. Shtrikman,  1963, A variational approach
    to the theory of the elastic behaviour of multiphase materials: Journal of
    the Mechanics and Physics of Solids, 11, 127--140,
    http://dx.doi.org/10.1016/0022-5096(63)90060-7.
.. [TAI.54.Poupon] Poupon, A., M. E. Loy, and M. P. Tixier,  1954, A
    contribution to electrical log interpretation in shaley sands: Trans. AIME,
    138--145.
.. [NGZ.51.Gassmann] Gassmann, F.,  1951, Über die Elastizität poröser Medien:
    Vier. der Natur. Gesellschaft Zürich,  1--23.
.. [TLA.88.raiga-clemenceau] Raiga-Clemenceau, J., J. P. Martin, and S.
    Nicoletis,  1988, The concept of acoustic formation factor for more
    accurate porosity determination from sonic transit time data: The Log
    Analyst, 29, 54--60,
    http://www.onepetro.org/mslib/app/Preview.do?paperNumber=SPWLA-1988-v29n1a4}.
.. [TLA.90.Krief] Krief, M., J. Garat, J. Stellingwerff, and J. Ventre,  1990,
    A petrophysical interpretation using the velocities of P and S waves
    (full-waveform sonic): The Log Analyst, 31, 355--369,
    http://www.onepetro.org/mslib/app/Preview.do?paperNumber=SPWLA-1990-v31n6a2&societyCode=SPWLA.
.. [GEO.07.Carcione] Carcione, J. M., B. Ursin, and J. I. Nordskag, 2007,
    Cross-property relations between electrical conductivity and the
    seismic velocity of rocks: Geophysics, 72,
    http://dx.doi.org/10.1190/1.2762224.
.. [GP.11.Aversana] Dell’Aversana, P., G. Bernasconi, F. Miotti, and D.
    Rovetta, 2011, Joint inversion of rock properties from sonic,
    resistivity and density well-log measurements: Geophysical Prospecting,
    59, 1144–1154, http://dx.doi.org/10.1111/j.1365-2478.2011.00996.x.
.. [GP.09.Chen] Chen, J., and T. A. Dickens,  2009, Effects of uncertainty in
    rock-physics models on reservoir parameter estimation using seismic
    amplitude variation with angle and controlled-source electromagnetics data:
    Geophysical Prospecting, 57, 61--74,
    http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2478.2008.00721.x/abstract.
..  [JSS.10.Patil] Patil, A., D. Huard, and C. J. Fonnesbeck,  2010, PyMC:
    Bayesian stochastic modelling in python: Journal of Statistical Software,
    35, 1--81, http://www.jstatsoft.org/v35/i04.


"""
_all__ = [#
        # 1.A Resistivity from porosity
        'rho_arit', 'rho_geom', 'rho_harm', 'rho_arch', 'rho_herm',
        'rho_glov', 'rho_self', 'rho_crim', 'rho_hsbs', 'rho_hslb',
        'rho_hsub', 'rho_poup', 'rho_hsub2',
        # 1.B Velocity from porosity
        'vp_arit', 'vp_geom', 'vp_harm', 'vp_raym', 'vp_hsbs', 'vp_hslb',
        'vp_hsub', 'vp_gass', 'vp_aff',
        # 2.A Porosity from resistivity
        'por_r_arit', 'por_r_geom', 'por_r_harm', 'por_r_arch', 'por_r_herm',
        'por_r_hsbs', 'por_r_hslb', 'por_r_hsub', 'por_r_hsub2', 'por_r_self',
        'por_r_crim', 'por_r_dems',
        # 2.B Porosity from velocity
        'por_v_arit', 'por_v_geom', 'por_v_harm', 'por_v_raym', 'por_v_hsbs',
        'por_v_hslb', 'por_v_hsub', 'por_v_gass', 'por_v_aff',
        # 3. Cross-relations via porosity
        'in2por2out',
        # 4. Cross-corelations directly
        'rho_faus', 'vp_faus',
        # 5. Relations for other parameters
        'rhof_cec', 'rhof_mol', 'm_e_folke', 'param_depth',
        # 6. Other necessary relations
        'den_bulk', 'por_dens', 'vp_modu', 'k_mu_hsbs', 'k_gass', 'k_mu_krie',
        ]

import numpy as np
from copy import deepcopy as dc
from scipy.optimize import brentq

#  1.A Resistivity from porosity                                    1A

def rho_arit(data):
    """Resistivities with arithmetic mean.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_arit(data)
    array([ 8.125,  7.75 ])

    """
    tdata = dc(data)

    try:
        rho_f = tdata['rho_f']
        rho_s = tdata['rho_s']
        por   = tdata['por']
    except NameError:
        raise

    return por*rho_f + (1.-por)*rho_s


def rho_geom(data):
    """Resistivities with geometric mean.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_geom(data)
    array([ 7.07106781,  6.59753955])

    """
    tdata = dc(data)

    try:
        rho_f = tdata['rho_f']
        rho_s = tdata['rho_s']
        por   = tdata['por']
    except NameError:
        raise

    return rho_f**por * rho_s**(1.-por)


def rho_harm(data):
    """Resistivities with harmonic mean.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_harm(data)
    array([ 5.71428571,  5.26315789])

    """
    tdata = dc(data)

    try:
        rho_f = tdata['rho_f']
        rho_s = tdata['rho_s']
        por   = tdata['por']
    except NameError:
        raise

    return (por/rho_f + (1.-por)/rho_s)**-1.


def rho_arch(data):
    """Resistivity with Archie.

    [TAI.42.Archie]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_f, por : scalar or vector
        - m_e : scalar or vector, optional; <2.>
        - a_f : scalar or vector, optional; <1.>
        - s_w : scalar or vector, optional; <1.>
        - n_e : scalar or vector, optional; <2.>
        - flag : int, optional; {<0>, 1, 2}
            1: If rho > rho_s = set rho to rho_s;
            2: If rho > rho_s = set rho to 0.

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_arch(data)
    array([ 40.        ,  27.77777778])

    """
    tdata = dc(data)

    try:
        rho_f = tdata['rho_f']
        por   = tdata['por']
    except NameError:
        raise

    m_e = tdata.get('m_e', np.array(2.))
    a_f = tdata.get('a_f', np.array(1.))
    s_w = tdata.get('s_w', np.array(1.))
    n_e = tdata.get('n_e', np.array(2.))
    flag = tdata.get('flag', np.array(0.))
    if flag > 0:
        try:
            rho_s = tdata['rho_s']
        except NameError:
            raise

    tpor = np.array(por, dtype=float, copy=True, ndmin=1)
    tpor[tpor == 0.] = np.NaN

    rho = a_f*rho_f*tpor**(-m_e)*s_w**(-n_e)
    if flag == 1:
        rho[rho > rho_s] = rho_s
    elif flag == 2:
        rho[rho > rho_s] = np.NaN

    return np.nan_to_num(rho)


def rho_herm(data):
    """Resistivity with Hermance.

    [GRL.79.Hermance]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector
        - m_e : scalar or vector, optional; <2.>

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_herm(data)
    array([ 8.42105263,  7.87401575])

    """
    tdata = dc(data)

    try:
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
        por   = tdata['por']
    except NameError:
        raise
    m_e = tdata.get('m_e', np.array(2.))

    return (por**m_e/rho_f + (1-por**m_e)/rho_s)**(-1)


def rho_glov(data):
    """Resistivity with Glover.

    The Glover et al. (2000) model [EPSL.00.Glover]_ is a special case of the
    Hermance model [GRL.79.Hermance]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector
            Resistivities, porosity.
        - m_e : scalar or vector, optional; <2.>
            Cementation exponent.
        - p_e : scalar or vector, optional; <-1. = m_e>
            If < 0, p_e is set to m_e.
        - flag : int, optional; {<0>, 1, 2}
            2: If rho > rho_s = set rho to np.NaN

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_glov(data)
    array([ 12.30769231,  11.76470588])

    """
    tdata = dc(data)

    try:
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
        por   = tdata['por']
    except NameError:
        raise
    m_e = tdata.get('m_e', np.array(2.))
    p_e = tdata.get('p_e', np.array(-1.))
    flag = tdata.get('flag', np.array(0.))

    if p_e < 0.:
        p_e = m_e
    rho = (por**m_e/rho_f + (1-por)**p_e/rho_s)**(-1)
    if flag == 2:
        rho[rho > rho_s] = np.NaN
    return rho


def rho_self(data):
    """Resistivity with the self-similar model.

    [GEO.81.Sen]_, [GEO.1983.Bussian]_. This model is in terms of porosity
    as a function of resistivity. Resistivity is found using a root-fining
    algorithm.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector
        - m_e : scalar or vector, optional; <2.>

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_self(data)
    array([ 6.88777642,  6.4       ])

    """
    tdata = dc(data)

    try:
        por = np.array(tdata['por'], dtype=float, copy=True, ndmin=1)
        rho_s = a_if_b_scal(tdata['rho_s'], por)
        rho_f = a_if_b_scal(tdata['rho_f'], por)
        # if 's_w' in tdata:
        #     s_w = a_if_b_scal(tdata['s_w'], por)
    except NameError:
        raise
    m_e = a_if_b_scal(tdata.get('m_e', np.array(2.)), por)
    # if 's_w' in tdata:
    #     n_e = a_if_b_scal(tdata.get('n_e', np.array(2.)), por)
    #     rho_f = rho_f/(s_w**n_e)

    rho_b = np.array([brentq(lambda x: por_r_self({'rho_f':rho_f[x_i],
                    'rho_s':rho_s[x_i], 'rho_b':x, 'm_e':m_e[x_i]}) - por[x_i],
                    min(rho_s[x_i], rho_f[x_i]), max(rho_s[x_i], rho_f[x_i]))
                    for x_i in range(len(por))])

    return rho_b


def rho_crim(data):
    """Resistivity with CRIM.

    [B.PER.96.Schon]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector
        - y_e : scalar or vector, optional; <1./2.>

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_crim(data)
    array([ 6.4       ,  5.91715976])

    """
    tdata = dc(data)

    try:
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
        por   = tdata['por']
    except NameError:
        raise
    y_e = tdata.get('y_e', np.array(1./2))

    return (por/rho_f**y_e + (1.-por)/rho_s**y_e)**(-1./y_e)


def rho_hsbs(data):
    """Resistivity with Hashin-Shtrikman bounds.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector

    Returns
    -------
    rho_l, rho_u : scalar or vector
        Resistivity, lower and upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_hsbs(data)
    (array([ 6.25   ,  5.78125]), array([ 7.        ,  6.53846154]))

    """
    tdata = dc(data)

    try:
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
        por   = tdata['por']
    except NameError:
        raise

    rho_u = (((1.-por)*rho_s/3. + por*rho_f*rho_s/(2.*rho_f+rho_s))**(-1.)
             -2./rho_s)**(-1.)
    rho_l = (((1.-por)*rho_f*rho_s/(rho_f+2*rho_s) + por*rho_f/3.)**(-1)
             -2./rho_f)**(-1.)

    if rho_s < rho_f:
        rho_t = dc(rho_u)
        rho_u = dc(rho_l)
        rho_l = dc(rho_t)

    return rho_l, rho_u


def rho_hslb(data):
    """Resistivity with Hashin-Shtrikman bounds - lower bound.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector

    Returns
    -------
    rho_l : scalar or vector
        Resistivity, lower bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_hslb(data)
    array([ 6.25   ,  5.78125])

    """
    rho_l, _ = rho_hsbs(data)
    return rho_l


def rho_hsub(data):
    """Resistivity with Hashin-Shtrikman bounds - upper bound.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_s, rho_f, por : scalar or vector

    Returns
    -------
    rho_u : scalar or vector
        Resistivity, upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.rho_hsub(data)
    array([ 7.        ,  6.53846154])

    """
    _, rho_u = rho_hsbs(data)
    return rho_u


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

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> data['vsh'] = np.array([0.85, 0.7])
    >>> vel2res.rho_poup(data)
    array([ 5.02433185,  5.15289877])

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


def rho_hsub2(data):
    """Resistivity with HS lower bound, Berryman.

    Resistivity with HS lower bound, Eq. 13, Berryman (1995).  NOT SURE ABOUT
    THIS, JUST TO REPRODUCE CARCIONE!  [B.AGU.95.Berryman]_, [GEO.07.Carcione]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_f, por : scalar or vector

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> data['vsh'] = np.array([0.85, 0.7])
    >>> vel2res.rho_hsub2(data)
    array([ 13.75,  11.25])

    """
    tdata = dc(data)

    try:
        por   = tdata['por']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    tpor = np.array(por, dtype=float, copy=True, ndmin=1)
    tpor[tpor == 0.] = np.NaN

    rho = (3.*rho_f/tpor - rho_f)/2.

    return np.nan_to_num(rho)

#  1.B Velocity from porosity                                       1B

def vp_arit(data):
    """P-wave velocity with arithmetic mean.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_s, vp_f, por : scalar or vector

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_arit(data)
    array([ 3.375,  3.25 ])

    """
    tdata = dc(data)

    try:
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
        por   = tdata['por']
    except NameError:
        raise

    return por*vp_f + (1.-por)*vp_s


def vp_geom(data):
    """P-wave velocity with geometric mean.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_s, vp_f, por : scalar or vector

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_geom(data)
    array([ 3.13016916,  2.98036443])

    """
    tdata = dc(data)

    try:
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
        por   = tdata['por']
    except NameError:
        raise

    return vp_f**por * vp_s**(1.-por)


def vp_harm(data):
    """P-wave velocity with harmonic mean (time-average eq.)

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_s, vp_f, por : scalar or vector

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_harm(data)
    array([ 2.82352941,  2.66666667])

    """
    tdata = dc(data)

    try:
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
        por   = tdata['por']
    except NameError:
        raise

    return (por/vp_f + (1.-por)/vp_s)**-1.


def vp_raym(data):
    """P-wave velocity with Raymer.

    [SPWLA.80.Raymer]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_s, vp_f, por : scalar or vector

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_raym(data)
    array([ 2.625,  2.41 ])

    """
    tdata = dc(data)

    try:
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
        por  = tdata['por']
    except NameError:
        raise

    return (1-por)**2*vp_s + por*vp_f


def vp_hsbs(data):
    """P-wave velocity with Hashin-Shtrikman bounds.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:
            - k_s, k_f, por : scalar or vector
            - mu_s : scalar or vector, optional; <0.>

    Returns
    -------
    vp_l, vp_u : scalar or vector
        P-wave velocity, lower and upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_hsbs(data)
    (array([ 1.7794873 ,  1.70230748]), array([ 3.38865497,  3.26380572]))

    """
    tdata = dc(data)

    try:
        por   = tdata['por']
        k_s   = tdata['k_s']
        k_f   = tdata['k_f']
    except NameError:
        raise
    mu_s = tdata.get('mu_s', np.array(0.))

    k_l, mu_l = k_mu_hsbs({'k_s':k_s, 'k_f':k_f, 'por':por, 'mu_s':0})
    k_u, mu_u = k_mu_hsbs({'k_s':k_s, 'k_f':k_f, 'por':por, 'mu_s':mu_s})
    den_ = den_bulk(tdata)
    vp_l = vp_modu({'den_b':den_, 'k_b':k_l, 'mu_b':mu_l})
    vp_u = vp_modu({'den_b':den_, 'k_b':k_u, 'mu_b':mu_u})

    return vp_l, vp_u


def vp_hslb(data):
    """P-wave velocity with Hashin-Shtrikman lower bound.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - k_s, k_f, por : scalar or vector
        - mu_s : scalar or vector, optional; <0.>

    Returns
    -------
    vp_l : scalar or vector
        P-wave velocity, lower bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_hslb(data)
    array([ 1.7794873 ,  1.70230748])

    """
    vp_l, _ = vp_hsbs(data)
    return vp_l


def vp_hsub(data):
    """P-wave velocity with Hashin-Shtrikman upper bound.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - k_s, k_f, por : scalar or vector
        - mu_s : scalar or vector, optional; <0.>

    Returns
    -------
    vp_u : scalar or vector
        P-wave velocity, upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_hsub(data)
    array([ 3.38865497,  3.26380572])

    """
    _, vp_u = vp_hsbs(data)
    return vp_u


def vp_gass(data):
    """P-wave velocity using Gassmann relations.

    [NGZ.51.Gassmann]_, [GEO.07.Carcione]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - k_s, k_f, por, mu_s : scalar or vector
        - a_k : scalar or vector, optional; <3.>

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_gass(data)
    array([ 2.75897016,  2.47602718])

    """
    tdata = dc(data)

    try:
        por   = tdata['por']
        k_s   = tdata['k_s']
        k_f   = tdata['k_f']
        mu_s  = tdata['mu_s']
    except NameError:
        raise
    a_k = tdata.get('a_k', np.array(3.))

    k_m, mu_m = k_mu_krie({'k_s':k_s, 'por':por, 'mu_s':mu_s, 'a_k':a_k})
    den_ = den_bulk(tdata)
    k_g = k_gass({'k_s':k_s, 'k_m':k_m, 'k_f':k_f, 'por':por})
    vp_b = vp_modu({'den_b':den_, 'k_b':k_g, 'mu_b':mu_m})

    return vp_b


def vp_aff(data):
    """P-wave velocity with AFF.

    Acoustic formation factor, [TLA.88.raiga-clemenceau]_, by default with
    m_e, v_s from [AAPGB.92.Issler]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - por : scalar or vector
        - vp_s : scalar or vector, optional; <1000./220>
        - m_e : scalar or vector, optional; <2.19>

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.vp_aff(data)
    array([ 2.13031663,  1.83157497])

    """
    tdata = dc(data)

    try:
        por = tdata['por']
    except NameError:
        raise
    vp_s = tdata.get('vp_s', np.array(1000./220))
    if 'x_e' in tdata: # Backwards compatibility
        m_e  = tdata.get('x_e', np.array(2.19))
    else:
        m_e  = tdata.get('m_e', np.array(2.19))

    return (1. - por)**m_e*vp_s

#  2.A Porosity from resistivity                                    2A

def por_r_arit(data):
    """Porosity with arithmetic mean from resistivities.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_arit(data)
    array([ 1. ,  0.8, -0. ])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    return (rho_b - rho_s)/(rho_f - rho_s)


def por_r_geom(data):
    """Porosity with geometric mean from resistivities.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_geom(data)
    array([ 1.        ,  0.66096405, -0.        ])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    return (np.log(rho_b) - np.log(rho_s))/(np.log(rho_f) - np.log(rho_s))


def por_r_harm(data):
    """Porosity with harmonic mean from resistivities.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_harm(data)
    array([ 1. ,  0.5, -0. ])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    return (rho_f/rho_b)*((rho_b) - rho_s)/ (rho_f - rho_s)


def por_r_arch(data):
    """Porosity with Archie from resistivities.

    [TAI.42.Archie]_, porosity is limited to 1.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_f : scalar or vector
        - m_e : scalar or vector, optional; <2.>
        - a_f : scalar or vector, optional; <1.>
        - s_w : scalar or vector, optional; <1.>
        - n_e : scalar or vector, optional; <1.>

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_arch(data)
    array([ 1.        ,  0.79056942,  0.5       ])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_f = tdata['rho_f']
    except NameError:
        raise
    m_e = tdata.get('m_e', np.array(2.))
    a_f = tdata.get('a_f', np.array(1.))
    s_w = tdata.get('s_w', np.array(1.))
    n_e = tdata.get('n_e', np.array(1.))

    rho_b = np.array(rho_b, dtype=float, copy=True, ndmin=1)

    trho_b = np.array(rho_b, dtype=float, copy=True, ndmin=1)
    trho_b[trho_b == 0.] = np.NaN

    por_ = ((a_f*rho_f) / (trho_b * s_w**n_e))**(1./m_e)
    por_[por_ > 1.0] = 1.0

    return np.nan_to_num(por_)


def por_r_herm(data):
    """Porosity with Hermance from resistivities.

    [GRL.79.Hermance]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector
        - m_e : scalar or vector, optional; <2.>

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_herm(data)
    array([ 1.        ,  0.70710678, -0.        ])

    """
    tdata = dc(data)

    m_e = tdata.get('m_e', np.array(2.))

    return por_r_harm(tdata)**(1./m_e)


def por_r_hsbs(data):
    """Porosity with Hashin-Shtrikman bounds from resistivity.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector

    Returns
    -------
    por_l, por_u : scalar or vector
        Porosity, lower and upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> lower, upper = vel2res.por_r_hsbs(data)
    >>> print(lower)
    [ 1.          0.66666667 -0.        ]
    >>> print(upper)
    [ 1.          0.57142857 -0.        ]

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    por_l = (rho_b - rho_s)*(2*rho_f + rho_s)/((rho_f - rho_s)*
                                               (2*rho_b + rho_s))
    por_u = ((rho_b - rho_s)/(rho_f - rho_s)*3.*rho_f / (2*rho_b + rho_f))

    if rho_s < rho_f:
        por_t = dc(por_u)
        por_u = dc(por_l)
        por_l = dc(por_t)

    return por_l, por_u


def por_r_hslb(data):
    """Porosity with HS lower bound from resistivity.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector

    Returns
    -------
    por_l : scalar or vector
        Porosity, lower bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_hslb(data)
    array([ 1.        ,  0.66666667, -0.        ])

    """
    por_l, _ = por_r_hsbs(data)
    return por_l


def por_r_hsub(data):
    """Porosity with HS upper bound from resistivity.

    [MPS.63.Hashin]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector

    Returns
    -------
    por_l : scalar or vector
        Porosity, upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_hsub(data)
    array([ 1.        ,  0.57142857, -0.        ])

    """
    _, por_u = por_r_hsbs(data)
    return por_u


def por_r_hsub2(data):
    """Porosity with HS lower bound from resistivity, Berryman.

    Porosity with HS lower bound, Eq. 13, Berryman (1995). This is the same as
    por_hsbs(rho_b, rho_s, rho_f) if rho_b == rho_f.  NOT SURE ABOUT THIS, JUST
    TO REPRODUCE CARCIONE. [B.AGU.95.Berryman], [GEO.07.Carcione]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_f : scalar or vector

    Returns
    -------
    por_l : scalar or vector
        Porosity, upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_hsub2(data)
    array([ 1.        ,  0.71428571,  0.33333333])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    return 3*rho_f/(rho_f + 2*rho_b)


def por_r_self(data):
    """Porosity with self-similar model from resistivities.

    [GEO.1983.Bussian]_, [GEO.81.Sen]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector
        - m_e : scalar or vector, optional; <2.>

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_self(data)
    array([ 1.        ,  0.63245553, -0.        ])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
    except NameError:
        raise
    m_e = tdata.get('m_e', np.array(2.))

    return ((rho_b - rho_s)/(rho_f - rho_s)) * (rho_f/rho_b)**(1./m_e)


def por_r_crim(data):
    """Porosity with CRIM from resistivities.

    [B.PER.96.Schon]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_s, rho_f : scalar or vector
        - y_e : scalar or vector, optional; <1./2.>

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_crim(data)
    array([ 1.        ,  0.58113883, -0.        ])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_s = tdata['rho_s']
        rho_f = tdata['rho_f']
    except NameError:
        raise
    y_e = tdata.get('y_e', np.array(1./2))

    por_ = ((rho_b**y_e-rho_s**y_e)/
            (rho_f**y_e-rho_s**y_e))*(rho_f/rho_b)**y_e

    return por_


def por_r_dems(data):
    """Porosity with DEM from resistivities.

    Differential effective medium DEM, [GEO.81.Sen]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_r_dems(data)
    array([ 1.        ,  0.73100443,  0.39685026])

    """
    tdata = dc(data)

    try:
        rho_b = tdata['rho_b']
        rho_f = tdata['rho_f']
    except NameError:
        raise

    trho_b = np.array(rho_b, dtype=float, copy=True, ndmin=1)

    por_ = (rho_f / trho_b )**(2./3)
    por_[por_ > 1.0] = 1.0

    return por_

#  2.B Porosity from velocity                                       2B

def por_v_arit(data):
    """Porosity with arithmetic mean from P-wave velocities.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, vp_s, vp_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_arit(data)
    array([ 1. ,  0.5, -0. ])

    """
    tdata = dc(data)

    try:
        vp_b = tdata['vp_b']
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
    except NameError:
        raise

    return (vp_b - vp_s)/(vp_f - vp_s)


def por_v_geom(data):
    """Porosity with geometric mean from P-wave velocities.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, vp_s, vp_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_geom(data)
    array([ 1.      ,  0.382017, -0.      ])

    """
    tdata = dc(data)

    try:
        vp_b = tdata['vp_b']
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
    except NameError:
        raise

    return (np.log(vp_b) - np.log(vp_s))/(np.log(vp_f) - np.log(vp_s))


def por_v_harm(data):
    """Porosity with harmonic mean from P-wave velocities.

    Also called time-average equation, or Wyllie equation.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, vp_s, vp_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_harm(data)
    array([ 1.        ,  0.27272727, -0.        ])

    """
    tdata = dc(data)

    try:
        vp_b = tdata['vp_b']
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
    except NameError:
        raise

    return (vp_f/vp_b)*(vp_b - vp_s)/(vp_f - vp_s)


def por_v_raym(data):
    """Porosity with Raymer from P-wave velocities.

    [SPWLA.80.Raymer]_

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, vp_s, vp_f : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_raym(data)
    array([ 0.625     ,  0.22287618,  0.        ])

    """
    tdata = dc(data)

    try:
        vp_b = tdata['vp_b']
        vp_s = tdata['vp_s']
        vp_f = tdata['vp_f']
    except NameError:
        raise

    tvp_b = np.array(vp_b, dtype=float, copy=True, ndmin=1)

    tvp_b[tvp_b < vp_f] = vp_f
    tvp_b[tvp_b > vp_s] = vp_s
    por_ = (-np.sqrt(4*vp_s*(tvp_b-vp_f)+vp_f**2)+2*vp_s-vp_f)/(2*vp_s)
    por_[por_ < 0.] = 0.

    return por_


def por_v_hsbs(data):
    """Porosity with HS bounds from P-wave velocities.

    Porosity with Hashin-Shtrikman bounds from P-wave velocities, using a
    root-finding algorithm, [MPS.63.Hashin]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, k_s, k_f, den_f, den_s : scalar or vector
        - mu_s : scalar or vector, optional; <0.>

    Returns
    -------
    por_l, por_u : scalar or vector
        Porosity, lower and upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> lower, upper = vel2res.por_v_hsbs(data)
    >>> print(lower)
    [ 1.          0.00797824  0.        ]
    >>> print(upper)
    [ 1.          0.50542094  0.        ]

    """
    tdata = dc(data)

    try:
        vp_b  = tdata['vp_b']
        k_s   = tdata['k_s']
        k_f   = tdata['k_f']
        den_f = tdata['den_f']
        den_s = tdata['den_s']
    except NameError:
        raise
    mu_s = dc(tdata.get('mu_s', np.array(0.)))

    tvp_b = np.array(vp_b, dtype=float, copy=True, ndmin=1)

    vp_s = vp_modu({'den_b':den_s, 'k_b':k_s, 'mu_b':mu_s})
    vp_f = vp_modu({'den_b':den_f, 'k_b':k_f, 'mu_b':0})
    tvp_b[tvp_b < vp_f] = vp_f
    tvp_b[tvp_b > vp_s] = vp_s

    por_u = np.array([brentq(lambda x: vp_hsub({'vp_b':tvp_b,
                     'k_s':k_s, 'k_f':k_f, 'mu_s':mu_s, 'den_f':den_f,
                     'den_s':den_s, 'por':x}) - tvp_b[x_i], 0., 1.)
                     for x_i in range(len(tvp_b))])

    vpm = max(vp_hslb({'vp_b':tvp_b, 'k_s':k_s, 'k_f':k_f, 'mu_s':mu_s,
              'den_f':den_f, 'den_s':den_s, 'por':np.linspace(0.,.1,10)}))
    vpm = vp_hslb({'vp_b':tvp_b, 'k_s':k_s, 'k_f':k_f, 'mu_s':mu_s,
              'den_f':den_f, 'den_s':den_s, 'por':0.})
    tvp_b[tvp_b > vpm] = vpm

    por_l = np.array([brentq(lambda x: vp_hslb({'vp_b':tvp_b,
                     'k_s':k_s, 'k_f':k_f, 'mu_s':mu_s, 'den_f':den_f,
                     'den_s':den_s, 'por':x}) - tvp_b[x_i],
                      0., 1.) for x_i in range(len(tvp_b))])

    return por_l, por_u


def por_v_hslb(data):
    """Porosity with HS lower bounds from P-wave velocities.

    Porosity with Hashin-Shtrikman lower bounds from P-wave velocities, using a
    root-finding algorithm, [MPS.63.Hashin]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, k_s, k_f, den_f, den_s : scalar or vector
        - mu_s : scalar or vector, optional; <0.>

    Returns
    -------
    por_l : scalar or vector
        Porosity, lower bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_hslb(data)
    array([ 1.        ,  0.00797824,  0.        ])

    """
    por_l, _ = por_v_hsbs(data)
    return por_l


def por_v_hsub(data):
    """Porosity with HS upper bounds from P-wave velocities.

    Porosity with Hashin-Shtrikman upper bounds from P-wave velocities, using a
    root-finding algorithm, [MPS.63.Hashin]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, k_s, k_f, den_f, den_s : scalar or vector
        - mu_s : scalar or vector, optional; <0.>

    Returns
    -------
    por_u : scalar or vector
        Porosity, upper bound.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_hsub(data)
    array([ 1.        ,  0.50542094,  0.        ])

    """
    _, por_u = por_v_hsbs(data)
    return por_u


def por_v_gass(data):
    """Porosity with Gassmann from P-wave velocities.

    Porosity with Gassmann from P-wave velocities, using a root-finding
    algorithm, [NGZ.51.Gassmann]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, k_s, k_f, mu_s, den_f, den_s: scalar or vector
        - a_k : scalar or vector, optional; <3.>

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_gass(data)
    array([ 1.        ,  0.25158234,  0.        ])

    """
    tdata = dc(data)

    try:
        vp_b = np.array(tdata['vp_b'], dtype=float, copy=True, ndmin=1)
        k_s   = a_if_b_scal(tdata['k_s'], vp_b)
        k_f   = a_if_b_scal(tdata['k_f'], vp_b)
        mu_s  = a_if_b_scal(tdata['mu_s'], vp_b)
        den_f  = a_if_b_scal(tdata['den_f'], vp_b)
        den_s  = a_if_b_scal(tdata['den_s'], vp_b)
    except NameError:
        raise
    a_k = a_if_b_scal(tdata.get('a_k', np.array(3.)), vp_b)

    vp_s = vp_modu({'den_b':den_s, 'k_b':k_s, 'mu_b':mu_s})
    vp_f = vp_modu({'den_b':den_f, 'k_b':k_f, 'mu_b':0})
    vp_b[vp_b < vp_f] = vp_f[vp_b < vp_f]
    vp_b[vp_b > vp_s] = vp_s[vp_b > vp_s]

    por = np.array([brentq(lambda x: vp_gass({'vp_b':vp_b,
                   'k_s':k_s[x_i], 'k_f':k_f[x_i], 'mu_s':mu_s[x_i],
                   'den_f':den_f[x_i], 'den_s':den_s[x_i], 'a_k':a_k[x_i],
                   'por':x}) - vp_b[x_i],
                    0., 1.) for x_i in range(len(vp_b))])

    return por


def por_v_aff(data):
    """Porosity with AFF from P-wave velocities.

    Acoustic formation factor, [TLA.88.raiga-clemenceau]_, by default with
    m_e, v_s from [AAPGB.92.Issler]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b : scalar or vector
        - vp_s : scalar or vector, optional; <1000./220>
        - m_e : scalar or vector, optional; <2.19>

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> vel2res.por_v_aff(data)
    array([ 0.36101049,  0.15725672,  0.        ])

    """
    tdata = dc(data)

    try:
        vp_b  = tdata['vp_b']
    except NameError:
        raise
    vp_s = tdata.get('vp_s', np.array(1000./220))
    if 'x_e' in tdata: # Backwards compatibility
        m_e  = tdata.get('x_e', np.array(2.19))
    else:
        m_e  = tdata.get('m_e', np.array(2.19))

    return 1. - (vp_b/vp_s)**(1./m_e)

#  3. Cross-relations via porosity                                   3

def in2por2out(data, in2por=0, por2out=0):
    """Cross-property vel->por->res or res->por->vel.

    `in2por` and `por2out` have to be provided either as parameter, or within
    data-dictionary.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - all the necessary parameters for transforms `in2por` and `por2out`
        - in2por : transform that yields porosity, optional
        - por2out : transform that uses porosity, optional

    in2por : transform that yields porosity, optional; <0>
        If <in2por> = -1, it assumes that <por2out> is a direct transform
        vel -> res.
    por2out : transform that uses porosity, optional; <0>


    Returns
    -------
    out : scalar or vector
        Result from in2por -> por2out.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> # Method one, transforms provided as parameters
    >>> in2por2out(data, por_v_raym, rho_crim)
    array([  3.78698225,   6.68705707,  10.        ])
    >>> # Method one, transforms provided within dictionary
    >>> data['in2por']  = por_v_raym
    >>> data['por2out'] = rho_crim
    >>> in2por2out(data)
    array([  3.78698225,   6.68705707,  10.        ])

    """
    tdata = dc(data)

    if in2por == 0:
        try:
            in2por = tdata['in2por']
        except NameError:
            raise
    if por2out == 0:
        try:
            por2out = tdata['por2out']
        except NameError:
            raise

    if in2por != -1:
        tdata['por'] = np.array(in2por(tdata))

    return por2out(tdata)

#   4. Cross-corelations directly                                    4

def rho_faus(data):
    """Resistivity with Faust.

    I had that wrong for all of the first year or so.  It is wrong in
    [GEO.07.Carcione]_, but corrected in [B.CUP.09.Mavko]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - vp_b, rho_f, depth : scalar or vector
            Velocity, resistivity, depth.
        - flag : int, optional; <0.>
            If flag = 1 : limit rho to rho_f, rho_s;
            if flag = 2: write NaN if rho exceeds rho_s or rho_f

    Returns
    -------
    rho : scalar or vector
        Resistivity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> data['depth'] = np.array([2.])
    >>> vel2res.rho_faus(data)
    array([  0.09903997,   3.76061702,  35.61417531])

    """
    tdata = dc(data)

    try:
        vp_b  = tdata['vp_b']
        depth  = tdata['depth']
        rho_f = tdata['rho_f']
    except NameError:
        raise
    flag = tdata.get('flag', np.array(0.))
    if flag > 0:
        try:
            rho_s = tdata['rho_s']
        except NameError:
            raise

    # Wrong equation (49) in Carcione et al (2007)
    #w rho = depth * rho_f * (2.2888/vp_b)**6.

    # Right equation in The Rock Physics Handbook, Mavko et al
    rho = rho_f * (vp_b/2.2888)**6. / depth

    if flag == 1:
        rho[rho < rho_f] = rho_f
        rho[rho > rho_s] = rho_s
    elif flag == 2:
        rho[rho < rho_f] = np.NaN
        rho[rho > rho_s] = np.NaN

    return np.nan_to_num(rho)


def vp_faus(data):
    """P-wave velocity with Faust.

    I had that wrong for all of the first year or so.  It is wrong in
    [GEO.07.Carcione]_, but corrected in [B.CUP.09.Mavko]_.

    See: Module documentation for a description of the parameters.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - rho_b, rho_f, depth : scalar or vector

    Returns
    -------
    vp : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data, nst=3)
    >>> data['depth'] = np.array([2.])
    >>> vel2res.vp_faus(data)
    array([ 2.56909114,  2.77843031,  3.236852  ])

    """
    tdata = dc(data)

    try:
        rho_b  = tdata['rho_b']
        rho_f = tdata['rho_f']
        depth  = tdata['depth']
    except NameError:
        raise

    #w return 2.2888*(depth * rho_f / rho_b)**(1./6)
    return 2.2888*(depth * rho_b / rho_f)**(1./6)

#   5. Relations for other parameters                                5

def rhof_cec(depth, t_0=4., grad=3./100):
    """Temperature dependent rho_f values.

    Temperature dependent rho_f values following CEC / dual water layer
    theory by [JPT.84.Clavier]_.

    Default to temperature gradient standard values typical for North Sea:
    4 + 0.03 * depth  [t_0 + grad * depth]

    Parameters
    ----------
    depth : scalar or vector
        Depth in meter.
    t_0 : scalar, optional; <4.>
        Temperature at depth 0 km, in degree Celsius.
    grad : scalar, optional; <0.03>
        Temperature gradient, in degree Celsius per m.

    Returns
    -------
    rho_f : scalar or vector
        Resistivity of the pore fluid.

    Examples
    --------
    >>> import vel2res
    >>> vel2res.rhof_cec(np.linspace(1000,4000,4))
    array([ 0.09927209,  0.04978154,  0.0348156 ,  0.02779453])

    """

    depth = np.array(depth, dtype=float, copy=True, ndmin=1)
    temp = t_0 + grad*depth
    rho_f = 1./(6.8*(1 + 0.0545*(temp-25) - 0.0001127*(temp-25)**2))

    return rho_f


def rhof_mol(depth, t_0=4., grad=3./100, mol=0.629711):
    """Temperature and molatility dependent rho_f values.

    Temperature and molatility dependent rho_f values after [GEO.92.Sen]_.

    Default to temperature gradient standard values typical for North Sea:
    4 + 0.03 * depth  [t_0 + grad * depth]

    Parameters
    ----------
    depth : scalar or vector
        Depth in meter.
    t_0 : scalar, optional; <4.>
        Temperature at depth 0 km, in degree Celsius.
    grad : scalar, optional; <0.03>
        Temperature gradient, in degree Celsius per m.
    mol : scalar, optional; <0.629711>
        Molality (default molatility: 0.629711 [35000 ppm NaCl])

    Returns
    -------
    rho_f : scalar or vector
        Resistivity of the pore fluid.

    Examples
    --------
    >>> import vel2res
    >>> vel2res.rhof_mol(np.linspace(1000,4000,4))
    array([ 0.14978062,  0.09812083,  0.07387391,  0.0598386 ])

    """

    depth = np.array(depth, dtype=float, copy=True, ndmin=1)
    temp = t_0 + grad*depth
    rho_f = 1./((5.6 + .27*temp - .00015*temp**2)*mol -
             ((2.36 + .099*temp)/(1. + .214*mol))*mol**(3./2))

    return rho_f


def m_e_folke(data):
    """Porosity corrected cementation-exponent after

    [SEG.10.Engelmark]_

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - m_e : scalar or vector
            Cementation exponent.
        - por : scalar or vector, optional; <data['in2por'](data)>
            If porosity is not provided, it is calculated with data['in2por']
            of the data. Data needs then 'in2por' and its parameters.

    Returns
    -------
    m_e : scalar or vector
        Porosity corrected cementation exponent.

    Examples
    --------
    >>> import vel2res
    >>> data = {}
    >>> data['por'] = np.array([0.3, 0.4])
    >>> data['m_e'] = 3.
    >>> vel2res.m_e_folke(data)
    array([ 2.7,  2.6])

    """
    tdata = dc(data)

    try:
        m_e  = tdata['m_e']
    except NameError:
        raise
    try:
        por  = tdata['por']
    except:
        try:
            por = tdata['in2por'](tdata)
        except NameError:
            raise

    return m_e - por


def param_depth(param, depth, grad, min_d):
    """Depth dependent parameters.

    Parameters
    ----------
    param : scalar or vector
        Parameter to apply depth trend (any unit).
    depth : scalar or vector
        Depths (any unit).
    grad : scalar
        Gradient. Units according to param and grad.
    min_d : scalar
        From which depth onwards the gradient applies. Unit according to depth.

    Returns
    -------
    newparam : scalar or vector
        Corrected parameter.

    Examples
    --------
    >>> import vel2res
    >>> vel2res.param_depth(np.linspace(1,4,4), np.linspace(1,4,4), 3., 2)
    array([  1.,   2.,   6.,  10.])

    """

    tdepth = np.array(depth, dtype=float, copy=True, ndmin=1)
    newparam = np.zeros(np.size(tdepth))
    newparam[tdepth < min_d] = param[tdepth < min_d]
    newparam[tdepth >= min_d] = param[tdepth >= min_d] + grad*(tdepth[tdepth
            >= min_d]-min_d)

    return newparam


#   6. Other necessary relations                                     6

def den_bulk(data):
    """Bulk density from fluid and grain densities.

    Arithmetic mean.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - den_s, den_f, por : scalar or vector

    Returns
    -------
    den_b : scalar or vector
        Bulk density.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.den_bulk(data)
    array([ 2.125,  2.05 ])

    """
    tdata = dc(data)

    try:
        por   = tdata['por']
        den_s = tdata['den_s']
        den_f = tdata['den_f']
    except NameError:
        raise

    return (1. - por) * den_s + por * den_f


def por_dens(data):
    """Porosity from fluid and grain densities.

    Arithmetic mean.

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - den_s, den_f, den_b : scalar or vector

    Returns
    -------
    por : scalar or vector
        Porosity.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['den_b'] = np.array([ 2.125,  2.05 ])
    >>> vel2res.por_dens(data)
    array([ 0.25,  0.3 ])

    """
    tdata = dc(data)

    try:
        den_b = tdata['den_b']
        den_s = tdata['den_s']
        den_f = tdata['den_f']
    except NameError:
        raise

    return (den_b - den_s)/(den_f - den_s)


def vp_modu(data=None, den_b=None, k_b=None, mu_b=None):
    """P-wave velocity from bulk and shear moduli.

    :Note: Example of how I could implement the use of input either a
        dictionary or the variables itself.

    Parameters
    ----------
    data : dict, optional; <None>
        If None, then the values have to be provided on their own.
        Containing the following entries:

        - den_b, k_b, mu_b : scalar or vector

    Returns
    -------
    v_p : scalar or vector
        P-wave velocity.

    Examples
    --------
    >>> import vel2res
    >>> vel2res.vp_modu({'den_b': 2.5, 'k_b': 20., 'mu_b': 15. })
    4.0
    >>> vel2res.vp_modu(den_b= 2.5, k_b= 20., mu_b= 15. )
    4.0

    """
    if data != None:
        tdata = dc(data)


        try:
            den_b = tdata['den_b']
            k_b   = tdata['k_b']
            mu_b  = tdata['mu_b']
        except NameError:
            raise

    return np.sqrt((k_b + 4.*mu_b/3.)/den_b)


def k_mu_hsbs(data):
    """Bulk and shear moduli with HS bounds.

    [MPS.63.Hashin]_, [GEO.07.Carcione]_

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - k_s, k_f, por : scalar
        - mu_s : scalar, optional; <0.>
            default mu_s = 0 => lower bound

    Returns
    -------
    k_m, mu_m : scalar or vector
        Bulk and shear moduli.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.k_mu_hsbs(data)
    (array([ 12.40133779,  11.06824615]), array([ 9.        ,  8.07692308]))

    """
    tdata = dc(data)

    try:
        k_s   = tdata['k_s']
        k_f   = tdata['k_f']
        por   = tdata['por']
    except NameError:
        raise
    mu_s = tdata.get('mu_s', np.array(0.))

    x_1 = mu_s/6.*((9.*k_s + 8.*mu_s)/(k_s + 2.*mu_s))

    if mu_s == 0.0:
        mu_m = 0.0
    else:
        x_1 = mu_s/6.*((9.*k_s + 8.*mu_s)/(k_s + 2.*mu_s))
        mu_m = ((1.-por)/(mu_s + x_1) + por/x_1)**-1 - x_1

    x_2 = 4.*mu_m/3
    k_m = (por/(k_f + x_2) + (1. - por)/(k_s + x_2))**-1 - x_2
    return k_m, mu_m


def k_gass(data):
    """Bulk modulus with Gassmann.

    [NGZ.51.Gassmann]_, [GEO.07.Carcione]_

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - k_s, k_f, k_m, por : scalar

    Returns
    -------
    k_g : scalar or vector
        Bulk modulus.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> data['k_m'] = k_mu_krie(data)
    >>> vel2res.k_gass(data)
    array([ 9.84719725,  8.23120603])

    """
    tdata = dc(data)

    try:
        k_s   = tdata['k_s']
        k_m   = tdata['k_m']
        k_f   = tdata['k_f']
        por   = tdata['por']
    except NameError:
        raise

    tpor = np.array(por, dtype=float, copy=True, ndmin=1)
    tk_m = np.array(k_m, dtype=float, copy=True, ndmin=1)

    k_g = np.zeros(tpor.shape)
    k_g[tpor == 0. ] = k_s
    k_g[tpor == 1. ] = k_f
    b_i = np.logical_and(tpor != 1., tpor != 0.)
    k_g[b_i] = ((k_s - tk_m[b_i] + tpor[b_i]*tk_m[b_i]*(k_s/k_f - 1.)) /
               (1. - tpor[b_i] - tk_m[b_i]/k_s + tpor[b_i] * k_s / k_f))

    return k_g


def k_mu_krie(data):
    """Bulk and shear moduli with Krief.

    [TLA.90.Krief]_, [B.CUP.09.Mavko]_, [GEO.07.Carcione]_

    In Mavko et al, Krief is defined differently than in Carcione et al.  If
    a_k is given as positive number, then Mavko is applied.  If a_k is given as
    negative number, then Carcione is applied (comparing to Krief reveals that
    Mavko is the correct one).

    Parameters
    ----------
    data : dict
        Containing the following entries:

        - k_s, mu_s, por : scalar or vector
            Moduli, porosity.
        - a_k : scalar or vector, optional; <3.>
            If positive, Mavko applied, if negative, Carcione.

    Returns
    -------
    k_m, mu_m : scalar or vector
        Bulk and shear moduli.

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> data['por'] = np.array([0.25, 0.3])
    >>> vel2res.k_mu_krie(data)
    (array([ 6.328125  ,  4.33675066]), array([ 4.74609375,  3.25256299]))

    """
    tdata = dc(data)

    try:
        k_s   = tdata['k_s']
        mu_s  = tdata['mu_s']
        por   = tdata['por']
    except NameError:
        raise
    a_k = tdata.get('a_k', np.array(3.))

    tpor = np.array(por, dtype=float, copy=True, ndmin=1)

    k_m = np.zeros(tpor.shape)
    mu_m = np.zeros(tpor.shape)
    a_exp = np.zeros(tpor.shape)
    b_i = (tpor != 1.)

    if a_k >= 0:
        a_exp[b_i] = np.array(a_k/(1. - tpor[b_i]))
    else:
        a_exp[b_i] = np.array(1. + -1*a_k/(1. - tpor[b_i]))

    k_m[tpor == 1.] = 0.
    mu_m[tpor == 1.] = 0.
    k_m[tpor == 0.] = k_s
    mu_m[tpor == 0.] = mu_s
    b_i = np.logical_and(tpor != 1., tpor != 0.)
    k_m[b_i] = k_s * (1. - tpor[b_i])**a_exp[b_i]
    mu_m[b_i] = (mu_s / k_s) * k_m[b_i]

    return k_m, mu_m


#   7. Functions that where in other files

def a_if_b_scal(b, a):
    """Check if b is scalar. If yes, make same size as a.

    If b is a vector, return b. If b is a scalar, it returns a vector of
    the same shape as a, filled with this scalar.

    Parameters
    ----------
    a, b : scalar or array
        Number or vectors, where b must be either a scalar or a vector of
        the same size as a.

    Returns
    -------
    out : array, float
        Adjusted b

    Examples
    --------
    >>> import vel2res
    >>> vector_a = np.array([1,2,3])
    >>> vector_b = 1.
    >>> vector_c = np.array([1., 10., 100.])
    >>> vel2res.a_if_b_scal(vector_b, vector_a)
    array([ 1.,  1.,  1.])
    >>> vel2res.a_if_b_scal(vector_c, vector_a)
    array([   1.,   10.,  100.])

    """
    if (np.size(b) != 1) and (np.size(b) != np.size(a)):
        ermsg = "<b> must be a scalar or of same size as <a>!"
        raise ValueError(ermsg)

    fa = np.array(a, dtype=float)
    out = np.ones_like(fa)
    out[:] = b
    return out


def carc_tab1(rock):
    """Carcione et al. 2007, Geophysics, Table 1.

    Rock physics parameter from Table 1 in [GEO.07.Carcione]_.

    Parameters
    ----------
    rock : string; {'shale', 'sand'}
        Either shale or sand.

    Returns
    -------
    data : dict
        Containing:

        - rho_s : array (resistivity of solid fraction)
        - rho_f :  array (resistivity of fluid fraction)
        - k_s :  array (bulk modulus of solid fraction)
        - k_f :  array (bulk modulus of fluid fraction)
        - mu_s :  array (shear modulus of solid fraction)
        - den_s :  array (density of solid fration)
        - den_f :  array (density of fluid fration)

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')

    """

    # Brine shale parameters from Table 1.
    if rock == 'shale':
        data = dict(rho_s = np.array(10.),
                    rho_f = np.array(2.5),
                    k_s = np.array(20.),
                    mu_s = np.array(15.),
                    den_s = np.array(2.5),
                    k_f = np.array(2.25),
                    den_f = np.array(1.),
                    )

    # Oil sand parameters from Table 1.
    elif rock == 'sand':
        data = dict(rho_s = np.array(1000.),
                    rho_f = np.array(100000.),
                    k_s = np.array(39.),
                    mu_s = np.array(40.),
                    den_s = np.array(2.65),
                    k_f = np.array(0.57),
                    den_f = np.array(0.7),
                    )

    return data


def carc_tab2(rock):
    """Carcione et al. 2007, Geophysics, Table 2.

    Rock physics parameter from Table 2 in [GEO.07.Carcione]_.

    Parameters
    ----------
    rock : string; {'shale', 'sand'}
        Either shale or sand.

    Returns
    -------
    data : dict
        Containing:

        - rho_s : array (resistivity of solid fraction)
        - rho_f :  array (resistivity of fluid fraction)
        - k_s :  array (bulk modulus of solid fraction)
        - k_f :  array (bulk modulus of fluid fraction)
        - mu_s :  array (shear modulus of solid fraction)
        - den_s :  array (density of solid fration)
        - den_f :  array (density of fluid fration)

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')

    """

    # Brine shale parameters from Table 2.
    if rock == 'shale':
        data = dict(rho_s = np.array(5.),
                    rho_f = np.array(1./15.),
                    k_s = np.array(25.),
                    mu_s = np.array(20.),
                    den_s = np.array(2.65),
                    k_f = np.array(2.25),
                    den_f = np.array(1.03),
                    )

    # Brine sand parameters from Table 2,
    elif rock == 'sand':
        data = dict(rho_s = np.array(100.),
                    rho_f = np.array(1./15.),
                    k_s = np.array(39.),
                    mu_s = np.array(40.),
                    den_s = np.array(2.65),
                    k_f = np.array(2.25),
                    den_f = np.array(1.03),
                    )

    return data


def carc_der(data, nst=101):
    """Derived values from Carcione et al. 2007.

    Rock physics parameter from Table 2 in [GEO.07.Carcione]_.
    Spacing in conductivites is better for linear plots!

    Parameters
    ----------
    data : dict
        As returned from `carc_tab1` or `carc_tab2`.
    nst : int, optional; <101>
        Number of samples.

    Returns
    -------
    data : dict
        Does not return, but puts in place:

        - rho_b : array (bulk resistivity)
        - rho_0 : array (bulk resistivity, with rho_0[0] = inf)
        - vp_s : array (P-wave velocity of solid fraction)
        - vp_f : array (P-wave velocity of fluid fraction)
        - vp_b : array (bulk P-wave velocity)

    Examples
    --------
    >>> import vel2res
    >>> data = vel2res.carc_tab1('shale')
    >>> vel2res.carc_der(data)

    """

    if data['rho_s'] > data['rho_f']:
        data['rho_b'] = 1./np.linspace(1./data['rho_f'], 1./data['rho_s'], nst)
    elif data['rho_s'] == data['rho_f']:
        data['rho_b'] = data['rho_s']
    else:
        data['rho_b'] = 1./np.linspace(1./data['rho_s'], 1./data['rho_f'], nst)

    data['rho_0'] = np.r_[np.inf, 1./np.linspace(0, 1./data['rho_f'], nst)[1:]]
    data['vp_s'] = vp_modu({'den_b': data['den_s'], 'k_b': data['k_s'],
                            'mu_b': data['mu_s']})
    data['vp_f'] = vp_modu({'den_b': data['den_f'], 'k_b': data['k_f'],
                            'mu_b': 0.})
    data['vp_b'] = np.linspace(data['vp_f'], data['vp_s'], nst)
