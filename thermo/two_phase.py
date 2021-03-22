# File to contain functions for calculating various two-phase parameters

import numpy as np

def mean_viscosity_McAdams(mu_l, mu_g, x):
    """Calculate mean viscosity in two-phase flow based on McAdams formula (10.42 in Carey2008)

    Args:
        mu_l (Pa*s): Liquid viscosity at saturation point
        mu_g (Pa*s): Gas viscosity at saturation point
        x (-): Vapour quality

    Returns:
        [type]: [description]
    """
    return 1 / (x / mu_g + (1-x) / mu_l )

def homogenous_void_fraction(x, rho_g, rho_l):
    """ Calculate void fraction based on homogenous assumption

    Args:
        x (-): Vapour quality
        rho_g (kg/m^3): Gas density at saturation
        rho_l (kg/m^3): Liquid density at saturation

    Returns:
        alpha (-): Void fraction based on homogenous assumption
    """
    return rho_l * x / ( rho_g * (1-x) + rho_l * x )

def mixture_density(alpha, rho_g, rho_l):
    """Calculate mixture density of a two-phase mixture

    Args:
        alpha (-): Void fraction
        rho_g (kg/m^3): Gas density
        rho_l (kg/m^3): Liquid density

    Returns:
        rho (-): Mixture density of two-phase flow
    """
    return rho_g * alpha + rho_l * (1-alpha)

def Nu_Kandlikar_NBD_CBD_dryout(args):
    """ Calculate Nusselt number in micro-channels, which are dominated by the nucleate boiling regime according to Kandlikar (eq. 12.164 in Cary2008)

    However, an adaptation is made, since the relations for dry-out do not fully extend to micro-channels.
    This is because the dryout_quality is calculated to be well above x=1 based on 12.77 and 12.92, which is impossible, of course.
    However, eventually the modelled nucleate boiling an convective film boiling in the channel tends to zero, resulting in predicted infinite channel lengths.

    As such, this function will be adapted to simply let Nusselt number that models convection in a dried-out channel to take over when it is the largest number.
    NOTE: the used conductivity for calculating the dry-out relation might be different, this must be taken into account.

    Args:
        args (-): Flow similarity parameters, two-phase parameters, etc...

    Returns:
        Nu (-): Returns the one of three Nusselt numbers: either Nucleate boiling, convective boiling or dryout number. Depending on which is the largest.
    """
    rho_ratio = args['rho_l']/args['rho_g'] # [-] Ratio of saturation densities
    x = args['x'] # Vapour quality
    
    # Slight adaptation, the Nusselt number of the flow as entirely as liquid is inputted instead of h_le
    Nu_NBD = ( 0.6683 * rho_ratio**0.1 * x**0.16 * (1 - x)**0.64  + \
            1058 * args['Bo']**0.7 * (1 - x)**0.8 ) * args['Nu_le'] # [-] Nusselt number of nucleate boiling regime
    Nu_CBD = (  1.1360 * rho_ratio**0.45 * x**0.72 * (1 - x)**0.08 + \
                667.2 * args['Bo']**0.7 * (1 - x)** 0.8 ) * args['Nu_le'] # Nusselt number of convective boiling regime

    # Check which of three Nusselt numbers are at a maximum
    Nu = np.maximum(Nu_NBD, Nu_CBD) # [-] Find the maximum of NBD and CBD
    Nu = np.maximum(Nu, args['Nu_dryout']) # [-] Of those two, check if the dry-out function is larger
    return Nu

def Nu_Kandlikar_NBD_dryout(args):
    """ TODO: Same as previous function, but no convective boiling
    Calculate Nusselt number in micro-channels, which are dominated by the nucleate boiling regime according to Kandlikar (eq. 12.164 in Cary2008)

    However, an adaptation is made, since the relations for dry-out do not fully extend to micro-channels.
    This is because the dryout_quality is calculated to be well above x=1 based on 12.77 and 12.92, which is impossible, of course.
    However, eventually the modelled nucleate boiling an convective film boiling in the channel tends to zero, resulting in predicted infinite channel lengths.

    As such, this function will be adapted to simply let Nusselt number that models convection in a dried-out channel to take over when it is the largest number.
    NOTE: the used conductivity for calculating the dry-out relation might be different, this must be taken into account.

    Args:
        args (-): Flow similarity parameters, two-phase parameters, etc...

    Returns:
        Nu (-): Returns the one of three Nusselt numbers: either Nucleate boiling, convective boiling or dryout number. Depending on which is the largest.
    """
    rho_ratio = args['rho_l']/args['rho_g'] # [-] Ratio of saturation densities
    x = args['x'] # Vapour quality
    
    # Slight adaptation, the Nusselt number of the flow as entirely as liquid is inputted instead of h_le
    Nu_NBD = ( 0.6683 * rho_ratio**0.1 * x**0.16 * (1 - x)**0.64  + \
            1058 * args['Bo']**0.7 * (1 - x)**0.8 ) * args['Nu_le'] # [-] Nusselt number of nucleate boiling regime

    print("Nu_NDB:")
    print(Nu_NBD)
    # Check which of three Nusselt numbers are at a maximum
    Nu = np.maximum(Nu_NBD, args['Nu_dryout']) # [-] Find the maximum of NBD and CBD
    return Nu

def Nu_DB_two_phase(args):
    """ Dittus-Boelter, reimplemented for two-phase flow. The formula is the same, but is important that Nu is calculated in relation with the liquid saturation conductivity kappa_l.
    This is done to maintan consistency with how the other two-phase Nusselt numbers are calculated (in relation to kappa_l)

    Args:
        args (-): Dictionary of flow similarity paramters

    Returns:
        Nu: Nusselt number, in relation to kappa_l
    """
    return 0.023 * args['Re']**0.8 * args['Pr']**0.4 * args['kappa'] / args['kappa_l']

def dryout_quality(p, m_dot, A_channel, D_hydr):
    """Calculate dry-out quality based on Equations 12.77 & and 12.92 for water in Carey2008

    Args:
        p (Pa): Pressure
        m_dot (kg/s): Mass flow
        A_channel (m^2): Channel cross-sectional area
        D_hydr (m): Hydraulic diameter

    Returns:
        x_crit/dryout: Dry-out/critcal quality
    """

    a = p * 1e-5 / 98 # [-] Pressure coefficient ( P/ 98 ) in formula 12.92 (P is pressure in bar)
    # First calculate dry-out quality at 8mm
    x_crit_8mm = (  # [-] Dryout quality at 8mm
        0.39 \
        + 1.57 * a
        - 2.04 * a**2
        + 0.68 * a**3 ) \
        * ( m_dot/A_channel / 1e3 )**-0.5
    # Now correct for diameter with 12.77
    x_crit = x_crit_8mm * ( 8e-3 / D_hydr )**0.15
    return x_crit # [-] Critical quality (dry-out)

def mean_conductivity(kappa_g, kappa_l, rho_l, rho_g, x):
    """Return the mean conductivity based on 12.129 of Carey2008

    Args:
        kappa_g (W/(m*K)): Gas conductivity
        kappa_l (W/(m*K)): Liquid conductivity
        rho_l (kg/m^3): Liquid saturation conductivity
        rho_g (kg/m^3): Gas saturation conductivity
        x (-): Vapour quality

    Returns:
        kappa_tp: Mean two-phase conductivity
    """
    return ( x * kappa_g * (rho_l / rho_g) + (1-x) * kappa_l ) /\
        ( x * (rho_l / rho_g) + (1-x) )

def mean_Prandtl(Pr_g, Pr_l, rho_l, rho_g, x):
    """Return the mean Prandtly based on 12.129 of Carey2008

    Args:
        Pr_g (-): Gas saturation Prandtl number
        Pr_l (-): Liquid saturation Prandtl number
        rho_l (kg/m^3): Liquid saturation density
        rho_g (kg/m^3): Gas saturation density
        x (-): Vapour quality

    Returns:
        Pr_tp: Mean two-phase conductivity
    """
    return ( x * Pr_g * (rho_l / rho_g) + (1-x) * Pr_l ) /\
        ( x * (rho_l / rho_g) + (1-x) )

def mean_viscosity(mu_g, mu_l, rho_l, rho_g, x):
    """Return the mean viscosity based on 12.129 of Carey2008

    Args:
        mu_g (-): Gas saturation viscosity
        mu_l (-): Liquid saturation viscosity
        rho_l (kg/m^3): Liquid saturation density
        rho_g (kg/m^3): Gas saturation density
        x (-): Vapour quality

    Returns:
        mu_tp: Mean two-phase conductivity
    """
    return ( x * mu_g * (rho_l / rho_g) + (1-x) * mu_l ) /\
        ( x * (rho_l / rho_g) + (1-x) )