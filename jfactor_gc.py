import numpy as np
from scipy.integrate import quad, dblquad

KPC_TO_CM = 3.086e21

class DensityProfile(object):
    ''' this class allows to define and correctly normalise the density profiles '''
    
    def __init__(self, name, functional_form, r_s, rho_sun, r_sun):
        ''' The normalisation is computed by rho_sun, r_sun.
            functional_form is the functional_form of the profile:
                - it does not contain the normalisation density rho_s
                - it does not contain r_s, but y := r/r_s
                - it depends solely on y.
        '''
        self.name = name
        self.r_s = r_s
        self._functional_form = functional_form
        self.set_normalization(rho_sun = rho_sun, r_sun = r_sun)

    def functional_form(self):
        return self._functional_form

    def full_form(self):
        return lambda r: self.rho_s * self._functional_form(r / self.r_s)

    def set_normalization(self, rho_sun, r_sun):
        self.r_sun = r_sun
        self.rho_s = rho_sun / self._functional_form(r_sun / self.r_s)

    def get_name(self):
        return self.name

    def __call__(self, r):
        return self.full_form().__call__(r)

    def is_optimized(self, other):
        ''' check if the profile is optimised for a certain ROI: same name and same parameters (related to the profile) '''
        if not isinstance(other, DensityProfile):
            raise TypeError("'is_optimized' not supported between instances of 'DensityProfile' and '%s'" % other.__class__.__name__)
        return self.name == other.name

    def __eq__(self, other):
        ''' check full equality between profile parameters '''
        if not isinstance(other, DensityProfile):
            raise TypeError("'==' not supported between instances of 'DensityProfile' and '%s'" % other.__class__.__name__)
        return self.eq_but_norm(other) and np.isclose(self.rho_s, other.rho_s, atol = 0., rtol = 1e-4)

    def eq_but_norm(self, other):
        ''' check if the profiles parameters, but rho_s (the normalization) are equal.
            If rho_s is equal it returns True as well.
            In case of True, we can simply rescale the J-factor on the basis of the new normalization (even if they are full equal)
        '''
        if not isinstance(other, DensityProfile):
            raise TypeError("'eq_but_norm' not supported between instances of 'DensityProfile' and '%s'" % other.__class__.__name__)
        return self.name == other.name and np.logical_and.reduce(np.isclose([self.r_s, self.r_sun], [other.r_s, other.r_sun], atol = 0., rtol = 1e-4))

    def __hash__(self):
        return hash((self.name, self.r_s, self.rho_s, self.r_sun))

    def __str__(self):
        return "%s(rho_s = %.4e GeV cm^-3, r_s = %.2e kpc)" % (self.name, self.rho_s, self.r_s)

    def get_parameters_items(self):
        return [("profile_rho_s", self.rho_s), ("profile_r_s", self.r_s)]

class PROFILES:
    class NFW(DensityProfile):
        def __init__(self, r_s, rho_sun = 0.4, r_sun = 8.5, **kwargs):
            self.gamma = kwargs.get("gamma", 1.0)
            functional_form = lambda y: np.power(y, -self.gamma) * np.power(1 + y, self.gamma-3.)
            super(PROFILES.NFW, self).__init__(name = "NFW", functional_form = functional_form, r_s = r_s, rho_sun = rho_sun, r_sun = r_sun)

        def is_optimized(self, other):
            test = super(PROFILES.NFW, self).is_optimized(other = other)
            return self.gamma == other.gamma if test else False

        def eq_but_norm(self, other):
            test = super(PROFILES.NFW, self).eq_but_norm(other = other)
            return self.gamma == other.gamma if test else False

        def __str__(self):
            return super(PROFILES.NFW, self).__str__() + "\b, gamma = %.2e)" % self.gamma

        def __hash__(self):
            return hash((self.name, self.r_s, self.r_sun, self.rho_s, self.gamma))

        def get_parameters_items(self):
            return super(PROFILES.NFW, self).get_parameters_items() + [("profile_gamma", self.gamma)]

    class Einasto(DensityProfile):
        def __init__(self, r_s, rho_sun = 0.4, r_sun = 8.5, **kwargs):
            self.alpha = kwargs.get("alpha", 0.17)
            functional_form = lambda y: np.exp(-2/self.alpha * (np.power(y, self.alpha) - 1))
            super(PROFILES.Einasto, self).__init__(name = "Einasto", functional_form = functional_form, r_s = r_s, rho_sun = rho_sun, r_sun = r_sun)

        def is_optimized(self, other):
            test = super(PROFILES.Einasto, self).is_optimized(other = other)
            return self.alpha == other.alpha if test else False

        def eq_but_norm(self, other):
            test = super(PROFILES.Einasto, self).eq_but_norm(other = other)
            return self.alpha == other.alpha if test else False

        def __str__(self):
            return super(PROFILES.Einasto, self).__str__() + "\b, alpha = %.2e)" % self.alpha

        def __hash__(self):
            return hash((self.name, self.r_s, self.r_sun, self.rho_s, self.alpha))

        def get_parameters_items(self):
            return super(PROFILES.Einasto, self).get_parameters_items() + [("profile_alpha", self.alpha)]

    class Burkert(DensityProfile):
        def __init__(self, r_s, rho_sun = 0.4, r_sun = 8.5, **kwargs):
            functional_form = lambda y: np.power( (1 + y) * (1 + np.power(y, 2)), -1)
            super(PROFILES.Burkert, self).__init__(name = "Burkert", functional_form = functional_form, r_s = r_s, rho_sun = rho_sun, r_sun = r_sun)

    class Isothermal(DensityProfile):
        def __init__(self, r_s, rho_sun = 0.4, r_sun = 8.5, **kwargs):
            functional_form = lambda y: np.power( 1 + np.power(y, 2), -1)
            super(PROFILES.Isothermal, self).__init__(name = "Isothermal", functional_form = functional_form, r_s = r_s, rho_sun = rho_sun, r_sun = r_sun)


def _jfactor_cone(ang_1, ang_2, profile_func, x_sun, integration_bound):
    ''' it computes the conic part of the J-factor:
            - ang_1: is the angle between the main axis of the cone and the slant height of the masked conic region (in radians)
            - ang_2: is the angle between the main axis of the cone and its slant height (in radians)
            - x_sun := r_sun/r_s
    '''
    if ang_1 >= ang_2:
        return 0.
    def coefft(t):
        ''' useful redefinition, t := cos(theta), where theta is the integration variable '''
        return np.power(x_sun, 2)*(1 - np.power(t,2))
    def func(y, t):
        ''' integrand function for the integral without singularity '''
        return np.power(profile_func(y),2) * y * np.power(np.power(y,2) - coefft(t), -0.5)
    def func_sing(t):
        ''' integrand function for the singular integral over y '''
        return lambda y: np.power(profile_func(y),2) * y * np.power(y + np.sqrt(coefft(t)), -0.5)
    def inte_func(t):
        ''' integration over r for the singular region. It will then be integrated over t '''
        return quad(func_sing(t), np.sqrt(coefft(t)), x_sun, weight = 'alg', wvar = (-0.5,0), epsrel = 1e-6, epsabs = 0)[0]
    return 2*np.pi*dblquad(func, np.cos(ang_2), np.cos(ang_1), lambda t: x_sun, lambda t: integration_bound, epsrel = 1e-6, epsabs = 0)[0] + 2*2*np.pi*quad(inte_func, np.cos(ang_2), np.cos(ang_1), epsrel = 1e-6, epsabs = 0)[0]

def _jfactor_mask(ang_1, ang_2, ang_lat, ang_long, profile_func, x_sun, integration_bound):
    ''' it computes the contribution to the J-factor given by a particular shape of the integration region
        sometimes used as a masked region:
            - ang_1: is the amplitude of the inner conic region (in radians)
            - ang_2: is the amplitude of the main conic region (in radians)
            - ang_lat: latitude angle covered by this region: |latitude| < ang_lat (in radians)
            - ang_long: minimum longitude covered by this region: ang_long < |longitude| < pi/2 (in radians)
            - x_sun := r_sun/r_s
    '''
    if ang_long >= ang_2 or ang_lat == 0: # the mask is zero both if it is outside the cone integration region or if the latitude has zero amplitude
        return 0.
    # useful definitions
    def coeffbl(b, l):
        ''' useful redefinition '''
        return np.power(x_sun, 2)*(1 - np.power(np.cos(b) * np.cos(l),2))
    def func_normal(b, l):
        ''' integrand function for the integral without singularity '''
        return lambda y: np.power(profile_func(y),2) * y * np.power(np.power(y,2) - coeffbl(b, l), -0.5)
    def func_sing(b, l):
        ''' integrand function for the singular integral over y '''
        return lambda y: np.power(profile_func(y),2) * y * np.power(y + np.sqrt(coeffbl(b, l)), -0.5)
    def func(l, b):
        ''' integrand function already integrated over r, b is the longitude and l is the latitude '''
        return np.cos(b) * (2.*quad(func_sing(b, l), np.sqrt(coeffbl(b, l)), x_sun, weight = 'alg', wvar = (-0.5,0), epsrel = 1e-6, epsabs = 0)[0] + quad(func_normal(b, l), x_sun, integration_bound, epsrel = 1e-6, epsabs = 0)[0])
    def b_prime(ang):
        if np.isclose(ang, np.pi/2., atol = 0., rtol = 1e-4):
            return np.pi/2.
        return np.arctan(np.sqrt(np.power(np.sin(ang),2) - np.power(np.sin(ang_long),2)/np.cos(ang)))
    def l_prime(ang, b):
        if np.isclose(ang, np.pi/2., atol = 0., rtol = 1e-4):
            return np.pi/2.
        return np.arcsin(np.cos(ang)*np.sqrt(np.power(np.tan(ang),2) - np.power(np.tan(b),2)))
    B_2 = np.amin([ang_lat, b_prime(ang_2)])
    B_3 = np.amin([ang_lat, ang_1, b_prime(ang_2)])
    if ang_long > ang_1: # condition 1
        if ang_lat > ang_1 and b_prime(ang_2) > ang_1: # condition 1 + 4
            return 4*dblquad(func, 0., B_2, lambda b: ang_long, lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0]
        else: # condition 1 only
            return 4*dblquad(func, 0., B_3, lambda b: ang_long, lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0]
    else: # condition 2
        B_1 = np.amin([ang_lat, b_prime(ang_1)])
        if ang_lat > ang_1 and b_prime(ang_2) > ang_1: # condition 2 + 4 (includes condition 3 as well)
            return 4*(dblquad(func, 0, B_1, lambda b: l_prime(ang_1, b), lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0] + dblquad(func, b_prime(ang_1), B_2, lambda b: ang_long, lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0])
        elif ang_lat > b_prime(ang_1) and not np.isclose(b_prime(ang_1), ang_lat, atol = 0., rtol = 1e-4): # condition 2 + 3, also use isclose, because in the case equality we prevent the computation of one integral (with a little approximation).
            return 4*(dblquad(func, 0, B_1, lambda b: l_prime(ang_1, b), lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0] + dblquad(func, b_prime(ang_1), B_3, lambda b: ang_long, lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0])
        else: # condition 2 only
            return 4*dblquad(func, 0, B_1, lambda b: l_prime(ang_1, b), lambda b: l_prime(ang_2, b), epsrel = 1e-6, epsabs = 0)[0]

def _jfactor_simplified_region(ang_2, ang_lat, profile_func, x_sun, ang_1, integration_bound):
    ''' it computes the J-factor for some particular values of the parameters, which make the computation more stable and faster:
            - ang_1: is the amplitude of the inner masked conic region (in radians)
            - ang_2: is the amplitude of the main conic region (in radians)
            - ang_lat: latitude angle covered by the mask: the mask covers: |latitude| < ang_lat (in radians)
            - x_sun := r_sun/r_s
    '''
    if ang_1 >= ang_2 or ang_lat >= ang_2:
        return 0.
    def coefft(t):
        ''' useful redefinition, t := cos(theta), where theta is the integration variable '''
        return np.power(x_sun, 2)*(1 - np.power(t,2))
    def phi_prime(t):
        ''' useful redefinition for integration over phi. '''
        return np.arctan(np.sqrt( 1/np.power(np.tan(ang_lat), 2) - (1 + 1/np.power(np.tan(ang_lat), 2)) * np.power(t, 2) ))
    def func(y, t):
        ''' integrand function for the integral without singularity '''
        return phi_prime(t) * np.power(profile_func(y),2) * y * np.power(np.power(y,2) - coefft(t), -0.5)
    def func_sing(t):
        ''' integrand function for the singular integral over y '''
        return lambda y: np.power(profile_func(y),2) * y * np.power(y + np.sqrt(coefft(t)), -0.5)
    def inte_func(t):
        ''' integration over r for the singular region. It will then be integrated over t '''
        return phi_prime(t) * quad(func_sing(t), np.sqrt(coefft(t)), x_sun, weight = 'alg', wvar = (-0.5,0), epsrel = 1e-6, epsabs = 0)[0]
    theta_prime = np.amax([ang_1, ang_lat])
    return 4 * (dblquad(func, np.cos(ang_2), np.cos(theta_prime), lambda t: x_sun, lambda t: integration_bound, epsrel = 1e-6, epsabs = 0)[0] + 2*quad(inte_func, np.cos(ang_2), np.cos(theta_prime), epsrel = 1e-6, epsabs = 0)[0])

def _jfactor_gc_with_mask(ang_2, profile, mask, ang_1=None, ang_lat=None, ang_long=None, integration_bound=np.inf):
    ''' if mask = False, then the mask is set to zero; while the circular mask can be included by setting ang_1 != 0 
        ang_lat and ang_long are the latitude and longitude of the mask according to the jfactor computation scheme '''
    if ang_1 >= ang_2:
        return 0.
    if mask and ang_lat != 0. and ang_long < ang_2:
        def b_prime(ang):
            if np.isclose(ang, np.pi/2., atol = 0., rtol = 1e-4):
                return np.pi/2.
            return np.arctan(np.sqrt(np.power(np.sin(ang),2) - np.power(np.sin(ang_long),2)/np.cos(ang)))
        if (ang_long == 0.) or (ang_long <= ang_1 and ang_lat <= b_prime(ang_1)):
            return KPC_TO_CM * np.power(profile.rho_s, 2)*profile.r_s * _jfactor_simplified_region(ang_1 = ang_1, ang_2 = ang_2, ang_lat = ang_lat, profile_func = profile.functional_form(), x_sun = profile.r_sun/profile.r_s, integration_bound = integration_bound)
        else:
            if ang_lat >= b_prime(ang_2): # compute the mask with the simpler formula, by inverting latitude and longitude, exploiting spherical symmetry of the density profiles
                jfactor_mask_value = _jfactor_simplified_region(ang_1 = ang_1, ang_2 = ang_2, ang_lat = ang_long, profile_func = profile.functional_form(), x_sun = profile.r_sun/profile.r_s, integration_bound = integration_bound)
            else:
                jfactor_mask_value = _jfactor_mask(ang_1 = ang_1, ang_2 = ang_2, ang_lat = ang_lat, ang_long = ang_long, profile_func = profile.functional_form(), x_sun = profile.r_sun/profile.r_s, integration_bound = integration_bound)
            return KPC_TO_CM * np.power(profile.rho_s, 2)*profile.r_s * (_jfactor_cone(ang_1 = ang_1, ang_2 = ang_2, profile_func = profile.functional_form(), x_sun = profile.r_sun/profile.r_s, integration_bound = integration_bound) - jfactor_mask_value)
    else: # no mask = only cone integral
        return KPC_TO_CM * np.power(profile.rho_s, 2)*profile.r_s * _jfactor_cone(ang_1 = ang_1, ang_2 = ang_2, profile_func = profile.functional_form(), x_sun = profile.r_sun/profile.r_s, integration_bound = integration_bound)

def jfactor_gc(ang_2, ang_lat, ang_long, profile, strategy="normal", mask=True, ang_1=0, r_max=np.inf):
    """ x_max = r_max/r_s
    y_max = sqrt(x_max**2 + x_sun**2 - 2 x_max x_sun t)
    """
    assert strategy in ["normal", "inverted"], "Strategy '%s' not allowed" % strategy
    integration_bound = np.inf
    if r_max < 1e100:
        integration_bound = np.sqrt(np.power(r_max/profile.r_s,2) + np.power(profile.r_sun/profile.r_s,2) - 2 * r_max/profile.r_s * profile.r_sun/profile.r_s)
    if strategy == "normal" or not mask:
        return _jfactor_gc_with_mask(ang_2=ang_2, ang_lat=ang_lat, ang_long=ang_long, profile=profile, mask=mask, ang_1=ang_1, integration_bound=integration_bound)
    else:
        # we can compute the mask of this strategy by using the general J-factor formula in _jfactor_gc_with_mask and using
        # a mask inverting latitude and longitude. In this way we can obtain the mask on the direction perpendicular
        # to the galactic plane. However, given the spherical symmetry of the density profile, this mask can be used
        # also for the galactic plane. Finally we subtrack this mask from the whole cone integral.
        jfactor_mask_value = _jfactor_gc_with_mask(ang_2=ang_2, ang_lat=ang_long, ang_long=ang_lat, profile=profile, mask=mask, ang_1=ang_1, integration_bound=integration_bound)
        return _jfactor_gc_with_mask(ang_2=ang_2, profile=profile, mask=False, integration_bound=integration_bound) - jfactor_mask_value
