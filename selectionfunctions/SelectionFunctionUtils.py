import numpy as np
from scipy import interpolate

class littlewoodpaley:

    def __init__(self, B = 2.0):

        self.B = B
        self.psi_spline = interpolate.splrep ( \
        np.arange (-1.01, 0.02, 0.01),
        np.array([  0.00000000e+00,   0.00000000e+00,   6.10726446e-26,
	        1.80473593e-14,   1.63146885e-10,   1.81011396e-08,
	        3.33941762e-07,   2.47115014e-06,   1.07501585e-05,
	        3.33635137e-05,   8.23638779e-05,   1.72785830e-04,
	        3.21411357e-04,   5.45573939e-04,   8.62196482e-04,
	        1.28711301e-03,   1.83464846e-03,   2.51740299e-03,
	        3.34618479e-03,   4.33004296e-03,   5.47636332e-03,
	        6.79099953e-03,   8.27842094e-03,   9.94186438e-03,
	        1.17834820e-02,   1.38044808e-02,   1.60052501e-02,
	        1.83854783e-02,   2.09442559e-02,   2.36801676e-02,
	        2.65913725e-02,   2.96756753e-02,   3.29305873e-02,
	        3.63533793e-02,   3.99411282e-02,   4.36907558e-02,
	        4.75990635e-02,   5.16627608e-02,   5.58784904e-02,
	        6.02428494e-02,   6.47524071e-02,   6.94037205e-02,
	        7.41933466e-02,   7.91178536e-02,   8.41738297e-02,
	        8.93578906e-02,   9.46666853e-02,   1.00096901e-01,
	        1.05645269e-01,   1.11308565e-01,   1.17083611e-01,
	        1.22967283e-01,   1.28956505e-01,   1.35048255e-01,
	        1.41239561e-01,   1.47527507e-01,   1.53909226e-01,
	        1.60381906e-01,   1.66942786e-01,   1.73589155e-01,
	        1.80318352e-01,   1.87127766e-01,   1.94014833e-01,
	        2.00977036e-01,   2.08011904e-01,   2.15117011e-01,
	        2.22289973e-01,   2.29528448e-01,   2.36830134e-01,
	        2.44192769e-01,   2.51614129e-01,   2.59092025e-01,
	        2.66624305e-01,   2.74208849e-01,   2.81843571e-01,
	        2.89526414e-01,   2.97255354e-01,   3.05028392e-01,
	        3.12843559e-01,   3.20698910e-01,   3.28592527e-01,
	        3.36522513e-01,   3.44486996e-01,   3.52484123e-01,
	        3.60512062e-01,   3.68568999e-01,   3.76653139e-01,
	        3.84762704e-01,   3.92895928e-01,   4.01051064e-01,
	        4.09226374e-01,   4.17420136e-01,   4.25630637e-01,
	        4.33856174e-01,   4.42095054e-01,   4.50345591e-01,
	        4.58606108e-01,   4.66874931e-01,   4.75150394e-01,
	        4.83430833e-01,   4.91714588e-01,   5.00000000e-01,
	        5.08285412e-01]))

    def psi (self, u):
        """Estimate the psi function.

        "Psi" is the name of a function defined in the article by Marinucci et al.
        (2008) that is used to build the actual needlet."""

        neg_u = np.clip (-np.abs (u), -1.0, 0.0)
        value = interpolate.splev (neg_u, self.psi_spline)

        if np.isscalar (u):
            if u > 0.0:
                return 1.0 - value
            else:
                return value
        else:
            u = np.array (u)  # Ensure that "u" is of the proper type
            return np.where (u > 0.0, 1 - value, value)

    def phi (self, t):
        """Estimate the phi function.

        "Phi" is the name of a function defined in the article by Marinucci et al.
        (2008) that is used to build the actual needlet."""

        # Ensure that "t" is of the correct type
        if not np.isscalar (t): t = np.array (t)
        val = np.clip (1 - 2*self.B/(self.B - 1) * (t - 1.0/self.B), -1.0, 1.0)
        return self.psi (val)

    def window_function (self, l, j):
        u = l * np.power(self.B,-j)
        return np.sqrt (np.clip (self.phi (u / self.B) - self.phi (u), 0.0, 5.0))

    def start(self, j):
        return int(np.floor(self.B**(j-1)))

    def end(self, j):
        return int(np.ceil(self.B**(j+1)))

class chisquare:

    def __init__(self, j, p = 1.0, B = 2.0, F = 1e-6, normalise = False):

        self.j = np.array([_j for _j in j if _j >= 0])
        self.p = p
        self.B = B
        self.F = F
        self.normalise = normalise
        self.compute_normalisation()
        self.compute_needlet_normalisation()

    def window_function(self, l, j):
        u = l*(l+1) / np.power(self.B,2.0*j)
        N = self.normalisation[l.astype(np.int)] if type(l) == np.ndarray else self.normalisation[int(l)]

        return N*np.power(u,self.p)*np.exp(-u)

    def compute_normalisation(self):

        self.lmax = self.end(max(self.j))
        self.normalisation = np.ones(self.lmax+1)

        if self.normalise == True:
            # Renormalise (Marinucci 2008) Equation 2
            jinf = np.arange(1000)
            for l in range(1,self.lmax+1):
                self.normalisation[l] = 1.0/np.sum(np.square(self.window_function(l,jinf)))

    def compute_needlet_normalisation(self):

        self.needlet_normalisaiton = np.ones(len(self.j)+1)

        for ineedlet, j in enumerate(self.j):
            if j==-1:
                self.needlet_normalisaiton[ineedlet]=1.0
                continue

            start = self.start(j)
            end = self.end(j)
            modes = np.arange(start, end + 1, dtype = 'float')
            window = self.window_function(modes,j)*(2.0*modes+1.0)/np.sqrt(4.0*np.pi)#*npix_needle)

            self.needlet_normalisaiton[ineedlet] = np.sum(window)

    def start(self, j):
        return 1

    def end(self, j):
        from scipy import special
        G = -self.p*special.lambertw(-np.power(self.F,1.0/self.p)/np.e,k=-1).real*np.power(self.B,2.0*j)
        return int(np.ceil(0.5*(-1.0+np.sqrt(1.0+4.0*G))))
