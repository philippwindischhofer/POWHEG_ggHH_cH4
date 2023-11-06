import numpy as np
import re
import math
import operator
import itertools
from scipy import interpolate
import scipy
import scipy.optimize
from math import sqrt
import random
import os, time

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import pylab as pl
#from mpl_toolkits.mplot3d import Axes3D
#import statsmodels.api as sm
#from phasespace import *

def combinegrids(grid_temp, cHHH, ct, ctt, cg, cgg):

    # Grid exists, proceed
    if os.path.exists(grid_temp):
        return

    # Lock to prevent other parallel processes from writing to the grid
    # Ensure lock file exists and create unique symlink to lock file
    os.system("touch lock")
    lockname = grid_temp + '.%s.lock' % os.getpid()
    os.link("lock", lockname)

    # If more than 1 symlink to lock file
    if os.stat("lock").st_nlink > 2:
        # Wait for other instance to create grid
        while not os.path.exists(grid_temp):
            print("Waiting for " + str(grid_temp) + " to be created")
            time.sleep(5)
        # Cleanup our symlink and return
        os.system("rm -f " + lockname)
        return
    # else this instance should produce the grid

    # Produce grid

    # Grid numbering format
    np.set_printoptions(formatter={'float': '{:.18E}'.format})

    print("******************************************************************************************")
    print(" Combining grids for chhh = " + str(cHHH) + ", ct = " + str(ct) + ", ctt = " + str(ctt) +", cg = " + str(cg) + ", cgg = " + str(cgg))
    print("******************************************************************************************")

    # Build grid for give value of cHHH
    incr = re.split('(_)', grid_temp)
    incr = "".join(incr[:len(incr)-10])
    cHHH_grids = [  incr + '_+5.000000E-01_+4.782609E-01_+1.000000E+00_+6.875000E-01_+8.888889E-01.grid',
  incr + '_-1.000000E+00_-2.500000E+00_+2.857143E-01_+4.666667E-01_+1.818182E-01.grid',
  incr + '_+4.545455E-02_-5.263158E-02_-6.875000E-01_+2.941176E-01_+6.923077E-01.grid',
  incr + '_+1.250000E-01_+1.111111E-01_-3.333333E-01_+2.173913E-01_-8.000000E-01.grid',
  incr + '_-1.142857E+00_+4.444444E-01_+9.090909E-02_+7.272727E-01_+6.428571E-01.grid',
  incr + '_+8.333333E-01_+4.545455E-01_-4.285714E-01_+9.523810E-02_+5.263158E-01.grid',
  incr + '_+1.000000E+00_+1.428571E-01_-3.333333E-01_+5.217391E-01_-1.500000E+00.grid',
  incr + '_-5.263158E-02_-7.142857E-02_+2.083333E-01_+6.923077E-01_+9.000000E-01.grid',
  incr + '_-3.750000E-01_+4.705882E-01_-3.750000E-01_-6.666667E-01_+4.166667E-01.grid',
  incr + '_-6.470588E-01_-8.333333E-01_+5.000000E-01_-1.333333E+00_+5.000000E-01.grid',
  incr + '_-6.666667E-01_+1.666667E-01_+1.100000E+00_+5.000000E-01_+1.000000E+00.grid',
  incr + '_+2.000000E-01_+6.000000E-01_+8.695652E-02_-3.478261E-01_-1.333333E+00.grid',
  incr + '_-4.210526E-01_+4.444444E-01_-2.222222E-01_+3.809524E-01_-5.714286E-01.grid',
  incr + '_-1.176471E-01_+2.200000E+00_-7.692308E-02_-1.875000E-01_+5.555556E-01.grid',
  incr + '_+2.000000E-01_-8.571429E-01_-1.000000E+00_+3.125000E-01_-1.166667E+00.grid',
  incr + '_+3.000000E-01_+1.111111E-01_+1.285714E+00_+1.285714E+00_-4.615385E-01.grid',
  incr + '_-4.347826E-01_-8.000000E-01_+1.111111E-01_-6.315789E-01_+4.347826E-02.grid',
  incr + '_-1.142857E+00_-3.333333E-01_-5.000000E-01_-5.000000E-01_+4.117647E-01.grid',
  incr + '_+2.250000E+00_-6.666667E-01_+2.727273E-01_+3.571429E-01_-1.000000E+00.grid',
  incr + '_+6.111111E-01_+2.777778E-01_+1.111111E-01_-8.000000E-01_+2.272727E-01.grid',
  incr + '_+2.173913E-01_+3.000000E+00_-5.263158E-01_+4.761905E-02_-3.809524E-01.grid',
  incr + '_+4.545455E-01_+4.000000E-01_-1.500000E+00_+5.454545E-01_+6.428571E-01.grid',
  incr + '_+2.500000E-01_+1.111111E+00_-4.166667E-01_-4.444444E-01_+5.000000E-02.grid']

    amps  = []
    ME2s  = []
    dME2s = []

    for grid in cHHH_grids:
        amps.append(np.loadtxt(grid, unpack=True))
    print("Grids loaded.")

    C = np.array([[14641/279841., 1., 121/2116., 121/1024., 64/81., 
  121/529., 1331/24334., 11/46., 11/32., 8/9., 
  1331/16928., 968/4761., 121/1472., 44/207., 11/36., 
  14641/194672., 121/368., 1331/11776., 121/414., 
  14641/135424., 121/256., 1331/8192., 121/288.], [625/16., 
  4/49., 25/4., 49/225., 4/121., 25/14., 125/8., 5/7., 
  -2/15., 4/77., -35/12., 25/22., -7/6., 5/11., 
  -14/165., -175/24., -1/3., 49/90., -7/33., 49/36., 
  14/225., -343/3375., 98/2475.], [1/130321., 121/256., 
  1/174724., 25/139876., 81/169., -11/5776., -1/150898., 
  1/608., -5/544., -99/208., 5/135014., 9/4693., 
  -5/156332., -9/5434., 45/4862., -5/116603., 55/5168., 
  -25/120802., -45/4199., 25/104329., -275/4624., 
  125/108086., 225/3757.], [1/6561., 1/9., 1/5184., 
  25/33856., 16/25., -1/243., 1/5832., -1/216., -5/552., 
  4/15., 5/14904., -4/405., 5/13248., -1/90., -1/46., 
  5/16767., -5/621., 25/38088., -4/207., 25/42849., 
  -25/1587., 125/97336., -20/529.], [256/6561., 1/121., 
  1024/3969., 4096/5929., 81/196., 16/891., -512/5103., 
  -32/693., -64/847., 9/154., -1024/6237., 8/63., 
  2048/4851., -16/49., -288/539., 512/8019., 32/1089., 
  -2048/7623., 16/77., 1024/9801., 64/1331., -4096/9317., 
  288/847.], [625/14641., 9/49., 625/4356., 25/3969., 
  100/361., -75/847., 625/7986., -25/154., -5/147., 
  -30/133., 125/7623., 250/2299., 125/4158., 125/627., 
  50/1197., 250/27951., -10/539., 50/14553., 100/4389., 
  100/53361., -4/1029., 20/27783., 40/8379.], [1/2401., 
  1/9., 1/49., 144/529., 9/4., -1/147., 1/343., -1/21., 
  -4/23., 1/2., 12/1127., -3/98., 12/161., -3/14., 
  -18/23., 12/7889., -4/161., 144/3703., -18/161., 
  144/25921., -48/529., 1728/12167., -216/529.], [1/38416., 
  25/576., 1/70756., 81/61009., 81/100., 5/4704., 
  1/52136., 5/6384., -15/1976., 3/16., -9/48412., 9/1960.,
   -9/65702., 9/2660., -81/2470., -9/35672., -15/1456., 
  81/44954., -81/1820., 81/33124., 135/1352., -729/41743., 
  729/1690.], [4096/83521., 9/64., 9/289., 1/16., 25/144.,
   -24/289., -192/4913., 9/136., -3/32., -5/32., 16/289., 
  80/867., -3/68., -5/68., 5/48., -1024/14739., 2/17., 
  -4/51., -20/153., 256/2601., -1/6., 1/9., 
  5/27.], [625/1296., 1/4., 3025/10404., 1936/2601., 1/4.,
   25/72., 1375/3672., 55/204., 22/51., 1/4., 275/459., 
  25/72., 1210/2601., 55/204., 22/51., 125/162., 5/9., 
  440/459., 5/9., 100/81., 8/9., 704/459., 
  8/9.], [1/1296., 121/100., 1/81., 1/9., 1., 11/360., 
  -1/324., -11/90., -11/30., 11/10., -1/108., 1/36., 
  1/27., -1/9., -1/3., 1/432., 11/120., -1/36., 1/12., 
  1/144., 11/40., -1/12., 1/4.], [81/625., 4/529., 
  9/625., 64/13225., 16/9., 18/575., 27/625., 6/575., 
  -16/2645., -8/69., -72/2875., -12/25., -24/2875., 
  -4/25., 32/345., -216/2875., -48/2645., 192/13225., 
  32/115., 576/13225., 128/12167., -512/60835., 
  -256/1587.], [256/6561., 4/81., 1024/29241., 4096/159201.,
   16/49., -32/729., -512/13851., 64/1539., 128/3591., 
  8/63., -1024/32319., -64/567., 2048/68229., 128/1197., 
  256/2793., 512/15309., -64/1701., -2048/75411., 
  -128/1323., 1024/35721., -128/3969., -4096/175959., 
  -256/3087.], [14641/625., 1/169., 484/7225., 9/18496., 
  25/81., -121/325., -2662/2125., 22/1105., -3/1768., 
  -5/117., 363/3400., 121/45., -33/5780., -22/153., 
  5/408., -3993/2000., 33/1040., -99/10880., -11/48., 
  1089/6400., -9/3328., 27/34816., 5/256.], [1296/2401., 
  1., 36/1225., 1/256., 49/36., -36/49., -216/1715., 
  6/35., -1/16., 7/6., 9/196., -6/7., -3/280., 1/5., 
  -7/96., -135/686., 15/56., -15/896., 5/16., 225/3136., 
  -25/256., 25/4096., -175/1536.], [1/6561., 81/49., 
  1/900., 729/4900., 36/169., 1/63., 1/2430., 3/70., 
  243/490., -54/91., 1/210., -2/351., 9/700., -1/65., 
  -81/455., 1/567., 9/49., 27/490., -6/91., 1/49., 
  729/343., 2187/3430., -486/637.], [256/625., 1/81., 
  64/529., 14400/190969., 1/529., 16/225., 128/575., 
  8/207., 40/1311., 1/207., 384/2185., 16/575., 
  960/10051., 8/529., 120/10051., 768/2375., 16/285., 
  1152/8303., 48/2185., 2304/9025., 16/361., 17280/157757., 
  144/8303.], [1/81., 1/4., 64/441., 16/49., 49/289., 
  -1/18., 8/189., -4/21., -2/7., -7/34., 4/63., 7/153., 
  32/147., 8/51., 4/17., 1/54., -1/12., 2/21., 7/102., 
  1/36., -1/8., 1/7., 7/68.], [16/81., 9/121., 9/4., 
  2025/3136., 1., 4/33., -2/3., -9/22., 135/616., 
  -3/11., 5/14., -4/9., -135/112., 3/2., -45/56., 
  -20/189., -5/77., -75/392., 5/21., 25/441., 75/2156., 
  1125/10976., -25/196.], [625/104976., 1/81., 3025/104976.,
   484/2025., 25/484., 25/2916., 1375/104976., 55/2916., 
  -22/405., 5/198., -55/1458., 125/7128., -121/1458., 
  25/648., -1/9., -25/1458., -2/81., 44/405., -5/99., 
  4/81., 16/225., -352/1125., 8/55.], [81., 100/361., 
  225/529., 25/233289., 64/441., -90/19., 135/23., 
  -150/437., -50/9177., 80/399., 15/161., -24/7., 
  25/3703., -40/161., -40/10143., 9/7., -10/133., 5/3381.,
   -8/147., 1/49., -10/8379., 5/213003., 
  -8/9261.], [16/625., 9/4., 4/121., 900/14641., 81/196., 
  -6/25., 8/275., -3/11., -45/121., -27/28., 24/605., 
  18/175., 60/1331., 9/77., 135/847., 48/1375., -18/55., 
  72/1331., 54/385., 144/3025., -54/121., 1080/14641., 
  162/847.], [10000/6561., 25/144., 25/324., 1/81., 
  1/400., -125/243., 250/729., -25/216., 5/108., -1/48., 
  -100/729., 5/81., -5/162., 1/72., -1/180., -4000/6561., 
  50/243., 40/729., -2/81., 1600/6561., -20/243., 
  -16/729., 4/405.]])

    Cinv = np.linalg.inv(C)
    coeffs = np.array([ct**4,ctt**2,ct**2*cHHH**2,cg**2*cHHH**2,cgg**2,ctt*ct**2,ct**3*cHHH,
                       ctt*ct*cHHH,ctt*cg*cHHH,ctt*cgg,ct**2*cg*cHHH,ct**2*cgg,
                       ct*cHHH**2*cg,ct*cHHH*cgg,cg*cHHH*cgg,
                       ct**3*cg,ct*ctt*cg,ct*cg**2*cHHH,ct*cg*cgg,
                       ct**2*cg**2,ctt*cg**2,cg**3*cHHH,cg**2*cgg])

    # Check that the grids have the same values for s, t
    for amp in amps[0:]:
      for amp2 in amps[0:]:

        if not (np.array_equal(amp[0], amp2[0]) and np.array_equal(amp[1], amp2[1])):
            print("The virtual grids do not contain the same phase-space points!")
    for i, psp in enumerate(amps[0][0]):
        A = np.matmul(Cinv, np.transpose(np.array([ [amps[j][2][i] for j in range(len(cHHH_grids))] ])))
#        A = np.matmul(Cinv, np.array([[amps[0][2][i]], [amps[1][2][i]], [amps[2][2][i]]]))
        ME2 = np.matmul(coeffs, A)

        # Compute the uncertainties on the final PSP amplitude (uncorr. between PSPs, corr. between A coeffs)
        sigmaF = np.diag( np.array([ amps[j][3][i]**2 for j in range(len(cHHH_grids))] ))
        dA = np.matmul(Cinv, np.matmul(sigmaF, np.transpose(Cinv)))
        dME2 = np.matmul(coeffs, np.matmul(coeffs, dA))

        ME2s.append(ME2)
        dME2s.append(np.sqrt(dME2))

    np.savetxt(grid_temp, np.transpose([amps[0][0], amps[0][1], ME2s, dME2s]))

    print("Saved grid " + str(grid_temp))
    os.system("rm -f " + lockname)

class Bin:
    def __init__(self):
        self.n = 0
        self.y = 0.
        self.y2 = 0.
        self.e2 = 0.

    def add(self, y, e=0.):
        self.n += 1
        self.y += y
        self.y2 += y * y
        self.e2 += e * e

    def gety(self, sample=False):
        if self.n == 0:
            return float('nan')
        if sample:
            return random.gauss(*self.getres())
        return self.y / self.n

    def gete(self):
        if self.n > 1:
            return math.sqrt((self.y2 - self.y ** 2 / self.n) / (self.n - 1) + self.e2 / self.n ** 2)
        if self.n > 0:
            return math.sqrt(self.e2) / self.n
        return float('nan')

    def getres(self):
        return (self.gety(), self.gete())

    def __str__(self):
        return str(self.getres())

class Grid:
    def __init__(self, method, dim, xmin, xmax, nbin, data=[], binned=False):
        self.dim = dim
        self.xmin = np.array(xmin, dtype=float)
        self.xmax = np.array(xmax, dtype=float)
        self.nbin = np.array(nbin)
        self.binned = binned
        self.method = method

        self.dx = np.abs((self.xmin - self.xmax).astype(float) / self.nbin)

        if (not binned):
            self.nbin = self.nbin + 1

        self.data = {}
        self.deg = {}  # degree of interpolation poynomial
        self.pol = {}  # coefficients of interpolation polynomial
        for k in np.ndindex(tuple(self.nbin)):
            self.data[k] = Bin()

        if (binned):
            if (data):
                self.addPoint(data)
        else:
            # nneighbor=data.size/np.prod(self.nbin)
            nneighbor = 5
            nsamples = 1

            def linfunc(x, a, b, c):
                return a + x[0] * b + x[1] * c

            def linfunc1d(x, a, b):
                return a + x[0] * b

            def constfunc(x, a):
                return a

            def func2(x, a, b, c, d, e, f):
                return a + x[0] * b + x[1] * c + d * x[0] ** 2 + e * x[1] ** 2 + f * x[0] * x[1]

            for _ in range(nsamples):
                if nsamples > 1:
                    dat = np.array(
                        [np.append(d[:self.dim], [random.gauss(d[self.dim], d[self.dim + 1]), d[self.dim + 1]]) for
                         d in data])  # dat = data with y values sampled according to gauss(y,e)
                else:
                    dat = data
                for k in np.ndindex(tuple(self.nbin)):
                    x = self.xmin + k * self.dx
                    points = np.array(sorted(dat, key=lambda p: np.linalg.norm(p[:self.dim] - x))[:nneighbor])
                    X = points[:, 0:2]
                    #          X=sm.add_constant(X)
                    Y = points[:, 2]
                    E = 1 / points[:, 3] ** 2  # * np.linalg.norm(points[:,:self.dim]- x)**(-2)

                    # sigma=[ p[3] * np.linalg.norm(p[:self.dim]-x)**2 for p in points] # with absolue_sigma=False, each point will contribute with relative weights w=1/sigma**2, non-linear dependence will increase quadracically with the distance and can be taken into account in the weight
                    sigma = points[:, 3]
                    X = points[:, 0:2]
                    X = np.array([xx - x for xx in X]).T
                    # popt, pcov = scipy.optimize.curve_fit(func2,X,Y)
                    #          popt, pcov = scipy.optimize.curve_fit(func2,X,Y,sigma=points[:,3],absolute_sigma=True)
                    # popt, pcov = scipy.optimize.curve_fit(linfunc,X,Y)
                    # popt, pcov = scipy.optimize.curve_fit(linfunc,X,Y,sigma=points[:,3],absolute_sigma=True)
                    #          popt, pcov = scipy.optimize.curve_fit(constfunc,X,Y,sigma=points[:,3],absolute_sigma=True)
                    #          popt, pcov = scipy.optimize.curve_fit(func2,X,Y,sigma=points[:,3])
                    popt, pcov2 = scipy.optimize.curve_fit(linfunc, X, Y, sigma=sigma, absolute_sigma=True)
                    popt, pcov1 = scipy.optimize.curve_fit(linfunc, X, Y, sigma=sigma, absolute_sigma=False)
                    # popt, pcov = scipy.optimize.curve_fit(linfunc,X,Y,sigma=points[:,3],absolute_sigma=False)
                    #          print "fitval: ",linfunc(x,*popt)
                    self.data[k].add(popt[0], max(sqrt(pcov1[0, 0]), sqrt(
                        pcov2[0, 0])))  # corresponds to multiplying error estimate of chi-sq fit with max(1,chisq)
                    '''
                       model y = f(x) = c
                       y1   e1     y2   e2     popt   pcov1   pcov2    chisq
                       1    1e-5   2    1e-5   1.5    5e-11   0.25     5e9      # pcov2 = pcov1*chisq
                       0.9  1      1.1  1      1      0.5     0.0025   0.02
                       ==> use cov2 if chisq>>1, cov1 if chisq<1
                    '''
                    # print x,np.mean(Y),popt[0],sqrt(pcov1[0,0]),sqrt(pcov2[0,0])
                    # self.data[k].add(linfunc(x,*popt))

                    #          print "polyval:  ",polyval2d(x[0],x[1],fit,degree)
                    ##          smx=sm.add_constant(x)
                    #          smx=[1.,x[0], x[1]]
                    #          print x,smx
                    #          print "wls.pred: ",wls_model.predict(smx)
                    #
                    #          raw_input('')

    def polyfit2d(self, x, y, z, orders):
        ncols = np.prod(orders + 1)
        G = np.zeros((x.size, ncols))
        ij = itertools.product(range(orders[0] + 1), range(orders[1] + 1))
        for k, (i, j) in enumerate(ij):
            G[:, k] = x ** i * y ** j
        m, _, _, _ = np.linalg.lstsq(G, z, rcond=None)
        return m

    def polyval2d(self, x, y, m, orders):
        ij = itertools.product(range(orders[0] + 1), range(orders[1] + 1))
        z = np.zeros_like(x)
        for a, (i, j) in zip(m, ij):
            z += a * x ** i * y ** j
        return z

    def x(self, k):  # coordinate after variable transformation applied to s,t
        if (self.binned):
            return (self.xmin + (np.array(k) + 0.5) * self.dx)
        else:
            return (self.xmin + (np.array(k) + 0.) * self.dx)

    def k(self, x):  # index of data array
        #    if (not all(x>self.xmin and x < self.xmax)):
        #      raise Exception('xout of bounds: ',x)
        return tuple((np.minimum((x - self.xmin) / self.dx, self.nbin - 1)).astype(int))

    def gridpoints(self, sample=False, flag='k', extendToBorder=False, returnError=False):
        if flag == 'k':
            return [[k, d.gety(sample)] for k, d in self.data.items()]
        if flag == 'x':
            if (returnError):
                res = [[self.x(k), d.gety(sample), d.gete()] for k, d in self.data.items()]
            else:
                res = [[self.x(k), d.gety(sample)] for k, d in self.data.items()]
            if (extendToBorder):
                for d in range(self.dim):
                    nbin = self.nbin.copy()
                    nbin[d] = 1
                    for k in np.ndindex(tuple(nbin)):
                        x = self.x(k)
                        x[d] = self.xmin[d]
                        res += [[x.copy(), self.data[k].gety(sample)], ]
                        k = list(k)
                        k[d] = self.nbin[d] - 1
                        x[d] = self.xmax[d]
                        res += [[x.copy(), self.data[tuple(k)].gety(sample)], ]
                        #            print "b", self.xmin,self.xmax
                for corner in np.ndindex(tuple([2, ] * self.dim)):
                    k = corner * (self.nbin - 1)
                    x = np.array(corner, float)
                    res += [[x.copy(), self.data[tuple(k)].gety(sample)], ]
            return res
        if flag == 'plain':
            return [self.x(k).tolist() + [d.gety(sample), ] for k, d in self.data.items()]

    def printgrid(self):
        for k, d in self.data.items():
            print(k, d, d.n)

    def addPoint(self, data):
        if (type(data[-1]) is float or type(data[-1]) is np.float64):
            assert (len(data) == self.dim + 2)
            xx = data[0:self.dim]
            kk = self.k(xx)
            self.data[kk].add(*data[-2:])
        elif (type(data[-1]) is np.ndarray):
            for d in data:
                self.addPoint(d)

    def initInterpolation(self, nsamples=1):
        self.interpolators = []
        self.nsamples = nsamples
        for i in range(nsamples):
            temp = list(zip(*self.gridpoints(sample=(nsamples != 1), flag='x', extendToBorder=True)))
            self.interpolators.append(interpolate.CloughTocher2DInterpolator(list(temp[0]), temp[1], fill_value=0.))
            #      temp=zip(*self.gridpoints(sample=(nsamples!=1),flag='plain',extendToBorder=False))
            #      self.interpolators.append(interpolate.SmoothBivariateSpline(temp[0],temp[1],temp[2],bbox=[0.,1.,0.,1.]))

            # polynomial interpolation
            nneighbor = 2
            for k in self.data:
                kk = np.asarray(k)
                kmin = np.maximum(kk - np.full(self.dim, nneighbor, dtype=int), np.zeros(self.dim)).astype(int)
                kmax = np.minimum(kk + np.full(self.dim, nneighbor, dtype=int), self.nbin - 1).astype(int)
                degree = kmax - kmin
                self.deg[k] = degree

                xx = self.x(k)
                dat = [np.append(self.x(x + kmin), self.data[tuple(x + kmin)].gety(sample=(nsamples != 1))) for x in
                       np.ndindex(*(degree + 1))]
                x, y, z = np.asarray(list(zip(*dat)))
                #        z=np.asarray(zip(*dat)[2])
                self.pol[k, i] = self.polyfit2d(x, y, z, orders=degree)


                #        xx=self.x(k)
                #        dat=[ [x+kmin[0], y+kmin[1], self.data[(x+kmin[0],y+kmin[1])].gety(sample=(nsamples!=1)) ] for x,y in np.ndindex(*(degree+1))]
                #        print kk,degree,xx
                #        print dat
                #        x,y=np.asarray(zip(*dat)[0:2],dtype=int)
                #        z=np.asarray(zip(*dat)[2])
                #        self.pol[k] = polyfit2d(x,y,z,degree)
                #        for x,y in np.ndindex(*(degree+1))+kmin:
                #          print np.array([x,y])
                #        pol=np.Poly

    def interpolate(self, x, y, selectsample=-1):
        if (self.method == 2):
            k = self.k(tuple([x, y]))
            temp = [self.polyval2d(x, y, self.pol[k, i], self.deg[k]) for i in range(self.nsamples)]
            # return (self.polyval2d(x,y,self.pol[k],self.deg[k]),0.)
        else:
            temp = [interpolator(x, y) for interpolator in self.interpolators]
            # temp=[ interpolator(float(x[0]),float(x[1])) for interpolator in self.interpolators]
            #      return (np.mean(temp),np.std(temp))
        if (selectsample == -1):
            return (np.mean(temp), np.std(temp))
        else:
            return (temp[selectsample], 0.)

    def __call__(self, X, Y, selectsample=-1):
        temp = np.array([X, Y]).T
        s = np.shape(X)
        temp = np.array([self.interpolate(x, y, selectsample)[0] for x, y in
                         np.array([x for x in np.vstack([X.ravel(), Y.ravel()])]).T]).reshape(s)

        #    s= tuple(list(np.shape(X)) + [2,])
        #    temp= np.array([ self.interpolate(x,y) for x,y in  np.array([ x for x in np.vstack([X.ravel(), Y.ravel()]) ]).T ]).reshape(s)

        return temp

class CreateGrid:
    def __init__(self,selected_grid):
        self.selected_grid = selected_grid
        self.selected_grid_dirname = os.path.dirname(self.selected_grid)
        self.mHs = 125.**2
        self.method = 1  # 1: CloughTocher;  2: Polynomial
        self.flatten = False
        self.polydegree = 2
        self.tolerance = 1.e-8 # costh is nudged to 1 if within this tolerance

        # cdf: transformation x=f(beta) leading to uniform distributions of points in x
        self.cdfdata = np.loadtxt(os.path.join(self.selected_grid_dirname,'events.cdf'))
        self.cdf = interpolate.splrep(*self.cdfdata.T, k=1)

        # fig=pl.figure()
        # ax = fig.add_subplot(111)
        # ax.scatter(*cdfdata.T)
        # x_grid = np.arange(0, 1, 0.01)
        # y_grid=interpolate.splev(x_grid,cdf)
        # ax.plot(x_grid,y_grid)
        # pl.show()
        # exit(1)

        # dataInUniform=np.loadtxt(selected_grid,converters={0: lambda sbeta: interpolate.splev(float(sbeta),cdf)})

        x0, x1, y, e = np.loadtxt(self.selected_grid, unpack=True)
        #dataIn = np.array([x0, x1, y, e]).T

        x0Uniform = interpolate.splev(x0, self.cdf)
        #dataInUniform = np.array([x0Uniform, x1, y, e]).T

        if (self.flatten):
            poly = np.poly1d(np.polyfit(x0Uniform, y, self.polydegree, w=1. / e))
            yFlat = y / poly(x0Uniform)

            # xp = np.linspace(0, 1, 100)
            # fig1=plt.figure()
            # ax1=fig1.gca()
            # ax1.plot(x0Uniform,y,'.',xp,poly(xp),'--')

            # fig2=plt.figure()
            # ax2=fig2.gca()
            # ax2.axhline(y=1., xmin=0, xmax=1,ls='--')
            # ax2.plot(x0Uniform,yFlat,'.')
        else:
            poly = np.poly1d([1.])
            yFlat = y

        dataInUniformFlat = np.array([x0Uniform, x1, yFlat, e]).T

        # fig0=plt.figure()
        # ax0 = fig0.add_subplot(111, projection='3d')
        # ax0.scatter(*(dataInUniformFlat.T),s=5,color='blue')
        # for d in dataInUniformFlat:
        #  ax0.plot([d[0],d[0]],[d[1],d[1]],[d[2]-d[3],d[2]+d[3]],'_',color='black')
        # plt.show()


        self.interpolator = Grid(self.method, 2, (0, 0), (1, 1), (100, 30), data=dataInUniformFlat, binned=False)

        # a.addPoint(dataInUniformFlat)
        self.interpolator.initInterpolation(1)
        #self.interpolator.printgrid()

    def u(self, s, t):
        return 2 * self.mHs - s - t

    def beta(self, s):
        return math.sqrt(1. - 4. * self.mHs / s)

    def costh(self, s, t):
        res = (t - self.u(s, t)) / (s * self.beta(s))
        if (abs(res) - 1.) > 0. :
            if (abs(res) - 1.) < self.tolerance:
                if res < 0.:
                    res = -1.
                else:
                    res = 1.
            else:
                raise ValueError("Grid called with cos(theta) > 1, cos(theta) = {:.40e}".format(res))
        return res

    def GetAmplitude(self, s, t):
        b = self.beta(s)
        cTh = abs(self.costh(s, t))
        xUniform = interpolate.splev(b, self.cdf)
        #  print s,t
        #  print b,cTh
        #  print a.interpolate(xUniform,cTh)[0]
        return self.interpolator.interpolate(xUniform, cTh)[0]


        # print a.gridpoints(flag='x')
        #
        # fout=open("output","w")
        #
        ####
        #####pspoints=[PSpoint() for i in range(400)]
        #####pspoints=[[s,cosTh,w,a.interpolate(interpolate.splev(beta(s),cdf),abs(cosTh))[0]] for s,cosTh,w in [PSpoint() for _ in range(10)] ]
        #####pspoints=[w*BornHTL(s,s/4,pdf) for s,cosTh,w in [PSpoint() for _ in range(10000)] ]
        ####pspoints=[w*a.interpolate(interpolate.splev(beta(s),cdf),abs(cosTh))[0] for s,cosTh,w in [PSpoint() for _ in range(1000000)] ]
        ####print pspoints
        ####print np.mean(pspoints)
        ####print np.std(pspoints,ddof=1)/sqrt(len(pspoints))
        ####
        ####exit(1)
        ####
        #
        # print a.interpolate(0.9032,0.4)
        #
        ##random.shuffle(dataIn)
        # diff=[]
        # nstd=[]
        # sumdiff=0.
        #
        #
        # for d in dataIn[1500:]:
        #  xUniform=interpolate.splev(d[0],cdf)
        #  gridres=a.interpolate(xUniform,d[1])
        #  if(flatten):
        #    gridres=[gridres[0]*poly(xUniform),gridres[1]*poly(xUniform)]
        #  fout.write("s/mTs= %5.2f  x=(%5.2f,%5.2f)    amp: %7.4e +- %7.4e     grid:  %7.4e +- %7.4e   diff[%%]:  %5.2f    #stddev(amp): %5.2f   #stddev(grid): %5.2f   #stddev(gr+am): %5.2f\n" % ((125./173)**2*4./(1.-d[0]**2),xUniform,d[1],d[2],d[3],gridres[0],gridres[1], ((gridres[0]-d[2])/d[2])*100.,abs(gridres[0]-d[2])/d[3],abs(gridres[0]-d[2])/gridres[1],abs(gridres[0]-d[2])/sqrt(gridres[1]**2+d[3]**2)))
        #  print("s/mTs= %5.2f  x=(%5.2f,%5.2f)    amp: %7.4e +- %7.4e     grid:  %7.4e +- %7.4e   diff[%%]:  %5.2f    #stddev(amp): %5.2f   #stddev(grid): %5.2f   #stddev(gr+am): %5.2f  " % ((125./173)**2*4./(1.-d[0]**2),xUniform,d[1],d[2],d[3],gridres[0],gridres[1],
        #                                             ((gridres[0]-d[2])/d[2])*100.,abs(gridres[0]-d[2])/d[3],abs(gridres[0]-d[2])/gridres[1],abs(gridres[0]-d[2])/sqrt(gridres[1]**2+d[3]**2)))
        #  diff+= [((gridres[0]-d[2])/d[2]) ,]
        #  nstd+= [abs((gridres[0]-d[2])/sqrt(gridres[1]**2+d[3]**2)) ,]
        # print 'N points:          ', len(diff)
        # print 'sum diff:          ', np.sum(diff)
        # print 'sum abs diff:      ', np.sum(np.abs(diff))
        # print 'mean sgn diff [%]: ', np.mean(diff)*100.
        # print 'mean abs diff [%]: ', np.mean(np.abs(diff))*100.
        # print 'medn abs diff [%]: ', np.median(np.abs(diff))*100.
        # print
        # print 'mean nstd: ', np.mean(nstd)
        # print 'medn nstd: ', np.median(nstd)
        # print '68 %-ile:   ', np.percentile(nstd,68.3)
        # print '95 %-ile:   ', np.percentile(nstd,95.45)
        # print '99 %-ile:   ', np.percentile(nstd,99.73)
        ##  print 'input data: ',d,
        ##  print interpolate.splev(float(d[0]),cdf),
        ##  print a.interpolate(interpolate.splev(d[0],cdf),d[1]),
        ##  print a.interpolate(interpolate.splev(d[0],cdf),d[1])/d[2]
        #
        # xmax=1.0
        # fout.close()
        #
        # if(True): #plot
        #  fig3=plt.figure()
        #  ax3 = fig3.add_subplot(111, projection='3d')
        #
        #  x_grid = np.arange(0, xmax+.00001, 0.01)
        #  y_grid = np.arange(0, 1.00001, 0.01)
        #  x_grid, y_grid = np.meshgrid(x_grid, y_grid)
        #
        #  ax3.set_xlabel(r'f(s)')
        #  ax3.set_ylabel(u'| cos(\u03B8) |')
        #  #ax3.set_ylabel(u'| cos() |')
        #  ax3.set_zlabel(r'M^2')
        #  ax3.set_xlim(0.,xmax)
        ##  ax3.set_zlim(0.,6e-6)
        #
        #
        #
        ##plot interpolation
        #  #z1_grid = a.interpolators[0].ev(x_grid,y_grid)
        #  #z1_grid = a.interpolators[0](x_grid,y_grid)
        #  z1_grid = a(x_grid,y_grid,selectsample=0)
        #  ax3.plot_surface(x_grid, y_grid, z1_grid,color='yellow',alpha=0.4)
        #
        #  #z1_grid = a.interpolators[1](x_grid,y_grid)
        #  z1_grid = a(x_grid,y_grid,selectsample=1)
        #  ax3.plot_surface(x_grid, y_grid, z1_grid,color='orange',alpha=0.4)
        #
        ## input points
        #  ax3.scatter(*(dataInUniformFlat.T),s=5)
        ##  for d in dataInUniformFlat:
        ##    ax3.plot([d[0],d[0]],[d[1],d[1]],[d[2]-d[3],d[2]+d[3]],'_',color='black')
        #
        ##grid
        #  temp=[[x[0],x[1],y] for x,y in a.gridpoints(flag='x',extendToBorder=False,returnError=False) ]
        #  ax3.scatter(*(zip(*temp)),color='red',s=20)
        #  temp=[[x[0],x[1],y,e] for x,y,e in a.gridpoints(flag='x',extendToBorder=False,returnError=True )]
        ##  for d in temp:
        ##    ax3.plot([d[0],d[0]],[d[1],d[1]],[d[2]-d[3],d[2]+d[3]],'|',color='red')
        #  z1_grid=a(x_grid,y_grid)
        ##  temp=[[x,x,a.interpolate(x[0],x[1])[0]] for x,y in zip(x_grid,y_grid) ]
        ##  ax3.plot_surface(x_grid,y_grid,z1_grid,color='red',s=20)
        ##  ax3.scatter(*(zip(*temp)),color='orange',s=20)
        ##
        #  tempdata=[ [x0,x1,y-a.interpolate(x0,x1)[0]] for x0,x1,y,e in dataInUniformFlat ]
        #  fig4=plt.figure()
        #  ax4 = fig4.add_subplot(111, projection='3d')
        #  ax4.scatter(*(zip(*tempdata)))
        #  ax4.set_title("absolute difference")
        #  ax4.set_xlim(0.,xmax)
        #  ax4.set_zlim(-6e-6,6e-6)
        #
        #  tempdata=[ [x0,x1,(y-a.interpolate(x0,x1)[0])/y] for x0,x1,y,e in dataInUniformFlat ]
        #  fig5=plt.figure()
        #  ax5 = fig5.add_subplot(111, projection='3d')
        #  ax5.scatter(*(zip(*tempdata)))
        #  ax5.set_title("relative difference")
        #  ax5.set_xlim(0.,xmax)
        #  ax5.set_zlim(-1.5,1.5)
        #
        #  tempdata=[ [x0,x1,(y-a.interpolate(x0,x1)[0])/sqrt(e**2+a.interpolate(x0,x1)[1]**2)] for x0,x1,y,e in dataInUniformFlat ]
        #  fig6=plt.figure()
        #  ax6 = fig6.add_subplot(111, projection='3d')
        #  ax6.scatter(*(zip(*tempdata)))
        #  ax6.set_title("#std. dev.")
        #  ax6.set_xlim(0.,xmax)
        #
        ##  z1_grid *= poly(x_grid)
        ##
        ###  ax.plot_surface(x_grid, y_grid, z1_grid,color='yellow')
        ##  temp=[[list(a.x(x))[0],list(a.x(x))[1],y*poly(a.x(x)[0])] for x,y in a.gridpoints() ]
        ##  ax3.scatter(*(dataInUniform.T))
        ##  ax3.scatter(*(zip(*temp)),color='red')
        #  #print temp
        #  plt.show()
        #  ##with np.loadtxt(selected_grid) as data:
        #  ##  for d in data:
        #  #

