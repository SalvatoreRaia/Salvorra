from math import *
from random import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import zipf
from scipy import stats
from _functools import  reduce
import pandas as pd


# function x=randht(n, varargin)

# RANDHT generates n observations distributed as some continous heavy-
# tailed distribution. Options are power law, log-normal, stretched
# exponential, power law with cutoff, and exponential. Can specify lower
# cutoff, if desired.
#
#    Example:
#       x = randht(10000,'powerlaw',alpha);
#       x = randht(10000,'xmin',xmin,'powerlaw',alpha);
#       x = randht(10000,'cutoff',alpha, lambda);
#       x = randht(10000,'exponential',lambda);
#       x = randht(10000,'lognormal',mu,sigma);
#       x = randht(10000,'stretched',lambda,beta);
#
#    See also PLFIT, PLVAR, PLPVA
#
#    Source: http://www.santafe.edu/~aaronc/powerlaws/


# Version 1.0.2 (2008 April)
# Copyright (C) 2007 Aaron Clauset (Santa Fe Institute)

# Ported to python by Joel Ornstein (2011 August)
# (joel_ornstein@hmc.edu)

# Distributed under GPL 2.0
# http://www.gnu.org/copyleft/gpl.html
# RANDHT comes with ABSOLUTELY NO WARRANTY
#
# Notes:
#
def randht(n, *varargin):
    Type = '';
    xmin = 1;
    alpha = 2.5;
    beta = 1;
    Lambda = 1;
    mu = 1;
    sigma = 1;

    # parse command-line parameters; trap for bad input
    i = 0;
    while i < len(varargin):
        argok = 1;
        if type(varargin[i]) == str:
            if varargin[i] == 'xmin':
                xmin = varargin[i + 1]
                i = i + 1
            elif varargin[i] == 'powerlaw':
                Type = 'PL'
                alpha = varargin[i + 1]
                i = i + 1
            elif varargin[i] == 'cutoff':
                Type = 'PC';
                alpha = varargin[i + 1]
                Lambda = varargin[i + 2]
                i = i + 2
            elif varargin[i] == 'exponential':
                Type = 'EX'
                Lambda = varargin[i + 1]
                i = i + 1
            elif varargin[i] == 'lognormal':
                Type = 'LN';
                mu = varargin[i + 1]
                sigma = varargin[i + 2]
                i = i + 2
            elif varargin[i] == 'stretched':
                Type = 'ST'
                Lambda = varargin[i + 1]
                beta = varargin[i + 2]
                i = i + 2
            else:
                argok = 0

        if not argok:
            print
            '(RANDHT) Ignoring invalid argument #', i + 1

        i = i + 1

    if n < 1:
        print
        '(RANDHT) Error: invalid ''n'' argument; using default.\n'
        n = 10000;

    if xmin < 1:
        print
        '(RANDHT) Error: invalid ''xmin'' argument; using default.\n'
        xmin = 1;

    x = []
    if Type == 'EX':
        x = []
        for i in range(n):
            x.append(xmin - (1. / Lambda) * log(1 - random()))
    elif Type == 'LN':
        y = []
        for i in range(10 * n):
            y.append(exp(mu + sigma * normalvariate(0, 1)))

        while True:
            y = filter(lambda X: X >= xmin, y)
            q = len(y) - n;
            if q == 0: break

            if q > 0:
                r = range(len(y));
                shuffle(r)
                ytemp = []
                for j in range(len(y)):
                    if j not in r[0:q]:
                        ytemp.append(y[j])
                y = ytemp
                break
            if (q < 0):
                for j in range(10 * n):
                    y.append(exp(mu + sigma * normalvariate(0, 1)))

        x = y

    elif Type == 'ST':
        x = []
        for i in range(n):
            x.append(pow(pow(xmin, beta) - (1. / Lambda) * log(1. - random()), (1. / beta)))
    elif Type == 'PC':

        x = []
        y = []
        for i in range(10 * n):
            y.append(xmin - (1. / Lambda) * log(1. - random()))
        while True:
            ytemp = []
            for i in range(10 * n):
                if random() < pow(y[i] / float(xmin), -alpha): ytemp.append(y[i])
            y = ytemp
            x = x + y
            q = len(x) - n
            if q == 0: break;

            if (q > 0):
                r = range(len(x))
                shuffle(r)

                xtemp = []
                for j in range(len(x)):
                    if j not in r[0:q]:
                        xtemp.append(x[j])
                x = xtemp
                break;

            if (q < 0):
                y = []
                for j in range(10 * n):
                    y.append(xmin - (1. / Lambda) * log(1. - random()))


    else:
        x = []
        for i in range(n):
            x.append(xmin * pow(1. - random(), -1. / (alpha - 1.)))

    return x


from math import *


# function [alpha, xmin, L]=plfit(x, varargin)
# PLFIT fits a power-law distributional model to data.
#    Source: http://www.santafe.edu/~aaronc/powerlaws/
#
#    PLFIT(x) estimates x_min and alpha according to the goodness-of-fit
#    based method described in Clauset, Shalizi, Newman (2007). x is a
#    vector of observations of some quantity to which we wish to fit the
#    power-law distribution p(x) ~  x^-alpha for x >= xmin.
#    PLFIT automatically detects whether x is composed of real or integer
#    values, and applies the appropriate method. For discrete data, if
#    min(x) > 1000, PLFIT uses the continuous approximation, which is
#    a reliable in this regime.
#
#    The fitting procedure works as follows:
#    1) For each possible choice of x_min, we estimate alpha via the
#       method of maximum likelihood, and calculate the Kolmogorov-Smirnov
#       goodness-of-fit statistic D.
#    2) We then select as our estimate of x_min, the value that gives the
#       minimum value D over all values of x_min.
#
#    Note that this procedure gives no estimate of the uncertainty of the
#    fitted parameters, nor of the validity of the fit.
#
#    Example:
#       x = [500,150,90,81,75,75,70,65,60,58,49,47,40]
#       [alpha, xmin, L] = plfit(x)
#   or  a = plfit(x)
#
#    The output 'alpha' is the maximum likelihood estimate of the scaling
#    exponent, 'xmin' is the estimate of the lower bound of the power-law
#    behavior, and L is the log-likelihood of the data x>=xmin under the
#    fitted power law.
#
#    For more information, try 'type plfit'
#
#    See also PLVAR, PLPVA

# Version 1.0.10 (2010 January)
# Copyright (C) 2008-2011 Aaron Clauset (Santa Fe Institute)

# Ported to Python by Joel Ornstein (2011 July)
# (joel_ornstein@hmc.edu)

# Distributed under GPL 2.0
# http://www.gnu.org/copyleft/gpl.html
# PLFIT comes with ABSOLUTELY NO WARRANTY
#
#
# The 'zeta' helper function is modified from the open-source library 'mpmath'
#   mpmath: a Python library for arbitrary-precision floating-point arithmetic
#   http://code.google.com/p/mpmath/
#   version 0.17 (February 2011) by Fredrik Johansson and others
#

# Notes:
#
# 1. In order to implement the integer-based methods in Matlab, the numeric
#    maximization of the log-likelihood function was used. This requires
#    that we specify the range of scaling parameters considered. We set
#    this range to be 1.50 to 3.50 at 0.01 intervals by default.
#    This range can be set by the user like so,
#
#       a = plfit(x,'range',[1.50,3.50,0.01])
#
# 2. PLFIT can be told to limit the range of values considered as estimates
#    for xmin in three ways. First, it can be instructed to sample these
#    possible values like so,
#
#       a = plfit(x,'sample',100)
#
#    which uses 100 uniformly distributed values on the sorted list of
#    unique values in the data set. Second, it can simply omit all
#    candidates above a hard limit, like so
#
#       a = plfit(x,'limit',3.4)
#
#    Finally, it can be forced to use a fixed value, like so
#
#       a = plfit(x,'xmin',3.4)
#
#    In the case of discrete data, it rounds the limit to the nearest
#    integer.
#
# 3. When the input sample size is small (e.g., < 100), the continuous
#    estimator is slightly biased (toward larger values of alpha). To
#    explicitly use an experimental finite-size correction, call PLFIT like
#    so
#
#       a = plfit(x,'finite')
#
#    which does a small-size correction to alpha.
#
# 4. For continuous data, PLFIT can return erroneously large estimates of
#    alpha when xmin is so large that the number of obs x >= xmin is very
#    small. To prevent this, we can truncate the search over xmin values
#    before the finite-size bias becomes significant by calling PLFIT as
#
#       a = plfit(x,'nosmall')
#
#    which skips values xmin with finite size bias > 0.1.

def plfit(x, *varargin):
    vec = []
    sample = []
    xminx = []
    limit = []
    finite = False
    nosmall = False
    nowarn = False

    # parse command-line parameters trap for bad input
    i = 0
    while i < len(varargin):
        argok = 1
        if type(varargin[i]) == str:
            if varargin[i] == 'range':
                Range = varargin[i + 1]
                if Range[1] > Range[0]:
                    argok = 0
                    vec = []
                try:
                    vec = map(lambda X: X * float(Range[2]) + Range[0], \
                              range(int((Range[1] - Range[0]) / Range[2])))


                except:
                    argok = 0
                    vec = []

                if Range[0] >= Range[1]:
                    argok = 0
                    vec = []
                    i -= 1

                i += 1


            elif varargin[i] == 'sample':
                sample = varargin[i + 1]
                i = i + 1
            elif varargin[i] == 'limit':
                limit = varargin[i + 1]
                i = i + 1
            elif varargin[i] == 'xmin':
                xminx = varargin[i + 1]
                i = i + 1
            elif varargin[i] == 'finite':
                finite = True
            elif varargin[i] == 'nowarn':
                nowarn = True
            elif varargin[i] == 'nosmall':
                nosmall = True
            else:
                argok = 0

        if not argok:
            print
            '(PLFIT) Ignoring invalid argument #', i + 1

        i = i + 1

    if vec != [] and (type(vec) != list or min(vec) <= 1):
        print
        '(PLFIT) Error: ''range'' argument must contain a vector or minimum <= 1. using default.\n'

        vec = []

    if sample != [] and sample < 2:
        print
        '(PLFIT) Error: ''sample'' argument must be a positive integer > 1. using default.\n'
        sample = []

    if limit != [] and limit < min(x):
        print
        '(PLFIT) Error: ''limit'' argument must be a positive value >= 1. using default.\n'
        limit = []

    if xminx != [] and xminx >= max(x):
        print
        '(PLFIT) Error: ''xmin'' argument must be a positive value < max(x). using default behavior.\n'
        xminx = []

    # select method (discrete or continuous) for fitting
    if reduce(lambda X, Y: X == True and floor(Y) == float(Y), x, True):
        f_dattype = 'INTS'
    elif reduce(lambda X, Y: X == True and (type(Y) == int or type(Y) == float ), x, True):
        f_dattype = 'REAL'
    else:
        f_dattype = 'UNKN'

    if f_dattype == 'INTS' and min(x) > 1000 and len(x) > 100:
        f_dattype = 'REAL'

    # estimate xmin and alpha, accordingly

    if f_dattype == 'REAL':
        xmins = unique(x)
        xmins.sort()
        xmins = xmins[0:-1]
        if xminx != []:
            xmins = [min(filter(lambda X: X >= xminx, xmins))]

        if limit != []:
            xmins = filter(lambda X: X <= limit, xmins)
            if xmins == []: xmins = [min(x)]

        if sample != []:
            step = float(len(xmins)) / (sample - 1)
            index_curr = 0
            new_xmins = []
            for i in range(0, sample):
                if round(index_curr) == len(xmins): index_curr -= 1
                new_xmins.append(xmins[int(round(index_curr))])
                index_curr += step
            xmins = unique(new_xmins)
            xmins.sort()

        dat = []
        z = sorted(x)

        for xm in range(0, len(xmins)):
            xmin = xmins[xm]
            z = filter(lambda X: X >= xmin, z)

            n = len(z)
            # estimate alpha using direct MLE

            a = float(n) / sum(map(lambda X: log(float(X) / xmin), z))
            if nosmall:
                if (a - 1) / sqrt(n) > 0.1 and dat != []:
                    xm = len(xmins) + 1
                    break

            # compute KS statistic
            # cx   = map(lambda X:float(X)/n,range(0,n))
            cf = map(lambda X: 1 - pow((float(xmin) / X), a), z)
            dat.append(max(map(lambda X: abs(cf[X] - float(X) / n), range(0, n))))
        D = min(dat)
        xmin = xmins[dat.index(D)]
        z = filter(lambda X: X >= xmin, x)
        z.sort()
        n = len(z)
        alpha = 1 + n / sum(map(lambda X: log(float(X) / xmin), z))
        if finite: alpha = alpha * float(n - 1) / n + 1. / n  # finite-size correction
        if n < 50 and not finite and not nowarn:
            print
            '(PLFIT) Warning: finite-size bias may be present.\n'

        L = n * log((alpha - 1) / xmin) - alpha * sum(map(lambda X: log(float(X) / xmin), z))
    elif f_dattype == 'INTS':

        x = map(int, x)
        if vec == []:
            for X in range(150, 351):
                vec.append(X / 100.)  # covers range of most practical
                # scaling parameters
        zvec = map(zeta, vec)

        print(x, unique(x))
        xmins = unique(x)
        xmins.sort()
        xmins = xmins[0:-1]
        if xminx != []:
            xmins = [min(filter(lambda X: X >= xminx, xmins))]

        if limit != []:
            limit = round(limit)
            xmins = filter(lambda X: X <= limit, xmins)
            if xmins == []: xmins = [min(x)]

        if sample != []:
            step = float(len(xmins)) / (sample - 1)
            index_curr = 0
            new_xmins = []
            for i in range(0, sample):
                if round(index_curr) == len(xmins): index_curr -= 1
                new_xmins.append(xmins[int(round(index_curr))])
                index_curr += step
            xmins = unique(new_xmins)
            xmins.sort()

        if xmins == []:
            print
            '(PLFIT) Error: x must contain at least two unique values.\n'
            alpha = 'Not a Number'
            xmin = x[0]
            D = 'Not a Number'
            return [alpha, xmin, D]

        xmax = max(x)

        z = x
        z.sort()
        datA = []
        datB = []

        for xm in range(0, len(xmins)):
            xmin = xmins[xm]
            z = filter(lambda X: X >= xmin, z)
            n = len(z)
            # estimate alpha via direct maximization of likelihood function

            # force iterative calculation
            L = []
            slogz = sum(map(log, z))
            xminvec = map(float, range(1, xmin))
            for k in range(0, len(vec)):
                L.append(-vec[k] * float(slogz) - float(n) * log(
                    float(zvec[k]) - sum(map(lambda X: pow(float(X), -vec[k]), xminvec))))

            I = L.index(max(L))
            # compute KS statistic
            fit = reduce(lambda X, Y: X + [Y + X[-1]], \
                         (map(lambda X: pow(X, -vec[I]) / (
                                     float(zvec[I]) - sum(map(lambda X: pow(X, -vec[I]), map(float, range(1, xmin))))),
                              range(xmin, xmax + 1))), [0])[1:]
            cdi = []
            for XM in range(xmin, xmax + 1):
                cdi.append(len(filter(lambda X: floor(X) <= XM, z)) / float(n))

            datA.append(max(map(lambda X: abs(fit[X] - cdi[X]), range(0, xmax - xmin + 1))))
            datB.append(vec[I])
        # select the index for the minimum value of D
        I = datA.index(min(datA))
        xmin = xmins[I]
        z = filter(lambda X: X >= xmin, x)
        n = len(z)
        alpha = datB[I]
        if finite: alpha = alpha * (n - 1.) / n + 1. / n  # finite-size correction
        if n < 50 and not finite and not nowarn:
            print
            '(PLFIT) Warning: finite-size bias may be present.\n'

        L = -alpha * sum(map(log, z)) - n * log(zvec[vec.index(max(filter(lambda X: X <= alpha, vec)))] - \
                                                sum(map(lambda X: pow(X, -alpha), range(1, xmin))))
    else:
        print
        '(PLFIT) Error: x must contain only reals or only integers.\n'
        alpha = []
        xmin = []
        L = []

    return alpha, xmin, L


# helper functions (unique and zeta)


def unique(seq):
    # not order preserving
    set = {}
    map(set.__setitem__, seq, [])
    return set.keys()


def _polyval(coeffs, x):
    p = coeffs[0]
    for c in coeffs[1:]:
        p = c + x * p
    return p


_zeta_int = [ \
    -0.5,
    0.0,
    1.6449340668482264365, 1.2020569031595942854, 1.0823232337111381915,
    1.0369277551433699263, 1.0173430619844491397, 1.0083492773819228268,
    1.0040773561979443394, 1.0020083928260822144, 1.0009945751278180853,
    1.0004941886041194646, 1.0002460865533080483, 1.0001227133475784891,
    1.0000612481350587048, 1.0000305882363070205, 1.0000152822594086519,
    1.0000076371976378998, 1.0000038172932649998, 1.0000019082127165539,
    1.0000009539620338728, 1.0000004769329867878, 1.0000002384505027277,
    1.0000001192199259653, 1.0000000596081890513, 1.0000000298035035147,
    1.0000000149015548284]

_zeta_P = [-3.50000000087575873, -0.701274355654678147,
           -0.0672313458590012612, -0.00398731457954257841,
           -0.000160948723019303141, -4.67633010038383371e-6,
           -1.02078104417700585e-7, -1.68030037095896287e-9,
           -1.85231868742346722e-11][::-1]

_zeta_Q = [1.00000000000000000, -0.936552848762465319,
           -0.0588835413263763741, -0.00441498861482948666,
           -0.000143416758067432622, -5.10691659585090782e-6,
           -9.58813053268913799e-8, -1.72963791443181972e-9,
           -1.83527919681474132e-11][::-1]

_zeta_1 = [3.03768838606128127e-10, -1.21924525236601262e-8,
           2.01201845887608893e-7, -1.53917240683468381e-6,
           -5.09890411005967954e-7, 0.000122464707271619326,
           -0.000905721539353130232, -0.00239315326074843037,
           0.084239750013159168, 0.418938517907442414, 0.500000001921884009]

_zeta_0 = [-3.46092485016748794e-10, -6.42610089468292485e-9,
           1.76409071536679773e-7, -1.47141263991560698e-6, -6.38880222546167613e-7,
           0.000122641099800668209, -0.000905894913516772796, -0.00239303348507992713,
           0.0842396947501199816, 0.418938533204660256, 0.500000000000000052]


def zeta(s):
    """
    Riemann zeta function, real argument
    """
    if not isinstance(s, (float, int)):
        try:
            s = float(s)
        except (ValueError, TypeError):
            try:
                s = complex(s)
                if not s.imag:
                    return complex(zeta(s.real))
            except (ValueError, TypeError):
                pass
            raise NotImplementedError
    if s == 1:
        raise ValueError("zeta(1) pole")
    if s >= 27:
        return 1.0 + 2.0 ** (-s) + 3.0 ** (-s)
    n = int(s)
    if n == s:
        if n >= 0:
            return _zeta_int[n]
        if not (n % 2):
            return 0.0
    if s <= 0.0:
        return 0
    if s <= 2.0:
        if s <= 1.0:
            return _polyval(_zeta_0, s) / (s - 1)
        return _polyval(_zeta_1, s) / (s - 1)
    z = _polyval(_zeta_P, s) / _polyval(_zeta_Q, s)
    return 1.0 + 2.0 ** (-s) + 3.0 ** (-s) + 4.0 ** (-s) * z


Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1
print(N_counties)

#choose the county by a number between [1,N_counties], (random choice or user setting)
c = np.random.choice(np.arange(1, N_counties + 1 ))
c = 987
county = Adata[1][c]

#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word = int(word.iloc[c, 1])
print(word)

#IMPORTING DATA
filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'('+str(c)+')_'+str(word)+"w.csv"
dat = pd.read_csv(filename)
print('dat',dat)

'''
alpha = 2
#alpha = 1 + 1/alpha
no_rand_zipf = 5000000
y = randht(no_rand_zipf,'powerlaw', alpha)
#y = y.astype(np.int64)
#y = [ int(i) for i in y ]
y = np.asarray(y, dtype=np.float)
print(y)

random_zipf = zipf.pmf(np.arange(1, word+1), a = alpha)
#random_zipf = -np.sort(-random_zipf)
#random_zipf = random_zipf/np.sum(random_zipf)
print(random_zipf)
#random_zipf = random_zipf/np.sum(random_zipf)
#print(random_zipf)

print('MAAAAAX', np.max(y))
binss = np.arange(1,np.max(y))
#y = np.histogram(y, bins = binss)
y = np.array(y)
y0 = np.array(y[0])
y1 = np.array(y[1])

#print(y, y0, y0[35], y1)
a = np.array(np.nonzero(y0))
print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', a[0])
a = a[0]


# GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

print(int(len(y0[a]))+1, a)
for i in range (0, int(len(a))) :
    print(y0[i])
y0 = -np.sort(-y0)
print('perché sti zeri?',y0[a[130]])
Y = y0[a]
Y = Y[np.nonzero(Y)[0]]
#Y = Y/np.sum(Y)


x_min = 1
n = int(len(Y))
#x = np.arange(0,n+1)
sum_ln_x_i = 0

for i in range (x_min, n) :
    sum_ln_x_i = sum_ln_x_i + np.int64(np.log(Y[i]))

print(np.log(1/(x_min - 0.5)), sum_ln_x_i, (sum_ln_x_i - np.log(1/(x_min - 0.5))  )**-1, n, n*(sum_ln_x_i - np.log(1/(x_min - 0.5))  )**-1)
est_alpha = 1 + n*(sum_ln_x_i - np.log(1/(x_min - 0.5))  )**-1
print('ESTIMATED ALPHA',est_alpha)




# GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
x = np.arange(1, int(len(Y))+1)

print(Y, len(Y), len(x))
ax1.scatter(x, Y, label = 'sampled data', color='red', s = 1)

x1 = np.arange(1, len(random_zipf) + 1)
print(len(x1), len(random_zipf))
ax1.scatter(x1,random_zipf, label="Zipf's sampled data", s=1, color = 'red')

ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
ax1.set_xticks(np.arange(1, int(len(Y))+1), minor=False)

a = 20
b = int(len(Y))

x = np.arange(a,b)
print(len(x), len(Y[a:b]))
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(np.log10(x), np.log10(Y[a:b] ))
print('scipy sampled data, linear regression',slope1, intercept1)

x3 = np.arange(1,  int(len(Y))+1)
y3 = x3**slope1 * 10**intercept1
ax1.plot(x3,y3 , linewidth=1, color='purple', label = "LR on Zipf's sampled data")

Y = np.asarray(Y,dtype=float)


#Setting log-log scale

ax1.set_yscale('log')
ax1.set_xscale('log')
'''
X = np.arange(1, 10000)
#alpha, xmin, L
a = plfit(X)

#plt.show()