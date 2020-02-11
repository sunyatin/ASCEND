#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: rtournebize
"""

# Usage:
# python3 LRT_1exp_2exp_lomax.py file_with_two_columns_X_Y

from scipy.optimize import curve_fit, OptimizeWarning
import numpy as np
import sys
import warnings
from scipy.stats import chi2, f
import matplotlib.pyplot as plt

if len(sys.argv)-1 != 1:
    sys.exit('You should provide 1 argument.')

FILE = sys.argv[1]
outputPrefix = FILE

print(FILE)

###############################################################################

XY = np.genfromtxt(FILE, dtype=float, usecols = (0, 1))
x = XY[:,0]
y = XY[:,1]

np.seterr(all='ignore')
warnings.simplefilter("error", OptimizeWarning)

###############################################################################
###############################################################################
###############################################################################

def residuals_jacquelin_exp1d(x, y):

    n = len(x)

    S = [0]
    for k in np.arange(1, n):
        S += [ S[-1] + 1/2 * (y[k]+y[k-1])*(x[k]-x[k-1]) ]
              
    S = np.asarray(S)
    x = np.asarray(x)
    y = np.asarray(y)

    M1 = [
        [
            sum( (x-x[0])**2 ),
            sum( (x-x[0])*S )
            
        ],
        [
            sum( (x-x[0])*S ),
            sum( S**2 )
        ]]
        
    M2 = [
          sum( (y-y[0])*(x-x[0]) ),
          sum( (y-y[0])*S )
          ]
        
    M1 = np.array(M1)
    M2 = np.array(M2)
    
    K = np.dot(np.linalg.inv(M1), M2)
    c = K[1]
    
    N1 = [
        [
            n,
            sum( np.exp(c*x) )
            
        ],
        [
            sum( np.exp(c*x) ),
            sum( np.exp(2*c*x) )
        ]]
        
    N2 = [
          sum( y ),
          sum( y*np.exp(c*x) )
          ]
  
    AB = np.dot(np.linalg.inv(N1), N2)
          
    a = AB[0]
    b = AB[1]
  
    def exp(x, a, b, c):
        return a + b*np.exp(c*x)
        
    if np.isnan(a):
        return None

    try:
        popt, pcov = curve_fit(exp, x, y, maxfev = 5000, p0 = [a, b, c])
        yfit = exp(x, popt[0], popt[1], popt[2])
        ssq = ((yfit-y)**2).sum()
        return [[popt[2]*-0.5*100], yfit, [ssq]]
    except:
        return None

###############################################################################
def residuals_jacquelin_exp2d(x, y):

    n = len(x)

    S = [0]
    for k in np.arange(1, n):
        S += [ S[-1] + 1/2 * (y[k]+y[k-1])*(x[k]-x[k-1]) ]

    SS = [0]
    for k in np.arange(1, n):
        SS += [ SS[-1] + 1/2 * (S[k]+S[k-1])*(x[k]-x[k-1]) ]
    
    S = np.asarray(S)
    SS = np.asarray(SS)
    x = np.asarray(x)
    y = np.asarray(y)

    M1 = [
        [
            sum(SS**2),
            sum(SS*S),
            sum(SS*x**2),
            sum(SS*x),
            sum(SS)
        ],
        [
            sum(SS*S),
            sum(S**2),
            sum(S*x**2),
            sum(S*x),
            sum(S)
        ],
        [
            sum(SS*x**2),
            sum(S*x**2),
            sum(x**4),
            sum(x**3),
            sum(x**2)
        ],
        [
            sum(SS*x),
            sum(S*x),
            sum(x**3),
            sum(x**2),
            sum(x)
        ],
        [
            sum(SS),
            sum(S),
            sum(x**2),
            sum(x),
            n
        ]]
        
    M1 = np.array(M1)
    
    M2 = [
        sum(SS*y),
        sum(S*y),
        sum(x**2*y),
        sum(x*y),
        sum(y)
    ]
    
    M2 = np.array(M2)
    
    K = np.dot(np.linalg.inv(M1), M2)
    
    p = 1/2 * (K[1] + np.sqrt(K[1]**2 + 4*K[0]))
    q = 1/2 * (K[1] - np.sqrt(K[1]**2 + 4*K[0]))
    
    N1 = [
        [
            n,
            sum(np.exp(p*x)),
            sum(np.exp(q*x))
        ],
        [
            sum(np.exp(p*x)),
            sum(np.exp(2*p*x)),
            sum(np.exp((p+q)*x))
        ],
        [
            sum(np.exp(q*x)),
            sum(np.exp((p+q)*x)),
            sum(np.exp(2*q*x))
        ]
    ]
    
    N1 = np.array(N1)
    
    N2 = np.array([sum(y), sum(y*np.exp(p*x)), sum(y*np.exp(q*x))])
    
    H = np.dot(np.linalg.inv(N1), N2)
    
    a = H[0]
    b = H[1]
    c = H[2]

    def exp2(x, a, b, c, p, q):
        return a + b*np.exp(p*x) + c*np.exp(q*x)
        
    if np.isnan(a):
        return None
    
    try:
        popt, pcov = curve_fit(exp2, x, y, maxfev = 5000, p0 = [a, b, c, p, q])
        yfit = exp2(x, popt[0], popt[1], popt[2], popt[3], popt[4])
        ssq = ((yfit-y)**2).sum()
        return [[popt[3]*-0.5*100, popt[4]*-0.5*100], yfit, [ssq]]
    except:
        return None

###############################################################################
def residuals_lomax(x, y):
    
    def lomax(x, alpha, lamb):
        return alpha*lamb**alpha / ( (x+lamb)**(alpha+1)  )

    try:
        popt, pcov = curve_fit(lomax, x, y, maxfev = 5000)
        yfit = lomax(x, popt[0], popt[1])
        ssq = ((yfit-y)**2).sum()
        return [[popt[1]*0.5*100], yfit, [ssq]]
    except:
        return None


###############################################################################

if len(x) != len(y):
    sys.exit('the x and y vectors must have same length')
    
x = np.array(x)
y = np.array(y)
N = len(x)

plt.figure(figsize=(15, 10))
plt.scatter(x, y, color='black')

np_exp1 = 3
np_exp2 = 5
np_lomax = 2

ssq_exp1 = residuals_jacquelin_exp1d(x, y)
ssq_exp2 = residuals_jacquelin_exp2d(x, y)
ssq_lomax = residuals_lomax(x, y)

dont_do = False

output = open(outputPrefix+'.log', 'w')

if ssq_exp1 != None:
    T_1 = int(ssq_exp1[0][0])
    plt.plot(x, ssq_exp1[1], linestyle='-', color='orange', label='1-term exponential | Tf='+str(T_1))
    ssq_exp1 = ssq_exp1[2][0]
    output.write('1-term exponential, Sum Squared Residuals: '+str(ssq_exp1)+'\n')
    logL_exp1 = -N/2 * ( np.log(2*np.pi) + 1 - np.log(N) + np.log(ssq_exp1) )
else:
    dont_do = True
    
if ssq_exp2 != None:
    T_2 = [int(x) for x in ssq_exp2[0]]
    plt.plot(x, ssq_exp2[1], linestyle='--', color='red', label='2-term exponential | Tf='+', '.join([str(x) for x in T_2]))
    ssq_exp2 = ssq_exp2[2][0]
    output.write('2-term exponential, Sum Squared Residuals: '+str(ssq_exp2)+'\n')
    logL_exp2 = -N/2 * ( np.log(2*np.pi) + 1 - np.log(N) + np.log(ssq_exp2) )
else:
    dont_do = True

if ssq_lomax != None:
    T_lomax = int(ssq_lomax[0][0])
    plt.plot(x, ssq_lomax[1], linestyle=':', color='blue', label='Lomax function | Tf='+str(T_lomax))
    ssq_lomax = ssq_lomax[2][0]
    output.write('Lomax distribution, Sum Squared Residuals: '+str(ssq_lomax)+'\n')
else:
    dont_do = True
    
output.write('\n')
    
# I think this is not motivated
#logL_lomax = -N/2 * ( np.log(2*np.pi) + 1 - np.log(N) + np.log(ssq_lomax) )

# vs lomax
#LRT = -2*logL_exp1 + 2*logL_lomax
#p  = 1 - chi2.cdf(LRT, np_lomax - np_exp1)

# vs 2-exp
    
if dont_do == False:

    LRT = -2*logL_exp1 + 2*logL_exp2
    p  = 1 - chi2.cdf(LRT, np_exp2 - np_exp1)
    
    # F-ratio
    # np2 must be greater than np1
    # ie model1 is restricted and model2 is unrestricted nested model
    
    F_ratio = ((ssq_exp1 - ssq_exp2) / (np_exp2 - np_exp1)) / ((ssq_exp2 / (N - np_exp2)))
    p_F = 1 - f.cdf(F_ratio, np_exp2 - np_exp1, N - np_exp2)
    
    plt.suptitle('LRT: '+str(p)+'\nF-ratio: '+str(p_F))

    output.write('log-likelihood 1-term exponential: '+str(logL_exp1)+' | Tf='+str(T_1)+'\n')
    output.write('log-likelihood 2-term exponential: '+str(logL_exp2)+'| Tf='+', '.join([str(x) for x in T_2])+'\n')
    output.write('\n')
    output.write('LRT: '+str(LRT)+' | P-value: '+str(p)+'\n')
    output.write('F-ratio: '+str(F_ratio)+' | P-value: '+str(p_F)+' | Tf='+str(T_lomax)+'\n')
    output.close()

# plots
plt.legend()
#plt.show()
plt.savefig(outputPrefix+'.png')




