#!/usr/bin/python -Wignore


# Example 'harmonic_wave_function.py 5' 

import sys
import sympy
import unicodedata
import numpy as np
from math import pi
from matplotlib import pyplot as plt

def harmonic(x,n):
    func = 0
    z = x * pi
    amp = np.true_divide(2,(n + 1))
    for i in xrange(1,n+1,2) :
        sign = (-1)**np.true_divide((i-1),2)
        func = func + amp*sign*np.sin(i * z)
    return func



if __name__ == "__main__":

    x = np.linspace(0,1,500) # pi multiples [0,1]*pi
    plt.plot(x, harmonic(x,int(sys.argv[1])))
    plt.xlabel(unicodedata.lookup("GREEK SMALL LETTER PI"))
    plt.ylabel(unicodedata.lookup("GREEK SMALL LETTER PSI"),rotation=0)
    plt.axhline(y=0, color='grey', linestyle='dashed')
    plt.show()
