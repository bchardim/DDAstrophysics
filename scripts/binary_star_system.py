#!/usr/bin/python -Wignore


# Force true division
from __future__ import division
import sys
import sympy
import unicodedata
import numpy as np
from math import pi
from astropy import constants as const

def harmonic(x,n):
    func = 0
    z = x * pi
    amp = 2/(n + 1)
    for i in xrange(1,n+1,2) :
        sign = (-1)**((i-1)/2)
        func = func + amp*sign*np.sin(i * z)
    return func



if __name__ == "__main__":

    # Define unit conversion
    day_to_sec = 86400
    degree_to_radian = 0.01745329252
    meter_to_au = 6.68459e-12


    # Ask for input data
    print "######################################"
    print "# Output file path to save results"
    print "######################################"
    outfile = raw_input("Path for output file: ")
    print "######################################"
    print "# Enter the data for Star #1"
    print "######################################"
    m1 = float(input("Mass (solar masses): "))
    R1 = float(input("Radius (solar radii): "))
    T1 = float(input("Effective Temperature (K): "))
    print "######################################"
    print "# Enter the data for Star #2"
    print "######################################"
    m2 = float(input("Mass (solar masses): "))
    R2 = float(input("Radius (solar radii): "))
    T2 = float(input("Effective Temperature (K): "))
    print "######################################"
    print "# Enter the desired orbital parameters"
    print "######################################"
    P = float(input("Orbital Period (days): "))
    e = float(input("Orbital Eccentricity: "))
    i = float(input("Orbital Inclination (deg): "))
    phi = float(input("Orientation of Periastron (deg): "))
    print "######################################"
    print "# Enter center of mass velocity vector"
    print "######################################"
    vcmxp = float(input("v_x' (km/s): "))
    vcmyp = float(input("v_y' (km/s): "))
    vcmzp = float(input("v_z' (km/s): "))

    # Convert values to conventional SI units and radians
    m1  = m1*const.M_sun.value
    m2  = m2*const.M_sun.value
    G   = const.G.value
    M   = m1 + m2 # Total mass of system
    R1  = R1*const.R_sun.value
    R2  = R2*const.R_sun.value
    P   = P*day_to_sec
    i   = i*degree_to_radian
    phi = phi*degree_to_radian
    vcmxp = vcmxp*1000
    vcmyp = vcmyp*1000
    vcmzp = vcmzp*1000
 
    # Compute the semimajor axes of the orbits and reduced mass
    mu = m1*m2/M                            # Reduced mass
    a = (((P**2)*G*M)/(4*pi**2))**(1/3)     # Kepler's 3rd Law
    a1 = (mu/m1)*a                          # Semimajor axis for m1
    a2 = (mu/m2)*a                          # Semimajor axis for m2
    print ""
    print "-----------------------------------------"
    print "- Semimajor axis results:                "
    print "-----------------------------------------"
    print "The semimajor axis of the reduced mass is " + str(a*meter_to_au) + " AU"
    print "a1 =  " + str(a1*meter_to_au) + " AU"
    print "a2 =  " + str(a2*meter_to_au) + " AU"

    # Check that stars are not in contact
    if (a < R1 + R2):
        print ""
        print "!!! Your two stars are in contact!!!"
        print "R1 + R2 :"  + str((R1+R2)*meter_to_au) + " AU > a"
        print "Ending calculation"
        sys.exit(0)


    # Compute the luminosity of each star and the total luminosity
    L1 = (4*pi)*(R1**2)*(const.sigma_sb.value)*(T1)**4    # Stefan-Boltzmann law for m1
    L2 = (4*pi)*(R2**2)*(const.sigma_sb.value)*(T2)**4    # Smame S-B law for m2
    L  = L1 + L2                                          # Total uneclipsed luminosity

