#!/usr/bin/python -Wignore


# Force true division
from __future__ import division
import sys
import sympy
import unicodedata
import numpy as np
from math import pi
from astropy import constants as const
import matplotlib.pyplot as plt

def Flux(r_prime,R,T):
    #
    # This function calculates the flux for a specific value
    # of the distance from the center of the stellar disk, r_prime.
    #
    # Van Hamme model (solar-like)
    x = 0.648
    y = 0.207
    mu = np.cos(r_prime/(R*(pi/2)))
    Limb_Darkening = 1 - x*(1 - mu) - y*mu*(np.log(mu))
    flux = (const.sigma_sb.value)*(T**4)*(Limb_Darkening)   #Eq. (3.18)
    return flux



def Transformation (x, y, i):
    #
    # This function performs the coordinate transformation between
    # the orbital plane coordinates (x,y) and the plane of the sky
    # coordinates (xp, yp, zp), based on the angle of inclination, i.
    #

    xprime =  x*np.sin(i)          # Eq. (J.1)
    yprime =  y                    # Eq. (J.2)
    zprime = -x*np.cos(i)          # Eq. (J.3)
    return xprime, yprime, zprime


def Eclipse(x1, y1, R1, T1, x2, y2, R2, T2, i):
    #
    #  This funtion computes the change in observed luminosity due to
    #  an eclipse.
    #
    
    # Number of steps for r', theta' integrations
    Nr = 100
    Ntheta = 500
    dtheta_prime = pi/Ntheta
    dS = 0

    # Perform coordinate transformations
    x1p, y1p, z1p = Transformation(x1,y1,i)
    x2p, y2p, z2p = Transformation(x2,y2,i)

    # Determine which star is in front (f) and which star is in back (b)
    if x1p > 0:
        xfp, yfp, zfp, Rf, Tf = x1p, y1p, z1p, R1, T1
        xbp, ybp, zbp, Rb, Tb = x2p, y2p, z2p, R2, T2
    else:
        xfp, yfp, zfp, Rf, Tf = x2p, y2p, z2p, R2, T2
        xbp, ybp, zbp, Rb, Tb = x1p, y1p, z1p, R1, T1

    # Are the two stars close enough for an eclipse?
    d = ((yfp - ybp)**2 + (zfp - zbp)**2)**(1/2)           # Eq. (J.4)
    if (d <= Rf + Rb):
        # Find the angle between y' and the projected line between the
        # centers of the stars in the (y',z') plane.  The polar coordinate
        # integration will be centered on this line to take advantage
        # of spherical symmetry.
        theta0_prime = np.arctan2((zfp - zbp), (yfp - ybp))  # Eq. (J.5)
   
        # Determine the starting radius for the integration
        if (d < Rb - Rf):
            r_prime = d + Rf    #Foreground star disk entirely
            r_stop  = d - Rf    #inside background star disk
            if (r_stop < 0): 
                r_stop = 0
        else:
            r_prime = Rb
            r_stop  = 0
        dr_prime = r_prime/Nr

        # The surface integration loop
        while True:
            # Determine the limits of the angular integration for the current r_prime
            theta_prime = theta0_prime

            while True:
                yp_dA = r_prime*np.cos(theta_prime + dtheta_prime) + ybp
                zp_dA = r_prime*np.sin(theta_prime + dtheta_prime) + zbp
                if (((yp_dA - yfp)**2 + (zp_dA - zfp)**2)**(1/2) > Rf):
                    break
                theta_prime = theta_prime + dtheta_prime
                if ((theta_prime - theta0_prime) > pi):
                    break
            
            # Add the luminosity change for differential area  (Eq. J.7)
            dS = dS + 2*Flux(r_prime - dr_prime/2, Rb, Tb)*(r_prime - dr_prime/2)*dr_prime*(theta_prime - theta0_prime)

            # Check to see that there is no remaining overlap or if center of disk has been reached
            r_prime = r_prime - dr_prime
            if (r_prime < r_stop):
                break

    # Return calculated values
    return dS, y1p, z1p, y2p, z2p


if __name__ == "__main__":

    # Define unit conversion and constants
    day_to_sec = 86400
    degree_to_radian = 0.01745329252
    meter_to_au = 6.68459e-12
    Mbol_Sun = 4.74
    Mbol_max = -99999

    # Prepare plot
    plt.figure()

    # Ask for input data
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

    # Determine the energy produced by the uneclipsed, projected disks
    Nr = 100 # Surface rings steps
    drS1 = R1/Nr
    drS2 = R2/Nr
    rS1  = 0
    rS2  = 0
    S1   = 0
    S2   = 0
    for index in range(0,Nr+1):    # Numerical integration loop (Eq. J.7)
        rS1 = rS1 + drS1
        rS2 = rS2 + drS2
        S1  = S1 + 2*pi*Flux(rS1 - drS1/2, R1, T1)*(rS1 - drS1/2)*drS1
        S2  = S2 + 2*pi*Flux(rS2 - drS2/2, R2, T2)*(rS2 - drS2/2)*drS2
    S = S1 + S2

    # Initial orbit loop conditions
    N  = 1000       # Time steps per orbit
    t      = 0
    theta  = 0
    dt     = P/N                             # time step
    L_ang  = mu*((G*M*a*(1 - e**2)))**(1/2)  # L, Eq. (2.30)
    dAdt   = L_ang/(2*mu)                    # 2nd Law, Eq. (2.32)
  

    # Reduced mass orbit loop
    print ""
    print "------------------------------------------------"
    print "- Print results with radial velocities in km/s  "
    print "------------------------------------------------"
    print ""
    print "t/P     v1r(km/s)     v2r(km/s)     Mbol    dS(W)"
    print ""
    while (t < P):
        r   =  a*(1 - (e**2))/(1 + e*np.cos(theta))  # position, Eq. (2.3)
        v   =  (G*M*(2/r - 1/a))**(1/2)            # velocity, Eq. (2.36)
        vr  = -v*np.sin(i)*np.sin(theta + phi)     # radial velocity
        v1r =  (mu/m1)*vr                          # Eq. (2.23)
        v2r = -(mu/m2)*vr                          # Eq. (2.24)

        # Determine (x,y) positions of centers of stars
        x   = r*np.cos(theta + phi)                # reduced mass
        y   = r*np.sin(theta + phi)                # reduced mass
        x1  =  (mu/m1)*x
        y1  =  (mu/m1)*y
        x2  = -(mu/m2)*x
        y2  = -(mu/m2)*y

        # Calculate Mbol of the system (including eclipse effects)
        dS, y1p, z1p, y2p, z2p = Eclipse(x1, y1, R1, T1, x2, y2, R2, T2, i)
        Lt = L*(1 - dS/S)
        Mbol = Mbol_Sun - (5/2)*np.log10(Lt/const.L_sun.value)   # Mbol, Eq. (3.8)
        if (Mbol > Mbol_max):
            Mbol_max = Mbol
            t_max = t

        # Print results
        print ("%s, %s, %s, %s, %s"  %(t/P, (v1r + vcmxp)/1000, (v2r + vcmxp)/1000, Mbol, Lt*dS/S)) 
        dtheta = (2*dAdt/r**2)*dt                   # Eq. (2.31)
        theta  = theta + dtheta
        t      = t + dt
   
        # Fill light curve plot
        plt.plot(t/P, Mbol, 'bo')

        # Finish loop
        if (t > P):
            break

    # Plot light curve
    plt.title('light curve')
    plt.ylim(plt.ylim()[::-1])
    plt.xlabel('t/P')
    plt.ylabel('Mbol');
    plt.show()

