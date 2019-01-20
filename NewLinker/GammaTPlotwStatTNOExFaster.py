#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:41:05 2017

The library used to determine the size and angles of the cone for each detection
"""
import time
import numpy as np
from math import atan2,cos,sin,sqrt, radians, degrees, asin, cos, pi
import math
from datetime import datetime
from astropy.time import Time

d = 35 # d is the hypothetical distance from earth to TNO at the starting
        # observation, it is set to such that the TNO has the greatest apparent
        # motion in the sky(to exaggerate the motion, we could further reduce the distance)
tau = 2*pi
auToM = 1.5 * 10**11
G = 6.67 * 10**(-11)
M = 2 * 10**30
def maxVel(d):
    return math.sqrt(2*G * M / (d * auToM))
def maxArcsecPerDay(d):
    return maxVel(d) / (d * auToM) * 24 * 3600 / math.pi * 180 * 3600
def maxDegPerDay(d):
    return maxVel(d) / (d * auToM) * 24 * 3600 / math.pi * 180
## distance in AU
## time in Days
def apparentMotion(distance, time):
    return maxArcsecPerDay(distance)*time

def predict(gammaS, betaS, deltaT):
    """
    Para: Gamma and beta at starting date, and time difference(in days)
    Return: New Gamma and New Beta(could be new lambda and new beta,
    see ******************************)
    Note: the Lambda Earth is approximated using the period of Earth, 
    the distance is in units of AU, the coordinates system is built with
    sun as the center, earth starts at (0,-1) and ecliptic plane as the 
    X-Y plane
    Do not omit the global variables above!
    """
    
    #Lambda Earth(in fact, the difference in lambda earth) 
    #is "approximated" using the Earth Period
    lambdaE = (deltaT/365.25)*2*pi
    
    #result[0] = new gamma, result[1] = new beta
    result = np.array([0.0,0.0])
    
    dp = d*cos(betaS)
    Xt = dp*sin(gammaS)
    Yt = -1-dp*cos(gammaS)
    Zt = d*sin(betaS)
    Xe = sin(lambdaE)
    Ye = -cos(lambdaE)
    #Ze = 0
    
    # calculation for new Gamma
    # vector(Xt-Xe,Yt-Ye)
    gammaN = atan2(Yt-Ye, Xt-Xe) - lambdaE + pi/2
    #gammaN = atan2(Yt-Ye, Xt-Xe) + pi/2 *************************For beta against lambda just change 
    #the line above to this code**************************
    if gammaN < 0:
        gammaN = gammaN + 4*pi
    gammaN = gammaN %(2*pi)
    
    result[0] = gammaN
    
    #calculation for new Beta
    proj = sqrt((Xt-Xe)**2 + (Yt-Ye)**2)
    
    result[1] = atan2(Zt, proj)
    
    return result

#find the arc of the circle that the paths subtend
def calcCone(lookAhead, ra, dec, nite):
    #convert to gamma beta
    (betaS, gammaS) = equatorialToEcliptic(ra, dec, nite)
    #find new gamma new beta by motion of earth
    angle1 = []
    angle2 = []
    radius = []
    #calculate the cone that covers possible motions at the lookahead nite
    # and at specific intervals
    # all other possible detections at other nights should fall in between
    for y in range(lookAhead):
        x = y + 1
        gammaN, betaN = predict(np.radians(gammaS), np.radians(betaS), x)
        (raN, decN) = eclipticToEquatorial(np.degrees(betaN), 
                        np.degrees(gammaN), nite+x)
        #calculate the cone
        if raN > 180:
            raN = raN -360
        deltaRA = raN - ra
        deltaDEC = decN - dec
        #distance traveled by TNO
        distance = d
        tnoTrav = apparentMotion(distance, x) / 3600
        totDist = math.sqrt(deltaRA**2 + deltaDEC**2) + tnoTrav
        
        angle = np.arctan2(deltaDEC, deltaRA)
        #initialize the cone to be the entire circle
        deltaAng = np.pi - 0.01
        if(totDist-tnoTrav > tnoTrav):
            deltaAng = np.arcsin(tnoTrav/(totDist-tnoTrav))
            deltaAng = deltaAng/np.cos(np.deg2rad(dec/2))
        #between 0 and 2pi
        # a bunch of math that may or may not make sense
        angle1temp = ((angle+deltaAng)+2*pi) %(2*pi)  
        angle2temp = ((angle-deltaAng)+2*pi) %(2*pi) - abs(np.sin(np.deg2rad(dec/3))) 
        angle1.append(angle1temp)
        angle2.append(angle2temp)
        radius.append(totDist)

    return np.array(angle1), np.array(angle2), np.array(radius)

#converts from beta gamma to ra and dec
def eclipticToEquatorial(beta, gamma, MJD):  
    beta = radians(beta)  
    gamma = radians(gamma)        

    tilt = radians(23.43703)     
    # uses MJD = 56192.447049 as earth longitude = 0  
    earthPos =  ((MJD - 56192.447049) / 365.25 * tau) % tau  
    longitude = (earthPos + gamma) % tau      
    dec = asin(sin(beta) * cos(tilt) + cos(beta) * sin(tilt) * sin(longitude)) 

    ra = atan2(cos(beta) * sin(longitude) * cos (tilt) - sin(beta) * sin(tilt), 
               cos(beta) * cos(longitude))     
    return (degrees((ra + tau) % tau), degrees(dec))   

#converts from ra and dec to beta and gamma
def equatorialToEcliptic(ra, dec, MJD):       
    tilt = radians(23.43703)    
    ra = radians(ra)    
    dec = radians(dec)              
    # takes in ra and dec in degrees  
    def Beta(ra, dec):  
        return asin(sin(dec) * cos(tilt) - cos(dec) * sin(tilt) * sin(ra))  
    # get beta of detection  
    beta = Beta(ra, dec)     
    # longitude of detection (geocentric)   
    longitude = degrees(atan2((cos(dec)*sin(ra)*cos(tilt)+    
                                   sin(dec)*sin(tilt))/cos(beta),    
                                    cos(dec)*cos(ra)/cos(beta)))                
    longitude = (longitude + 360) % 360  
    # finding gamma            
    # uses MJD = 56192.447049 as earth longitude = 0    
    earthPos =  ((MJD - 56192.447049) / 365.25 * 360) % 360     
    gamma = (longitude - earthPos + 360) % 360   
    return (degrees(beta), gamma)   

#takes date in float form MJD and converts to string form YYYY/MM/DDThh:mm:ss.ssssss
def mjdToDate(mjd):
    dt = Time(mjd, format='mjd').datetime
    str = dt.isoformat()
    str = str.replace('-','/')
    return str
