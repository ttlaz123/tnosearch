import numpy as np

def sphdist(long1, lat1, long2, lat2):
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)

    dphi = abs(phi1-phi2)
    dl = np.radians(abs(lat1-lat2))

    a = np.sin(dphi/2) **2 
    b = np.sin(dl/2) **2
    c = np.cos(phi1)*np.cos(phi2)

    x = a + b*c
    d = 2*np.arctan2(np.sqrt(x), np.sqrt(1-a))
    return d
 
def eudist(long1, lat1, long2, lat2):
    a = np.radians(long1-long2)**2
    b = np.radians(lat1-lat2)**2
    return np.sqrt(a+b)

def prodist1(long1, lat1, long2, lat2):
    lat = np.radians((lat1+lat2)/2.)
    long1 = long1*np.cos(lat)
    long2 = long2*np.cos(lat)
    a = np.radians(long1-long2)**2
    b = np.radians(lat1-lat2)**2
    return np.sqrt(a+b)

def prodist2(long1, lat1, long2, lat2):
    lat = np.radians((lat1+lat2)/2.)
    long1 = long1*np.sqrt(np.cos(lat))
    long2 = long2*np.sqrt(np.cos(lat))
    a = np.radians(long1-long2)**2
    b = np.radians(lat1-lat2)**2
    return np.sqrt(a+b)


