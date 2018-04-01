#!/usr/bin/python

from pysextant import *
from math import pi

# sun pos from https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.htmlf
# sun elevation at feb 12 2017 13:20:00 in seattle
theta1 = 27.73*pi/180
posix1 = 1486934400
# sun elevation at feb 12 2017 10:20:00 in seattle 
theta2 = 23.06*pi/180
posix2 = 1486923600

(ll1, ll2) = est_lla(theta1, theta2, posix1, posix2)
print('In Seattle, sun at angle %5.2f deg at POSIX time %i' % (theta1*180/pi, posix1))
print('In Seattle, sun at angle %5.2f deg at POSIX time %i' % (theta2*180/pi, posix2))
print('Seatle     lat = % 6.2f, lon = % 6.2f' % (47.0+36.0/60, -122.0-19.0/60))
print('Solution 1 lat = % 6.2f, lon = % 6.2f' % (ll1[0,0]*180/pi, ll1[1,0]*180/pi))
print('Solution 2 lat = % 6.2f, lon = % 6.2f\n' % (ll2[0,0]*180/pi, ll2[1,0]*180/pi))

# sun elevation at feb 12 2017 13:20:00 in beijing
theta1 = 35.16*pi/180
posix1 = 1486876800
# sun elevation at feb 12 2017 10:20:00 in beijing 
theta2 = 28.67*pi/180
posix2 = 1486866000

(ll1, ll2) = est_lla(theta1, theta2, posix1, posix2)
print('In Beijing, sun at angle %5.2f deg at POSIX time %i' % (theta1*180/pi, posix1))
print('In Beijing, sun at angle %5.2f deg at POSIX time %i' % (theta2*180/pi, posix2))
print('Beijing    lat = % 6.2f, lon = % 6.2f' % (39+55.0/60, +116.0+25.0/60))
print('Solution 1 lat = % 6.2f, lon = % 6.2f' % (ll1[0,0]*180/pi, ll1[1,0]*180/pi))
print('Solution 2 lat = % 6.2f, lon = % 6.2f\n' % (ll2[0,0]*180/pi, ll2[1,0]*180/pi))

# sun elevation at mar 31 2018 17:00:00 in cape town
theta1 = 20.06*pi/180
posix1 = 1522508400
# sun elevation at mar 31 2018 12:00:00 in cape town
theta2 = 50.05*pi/180
posix2 = 1522490400

(ll1, ll2) = est_lla(theta1, theta2, posix1, posix2)
print('In Cape Town, sun at angle %5.2f deg at POSIX time %i' % (theta1*180/pi, posix1))
print('In Cape Town, sun at angle %5.2f deg at POSIX time %i' % (theta2*180/pi, posix2))
print('Cape Town  lat = % 6.2f, lon = % 6.2f' % (-33-92.49/60, 18+42.41/60))
print('Solution 1 lat = % 6.2f, lon = % 6.2f' % (ll1[0,0]*180/pi, ll1[1,0]*180/pi))
print('Solution 2 lat = % 6.2f, lon = % 6.2f\n' % (ll2[0,0]*180/pi, ll2[1,0]*180/pi))