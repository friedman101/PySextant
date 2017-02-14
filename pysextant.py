
import numpy as np
from math import *

def _posix2julian(posix):
	# posix time to julian date
	jd = posix/86400.0+2440587.5
	return jd

def _jd2T(jd):
	# centuries from J2000
	T = (jd - 2451545)/36525
	return T

def _R3(phi):
	R3 = np.matrix([ \
		[cos(phi), sin(phi), 0], \
		[-sin(phi), cos(phi), 0], \
		[0, 0, 1]\
		])
	return R3

def _sun_eci(posix):
	# source: http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4140/attde.pdf
	jd = _posix2julian(posix)
	t_ut1 = _jd2T(jd)
	lambda_sun = 280.4606184+36000.77005361*t_ut1
	m_sun = 357.5277233 + 35999.05034*t_ut1
	lambda_ecliptic = lambda_sun +  1.914666471+sin(m_sun*pi/180) + \
	0.918994643*sin(2*m_sun*pi/180)
	eps = 23.439291 - 0.0130042*2*t_ut1
	sun_eci = np.matrix([[cos(lambda_ecliptic*pi/180)], \
		[cos(eps*pi/180) * sin(lambda_ecliptic*pi/180)], \
		[sin(eps*pi/180) * sin(lambda_ecliptic*pi/180)]])
	return sun_eci

def _eci2ecef(posix):
	# source: http://aa.usno.navy.mil/publications/docs/Circular_179.pdf
	jd_utc = _posix2julian(posix)
	jd_tt = _posix2julian(posix + 32.184)
	DU = jd_utc - 2451545.0
	theta = 0.7790572732640 + 1.00273781191135448*DU
	T = _jd2T(jd_tt)
	GMST =  86400*theta+ (0.014506 + 4612.156534*T+ 1.3915817*T**2-0.00000044*T**3-0.000029956*T**4-0.0000000368*T**5)/15
	GMST_rad = GMST*2*pi/86400
	R = _R3(GMST_rad)
	return R

def _ecef2lla(ecef):
	x = ecef[0,0]
	y = ecef[1,0]
	z = ecef[2,0]
	a = 6378137
	e = 8.1819190842622e-2

	# calculations:
	b   = sqrt(a**2*(1-e**2));
	ep  = sqrt((a**2-b**2)/b**2);
	p   = sqrt(x**2+y**2);
	th  = atan2(a*z,b*p);
	lon = atan2(y,x);
	lat = atan2((z+ep**2*b*sin(th)**3),(p-e**2*a*cos(th)**3));
	N   = a/sqrt(1-e**2*sin(lat)**2);
	alt = p/cos(lat)-N;

	# return lon in range [0,2*pi)
	#lon = lon % 2*pi;

	k=abs(x)<1 and abs(y)<1

	if k:
		alt = abs(z)-b;

	lla = np.matrix([ [lat] , [lon], [alt] ])
	return lla

def est_lla(theta1, theta2, posix_1, posix_2):
	r_sun_eci = _sun_eci(posix_1)
	R_eci2ecef = _eci2ecef(posix_1)
	r_sun1 = R_eci2ecef*r_sun_eci;

	r_sun_eci = _sun_eci(posix_2)
	R_eci2ecef = _eci2ecef(posix_2)
	r_sun2 = R_eci2ecef*r_sun_eci

	r_me = _est_ecef(theta1, theta2, r_sun1, r_sun2)

	a = 6378137
	r_ecef1 = np.transpose(r_me[0,:]*a)
	r_ecef2 = np.transpose(r_me[1,:]*a)

	lla1 = _ecef2lla(r_ecef1)
	lla2 = _ecef2lla(r_ecef2)

	ll1 = lla1[ [0,1],:]
	ll2 = lla2[ [0,1],:]

	return ll1, ll2


def _est_ecef(theta1, theta2, r_sun1, r_sun2):
	r_sun1x = r_sun1[0,0]
	r_sun1y = r_sun1[1,0]
	r_sun1z = r_sun1[2,0]
	r_sun2x = r_sun2[0,0]
	r_sun2y = r_sun2[1,0]
	r_sun2z = r_sun2[2,0]

	theta1 = cos(pi/2 - theta1)
	theta2 = cos(pi/2 - theta2)

	# dump from solve_eqs.py
	r_me = \
	np.matrix([((-r_sun1x*r_sun1y*r_sun2y*theta2 - r_sun1x*r_sun1z*r_sun2z*theta2 + r_sun1x*r_sun2y**2*theta1 + r_sun1x*r_sun2z**2*theta1 + r_sun1y**2*r_sun2x*theta2 - r_sun1y*r_sun2x*r_sun2y*theta1 - r_sun1y*r_sun2z*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2) + r_sun1z**2*r_sun2x*theta2 - r_sun1z*r_sun2x*r_sun2z*theta1 + r_sun1z*r_sun2y*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2))/(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2), (r_sun1x**2*r_sun2y*theta2 - r_sun1x*r_sun1y*r_sun2x*theta2 - r_sun1x*r_sun2x*r_sun2y*theta1 + r_sun1x*r_sun2z*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2) - r_sun1y*r_sun1z*r_sun2z*theta2 + r_sun1y*r_sun2x**2*theta1 + r_sun1y*r_sun2z**2*theta1 + r_sun1z**2*r_sun2y*theta2 - r_sun1z*r_sun2x*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2) - r_sun1z*r_sun2y*r_sun2z*theta1)/(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2), (r_sun1x**2*r_sun2z*theta2 - r_sun1x*r_sun1z*r_sun2x*theta2 - r_sun1x*r_sun2x*r_sun2z*theta1 + r_sun1y**2*r_sun2z*theta2 - r_sun1y*r_sun1z*r_sun2y*theta2 - r_sun1y*r_sun2y*r_sun2z*theta1 + r_sun1z*r_sun2x**2*theta1 + r_sun1z*r_sun2y**2*theta1 + (-r_sun1x*r_sun2y + r_sun1y*r_sun2x)*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2))/(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2)), ((-r_sun1x*r_sun1y*r_sun2y*theta2 - r_sun1x*r_sun1z*r_sun2z*theta2 + r_sun1x*r_sun2y**2*theta1 + r_sun1x*r_sun2z**2*theta1 + r_sun1y**2*r_sun2x*theta2 - r_sun1y*r_sun2x*r_sun2y*theta1 + r_sun1y*r_sun2z*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2) + r_sun1z**2*r_sun2x*theta2 - r_sun1z*r_sun2x*r_sun2z*theta1 - r_sun1z*r_sun2y*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2))/(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2), (r_sun1x**2*r_sun2y*theta2 - r_sun1x*r_sun1y*r_sun2x*theta2 - r_sun1x*r_sun2x*r_sun2y*theta1 - r_sun1x*r_sun2z*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2) - r_sun1y*r_sun1z*r_sun2z*theta2 + r_sun1y*r_sun2x**2*theta1 + r_sun1y*r_sun2z**2*theta1 + r_sun1z**2*r_sun2y*theta2 + r_sun1z*r_sun2x*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2) - r_sun1z*r_sun2y*r_sun2z*theta1)/(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2), (r_sun1x**2*r_sun2z*theta2 - r_sun1x*r_sun1z*r_sun2x*theta2 - r_sun1x*r_sun2x*r_sun2z*theta1 + r_sun1y**2*r_sun2z*theta2 - r_sun1y*r_sun1z*r_sun2y*theta2 - r_sun1y*r_sun2y*r_sun2z*theta1 + r_sun1z*r_sun2x**2*theta1 + r_sun1z*r_sun2y**2*theta1 + (r_sun1x*r_sun2y - r_sun1y*r_sun2x)*sqrt(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - r_sun1x**2*theta2**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + 2*r_sun1x*r_sun2x*theta1*theta2 + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - r_sun1y**2*theta2**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + 2*r_sun1y*r_sun2y*theta1*theta2 + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2 - r_sun1z**2*theta2**2 + 2*r_sun1z*r_sun2z*theta1*theta2 - r_sun2x**2*theta1**2 - r_sun2y**2*theta1**2 - r_sun2z**2*theta1**2))/(r_sun1x**2*r_sun2y**2 + r_sun1x**2*r_sun2z**2 - 2*r_sun1x*r_sun1y*r_sun2x*r_sun2y - 2*r_sun1x*r_sun1z*r_sun2x*r_sun2z + r_sun1y**2*r_sun2x**2 + r_sun1y**2*r_sun2z**2 - 2*r_sun1y*r_sun1z*r_sun2y*r_sun2z + r_sun1z**2*r_sun2x**2 + r_sun1z**2*r_sun2y**2))])

	return r_me
