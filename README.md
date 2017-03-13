# PySextant

simple python tool to calculate latitude and longitude based on two observations of the elevation of the sun

## Math

Given two measurements of the angle between the horizon and the sun and two different times of day. Let's call the angle `theta_1` and `theta_2` and the times `tau_1` and `tau_2`.

The zenith angle, the compliement of `theta`, is the angle between the sky's zenith and the sun. These angles are:

```
phi_1 = pi/2 - theta_1
phi_2 = pi/2 - theta_2
```

If we assume the sun and position vectors are normalized then these angles are the arccosine of the [dot product](https://en.wikipedia.org/wiki/Dot_product) of vector pointing to the sun and the vector pointing to your location.

Using the assumption that your movement on the earth's surface is insignificant between `tau_1` and `tau_2` this gives us

```
cos(pi/2-theta_1) = a1 = sun_1 . pos
cos(pi/2-theta_2) = a2 = sun_2 . pos
```

Fortunately `sun_1` and `sun_2` are predictable functions of time, which we know (`tau_1` and `tau_2`).

Now we just have to solve for `pos` as a function of `sun_1`, `sun_2`, `a1`, and `a2`. To do this we'll use the [vector triple product](https://en.wikipedia.org/wiki/Triple_product) which states

```
a X (b X c) = b(a dot c) - c(a dot b)
```

Let's use the above with the vectors of interest for this problem and introduce some new helper variables to make the problem more tractable

```
pos X (sun_1 X sun_2) = sun_1(pos dot sun_2) - sun_2(pos dot sun_1)
pos X (sun_1 X sun_2) = pos X b = -c = sun_1*a2 - sun_2*a1
b X pos = c
```

Figuring out pos is the last step here. With some understanding of cross products and assuming each vector involved here is normalized we can arrive at

```
pos = (c X b)/(b.b) + t*b
pos = d + t*b;
```

Now we can make use of the fact that `||pos||=1` which means `pos . pos = 1`.

```
(d+tb).(d+tb) = 1
t^2(b.b) + 2t(b.d) + d.d - 1 = 0
```

And we can solve for t

```
t = [-2(b.d) +- sqrt(4(b.d)^2 - 4(b.b)*(d.d-1))]/(2b.b)
```

And to recap

```
b = sun_1 X sun_2
a1 = cos(pi/2 - theta_1)
a2 = cos(pi/2 - theta_2)
d = (c X b)/(b.b) = (sun_2*a1 - sun_1*a2) X (sun_1 X sun_2) / (b.b)
pos = d + t*b
```

You'll have two solutions for `t` and thus two solutions for `pos`. This will usually mean you need to know which hemisphere you're in.

There is a bit more that this utility does. The process is really more like this

1) calculate sun position in inertial coordinates at `tau_1`
2) calculate sun position in inertial coordinates at `tau_2`
3) rotate sun positions into an earth fixed frame
4) solve for `pos` using the math above
5) convert `pos` from cartesian into spherical coordinates (latitude, longitude)


## Running

a simple script is provided which shows how to use this module

```
git clone https://github.com/friedman101/PySextant.git
cd PySextant
./test.py
```