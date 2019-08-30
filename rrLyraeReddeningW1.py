# get current dir
import os
#here = os.path.dirname(__file__)
# use matplotlib without X server
import matplotlib
matplotlib.use('Agg')

import numpy as np
import scipy.interpolate as spi
import h5py
from astropy import units as u
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarWebQuery
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from matplotlib import *
import matplotlib
import time
from scipy.optimize import leastsq
from scipy.stats import norm
import matplotlib.mlab as mlab

# Module global vars
hdf5file ='stilism_cube.h5'
headers = None
cube = None
cubeXErr = None
cubeYErrMin = None
cubeYErrMax = None
axes = None
min_axes = None
max_axes = None

# For map, used by S. Ferron
step = None
hw = None
points = None
s = None

# Load HDF5 file
def init(hdf5file):
    """Load hdf5, calculate axes values corresponding to data.

    Args:
        hdf5file (str): full path for STILISM HDF5 file.

    Returns:
        dict: headers contains in HDF5 file.
        :func:`np.array`: 3D array which contains E(B-V).
        tuple: (x, y, z) where x,y,z contains array of axes
            corresponding to cube values.
        array: value min for x, y, z axes.
        array: value max for x, y, z axes.

    """
    # read hdf5 file
    with h5py.File(hdf5file, 'r') as hf:
        cube = hf['stilism/cube_datas'][:]
        cubeXerr = hf['stilism/cube_err_distance'][:]
        cubeYerrMin = hf['stilism/cube_err_magnitudemin'][:]
        cubeYerrMax = hf['stilism/cube_err_magnitudemax'][:]
        dc = hf['stilism/cube_datas']

        # Less method call are done with this version:
        headers = {k: v for k, v in dc.attrs.items()}

    sun_position = headers["sun_position"]
    gridstep_values = headers["gridstep_values"]

    # Calculate axes for cube value, with sun at position (0, 0, 0)
    min_axes = -1 * sun_position * gridstep_values
    max_axes = np.abs(min_axes)
    axes = (
        np.linspace(min_axes[0], max_axes[0], cube.shape[0]),
        np.linspace(min_axes[1], max_axes[1], cube.shape[1]),
        np.linspace(min_axes[2], max_axes[2], cube.shape[2])
    )

    # S. Ferron variable for map
    step = np.array(headers["gridstep_values"])
    hw = (np.copy(cube.shape) - 1) / 2.
    points = (
        np.arange(0, cube.shape[0]),
        np.arange(0, cube.shape[1]),
        np.arange(0, cube.shape[2])
    )
    s = hw * step

    return (headers, cube, cubeXerr, cubeYerrMin, cubeYerrMax,
        axes, min_axes, max_axes,
        step, hw, points, s)


# Initialisation
headers, cube, cubeXErr, cubeYErrMin, cubeYErrMax, axes, min_axes, max_axes, step, hw, points, s = init(hdf5file)


def reddening(vlong, ulong, vlat, ulat, frame, step_pc=5):
    """Calculate Reddening versus distance from Sun.

    Args:
        vlong (str or double): Longitude value.
        ulong (str): Longitude unit used in :class:`SkyCoord`.
        vlat (str or double): Latitude value.
        ulat (str): Latitude unit used in :class:`SkyCoord`.
        frame (str): Galactic, icrs ... values supported by :class:`SkyCoord`.

    Kwargs:
        step_pc (int): Incremental distance in parsec

    Returns:
        array: Parsec values.
        array: E(B-V) value obtain with integral of linear extrapolation.

    """
    # Calculate the position for 1pc
    sc = SkyCoord(
        vlong,
        vlat,
        distance=1 * u.pc,
        unit=(ulong, ulat),
        frame=frame
    )
    coords_xyz = sc.transform_to('galactic').represent_as('cartesian').get_xyz().value

    # Find the number of parsec I can calculate before go out the cube
    # (exclude divide by 0)
    not0 = np.where(coords_xyz != 0)
    max_pc = np.amin(
        np.abs(
            np.take(max_axes, not0) / np.take(coords_xyz, not0)
        )
    )

    # Calculate all coordinates to interpolate (use step_pc)
    distances = np.arange(0, max_pc, step_pc)
    sc = SkyCoord(
        vlong,
        vlat,
        distance=distances,
        unit=(ulong, ulat, 'pc'),
        frame=frame
    )
    sc = sc.transform_to('galactic').represent_as('cartesian')
    coords_xyz = np.array([coord.get_xyz().value for coord in sc])

    # linear interpolation with coordinates
    interpolation = spi.interpn(
        axes,
        cube,
        coords_xyz,
        method='linear'
    )
    xvalues = np.arange(0, len(interpolation) * step_pc, step_pc)
    yvalues = np.cumsum(interpolation) * step_pc

    # errors
    xerrors = spi.interpn(
        axes,
        cubeXErr,
        coords_xyz,
        method='linear'
    )
    yerrorsMin = spi.interpn(
        axes,
        cubeYErrMin,
        coords_xyz,
        method='linear'
    )
    yerrorsMax = spi.interpn(
        axes,
        cubeYErrMax,
        coords_xyz,
        method='linear'
    )

    return (
        xvalues,
        np.around(yvalues, decimals=3),
        np.around(xerrors, decimals=0),
        np.around(yerrorsMin, decimals=3),
        np.around(yerrorsMax, decimals=3)
    )


def cube_cut(vlong, ulong, vlat, ulat, frame, vdist, udist, vnlong, unlong, vnlat, unlat):
    """Calculate map cut of cube.

    Args:
        vlong (str or double): Longitude value.
        ulong (str): Longitude unit used in :class:`SkyCoord`.
        vlat (str or double): Latitude value.
        ulat (str): Latitude unit used in :class:`SkyCoord`.
        frame (str): Galactic, icrs, etc
            (values supported by :class:`SkyCoord`).
        vdist (str or double): Distance.
        udist (str): Distance unit used in Skycoord.
        vnlong (str or double): Normal longitude value.
        unlong (str or double): Longitude unit used in :class:`SkyCoord`
            for normal.
        vnlat (str or double): Normal latitude value.
        unlat (str): Latitude unit used in :class:`SkyCoord` for normal.

    Returns:
        img (str):plot of map

    """
    #  Transforming the reference position point into cartesian coordinates:
    sc = SkyCoord(
        vlong,
        vlat,
        distance=vdist,
        unit=(ulong, ulat, udist),
        frame=frame
    )
    r = sc.represent_as('cartesian').get_xyz().value

    #  Getting the normal to the plane in cartesian coordinates:
    wp = SkyCoord(
        vnlong,
        vnlat,
        unit=(unlong, unlat),
        frame=frame
    )
    w = wp.represent_as('cartesian').get_xyz().value

    #  Creating a direct base (u, v, w) using the normal vector:
    lon = np.radians(wp.l.degree)
    u = [np.sin(lon), -np.cos(lon), 0]
    u_latitude = 0 # arcsin(0) = 0 because u[2] = 0
    u_longitude = np.degrees(np.arctan2(u[1], u[0]))

    #  The last vector is just the results of the vector product of w and u
    #  (the order is important to keep the referential oriented as expected):
    v = np.cross(w, u)
    v_latitude = np.degrees(np.arcsin(v[2]))
    v_longitude = np.degrees(np.arctan2(v[1], v[0]))

    # maximum extension maximal of the slice
    #  Taking the two norm of the vector giving the sun position
    #  (if I have understand well the meaning of hw and step):
    f = np.linalg.norm(2 * hw * step, 2)
    #  Then we construct an array giving the cube border:
    c = np.arange(-f, f, step.min())

    #  Transforming this array into a mesh grid (it consist essentially of a
    #  coordinate matrix):
    cu, cv = np.meshgrid(c, c)

    #  The outer product of two vector u_1 and u_2 is equivalent to do the
    #  matrix multiplication u_1 u_2^T (T:: transpose). Here, the
    #  :func:`numpy.outer` function will flatten the matrix of size
    #  (amin(step), amin(step)) to have a vector of dimension $amin(step)^2$.
    #  Then it will multiple, following the matrix multiplication rule all
    #  coordinates of the argument to get in return a matrix of size
    #  (amin(step)^2, len(u)) (same for v and $r+s$ which should be of the same
    #  length).
    #  You can also search for the definition of the kroenecker product as the
    #  outer product is supposed to be a special case of it.
    #  This should give us a matrix to transform the cube coordinates into the
    #  plane coordinate, if I have understand everything correctly:
    cz = np.outer(cu, u) + np.outer(cv, v) + np.outer(np.ones_like(cu), r + s)

    z = np.ones([c.size * c.size])
    z[:] = None

    ix = np.floor((cz[:, 0]) / step[0])
    iy = np.floor((cz[:, 1]) / step[1])
    iz = np.floor((cz[:, 2]) / step[2])
    w = np.where(
        (ix >= 0) & (ix <= cube.shape[0] - 1) &
        (iy >= 0) & (iy <= cube.shape[1] - 1) &
        (iz >= 0) & (iz <= cube.shape[2] - 1)
    )

    z[w] = spi.interpn(
        points,
        cube,
        (ix[w], iy[w], iz[w]),
        method='linear',
        fill_value=1
    )
    z = np.reshape(z, [c.size, c.size])

    f = np.take(
        np.reshape(cu, [c.size * c.size]),
        w
    )
    wx = np.squeeze(
        np.where(
            (c >= np.amin(f)) & (c <= np.amax(f))
        )
    )
    f = np.take(
        np.reshape(cv, [c.size * c.size]),
        w
    )
    wy = np.squeeze(
        np.where(
            (c >= np.amin(f)) & (c <= np.amax(f))
        )
    )

    smap = z[np.ix_(wy, wx)]
    
    if v_latitude < 0:
        v_latitude = -v_latitude
        v_longitude = (v_longitude + 180)%360
        smap = smap[::-1]

    result = {}
    
    # x,y values for axes, log(z) values for a better representation of contour map
    addforLog0 = 1e-7
    result["addforlog0"] = addforLog0
    result["xJsTab"] = np.array2string(c[wx][::2], max_line_width=80, separator=',', precision=2, threshold=np.nan).replace("nan", "NaN")
    result["xTitle"] = "'Left to right towards ⇒ (l=%.1f°,b=%.1f°)'"%(u_longitude, u_latitude)
    result["yJsTab"] = np.array2string(c[wy][::2], max_line_width=80, separator=',', precision=2, threshold=np.nan).replace("nan", "NaN")
    result["yTitle"] = "'Bottom to top towards ⇒ (l=%.1f°,b=%.1f°)'"%(v_longitude, v_latitude)
    logsmap = np.log(smap[::2,::2] + addforLog0)
    result["zJsTab"] = np.array2string(logsmap, max_line_width=80, separator=',', precision=2, threshold=np.nan).replace("nan", "NaN")
    
    # Calculate color bar values (5 values between [min, max])
    # Use log because log(z), and corresponding true value
    logScaletmp = np.linspace(np.nanmin(logsmap), np.nanmax(logsmap), 5)
    result["logScale"] = np.array2string(logScaletmp, separator=",")
    scaletmp = [format(np.exp(v) - addforLog0, '.2e').replace("e-0", "e-") for v in logScaletmp]
    result["scale"] = "[" + ", ".join(["'%s'"%v for v in scaletmp]) + "]"
    
    result["title"] = "'Origin (l=%.1f°,b=%.1f°), distance=%.1fpc<br />Normal to the plane (l=%.1f°,b=%.1f°)'"%(
            sc.l.degree, sc.b.degree,  sc.distance.value,
            wp.l.degree, wp.b.degree)
    
    return result

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


np.set_printoptions(suppress=True,formatter={'float_kind':'{:f}'.format})

gaiaid=[]
parallax=[]
parallaxError=[]
period=[]
w1mean=[]
w1error=[]
w1ampl=[]
w1chi2=[]
w1sig=[]
w1_t0=[]
w2mean=[]
w2err=[]
w2ampl=[]
w2sig=[]
w2_t0=[]
feh=[]
Gmag=[]
chi2=[]
n_good_obs=[]
UWE=[]
Wise_qflag=[]
gaia_id=[]
RA=[]
DEC=[]


EBV=[]
newDistance=[]
catalog=[]
flag=[]

file = open('W1_fehUWE2.dat')
for line in file:
    if '#' not in line:
        line=line.strip()
        line=line.split()
        line=np.asarray(line)
        gaiaid.append(line[0])
        parallax.append(line[1])
        parallaxError.append(line[2])
        period.append(line[3])
        w1mean.append(line[4])
        w1error.append(line[5])
        w1ampl.append(line[6])
        w1chi2.append(line[7])
        w1sig.append(line[8])
        w1_t0.append(line[9])
        w2mean.append(line[10])
        w2err.append(line[11])
        w2ampl.append(line[12])
        w2chi2.append(line[13])
        w2sig.append(line[14])
        w2_t0.append(line[15])
        feh.append(line[16])
        efeh.append(line[17])
        Gmag.append(line[18])
        chi2.append(line[19])
        n_good_obs.append(line[20])
        UWE.append(line[21])
        Wise_qflag.append(line[21])
        gaia_id.append(line[22])
        RA.append(line[23])
        DEC.append(line[24])


gaiaid=np.array(gaiaid,dtype=float)
parallax=np.array(parallax,dtype=float)
parallaxError=np.array(parallaxError,dtype=float)
period=np.array(period,dtype=float)
w1mean=np.array(w1mean,dtype=float)
w1error=np.array(w1error,dtype=float)
w1ampl=np.array(w1ampl, dtype=float)
w1chi2=np.array(w1chi2, dtype=float)
w1sig=np.array(w1sig, dtype=float)
w1_t0=np.array(w1_t0, dtype=float)
w2mean=np.array(w2mean, dtype=float)
w2err=np.array(w2err, dtype=float)
w2ampl=np.array(w2ampl, dtype=float)
w2sig=np.array(w2sig, dtype=float)
w2_t0=np.array(w2_t0, dtype=float)
feh=np.array(feh, dtype=float)
Gmag=np.array(Gmag, dtype=float)
chi2=np.array(chi2, dtype=float)
n_good_obs=np.array(n_good_obs, dtype=float)
UWE=np.array(UWE, dtype=float)
Wise_qflag=np.array(Wise_qflag, dtype=str)
gaia_id=np.array(gaia_id, dtype=float)
RA=np.array(RA, dtype=float)
DEC=np.array(DEC, dtype=float)


absW1=−2.6*(np.log10(period)+0.27)+0.11*(feh+1.3)−0.42

distanceExp=(w1mean-absW1+5)/5

distance=10**distanceExp #in pc


file = open('ReddeningW1.dat','w')
file.write('#gaia_id parallax parallax_err period w1mean w1err w1ampl w1chi2 w1sig w1_t0 w2mean w2err w2ampl w2chi2 w2sig w2_t0 feh efeh  Gmag chi^2  n_good_obs UWE  Wise_qflag  gaia_id RA  DEC EBV catalog'+'\n')
file.write('#1        2        3            4       5       6     7      8     9     10    11     12    13    14      15    16   17   18    19    20      21      22   23          24     25  26    27 28'+'\n')
file.close()

for index in range(0,20):
    try:
        time.sleep(1)
        coords=SkyCoord(RA[index]*u.deg,Dec[index]*u.deg,distance=distance[index]*u.pc,frame='icrs')
        reddeningValue=bayestar(coords,mode='median')
#        if np.isnan(reddeningValue)==True:
#            catalogFlag=1
#            reddeningList=reddening(RA[index],'deg',Dec[index],'deg','icrs')
#            nearestDistanceIndex=find_nearest(reddeningList[0],distance[index])
#            nearestDistance=reddeningList[0][nearestDistanceIndex]
#            reddeningValue=reddeningList[1][nearestDistanceIndex]
#            distanceDifference=distance[index]-nearestDistance

 #       AW1=0.19*reddeningValue

#        newDistance=10**(0.2*(w1mean[index]-absmag[index]+5-AW1))
        #coords=SkyCoord(RA[index]*u.deg,Dec[index]*u.deg,distance=newDistance*u.pc,frame='icrs')
        #reddeningValue=bayestar(coords,mode='median')
        if np.isnan(reddeningValue)==True:
            catalogFlag="Stilism"
            reddeningList=reddening(RA[index],'deg',Dec[index],'deg','icrs')
            nearestDistanceIndex=find_nearest(reddeningList[0],distance[index])
            nearestDistance=reddeningList[0][nearestDistanceIndex]
            reddeningValue=reddeningList[1][nearestDistanceIndex]
            distanceDifference=newDistance-nearestDistance
            print("Stilism")
        else:
            catalogFlag="Green"
            #coordsMax=SkyCoord(RA[index]*u.deg,Dec[index]*u.deg,distance=(newDistance+100)*u.pc,frame='icrs')
            #coordsMin=SkyCoord(RA[index]*u.deg,Dec[index]*u.deg,distance=(newDistance-100)*u.pc,frame='icrs')
            #reddeningValueMax=bayestar(coordsMax,mode='median')
            #reddeningValueMin=bayestar(coordsMin,mode='median')
            print("Green")

        file = open('ReddeningW1.dat','a')
        file.write("{:.0f}".format(float(gaiaid[index]))+' '+str(reddeningValue)+' '+str(catalogFlag)+'\n')
        file.close()

    except:
        file = open('ReddeningVersion4.dat','a')
        file.write(str(gaiaid[index])+' ERROR \n')
        file.close()

