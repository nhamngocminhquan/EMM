import alphashape
from shapely.geometry import Point
from astropy.io import fits
from astropy.time import TimeDelta, Time
from astropy import units as u
from descartes import PolygonPatch
import emm_api as ep
import matplotlib
from matplotlib import cm
from matplotlib import patches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import pickle
import scipy.interpolate as spint
from scipy.spatial import Delaunay
import sys
plt.style.use('dark_background')

PI = 180  # 3.141592653589793238

# Get files
file_dir = "emm_data/emirs"
file_character = "emr"
fits_files = [filename for filename in listdir(file_dir) if (isfile(join(file_dir, filename))
                                                             and file_character in filename)]

# print(fits_files)
# print(img_files)

offset = 0 # 435780
initial = 666656480
mars_day_secs = 88619.616
num = 12
xx = np.linspace(-PI, PI, num=num, endpoint=False).tolist()
y = 0
zz = []

endval = len(fits_files)
for index in np.linspace(0, endval, num=endval, endpoint=False):# [26, 27]: # range(len(fits_files)):
    index = int(index)

    # for index in range(len(fits_files)):
    # Spectro view
    hdul = fits.open(join(file_dir, fits_files[index]))
    hdul.verify('fix')
    hdr = hdul[1].header
    # print(repr(hdul[1].header))
    data = hdul[1].data
    utc = data[0][9]
    utc2 = data[-1][9]
    second = data[0][0] + offset - initial + data[0][1] * 1.0 / 65536
    hdul.close()
    # print(data)

    # Field 3 is calibrated radiance
    # print(data[5][3])

    # Field 7 is max brightness temp
    # Fields 21 and 22 are latitudes and longitudes of IFOV
    dstgui = []
    for i in range(len(data)):
        if (data[i][7] > 0.001 and data[i][7] < 320) and (abs(data[i][21]) > 0.001 or abs(data[i][22]) > 0.001) and data[i][14]: #
            dstgui.append(1)
        else:
            dstgui.append(0)

    max_temps = [data[i][7] for i in range(len(data)) if dstgui[i] == 1]
    lats = [data[i][21] for i in range(len(data)) if dstgui[i] == 1]
    longs = [data[i][22] for i in range(len(data)) if dstgui[i] == 1]
    # print(max_temps)

    # Surround OG map with 6 copies for boundary condition
    p_lats = lats + lats + lats + (-np.array(lats) + PI).tolist() + (-np.array(lats) + PI).tolist() + \
             (-np.array(lats) - PI).tolist() + (-np.array(lats) - PI).tolist()
    p_longs = longs + (np.array(longs) + 2 * PI).tolist() + (np.array(longs) - 2 * PI).tolist() + \
              (-np.array(longs) + PI).tolist() + (-np.array(longs) - PI).tolist() + \
              (-np.array(longs) + PI).tolist() + (-np.array(longs) - PI).tolist()
    p_max_temps = max_temps + max_temps + max_temps + max_temps + max_temps + max_temps + max_temps

    # Finding boundary of data
    if PI == 180:
        alpha = 0.025
    else:
        alpha = 0.025 * 180 / PI
    points = np.array(list(zip(
        [p for p in p_longs if (p < 3 * PI) and (p > -3 * PI)],
        [p for p in p_lats if (p < 3 * PI / 2) and (p > -3 * PI / 2)]
    )))
    ashape = alphashape.alphashape(points, alpha)

    # Initial interpolation
    xy_bound = [-PI, PI, (-PI / 2), (PI / 2)]
    XX, YY = np.mgrid[xy_bound[0]:xy_bound[1]:30j, xy_bound[2]:xy_bound[3]:2j]
    ZZ = spint.griddata(
        list(zip(p_longs, p_lats)),
        p_max_temps,
        (XX, YY), method='cubic'
    )

    add_p_longs = XX[:, 0].tolist() + XX[:, 1].tolist() + \
                  (np.array(XX[:, 0]) - 2 * PI).tolist() + (np.array(XX[:, 0]) - 2 * PI).tolist() + \
                  (np.array(XX[:, 0]) + 2 * PI).tolist() + (np.array(XX[:, 0]) + 2 * PI).tolist()
    add_p_lats = YY[:, 0].tolist() + YY[:, 1].tolist() + \
                 YY[:, 0].tolist() + YY[:, 1].tolist() + \
                 YY[:, 0].tolist() + YY[:, 1].tolist()
    # print(add_p_lats, add_p_longs)
    add_p_max_temps = ([sum(ZZ[:, 0]) / len(ZZ[:, 0])] * len(ZZ[:, 0])) + (
                [sum(ZZ[:, 1]) / len(ZZ[:, 1])] * len(ZZ[:, 1])) + \
                      ([sum(ZZ[:, 0]) / len(ZZ[:, 0])] * len(ZZ[:, 0])) + (
                                  [sum(ZZ[:, 1]) / len(ZZ[:, 1])] * len(ZZ[:, 1])) + \
                      ([sum(ZZ[:, 0]) / len(ZZ[:, 0])] * len(ZZ[:, 0])) + (
                                  [sum(ZZ[:, 1]) / len(ZZ[:, 1])] * len(ZZ[:, 1]))

    # print(len(list(zip([*p_longs, *add_p_longs], [*p_lats, *add_p_lats]))))
    # print(len([*p_max_temps, *add_p_max_temps]))

    # Re-interpolation with approximated poles
    xy_bound = [-PI, PI, y, y]
    XX, YY = np.mgrid[xy_bound[0]:xy_bound[1]:(PI * 2 / num), xy_bound[2]:xy_bound[3]:1j]
    ZZ = spint.griddata(
        list(zip(p_longs + add_p_longs, p_lats + add_p_lats)),
        p_max_temps + add_p_max_temps,
        (XX, YY), method='cubic'
    )
    ZZ = np.ravel(ZZ).tolist()
    for i in range(num):
        inshape = False
        for poly in ashape.geoms:
            if poly.contains(Point(xx[i], y)):
                inshape = True

        if not inshape:
            ZZ[i] = np.nan

    ZZ.append(float(second))
    zz.append(ZZ)
    print(index)
    # if float(second) > 15 * mars_day_secs:
    #     print(ZZ)

    # ax = fig.add_subplot(111)
    # mn = 100  # min(max_temps)
    # mx = 300  # max(max_temps)
    # ax.imshow(ZZ.T, extent=(xy_bound[0], xy_bound[1], xy_bound[2], xy_bound[3]), origin='lower',
    #           cmap='plasma', vmin=mn, vmax=mx)
    # ax.scatter(p_longs, p_lats, c=p_max_temps, cmap='plasma', marker='o',
    #            vmin=mn, vmax=mx, edgecolor='white', s=6, linewidth=0.5)
    # ax.set_aspect(1)
    # ax.set_title("EMIRS max brightness temperature" + "\n" +
    #              "Equirectangular projection with cubic interpolation" + "\n" + utc + " -\n" + utc2, pad=12)
    # ax.set_xlim(xy_bound[0], xy_bound[1])
    # ax.set_ylim(xy_bound[2], xy_bound[3])
    # ax.set_xlabel("Longitude (deg)")
    # ax.set_ylabel("Latitude (deg)")
    # norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
    # cb = fig.colorbar(cm.ScalarMappable(cmap=cm.plasma, norm=norm), ax=ax, shrink=0.522)
    # cb.set_label("Maximum brightness temp (K)", labelpad=10)
    # rect = patches.Rectangle((-PI, -PI / 2), 2 * PI, PI, linewidth=1,
    #                          linestyle='dotted', edgecolor='w', facecolor='none')
    # ax.add_patch(rect)
    # ax.add_patch(PolygonPatch(ashape, alpha=.4, linewidth=.5,
    #                           linestyle='-', edgecolor='w', facecolor='none'))

# zz = list(zip(*zz))
# print(zz)
# print(zz[:, 1])
# print(zz[num - 1])
# fig, axes = plt.subplots(nrows=round(np.sqrt(num)), ncols=int(np.ceil(num/np.sqrt(num)) + 0.01))
max_sols = int(np.ceil(zz[endval - 1][num] / mars_day_secs) + 0.01)
col_count = 4
norm = matplotlib.colors.BoundaryNorm(range(max_sols + 1), cm.plasma.N, extend='neither')
fig, axes = plt.subplots(nrows=round(np.sqrt(num)), ncols=int(np.ceil(num/np.sqrt(num)) + 0.01),
                         sharex=True, sharey=True)
fig.canvas.manager.full_screen_toggle()
for i in range(num):
    buffr = np.floor(i / col_count)
    ax = axes[round(buffr), round(i - buffr * col_count)] # axes[round(i - buffr * round(np.sqrt(num))), round(buffr)]
    sols = 0
    x = []
    max_temps = []
    for j in range(len(zz)):
        sols_passed = int(zz[j][num] / mars_day_secs)
        if not sols_passed == sols:
            # ax.plot(np.array(x), np.array(max_temps))
            ax.scatter(np.array(x), np.array(max_temps), c=[sols_passed] * len(x),
                       cmap='plasma', norm=norm) # vmin=0, vmax=max_sols)
            x = []
            max_temps = []
            sols = sols + 1

        x.append((zz[j][num] - mars_day_secs * sols_passed) / 10000)
        max_temps.append(zz[j][i])

    mn = 120
    mx = 280

    # ax.set_aspect(1)
    # plt.gca().set_aspect('equal')
    # ax.set_title("EMIRS max brightness temperature" + "\n" +
    #              "Equirectangular projection with cubic interpolation" + "\n" + utc + " -\n" + utc2, pad=12)
    ax.set_xlim(0, np.ceil(mars_day_secs) * 1.0 / 10000)
    ax.set_ylim(mn, mx)
    ax.set_title("{}°".format(xx[i]))
    # ax.set_xlabel("Seconds in day (s)")
    # ax.set_ylabel("Max brightness temp (K)")
    # norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
    # cb = fig.colorbar(cm.ScalarMappable(cmap=cm.plasma, norm=norm), ax=ax, shrink=0.522)
    # cb.set_label("Maximum brightness temp (K)", labelpad=10)

# norm = matplotlib.colors.Normalize(vmin=0, vmax=max_sols)
cb = fig.colorbar(cm.ScalarMappable(cmap=cm.plasma, norm=norm), ax=axes.ravel().tolist())#, ticks=range(max_sols))#, shrink=0.5)
cb.set_label("Mars days", labelpad=10)
ax = fig.add_subplot(111, frameon=False)
ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
ax.grid(False)
ax.set_xlabel("Time of Mars day (10^5 seconds)\n(initial day starts with first EMIRS observation)")
ax.set_ylabel("Maximum brightness temperature (K)")
ax.xaxis.set_label_coords(0.5, -0.055)
ax.yaxis.set_label_coords(-0.045, 0.5)
fig.suptitle("Daily variation in maximum brightness temperature \n at {}° latitude and various longitudes".format(y))

c_fig = plt.gcf()
plt.show()  # block=True
# c_fig.set_size_inches((11, 8.5), forward=False)
# c_fig.savefig("export/10_{}".format(index), dpi=300)
c_fig.savefig("export/12_{}_1".format(3), dpi=300, transparent=True)
# plt.clf()
# plt.close(fig)
