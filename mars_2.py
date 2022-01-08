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
index = 2
# print(fits_files)
# print(img_files)

for index in np.linspace(22, 42, endpoint=False): #[26, 27]: #range(len(fits_files)):
    index = int(index)
    fig = plt.figure()
    fig.canvas.manager.full_screen_toggle()

    # for index in range(len(fits_files)):
    # Spectro view
    hdul = fits.open(join(file_dir, fits_files[index]))
    hdul.verify('fix')
    hdr = hdul[1].header
    # print(repr(hdul[1].header))
    data = hdul[1].data
    utc = data[0][9]
    utc2 = data[-1][9]
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
    xy_bound = [-PI * 1.5, PI * 1.5, (-PI / 2) * 1.5, (PI / 2) * 1.5]
    XX, YY = np.mgrid[xy_bound[0]:xy_bound[1]:1080j, xy_bound[2]:xy_bound[3]:540j]
    ZZ = spint.griddata(
        list(zip(p_longs + add_p_longs, p_lats + add_p_lats)),
        p_max_temps + add_p_max_temps,
        (XX, YY), method='cubic'
    )

    ax = fig.add_subplot(111)
    mn = 100  # min(max_temps)
    mx = 300  # max(max_temps)
    ax.imshow(ZZ.T, extent=(xy_bound[0], xy_bound[1], xy_bound[2], xy_bound[3]), origin='lower',
              cmap='plasma', vmin=mn, vmax=mx)
    ax.scatter(p_longs, p_lats, c=p_max_temps, cmap='plasma', marker='o',
               vmin=mn, vmax=mx, edgecolor='white', s=6, linewidth=0.5)
    ax.set_aspect(1)
    ax.set_title("EMIRS max brightness temperature" + "\n" +
                 "Equirectangular projection with cubic interpolation" + "\n" + utc + " -\n" + utc2, pad=12)
    ax.set_xlim(xy_bound[0], xy_bound[1])
    ax.set_ylim(xy_bound[2], xy_bound[3])
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
    cb = fig.colorbar(cm.ScalarMappable(cmap=cm.plasma, norm=norm), ax=ax, shrink=0.522)
    cb.set_label("Maximum brightness temp (K)", labelpad=10)
    rect = patches.Rectangle((-PI, -PI / 2), 2 * PI, PI, linewidth=1,
                             linestyle='dotted', edgecolor='w', facecolor='none')
    ax.add_patch(rect)
    ax.add_patch(PolygonPatch(ashape, alpha=.4, linewidth=.5,
                              linestyle='-', edgecolor='w', facecolor='none'))

    c_fig = plt.gcf()
    # plt.show()  # block=True
    c_fig.set_size_inches((11, 8.5), forward=False)
    # c_fig.savefig("export/10_{}".format(index), dpi=300)
    c_fig.savefig("export/12_{}".format(index), dpi=300, transparent=True)
    plt.clf()
    plt.close(fig)
