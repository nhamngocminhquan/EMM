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
from matplotlib.patches import PathPatch
from matplotlib.path import Path
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

PI = 180  # 3.141592653589793238


# Alpha convex hull code
# https://stackoverflow.com/a/50159452
def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it is not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    # if ax is not None:
    #     for t in tri.simplices:
    #         ax.plot([points[t[0]][0], points[t[1]][0]], [points[t[0]][1], points[t[1]][1]], c='w')
    #         ax.plot([points[t[2]][0], points[t[1]][0]], [points[t[2]][1], points[t[1]][1]], c='w')
    #         ax.plot([points[t[0]][0], points[t[2]][0]], [points[t[0]][1], points[t[2]][1]], c='w')
    #
    # edges = set()
    # for t in tri.simplices:
    #     a = points[t[0]]
    #     b = points[t[1]]
    #     c = points[t[2]]

    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges


# Get files
file_dir = "emm_data/0224"
file_character = "emr"
fits_files = [filename for filename in listdir(file_dir) if (isfile(join(file_dir, filename))
                                                             and file_character in filename)]
file_character = "xos1_f635_r_v03-01.fits"
img_files = [filename for filename in listdir(file_dir) if (isfile(join(file_dir, filename))
                                                            and file_character in filename)]
index = 2
# print(fits_files)
# print(img_files)

for index in range(len(fits_files)):
    # Time in sclk using conversion diverges 5 days 1 hour from labeled time
    # offset = 435600
    # for file in fits_files:
    #     hdul = fits.open(join(file_dir, file))
    #     hdul.verify('fix')
    #     hdr = hdul[1].header
    #     data = hdul[1].data
    #     hdul.close()
    #     # print(data[0][9])
    #     # print(data[hdr['NAXIS2'] - 1][9])
    #     print((Time(2000, format='jyear') + TimeDelta((data[0][0] + offset + data[0][1] * 1.0 / 65536) * u.s)).iso)
    #     print((Time(2000, format='jyear') + TimeDelta((data[hdr['NAXIS2'] - 1][0] + offset +
    #                                                    data[hdr['NAXIS2'] - 1][1] * 1.0 / 65536)*u.s)).iso)

    # Image view
    hdul = fits.open(join(file_dir, img_files[index]))
    # hdul = fits.open("emm_data/0224/emm_exi_l2a_20210224T095615_0011_xos1_f635_r_v03-01.fits")
    hdul.verify('fix')
    utc = hdul[0].header['DATE-OBS']
    utc = utc.partition('T')[0] + " " + utc.partition('T')[2]
    # print(repr(hdr))

    # 1: SCI img, 5: PHA, 7: LAT, 8: LONG, 9: HEIGHT, 10: INA
    interest = [1, 5, 7, 8, 9, 10]
    img = [hdul[i].data for i in interest]
    # print(img[4])
    # Fixing data
    hdr = hdul[1].header
    mn = hdr['DATAMIN']
    mx = hdr['DATAMAX']
    w = 2048
    h = 1536
    rg = mx - mn
    mask = np.zeros((h, w), dtype=bool)
    max_altitude = 10
    for i in range(h):
        for j in range(w):
            if not np.isnan(img[0][i][j]):
                if img[4][i][j] < max_altitude:
                    continue
                    # img[0][i][j] = (img[0][i][j] - mn) * 1.0 / rg
                else:
                    for k in range(len(interest)):
                        img[k][i][j] = np.nan

                # mask[i][j] = 1
            # else:
            # img[0][i][j] = 0

    hdul.close()

    fig = plt.figure()
    fig.canvas.manager.full_screen_toggle()
    radius = 1
    azim = 70
    elev = 12

    # Plot help
    # https://stackoverflow.com/a/42927880
    # ax = fig.add_subplot(1, 2, 1, projection='3d')
    # imgs = np.array(img[0])
    # lats = np.array(np.deg2rad(img[2]))
    # longs = np.array(np.deg2rad(img[3]))
    # norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
    # sf = ax.plot_surface(radius * np.cos(lats) * np.cos(longs),
    #                      radius * np.cos(lats) * np.sin(longs),
    #                      radius * np.sin(lats),
    #                      rstride=8, cstride=8,
    #                      facecolors=cm.gray(norm(imgs)))
    # facecolors=cm.gray(imgs))

    # ax.set_box_aspect((1, 1, 1))
    # ax.set_title("EXI image @ f635 on unit sphere" + "\n" + utc, pad=-10)
    # ax.set_xticks(np.arange(-radius, radius, 0.5))
    # ax.set_yticks(np.arange(-radius, radius, 0.5))
    # ax.set_zticks(np.arange(-radius, radius, 0.5))
    # ax.set_xlim(-radius, radius)
    # ax.set_ylim(-radius, radius)
    # ax.set_zlim(-radius, radius)
    # ax.view_init(azim=azim, elev=elev)
    # cb = fig.colorbar(cm.ScalarMappable(cmap=cm.gray, norm=norm), ax=ax, shrink=0.5)
    # cb.set_label("Calibrated DN", labelpad=10)
    # fig.colorbar(cm.ScalarMappable(cmap=cm.gray, norm=imgs).set_array([]), ax=ax, shrink=0.5)

    # for index in range(len(fits_files)):
    # Spectro view
    hdul = fits.open(join(file_dir, fits_files[index]))
    # hdul = fits.open("emm_data/0224/emm_emr_l2_20210224t095920_0011_r_v00-01.fits")
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
        if data[i][7] > 0.001 and (abs(data[i][21]) > 0.001 or abs(data[i][22]) > 0.001) and data[i][14]:
            dstgui.append(1)
        else:
            dstgui.append(0)

    # max_temps = [reading[7] for reading in data]
    # lats = [np.deg2rad(reading[21]) for reading in data]
    # longs = [np.deg2rad(reading[22]) for reading in data]
    max_temps = [data[i][7] for i in range(len(data)) if dstgui[i] == 1]
    lats = [data[i][21] for i in range(len(data)) if dstgui[i] == 1]
    longs = [data[i][22] for i in range(len(data)) if dstgui[i] == 1]

    # Surround OG map with 6 copies for boundary condition
    p_lats = lats + lats + lats + (-np.array(lats) + PI).tolist() + (-np.array(lats) + PI).tolist() + \
             (-np.array(lats) - PI).tolist() + (-np.array(lats) - PI).tolist()
    p_longs = longs + (np.array(longs) + 2 * PI).tolist() + (np.array(longs) - 2 * PI).tolist() + \
              (-np.array(longs) + PI).tolist() + (-np.array(longs) - PI).tolist() + \
              (-np.array(longs) + PI).tolist() + (-np.array(longs) - PI).tolist()
    p_max_temps = max_temps + max_temps + max_temps + max_temps + max_temps + max_temps + max_temps

    ax = fig.add_subplot(1, 2, 1)

    # Initial interpolation
    xy_bound = [-PI, PI, (-PI / 2), (PI / 2)]
    XX, YY = np.mgrid[xy_bound[0]:xy_bound[1]:30j, xy_bound[2]:xy_bound[3]:2j]
    ZZ = spint.griddata(
        list(zip(p_longs, p_lats)),
        p_max_temps,
        (XX, YY), method='cubic'
    )
    # print(ZZ)
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
    mn = 100  # min(max_temps)
    mx = 300  # max(max_temps)
    ax.imshow(ZZ.T, extent=(xy_bound[0], xy_bound[1], xy_bound[2], xy_bound[3]), origin='lower',
              cmap='plasma', vmin=mn, vmax=mx)
    ax.scatter(p_longs, p_lats, c=p_max_temps, cmap='plasma', marker='o',
               vmin=mn, vmax=mx, edgecolor='white', s=8, linewidth=0.5)
    ax.set_aspect(1)
    ax.set_title("EMIRS max brightness temperature" + "\n" +
                 "Equirectangular projection with cubic interpolation" + "\n" + utc + " -\n" + utc2, pad=12)
    ax.set_xlim(xy_bound[0], xy_bound[1])
    ax.set_ylim(xy_bound[2], xy_bound[3])
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    rect = patches.Rectangle((-PI, -PI / 2), 2 * PI, PI, linewidth=1,
                             linestyle='dotted', edgecolor='w', facecolor='none')
    ax.add_patch(rect)

    if PI == 180:
        alpha = 0.025
    else:
        alpha = 0.025 * 180 / PI
    # for i in range(len(longs)):
    #     if longs[i] < 0:
    #         longs[i] = longs[i] + 2 * PI

    pp_longs = [p for p in p_longs if (p < 3 * PI) and (p > -3 * PI)]
    pp_lats = [p for p in p_lats if (p < 3 * PI / 2) and (p > -3 * PI / 2)]
    points = np.array(list(zip([float(i) for i in pp_longs], [float(i) for i in pp_lats])))
    # alpha_shape(points, ax)
    # ashape = alpha_shape(points, alpha=40, only_outer=True)
    # print(ashape)
    # for i, j in ashape:
    #     ax.plot(points[[i, j], 0], points[[i, j], 1])

    # for i in range(len(longs)):
    #     if longs[i] > PI:
    #         longs[i] = longs[i] - 2 * PI

    # print(points)
    # ashape = alpha_shape(points, alpha=0.25, only_outer=True)
    # print(ashape)

    # np.random.seed(0)
    # x = 100.0 * np.random.rand(2000)
    # y = 100.0 * np.random.rand(2000) - 50.0
    # inside = ((x ** 2 + y ** 2 > 1.0) & ((x - 3) ** 2 + y ** 2 > 1.0) & ((x - 1.5) ** 2 + y ** 2 > 0.09))
    # points = np.vstack([x[inside], y[inside]]).T
    # print(points)

    # Computing the alpha shape
    # edges = alpha_shape(points, alpha=0.25, only_outer=True)
    # print(edges)
    # for i, j in edges:
    #     ax.plot(points[[i, j], 0], points[[i, j], 1])

    ashape = alphashape.alphashape(points, alpha)
    # poly = shapely.geometry.MultiPolygon(ashape)
    # print(ashape.propeties)
    # print(len(ashape.geoms))
    # patch = PathPatch(
    #     Path.make_compound_path(
    #         Path(np.asarray(ashape.geoms.exterior.coords)[:, :2]),
    #         *[Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors]
    #     )
    # )
    # for point in points:
    #     print(point, ashape.contains(point))
    # mycoordslist = [poly.exterior.coords for poly in ashape.geoms]
    # mycoordslist = [poly.exterior.coords for poly in ashape.geoms]
    # ax.plot(mycoordslist)
    # ashape
    # print(ashape.geoms)
    # print(ashape.convex_hull)
    # ax.plot(*ashape.exterior.xy)
    # print(ashape.convex_hull)
    # ax.add_patch(patches.Polygon(ashape))
    ax.add_patch(PolygonPatch(ashape, alpha=.4, linewidth=0.5,
                 linestyle='-', edgecolor='w', facecolor='none'))

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    res = 1
    xy_bound = [-PI, PI, (-PI / 2), (PI / 2)]
    XX, YY = np.mgrid[xy_bound[0]:xy_bound[1]:(180j * res), xy_bound[2]:xy_bound[3]:(90j * res)]
    ZZ = spint.griddata(
        list(zip(p_longs + add_p_longs, p_lats + add_p_lats)),
        p_max_temps + add_p_max_temps,
        (XX, YY), method='cubic'
    )
    for i in range(len(ZZ[:, 0])):
        for j in range(len(ZZ[0, :])):
            inshape = False
            for poly in ashape.geoms:
                if poly.contains(Point(XX[i, j], YY[i, j])):
                    inshape = True

            print(i, j)
            if not inshape:
                ZZ[i, j] = np.nan

    lats = np.deg2rad(lats)
    longs = np.deg2rad(longs)
    if PI == 180:
        XX = np.deg2rad(XX)
        YY = np.deg2rad(YY)

    norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
    sf = ax.plot_surface(radius * np.cos(YY) * np.cos(XX),
                         radius * np.cos(YY) * np.sin(XX),
                         radius * np.sin(YY),
                         rstride=2, cstride=2,
                         edgecolors=None,
                         facecolors=cm.plasma(norm(ZZ)))

    sc = ax.scatter(np.array(radius * np.cos(lats) * np.cos(longs)),
                    np.array(radius * np.cos(lats) * np.sin(longs)),
                    np.array(radius * np.sin(lats)),
                    c=max_temps, cmap='plasma', marker='o',
                    edgecolor='white', s=1, linewidth=0.5, vmin=mn, vmax=mx)

    ax.set_box_aspect((1, 1, 1))
    ax.set_title("EMIRS max brightness temperature\non unit sphere" + "\n" + utc + " -\n" + utc2, pad=-20)
    ax.set_xticks(np.arange(-radius, radius, 0.5))
    ax.set_yticks(np.arange(-radius, radius, 0.5))
    ax.set_zticks(np.arange(-radius, radius, 0.5))
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_zlim(-radius, radius)
    ax.view_init(azim=azim, elev=elev)
    cb = fig.colorbar(sc, ax=ax, shrink=0.5)
    cb.set_label("Max brightness temperature (K)", labelpad=10)

    # ax = fig.add_subplot(1, 1, 1)
    # max_temps = np.array([reading[7] for reading in data])
    # lats = np.array([reading[21] for reading in data])
    # longs = np.array([reading[22] for reading in data])
    # ax.scatter(np.array(lats), np.array(longs), c=max_temps, cmap='plasma', alpha=1, s=np.array(dstgui)*10)
    # # for i in range(hdr['NAXIS2']):
    # #     ax.scatter(data[i][21], data[i][22], alpha=0.5, s=1.5)

    c_fig = plt.gcf()
    plt.show()  # block=True
    c_fig.set_size_inches((11, 8.5), forward=False)
    c_fig.savefig("9_{}".format(index), dpi=300)
    plt.clf()

# Find time
# Time in utc field diverges 3 mins from labeled time
# Time in sclk using conversion diverges 5 days 1 hour from labeled time
# https://stackoverflow.com/a/46493066
# print((Time(2000, format='jyear') + TimeDelta((data[0][0] + data[0][1] * 1.0 / 65536)*u.s)).iso)
# print((Time(2000, format='jyear') + TimeDelta((data[hdr['NAXIS2'] - 1][0] + data[hdr['NAXIS2'] - 1][1] * 1.0 / 65536)*u.s)).iso)
# print(data[0][9])
# print(data[hdr['NAXIS2'] - 1][9])

# Understanding FOV and IFOV
# https://www.thermal-imaging-blog.com/thermography-terms-explained-fov-ifov-ifovmeasurement-on-your-infrared-camera/
