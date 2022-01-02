import matplotlib
from astropy.io import fits
from astropy.time import TimeDelta, Time
from astropy import units as u
import emm_api as ep
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm
import numpy as np
from os import listdir
from os.path import isfile, join
import pickle5 as pickle

# ep.get_sdc_token(user_name='QuanNham', password='irQPmFciRTa6ajP?')
# print(ep.USER_TOKEN)

# #ep.emm_search_and_download("emm_data", instrument='emr', level="l2", latest='true', start_date='2021-02-15', end_date='2021-11-01')
# em_q = ep.query_science_files(instrument="emr", latest='true', level='l2', start_date='2021-02-20', end_date='2021-03-16')
# ex_q = ep.query_science_files(instrument="exi", latest='true', level='l2a', start_date='2021-02-20', end_date='2021-03-16')
# for query in ex_q:
#     if query['descriptor'] != 'f635':
#         # print(query)
#         ex_q.remove(query)

# print(em_q)
# print(ex_q)

# Saving query
# https://stackoverflow.com/a/31891445
filename = 'emir_query.pk'
# with open(filename, 'wb') as fi:
#     pickle.dump(em_q, fi)
with open(filename, 'rb') as fi:
    em_q = pickle.load(fi)
    # print(em_q)

filename = 'exi_query.pk'
# with open(filename, 'wb') as fi:
#     pickle.dump(ex_q, fi)
with open(filename, 'rb') as fi:
    ex_q = pickle.load(fi)
    # print(ex_q)

i = 0
while i < len(ex_q):
    if ex_q[i]['descriptor'] != 'f635':
        # print(ex_q[i])
        ex_q.pop(i)
        try:
            while ex_q[i]['descriptor'] != 'f635':
                # print(i, len(ex_q))
                # print(ex_q[i])
                ex_q.pop(i)
                if i == len(ex_q):
                    break
        except Exception:
            continue
    i = i + 1
# [print(query['descriptor']) for query in ex_q]

# print(em_q[0])
# # print(em_q[0]['file_name'])
# # print(em_q[0]['timetag'].rsplit('+')[0] + '.000')
# print(ex_q[0])
# t1 = Time(em_q[0]['timetag'].partition('+')[0] + '.000', format='iso')
# t2 = Time(ex_q[0]['timetag'].partition('+')[0] + '.000', format='iso')
# print(t1.to_value('unix'))
# print((t1 - t2).to_value('hr'))

em_unixes = []
ex_unixes = []

for i in range(len(em_q)):
    em_unixes.append(((Time(em_q[i]['timetag'].partition('+')[0] + '.000', format='iso').to_value('unix')), i))
    # print(em_unixes[-1])

for i in range(len(ex_q)):
    ex_unixes.append(((Time(ex_q[i]['timetag'].partition('+')[0] + '.000', format='iso').to_value('unix')), i))
    # print(ex_unixes[-1])

em_unixes.sort()
ex_unixes.sort()

# 24-02 and 27-02 have 5 good measurements each
pairs = []
j = 0
i = 0
dist = 1015
while i < len(em_unixes):
    while (abs(ex_unixes[j][0] - em_unixes[i][0]) < dist or ex_unixes[j][0] < em_unixes[i][0]) and j < len(ex_unixes):
        if abs(ex_unixes[j][0] - em_unixes[i][0]) < dist:
            pairs.append((i, j))
            # print(em_unixes[i][1], ex_unixes[j][1],
            #     em_q[em_unixes[i][1]]['timetag'].partition('+')[0] +
            #         ' ' + ex_q[ex_unixes[j][1]]['timetag'].partition('+')[0])
            # print(em_unixes[i][1], ex_unixes[j][1], ex_unixes[j][0] - em_unixes[i][0])
        j = j + 1
    i = i + 1
    if j == len(ex_unixes):
        break

# 23, 126 and 23, 83 is good
# 2021-02-25 20:03:27 2021-02-25 20:00:17
# 2021-02-25 20:03:27 2021-02-25 20:20:22
# ep.get_sdc_token(user_name='QuanNham', password='irQPmFciRTa6ajP?')
# ep.download_from_query([em_q[23]], 'emm_data')
# print("\n\n")
# print(em_q[23]['timetag'].partition('+')[0] + ' ' + ex_q[126]['timetag'].partition('+')[0])
# print(em_q[23]['timetag'].partition('+')[0] + ' ' + ex_q[83]['timetag'].partition('+')[0])

# print(pairs)
# [print(em_q[em_unixes[i][1]]['timetag'].partition('+')[0] + ' ' + ex_q[ex_unixes[j][1]]['timetag'].partition('+')[0]) for i, j in pairs]

# Get files
file_dir = "emm_data/0224"
file_character = "emr"
fits_files = [filename for filename in listdir(file_dir) if (isfile(join(file_dir, filename))
                                                             and file_character in filename)]
file_character = "exi"
img_files = [filename for filename in listdir(file_dir) if (isfile(join(file_dir, filename))
                                                            and file_character in filename)]
index = 2
# print(fits_files)
# print(img_files)

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


# hdul = fits.open("emm_emr_l2_20210305t094430_0017_r_v00-01.fits")
# hdul = fits.open("emm_data/emm_emr_l2_20210308t045230_0019_r_v00-01.fits")
# hdul = fits.open("emm_data/emm_emr_l2_20210306t192330_0018_r_v00-01.fits")
# hdul = fits.open("emm_data/emm_emr_l2_20210311t000756_0021_r_v00-01.fits")

# Image view
hdul = fits.open(join(file_dir, img_files[index]))
# hdul = fits.open("emm_data/0224/emm_exi_l2a_20210224T095615_0011_xos1_f635_r_v03-01.fits")
hdul.verify('fix')
utc = hdul[0].header['DATE-OBS']
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
            if (img[4][i][j] < max_altitude):
                continue
                # img[0][i][j] = (img[0][i][j] - mn) * 1.0 / rg
            else:
                for k in range(len(interest)):
                    img[k][i][j] = np.nan

            # mask[i][j] = 1
        # else:
        # img[0][i][j] = 0

hdul.close()
# mn = min(img[0])
# mx = max(img[0])
# range = mx - mn
# img[0] = [(i - mn) * 1.0 / range for i in img[0]]
# print(img[0])

fig = plt.figure()
fig.canvas.manager.full_screen_toggle()
radius = 1
azim = 105
elev = 10

# Plot help
# https://stackoverflow.com/a/42927880
ax = fig.add_subplot(1, 2, 1, projection='3d')
imgs = np.array(img[0])
lats = np.array(np.deg2rad(img[2]))
longs = np.array(np.deg2rad(img[3]))
norm = matplotlib.colors.Normalize(vmin=mn, vmax=mx)
sf = ax.plot_surface(radius * np.cos(lats) * np.cos(longs),
                     radius * np.cos(lats) * np.sin(longs),
                     radius * np.sin(lats),
                     rstride=1, cstride=1,
                     facecolors=cm.gray(norm(imgs)))
                     # facecolors=cm.gray(imgs))

ax.set_box_aspect((1, 1, 1))
ax.set_title("EXI image @ f635 on unit sphere" + "\n" + utc, pad=-10)
ax.set_xticks(np.arange(-radius, radius, 0.5))
ax.set_yticks(np.arange(-radius, radius, 0.5))
ax.set_zticks(np.arange(-radius, radius, 0.5))
ax.set_xlim(-radius, radius)
ax.set_ylim(-radius, radius)
ax.set_zlim(-radius, radius)
ax.view_init(azim=azim, elev=elev)
cb = fig.colorbar(cm.ScalarMappable(cmap=cm.gray, norm=norm), ax=ax, shrink=0.5)
cb.set_label("Calibrated DN", labelpad=10)
# fig.colorbar(cm.ScalarMappable(cmap=cm.gray, norm=imgs).set_array([]), ax=ax, shrink=0.5)

# axImg = fig.add_subplot(1, 2, 1)
# XX, YY = np.meshgrid(np.linspace(0, w - 1, w), np.linspace(0, h - 1, h))
# axImg.imshow(img[0], cmap='gray', vmin=mn, vmax=mx)
# axImg.imshow(img[4], cmap='gray')
# axImg.imshow(np.sin(img[2]), cmap='gray')
# print(min([img[4][i][j] for i in range(h) for j in range(w) if mask[i][j] == 1]))

# Spectro view
hdul = fits.open(join(file_dir, fits_files[index]))
# hdul = fits.open("emm_data/0224/emm_emr_l2_20210224t095920_0011_r_v00-01.fits")
hdul.verify('fix')
hdr = hdul[1].header
# print(repr(hdr))
data = hdul[1].data
hdul.close()
# print(data)

# Field 3 is calibrated radiance
# print(data[5][3])

# Field 7 is max brightness temp
# Fields 21 and 22 are latitudes and longitudes of IFOV
ax = fig.add_subplot(1, 2, 2, projection='3d')
dstgui = []
for i in range(len(data)):
    if data[i][7] > 0.001 and (abs(data[i][21]) > 0.001 or abs(data[i][22]) > 0.001):
        dstgui.append(1)
    else:
        dstgui.append(0)

# max_temps = [reading[7] for reading in data]
# lats = [np.deg2rad(reading[21]) for reading in data]
# longs = [np.deg2rad(reading[22]) for reading in data]
max_temps = [data[i][7] for i in range(len(data)) if dstgui[i] == 1]
lats = [np.deg2rad(data[i][21]) for i in range(len(data)) if dstgui[i] == 1]
longs = [np.deg2rad(data[i][22]) for i in range(len(data)) if dstgui[i] == 1]

sc = ax.scatter(np.array(radius * np.cos(lats) * np.cos(longs)),
                np.array(radius * np.cos(lats) * np.sin(longs)),
                np.array(radius * np.sin(lats)),
                c=max_temps, cmap='plasma', marker='o')

# ax.plot(np.array(radius * np.cos(lats) * np.cos(longs)),
#            np.array(radius * np.cos(lats) * np.sin(longs)),
#            np.array(radius * np.sin(lats)))

ax.set_box_aspect((1, 1, 1))
ax.set_title("EMIR max brightness on unit sphere" + "\n" + utc, pad=-10)
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

# plt.show() # block=True
plt.savefig("6_{}".format(index), dpi=200)

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
