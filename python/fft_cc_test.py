import numpy as np
import scipy.io
import pathlib
import matplotlib.pyplot as plt
from kgpy.plot import CubeSlicer
from kgpy.img.coalignment import image_coalignment as img_align
import scipy.optimize
import scipy.signal
import scipy.fft as fft
import time

file_path3 = pathlib.Path(__file__).parents[1] / 'moses_super.sav'
sav3 = scipy.io.readsav(file_path3)
p = sav3.super_plus
m = sav3.super_minus
z = sav3.super_zero
imgshp = p.shape

pz = p - z
mz = m - z
scale = .5
fig, ax = plt.subplots()
ax.imshow(pz, vmin=-scale, vmax=scale)
fig, ax = plt.subplots()
ax.imshow(mz, vmin=-scale, vmax=scale)

start = time.time()
pz_normalized = (pz - pz.mean(axis=1, keepdims=True)) / pz.std(axis=1, keepdims=True)
mz_normalized = (mz - mz.mean(axis=1, keepdims=True)) / mz.std(axis=1, keepdims=True)

f_pz = fft.fft(pz_normalized, imgshp[1] * 2 - 1, axis=-1)
f_mz = fft.fft(mz_normalized, imgshp[1] * 2 - 1, axis=-1)
cc = fft.ifft(f_pz * f_mz.conj(), imgshp[1] * 2 - 1, axis=-1) / imgshp[1]
cc = cc.real
cc = np.roll(cc, imgshp[1] - 1,axis=1)
tcor = cc.mean(axis=0)
print(time.time()-start)

start = time.time()
cc = []
for i in range(imgshp[0]):
    cc.append(img_align.normalized_cc(pz[i], mz[i], mode='full'))
cc = np.array(cc)
mcor = np.mean(cc, axis=0)
print(time.time()-start)


fig, ax = plt.subplots()
ax.plot(mcor)
ax.plot(tcor, color='r')

plt.show()
