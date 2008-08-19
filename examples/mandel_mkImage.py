#!/usr/bin/env python

# Before running this be sure to apply Travis Oliphant's patch to PIL:
# http://www.scipy.org/Cookbook/PIL

from numpy import asarray, concatenate, ogrid, uint8
from PIL import Image


import sys
sys.path.append('../mkufunc')
from fast_vectorize import fast_vectorize


from mandel_c import mandel


@fast_vectorize(int)
def red(i):
    if i == -1: return 0
    return (i * 5) % 256

@fast_vectorize(int)
def green(i):
    if i == -1: return 0
    return (i % 16) * 15

@fast_vectorize(int)
def blue(i):
    if i == -1: return 0
    return 255


w, h = 1200, 900

y, x = ogrid[-1.5:+1.5:h*1j, -2.75:+1.15:w*1j]

mand = mandel(x, y)

r = asarray(red(mand),   dtype=uint8).reshape(h, w, 1)
g = asarray(green(mand), dtype=uint8).reshape(h, w, 1)
b = asarray(blue(mand),  dtype=uint8).reshape(h, w, 1)

a = concatenate((r, g, b), axis=2).reshape(h, w, 3)

im = Image.fromarray(a)
im.save('mandel.png')
