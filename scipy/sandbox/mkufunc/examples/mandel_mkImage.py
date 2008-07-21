#!/usr/bin/env python

# Before running this be sure to apply Travis Oliphant's patch to PIL:
# http://www.scipy.org/Cookbook/PIL

import numpy
from PIL import Image

from mkufunc.api import mkufunc

from mandel_c import mandel


@mkufunc(int)
def color(i):
    return (i * 10) % 256
    
    n = i % 3
    if n == 0:
        c = (255, 127, 128)
    elif n == 1:
        c = (128, 255, 0)
    else:
        c = (0, 128, 255)

    return c[0] + (c[1] + c[2]*256)*256



w, h = 1024, 768

x, y = numpy.ogrid[-2.5:+1.5:w*1j, -1.5:+1.5:h*1j]

img = mandel(x, y)

print img.dtype
print img.shape

img = color(img)
img.dtype = numpy.uint8
img = img.reshape(h, w, 4)
    
print img.dtype
print img.shape

pilImage = Image.fromarray(img)
pilImage.save('mandel.png')
