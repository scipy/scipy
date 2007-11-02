import numpy as N
import volumeInput as V
import EdgeFilters as E
import edgeSegmenter as S

slice = V.getFilteredSlice('slice112.raw')

edges, groups = S.Canny(slice)

ShenCastanLow = 0.3
b             = 0.8
window        = 7
lowThreshold  = 220 + 2048
highThreshold = 600 + 2048

edges, groups = S.ShenCastan(slice)

edges, groups = E.ShenCastanEdges(ShenCastanLow, b, window, lowThreshold, highThreshold, slice)

