      subroutine mdu
     *     (ek,dmin, v,l, head,last,next, mark)
clll. optimize
c***********************************************************************
c  mdu -- update degrees of uneliminated vertices in ek
c***********************************************************************
      integer  ek, dmin,  v(1), l(1),  head(1), last(1), next(1),
     *   mark(1),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,
     *   blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag
      tag = mark(ek) - last(ek)
c
c----for each vertex vi in ek
      i = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 11
      do 10 ilp=1,ilpmax
        i = l(i)
        vi = v(i)
        if (last(vi))  1, 10, 8
c
c------if vi neither prototype nor duplicate vertex, then merge elements
c------to compute degree
   1      tag = tag + 1
          dvi = last(ek)
c
c--------for each vertex/element vs/es in element list of vi
          s = l(vi)
   2      s = l(s)
          if (s.eq.0)  go to 9
            vs = v(s)
            if (next(vs).lt.0)  go to 3
c
c----------if vs is uneliminated vertex, then tag and adjust degree
              mark(vs) = tag
              dvi = dvi + 1
              go to 5
c
c----------if es is active element, then expand
c------------check for outmatched vertex
   3          if (mark(es).lt.0)  go to 6
c
c------------for each vertex vb in es
              b = es
              blpmax = last(es)
              do 4 blp=1,blpmax
                b = l(b)
                vb = v(b)
c
c--------------if vb is untagged, then tag and adjust degree
                if (mark(vb).ge.tag)  go to 4
                  mark(vb) = tag
                  dvi = dvi + 1
   4            continue
c
   5        go to 2
c
c------else if vi is outmatched vertex, then adjust overlaps but do not
c------compute degree
   6      last(vi) = 0
          mark(es) = mark(es) - 1
   7      s = l(s)
          if (s.eq.0)  go to 10
            es = v(s)
            if (mark(es).lt.0)  mark(es) = mark(es) - 1
            go to 7
c
c------else if vi is prototype vertex, then calculate degree by
c------inclusion/exclusion and reset overlap count
   8      evi = last(vi)
          dvi = last(ek) + last(evi) + mark(evi)
          mark(evi) = 0
c
c------insert vi in appropriate degree list
   9    next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        if (next(vi).gt.0)  last(next(vi)) = vi
        if (dvi.lt.dmin)  dmin = dvi
c
  10    continue
c
  11  return
      end
