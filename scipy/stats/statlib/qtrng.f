      double precision function prtrng(q, v, r, ifault)

      implicit double precision (a-h, o-z)
c
c        Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2
c	 Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)
c
c        Evaluates the probability from 0 to q for a studentized
c        range having v degrees of freedom and r samples.
c
c        Uses subroutine ALNORM = algorithm AS66.
c
c        Arrays vw and qw store transient values used in the
c        quadrature summation.  Node spacing is controlled by
c        step.  pcutj and pcutk control truncation.
c        Minimum and maximum number of steps are controlled by
c        jmin, jmax, kmin and kmax.  Accuracy can be increased
c        by use of a finer grid - Increase sizes of arrays vw
c        and qw, and jmin, jmax, kmin, kmax and 1/step proportionally.
c
      double precision q,v,r,vw(30), qw(30), pcutj, pcutk, step,
     #    vmax,zero,fifth,half,one,two,cv1,cv2,cvmax,cv(4)
      double precision g, gmid, r1, c, h, v2, gstep, pk1, pk2, gk, pk
      double precision w0, pz, x, hj, ehj, pj
      data pcutj, pcutk, step, vmax /0.0000000001d0, 0.0000000001d0, 0.45d0,
     #   120.0d0/, zero, fifth, half, one, two /0.0d0, 0.2d0, 0.5d0,
     #   1.0d0, 2.0d0/, cv1, cv2, cvmax /0.193064705d0, 0.293525326d0,
     #   0.39894228d0/, cv(1), cv(2), cv(3), cv(4) /0.318309886d0,
     #   -0.268132716d-2, 0.347222222d-2, 0.833333333d-1/
      data jmin, jmax, kmin, kmax/2, 15, 7, 15/
c
c        Check initial values
c
      prtrng = zero
      ifault = 0
      if (v .lt. one .or. r .lt. two) ifault = 1
      if (q .le. zero .or. ifault .eq. 1) goto 99
c
c        Computing constants, local midpoint, adjusting steps.
c
      g = step * r ** (-fifth)
      gmid = half * log(r)
      r1 = r - one
      c = log(r * g * cvmax)
      if(v.gt.vmax) goto 20
c
      h = step * v ** (-half)
      v2 = v * half
      if (v .eq. one) c = cv1
      if (v .eq. two) c = cv2
      if (.not. (v .eq. one .or. v .eq. two)) c = sqrt(v2)
     #    * cv(1) / (one + ((cv(2) / v2 + cv(3)) / v2 + cv(4)) / v2)
      c = log(c * r * g * h)
c
c        Computing integral
c        Given a row k, the procedure starts at the midpoint and works
c        outward (index j) in calculating the probability at nodes
c        symmetric about the midpoint.  The rows (index k) are also
c        processed outwards symmetrically about the midpoint.  The
c        centre row is unpaired.
c
   20 gstep = g
      qw(1) = -one
      qw(jmax + 1) = -one
      pk1 = one
      pk2 = one
      do 28 k = 1, kmax
	gstep = gstep - g
   21   gstep = -gstep
	gk = gmid + gstep
	pk = zero
	if (pk2 .le. pcutk .and. k .gt. kmin) goto 26
	w0 = c - gk * gk * half
        pz = alnorm(gk, .true.)
        x = alnorm(gk - q, .true.) - pz
	if (x .gt. zero) pk = exp(w0 + r1 * log(x))
	if (v .gt. vmax) goto 26
c
	jump = -jmax
   22   jump = jump + jmax
	do 24 j = 1, jmax
	  jj = j + jump
	  if (qw(jj) .gt. zero) goto 23
	  hj = h * j
	  if (j .lt. jmax) qw(jj + 1) = -one
	  ehj = exp(hj)
	  qw(jj) = q * ehj
	  vw(jj) = v * (hj + half - ehj * ehj * half)
c
   23     pj = zero
          x = alnorm(gk - qw(jj), .true.) - pz
	  if (x .gt. zero) pj = exp(w0 + vw(jj) + r1 * log(x))
	  pk = pk + pj
	  if (pj .gt. pcutj) goto 24
	  if (jj .gt. jmin .or. k .gt. kmin) goto 25
   24   continue
   25   h = -h
	if (h .lt. zero) goto 22
c
   26   prtrng = prtrng + pk
	if (k .gt. kmin .and. pk .le. pcutk .and. pk1 .le. pcutk)goto 99
	pk2 = pk1
	pk1 = pk
	if (gstep .gt. zero) goto 21
   28 continue
c
   99 return
      end
c
c
      double precision function qtrng(p, v, r, ifault)
      implicit double precision (a-h, o-z)
c
c        Algorithm AS 190.1  Appl. Statist. (1983) Vol.32, No. 2
c
c        Approximates the quantile p for a studentized range
c        distribution having v degrees of freedom and r samples
c        for probability p, p.ge.0.90 .and. p.le.0.99.
c
c        Uses functions  alnorm, ppnd, prtrng and qtrng0 -
c        Algorithms AS 66, AS 111, AS 190 and AS 190.2
c
      double precision p, v, r, pcut, p75, p80, p90, p99, p995
      double precision p175, one, two, five
      double precision q1, p1, q2, p2, e1, e2
      double precision eps
      data jmax, pcut, p75, p80, p90, p99, p995, p175, one, two, five
     #  /8, 0.001d0, 0.75d0, 0.80d0, 0.90d0, 0.99d0, 0.995d0, 1.75d0,
     #  1.0d0, 2.0d0, 5.0d0/
      data eps/1.0d-4/
c
c        Check input parameters
c
      ifault = 0
      nfault = 0
      if (v .lt. one .or. r.lt. two) ifault = 1
      if (p .lt. p90 .or. p .gt. p99) ifault = 2
      if (ifault .ne. 0) goto 99
c
c        Obtain initial values
c
      q1 = qtrng0(p, v, r, nfault)
      if (nfault .ne. 0) goto 99
      p1 = prtrng(q1, v, r, nfault)
      if (nfault .ne. 0) goto 99
      qtrng = q1
      if (abs(p1-p) .lt. pcut) goto 99
      if (p1 .gt. p) p1 = p175 * p - p75 * p1
      if (p1 .lt. p) p2 = p + (p - p1) * (one - p) / (one - p1) * p75
      if (p2 .lt. p80) p2 = p80
      if (p2 .gt. p995) p2 = p995
      q2 = qtrng0(p2, v, r, nfault)
      if (nfault .ne. 0) goto 99
c
c        Refine approximation
c
      do 14 j = 2, jmax
	p2 = prtrng(q2, v, r, nfault)
	if (nfault .ne. 0) goto 99
	e1 = p1 - p
	e2 = p2 - p
	qtrng = (q1 + q2) / two
	d = e2 - e1
	if (abs(d) .gt. eps) qtrng = (e2 * q1 - e1 * q2) / d
	if(abs(e1) .lt. abs(e2)) goto 12
	q1 = q2
	p1 = p2
   12   if (abs(p1 - p) .lt. pcut * five) goto 99
	q2 = qtrng
   14 continue
c
   99 if (nfault .ne. 0) ifault = 9
      return
      end
c
c
      double precision function qtrng0(p, v, r, ifault)
      implicit double precision (a-h, o-z)
c
c        Algorithm AS 190.2  Appl. Statist. (1983) Vol.32, No.2
c
c        Calculates an initial quantile p for a studentized range
c        distribution having v degrees of freedom and r samples
c        for probability p, p.gt.0.80 .and. p.lt.0.995.
c
c        Uses function ppnd - Algorithm AS 111
c
      double precision p, v, r, q, t, vmax, half, one, four, c1, c2, c3
      double precision c4, c5
      data vmax, half, one, four, c1, c2, c3, c4, c5 / 120.0d0, 0.5d0,
     #  1.0d0, 4.0d0, 0.8843d0, 0.2368d0, 1.214d0, 1.208d0, 1.4142d0/
c
      t=ppnd(half + half * p,ifault)
      if (v .lt. vmax) t = t + (t * t* t + t) / v / four
      q = c1 - c2 * t
      if (v .lt. vmax) q = q - c3 / v + c4 * t / v
      qtrng0 = t * (q * log(r - one) + c5)
      return
      end
      
