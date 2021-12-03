import numpy as np

"""
Created on Fri Dec  3 01:43:13 2021

@author: michele
"""

"""

root = Box()
child = Box(parent, level, edge)

Create a new box,parent is the parent box, level is a value that could be
associated to the dimension of the box and ,  as in the diRECT algorithm
determine how likely a box is selected for further split finally  edge stores
the value of one particular coordinate definingthe boundary between this box
and its "partner"

see fig 2 in the citated paper to better understand how a box is created .
Main rule of thumb when i call a box with a parent, a new child is added in
the parent instance
"""


class Box:
    def __init__(
        self,
        parent=0,
        level=1,
        edge=0,
        children=[],
        splitdim=0,
        xvalues=[],
        fvalues=[],
        cindex=0,
    ):
        self.parent = parent  # default 0(it do not have a parent the root )
        self.level = level
        self.edge = edge  # default 0 root
        self.children = children
        self.splitdim = splitdim
# since i have the splitting dimension i only need a point to identify the box
        self.xvalues = xvalues

        self.fvalues = fvalues
# is false when is inserted in vector children (to avoid infinite recursion )
        if (parent == 0 or parent is False):
            self.cindex = cindex  # identifier of children
        else:
            self.cindex = len(parent.children) + 1  # since ive added a child
            parent.children.append(
                Box(
                    False,
                    level,
                    edge,
                    [],
                    splitdim,
                    xvalues,
                    fvalues,
                    self.cindex)
                    )
            self.children = []
            parent.children[-1].parent = parent  # assing parent to children


# polinomial interpolation
def qfit(xfm, xf0, xfp):
    xm, fm = xfm
    x0, f0 = xf0
    xp, fp = xfp
    assert xp > x0 and x0 > xm and np.isfinite(xm) and np.isfinite(xp)
    cm = fm / ((xm - x0) * (xm - xp))  # coefficients of Lagrange polynomial
    c0 = f0 / ((x0 - xm) * (x0 - xp))
    cp = fp / ((xp - xm) * (xp - x0))
    csum = cm + c0 + cp
    qvalue = (
        lambda x: cm * (x - x0) * (x - xp) + c0 * (x - xm) * (x - xp) + cp * (x - xm) * (x - x0)
    )
    if csum > 0:
        # It's convex
        xmin = (cm * (x0 + xp) + c0 * (xm + xp) + cp * (xm + x0)) / (2 * csum)
        return xmin, qvalue(xmin)
    if (fm == f0) and (fm == fp):
        return x0, f0
    assert (f0 >= fm) or (f0 >= fp)
    if (f0 > fm) and (fp > fm):
        xmin = xm - (x0 - xm)
    else:
        xmin = xp + (xp - x0)

    return xmin, qvalue(xmin)


def init_boxes(f, x0, splits, u, v):
    xstar = x0.copy()
    fstar = f(xstar)
    box = Box()  # root
    n = len(x0)
    if len(splits) != n:
        raise Exception("need one vector of splits per dimension")

    for i in range(n):
        xtmp = xstar.copy()
        xsplit = splits[i]  # vector of s for dimension i
        L = len(xsplit)
        if L < 3:
            raise Exception("there have to be at least 3 evaluations per coordinate")
        if sorted(xsplit) != xsplit:
            raise Exception("splits must be in increasing order")
        fsplit = []
        fmin = float("inf")
        idxmin = 0
        for l in range(L):
            xtmp[i] = xsplit[
                l
            ]  # replace i coordinate of x_tmp with i coordinate of level i
            ftmp = f(xtmp)
            if ftmp < fmin:
                fmin = ftmp
                idxmin = l
            fsplit.append(ftmp)
        box.splitdim = i
        box.xvalues = xsplit.copy()
        box.fvalues = fsplit
        if idxmin == 0:
            raise Exception("function was not finite at any evaluation point")
        # child box creation
        overhang_left = xsplit[0] > u[i]
        idx = int((xsplit[0]) <= u[i])  # 1 if xsplit[0] out bound(<u)
        level = box.level + 1
        while idx <= L:
            if idx == 0:
                # call box and update parent with this children
                Box(box, level, u[i])
            elif idx == L:
                if xsplit[-1] < v[i]:
                    Box(box, level, v[i])  # beyond edge
            else:
                if fsplit[idx] < fsplit[idx + 1]:
                    # additional point with golden rule
                    edge = xsplit[idx] + (xsplit[idx + 1] - xsplit[idx]) * (
                        (np.sqrt(5) - 1) / 2
                    )
                    Box(box, level, edge)
                    Box(box, level + 1, edge)
                else:
                    edge = xsplit[idx + 1] - (xsplit[idx + 1] - xsplit[idx]) * ((np.sqrt(5) - 1) / 2)
                    Box(box, level + 1, edge)
                    Box(box, level, edge)
            idx = idx + 1
            # update xstar
        if fmin < fstar:
            fstar = fmin
            xstar[i] = xsplit[idxmin]
        # pick the next current box for splitting
        if idxmin == 1:
            im = 1
            i0 = 2
            ip = 3
        elif idxmin == L:
            im = L - 2
            i0 = L - 1
            ip = L
        else:
            im = idxmin - 1
            i0 = idxmin
            ip = idxmin + 1
        xqmin, _ = (
            qfit(xsplit[im], fsplit[im]),
            (xsplit[i0], fsplit[i0]),
            (xsplit[ip], fsplit[ip]),
        )
        if xqmin < xsplit[im]:
            assert idxmin == 1
            cindex = 1
        elif xqmin > xsplit[ip]:
            assert idxmin == L
            cindex = len(box.children)
        elif xqmin < xsplit[i0]:
            cindex = 2 * idxmin - 2 + overhang_left
        else:
            cindex = 2 * idxmin - 1 + overhang_left
        box = box.children[cindex]
    return (xstar, fstar, box)
