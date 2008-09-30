# Copyright Anne M. Archibald 2008
# Released under the scipy license
import numpy as np
from heapq import heappush, heappop

def distance_p(x,y,p=2):
    if p==np.inf:
        return np.amax(np.abs(y-x),axis=-1)
    elif p==1:
        return np.sum(np.abs(y-x),axis=-1)
    else:
        return np.sum(np.abs(y-x)**p,axis=-1)
def distance(x,y,p=2):
    if p==np.inf or p==1:
        return distance_p(x,y,p)
    else:
        return distance_p(x,y,p)**(1./p)


class KDTree(object):
    """kd-tree for quick nearest-neighbor lookup

    This class provides an index into a set of k-dimensional points
    which can be used to rapidly look up the nearest neighbors of any
    point. 

    The algorithm used is described in Maneewongvatana and Mount 1999. 
    The general idea is that the kd-tree is a binary trie, each of whose
    nodes represents an axis-aligned hyperrectangle. Each node specifies
    an axis and splits the set of points based on whether their coordinate
    along that axis is greater than or less than a particular value. 

    During construction, the axis and splitting point are chosen by the 
    "sliding midpoint" rule, which ensures that the cells do not all
    become long and thin. 

    The tree can be queried for the r closest neighbors of any given point 
    (optionally returning only those within some maximum distance of the 
    point). It can also be queried, with a substantial gain in efficiency, 
    for the r approximate closest neighbors.

    For large dimensions (20 is already large) do not expect this to run 
    significantly faster than brute force. High-dimensional nearest-neighbor
    queries are a substantial open problem in computer science.
    """

    def __init__(self, data, leafsize=10):
        """Construct a kd-tree.

        Parameters:
        ===========

        data : array-like, shape (n,k)
            The data points to be indexed. This array is not copied, and
            so modifying this data will result in bogus results.
        leafsize : positive integer
            The number of points at which the algorithm switches over to
            brute-force.
        """
        self.data = np.asarray(data)
        self.n, self.k = np.shape(self.data)
        self.leafsize = int(leafsize)
        if self.leafsize<1:
            raise ValueError("leafsize must be at least 1")
        self.maxes = np.amax(self.data,axis=0)
        self.mins = np.amin(self.data,axis=0)

        self.tree = self.__build(np.arange(self.n), self.maxes, self.mins)

    class node(object):
        pass
    class leafnode(node):
        def __init__(self, idx):
            self.idx = idx
    class innernode(node):
        def __init__(self, split_dim, split, less, greater):
            self.split_dim = split_dim
            self.split = split
            self.less = less
            self.greater = greater
    
    def __build(self, idx, maxes, mins):
        if len(idx)<=self.leafsize:
            return KDTree.leafnode(idx)
        else:
            data = self.data[idx]
            #maxes = np.amax(data,axis=0)
            #mins = np.amin(data,axis=0)
            d = np.argmax(maxes-mins)
            maxval = maxes[d]
            minval = mins[d]
            if maxval==minval:
                # all points are identical; warn user?
                return KDTree.leafnode(idx)
            data = data[:,d]

            # sliding midpoint rule; see Maneewongvatana and Mount 1999
            # for arguments that this is a good idea.
            split = (maxval+minval)/2
            less_idx = np.nonzero(data<=split)[0]
            greater_idx = np.nonzero(data>split)[0]
            if len(less_idx)==0:
                split = np.amin(data)
                less_idx = np.nonzero(data<=split)[0]
                greater_idx = np.nonzero(data>split)[0]
            if len(greater_idx)==0:
                split = np.amax(data)
                less_idx = np.nonzero(data<split)[0]
                greater_idx = np.nonzero(data>=split)[0]
            if len(less_idx)==0:
                # _still_ zero? all must have the same value
                assert np.all(data==data[0]), "Troublesome data array: %s" % data
                split = data[0]
                less_idx = np.arange(len(data)-1)
                greater_idx = np.array([len(data)-1])

            lessmaxes = np.copy(maxes)
            lessmaxes[d] = split
            greatermins = np.copy(mins)
            greatermins[d] = split
            return KDTree.innernode(d, split, 
                    self.__build(idx[less_idx],lessmaxes,mins),
                    self.__build(idx[greater_idx],maxes,greatermins))

    def __query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        
        side_distances = [max(0,x[i]-self.maxes[i],self.mins[i]-x[i]) for i in range(self.k)]
        # priority queue for chasing nodes
        # entries are:
        #  minimum distance between the cell and the target
        #  distances between the nearest side of the cell and the target
        #  the head node of the cell
        q = [(distance_p(np.array(side_distances),0),
              tuple(side_distances),                   
              self.tree)]
        # priority queue for the nearest neighbors
        # furthest known neighbor first
        # entries are (-distance**p, i)
        neighbors = []

        if eps==0:
            epsfac=1
        elif p==np.inf:
            epsfac = 1/(1+eps)
        else:
            epsfac = 1/(1+eps)**p

        if p!=np.inf and distance_upper_bound!=np.inf:
            distance_upper_bound = distance_upper_bound**p

        while q:
            min_distance, side_distances, node = heappop(q)
            if isinstance(node, KDTree.leafnode):
                # brute-force
                data = self.data[node.idx]
                a = np.abs(data-x[np.newaxis,:])
                if p==np.inf:
                    ds = np.amax(a,axis=1)
                elif p==1:
                    ds = np.sum(a,axis=1)
                else:
                    ds = np.sum(a**p,axis=1)
                for i in range(len(ds)):
                    if ds[i]<distance_upper_bound:
                        if len(neighbors)==k:
                            heappop(neighbors)
                        heappush(neighbors, (-ds[i], node.idx[i]))
                        if len(neighbors)==k:
                            distance_upper_bound = -neighbors[0][0]
            else:
                # we don't push cells that are too far onto the queue at all,
                # but since the distance_upper_bound decreases, we might get 
                # here even if the cell's too far
                if min_distance>distance_upper_bound*epsfac:
                    # since this is the nearest cell, we're done, bail out
                    break
                # compute minimum distances to the children and push them on
                if x[node.split_dim]<node.split:
                    near, far = node.less, node.greater
                else:
                    near, far = node.greater, node.less

                # near child is at the same distance as the current node
                heappush(q,(min_distance, side_distances, near))

                # far child is further by an amount depending only
                # on the split value
                sd = list(side_distances)
                if p == np.inf:
                    min_distance = max(min_distance, abs(node.split-x[node.split_dim]))
                elif p == 1:
                    sd[node.split_dim] = np.abs(node.split-x[node.split_dim])
                    min_distance = min_distance - side_distances[node.split_dim] + sd[node.split_dim]
                else:
                    sd[node.split_dim] = np.abs(node.split-x[node.split_dim])**p
                    min_distance = min_distance - side_distances[node.split_dim] + sd[node.split_dim]

                # far child might be too far, if so, don't bother pushing it
                if min_distance<=distance_upper_bound*epsfac:
                    heappush(q,(min_distance, tuple(sd), far))

        if p==np.inf:
            return sorted([(-d,i) for (d,i) in neighbors])
        else:
            return sorted([((-d)**(1./p),i) for (d,i) in neighbors])

    def query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf):
        """query the kd-tree for nearest neighbors

        Parameters:
        ===========

        x : array-like, last dimension self.k
            An array of points to query.
        k : integer
            The number of nearest neighbors to return.
        eps : nonnegative float
            Return approximate nearest neighbors; the kth returned value 
            is guaranteed to be no further than (1+eps) times the 
            distance to the real kth nearest neighbor.
        p : float, 1<=p<=infinity
            Which Minkowski p-norm to use. 
            1 is the sum-of-absolute-values "Manhattan" distance
            2 is the usual Euclidean distance
            infinity is the maximum-coordinate-difference distance
        distance_upper_bound : nonnegative float
            Return only neighbors within this distance. This is used to prune
            tree searches, so if you are doing a series of nearest-neighbor
            queries, it may help to supply the distance to the nearest neighbor
            of the most recent point.

        Returns:
        ========
        
        d : array of floats
            The distances to the nearest neighbors. 
            If x has shape tuple+(self.k,), then d has shape tuple if 
            k is one, or tuple+(k,) if k is larger than one.  Missing 
            neighbors are indicated with infinite distances.  If k is None, 
            then d is an object array of shape tuple, containing lists 
            of distances. In either case the hits are sorted by distance 
            (nearest first).
        i : array of integers
            The locations of the neighbors in self.data. i is the same
            shape as d.
        """
        x = np.asarray(x)
        if np.shape(x)[-1] != self.k:
            raise ValueError("x must consist of vectors of length %d but has shape %s" % (self.k, np.shape(x)))
        if p<1:
            raise ValueError("Only p-norms with 1<=p<=infinity permitted")
        retshape = np.shape(x)[:-1]
        if retshape!=():
            if k>1:
                dd = np.empty(retshape+(k,),dtype=np.float)
                dd.fill(np.inf)
                ii = np.empty(retshape+(k,),dtype=np.int)
                ii.fill(self.n)
            elif k==1:
                dd = np.empty(retshape,dtype=np.float)
                dd.fill(np.inf)
                ii = np.empty(retshape,dtype=np.int)
                ii.fill(self.n)
            elif k is None:
                dd = np.empty(retshape,dtype=np.object)
                ii = np.empty(retshape,dtype=np.object)
            else:
                raise ValueError("Requested %s nearest neighbors; acceptable numbers are integers greater than or equal to one, or None")
            for c in np.ndindex(retshape):
                hits = self.__query(x[c], k=k, p=p, distance_upper_bound=distance_upper_bound)
                if k>1:
                    for j in range(len(hits)):
                        dd[c+(j,)], ii[c+(j,)] = hits[j]
                elif k==1:
                    if len(hits)>0:
                        dd[c], ii[c] = hits[0]
                    else:
                        dd[c] = np.inf
                        ii[c] = self.n
                elif k is None:
                    dd[c] = [d for (d,i) in hits]
                    ii[c] = [i for (d,i) in hits]
            return dd, ii
        else:
            hits = self.__query(x, k=k, p=p, distance_upper_bound=distance_upper_bound)
            if k==1:
                if len(hits)>0:
                    return hits[0]
                else:
                    return np.inf, self.n
            elif k>1:
                dd = np.empty(k,dtype=np.float)
                dd.fill(np.inf)
                ii = np.empty(k,dtype=np.int)
                ii.fill(self.n)
                for j in range(len(hits)):
                    dd[j], ii[j] = hits[j]
                return dd, ii
            elif k is None:
                return [d for (d,i) in hits], [i for (d,i) in hits]
            else:
                raise ValueError("Requested %s nearest neighbors; acceptable numbers are integers greater than or equal to one, or None")



