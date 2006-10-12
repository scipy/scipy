import numpy as N
import survival
import model
import tempfile, shutil

class DiscreteRV:

    """
    A simple little class for working with discrete random vectors.
    """

    def __init__(self, x, w=None):
        self.x = N.squeeze(x)
        if self.x.shape == ():
            self.x = N.array([self.x])
        self.n = self.x.shape[0]
        if w is None:
            w = N.ones(self.n, N.float64)
        else:
            if w.shape[0] != self.n:
                raise ValueError, 'incompatible shape for weights w'
            if N.any(N.less(w, 0)):
                raise ValueError, 'weights should be non-negative'
            self.w = w / w.sum()

    def mean(self, f=None):
        if f is None:
            fx = self.x
        else:
            fx = f(self.x)
        return (fx * self.w).sum()

    def cov(self):
        mu = self.moment()
        dx = self.x - N.multiply.outer(mu, self.x.shape[1])
        return N.dot(dx, N.transpose(dx))

class Observation(survival.RightCensored):

    def __getitem__(self, item):
        if self.namespace is not None:
            return self.namespace[item]
        else:
            return getattr(self, item)

    def __init__(self, time, delta, namespace=None):
        self.namespace = namespace
        survival.RightCensored.__init__(self, time, delta)

    def __call__(self, formula, time=None, **extra):
        return formula(namespace=self, time=time, **extra)

class ProportionalHazards(model.LikelihoodModel):

    def __init__(self, subjects, formula, time_dependent=False):
        self.subjects, self.formula = subjects, formula
        self.time_dependent = time_dependent
        self.initialize(self.subjects)
        

    def initialize(self, subjects):

        self.failures = {}
        for i in range(len(subjects)):
            s = subjects[i]
            if s.delta:
                if not self.failures.has_key(s.time):
                    self.failures[s.time] = [i]
                else:
                    self.failures[s.time].append(i)
        
        self.failure_times = self.failures.keys()
        self.failure_times.sort()

    def cache(self):
        if self.time_dependent:
            self.cachedir = tempfile.mkdtemp()

        self.design = {}
        self.risk = {}
        first = True
        
        for t in self.failures.keys():
            if self.time_dependent:
                d = N.array([s(self.formula, time=t)
                             for s in self.subjects]).astype('<f8')
                dshape = d.shape
                dfile = file(tempfile.mkstemp(dir=self.cachedir)[1], 'w')
                d.tofile(dfile)
                dfile.close()
                del(d)
                self.design[t] = N.memmap(dfile.name,
                                          dtype=N.dtype('<f8'),
                                          shape=dshape)
            elif first:
                d = N.array([s(self.formula, time=t)
                             for s in self.subjects]).astype(N.float64)
                self.design[t] = d
            else:
                self.design[t] = d
            self.risk[t] = N.compress([s.atrisk(t) for s in self.subjects],
                                      N.arange(self.design[t].shape[0]),axis=-1)
    def __del__(self):

        shutil.rmtree(self.cachedir, ignore_errors=True)

    def logL(self, b, ties='breslow'):

        logL = 0
        for t in self.failures.keys():
            fail = self.failures[t]
            d = len(fail)
            risk = self.risk[t]
            Zb = N.dot(self.design[t], b)

            logL += Zb[fail].sum()

            if ties == 'breslow':
                s = N.exp(Zb[risk]).sum()
                logL -= N.log(N.exp(Zb[risk]).sum()) * d
            elif ties == 'efron':
                s = N.exp(Zb[risk]).sum()
                r = N.exp(Zb[fail]).sum()
                for j in range(d):
                    logL -= N.log(s - j * r / d)
            elif ties == 'cox':
                raise NotImplementedError, 'Cox tie breaking method not implemented'
            else:
                raise NotImplementedError, 'tie breaking method not recognized'
        return logL

    def score(self, b, ties='breslow'):

        score = 0
        for t in self.failures.keys():
            fail = self.failures[t]
            d = len(fail)
            risk = self.risk[t]
            Z = self.design[t]

            score += Z[fail].sum()

            if ties == 'breslow':
                w = N.exp(N.dot(Z, b))
                rv = DiscreteRV(Z[risk], w=w[risk])
                score -= rv.mean() * d
            elif ties == 'efron':
                w = N.exp(N.dot(Z, b))
                score += Z[fail].sum()
                for j in range(d):
                    efron_w = w
                    efron_w[fail] -= i * w[fail] / d
                    rv = DiscreteRV(Z[risk], w=efron_w[risk])
                    score -= rv.mean()
            elif ties == 'cox':
                raise NotImplementedError, 'Cox tie breaking method not implemented'
            else:
                raise NotImplementedError, 'tie breaking method not recognized'
        # FIXME: score is an int. it has no shape
        # is it that we shouldn't be using an int above
        # or that we shouldn't be looking at shape here
        if score.shape == ():
            score = N.array([score])
        return score

    def information(self, b, ties='breslow'):

        info = 0
        score = 0
        for t in self.failures.keys():
            fail = self.failures[t]
            d = len(fail)
            risk = self.risk[t]
            Z = self.design[t]

            if ties == 'breslow':
                w = N.exp(N.dot(Z, b))
                rv = DiscreteRV(Z[risk], w=w[risk])
                info += rv.cov()
            elif ties == 'efron':
                w = N.exp(N.dot(Z, b))
                score += Z[fail].sum()
                for j in range(d):
                    efron_w = w
                    efron_w[fail] -= i * w[fail] / d
                    rv = DiscreteRV(Z[risk], w=efron_w[risk])
                    info += rv.cov()
            elif ties == 'cox':
                raise NotImplementedError, 'Cox tie breaking method not implemented'
            else:
                raise NotImplementedError, 'tie breaking method not recognized'
        return score

if __name__ == '__main__':
    import numpy.random as R
    n = 100
    X = N.array([0]*n + [1]*n)
    b = 0.4
    lin = 1. + b * X
    Y = R.standard_exponential((2*n,)) / lin
    delta = R.binomial(1, 0.9, size=(2*n,))

    subjects = [Observation(Y[i], delta[i]) for i in range(2*n)]
    for i in range(2*n):
        subjects[i].X = X[i]

    import formula as F
    x = F.Quantitative('X')
    f = F.Formula(x)

    c = ProportionalHazards(subjects, f)

    c.cache()
    c.newton([0.4])
