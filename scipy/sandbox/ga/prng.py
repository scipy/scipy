""" Package-global pseudo-random number generator.

This global is a transitional hack from the old code. Ideally, each run of a GA
should control its own in order to allow multiple concurrent runs. However, we
are transitioning from an older implementation that used a really global PRNG.
"""

from numpy.random import RandomState


class GAPRNG(RandomState):
    """ PRNG for the GA package.

    In addition to all of the functionality derived from RandomState, we also
    store the seed values that were used.
    """

    def seed(self, seed=None):
        """ Seed the generator.

        seed can be an integer, an array (or other sequence) of integers of any
        length, or None. If seed is None, then RandomState will try to read data
        from /dev/urandom (or the Windows analogue) if available or seed from
        the clock otherwise.
        """
        RandomState.seed(self, seed)
        self.initial_seed = seed

    def choice(self, seq):
        """ Randomly and uniformly select an item from a sequence.
        """
        i = self.randint(len(seq))
        return seq[i]


prng = GAPRNG()
