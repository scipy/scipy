import numpy as N

class survival_time:
    def __init__(self, time, delta):
        self.time, self.delta = time, delta

    def atrisk(self, time):
        raise NotImplementedError

class right_censored(survival_time):

    def atrisk(self, time):
        return N.less_equal.outer(time, self.time)

class left_censored(survival_time):

    def atrisk(self, time):
        return N.greater_equal.outer(time, self.time)
