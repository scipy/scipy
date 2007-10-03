import numpy as N

class SurvivalTime:
    def __init__(self, time, delta):
        self.time, self.delta = time, delta

    def atrisk(self, time):
        raise NotImplementedError

class RightCensored(SurvivalTime):

    def atrisk(self, time):
        return N.less_equal.outer(time, self.time)

class LeftCensored(SurvivalTime):

    def atrisk(self, time):
        return N.greater_equal.outer(time, self.time)
