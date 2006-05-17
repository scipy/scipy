import numpy as N

class SurvivalTime:
    pass

class RightCensored(SurvivalTime):

    def __init__(self, time, delta):
        self.time, self.delta = time, delta

    def atrisk(self, time):
        return N.less_equal.outer(time, self.time)

class LeftCensored(RightCensored):

    def atrisk(self, time):
        return N.greater_equal.outer(time, self.time)
