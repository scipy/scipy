"""The following module describes how to use 
the average function of the numpy class"""

# From the existing set of classes we import numpy
import numpy as np
# We can rename numpy as np. Once may use anything

# -----METHOD 1-----
"""What is Average?
average = [Sum of all the elements of array]/[Number of elements in array]"""

# Let us declare an array of integers as follows
data = [1,2,3,4]

# Let us store the average in a variable 'av'
av = np.average(data)

print('The average using Method 1 is:', av)

# -----METHOD 2 - Weights -----
"""What is Weighted Average?
Given data=[1,2,3,4] and weight=[0.1,0.2,0.3,0.4]
average = [(1*0.2)+(2*0.2)+(3*0.3)+(4*0.4)]/[0.1+0.2+0.3+0.4]
"""
# We are going to use the same data as in METHOD 1
# Here we are going to have weight as [0.1,0.2,0.3,0.4]

weight_array = [0.1,0.2,0.3,0.4]

"""SYNTAX for Weighted Average
Weighted_Average = numpy.average(<data_name>, weights = <weight_array>)
Remember, here, 'weights' is a Python Keyword"""

weighted_av = np.average(data, weights = weight_array)

print('The Weighted Average using Method 2 is: ', weighted_av)
