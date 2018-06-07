#!/usr/local/bin/python3
"""
RunningStat class and some functions for limiting and checking values.

The MIT License (MIT)

Copyright (c) 2018 Justin Ng

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import math
import numpy as np
from collections import deque
from itertools import islice

'''
This class is useful for keeping track of the value of a variable as well as the
mean, variance, and standard deviation. Stats are calculated iteratively so it is not
necessary to retain all past values. The stats come in two flavors: total and windowed.
Total mean, variance, and standard deviation are computed over all past values of the
variable (unless you clear()). Windowed calculates over a fixed number of past values set
by window. The median is calculated over a separate window defined by median_window.

References:
    http://www.taylortree.com/2010/11/running-variance.html
    https://www.johndcook.com/blog/standard_deviation/
'''
class RunningStat:
    def __init__(self, window, median_window = 30):
        self.window = window
        self.median_window = median_window
        # Keep a FIFO of the past values
        self.vals = deque((self.window) * [0.0], self.window)
        self.n = 0
        self._mean = 0.0
        self.square = 0.0
        self.win_mean = 0.0
        self.powsumavg = 0.0

    # Clears the total stats calculation but has no effect on windowed stats
    def clear(self):
        self.n = 0

    def push(self, val):
        val = float(val)
        self.n = self.n + 1
        if (self.n == 1):
            self._mean = val
            self.square = 0.0
        else:
            new_mean = self._mean + (val - self._mean)/self.n
            self.square = self.square + (val - self._mean)*(val - new_mean)
            # set up for next iteration
            self._mean = new_mean

        self.win_mean = self.win_mean + ((val - self.vals[self.window - 1]) / self.window)
        newamt = val
        oldamt = self.vals[self.window - 1]
        self.powsumavg = self.powsumavg + (((newamt * newamt) - (oldamt * oldamt)) / self.window)
        self.vals.appendleft(val)

    def numDataValues(self):
        return self.n

    def curVal(self):
        return self.vals[0]

    def mean(self):
        return self._mean

    def var(self):
        return  self.square/(self.n - 1) if (self.n > 1) else 0.0

    def stdDev(self):
        return math.sqrt(self.var())

    def winMean(self):
        return self.win_mean

    def winVar(self):
        winvar = (self.powsumavg * self.window - self.window * self.win_mean * self.win_mean) / self.window
        if (winvar < 0):
            winvar = 0
        return winvar

    def winStdDev(self):
        return math.sqrt(self.winVar())

    def winMedian(self):
        return np.median(list(islice(self.vals, 0, self.median_window)))

def limitValue(value, min_val = None, max_val = None):
    if ((max_val != None) and (value > max_val)):
        return max_val
    if ((min_val != None) and (value < min_val)):
        return min_val
    return value

def checkValue(value, min_val, max_val):
    if (value > max_val):
        return 1
    if (value < min_val):
        return -1
    return 0

