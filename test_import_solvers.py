#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 09:42:27 2018

@author: jonaswillmann
"""

import solvers
import numpy as np
import sys


def input_matrix(fname):

    try:
        mm = open(fname, "r")
    except FileNotFoundError:
        print("Input file {} not found".format(fname))
        sys.exit(1)
    else:
        matrix = []
        for line in open(fname).readlines():
            matrix.append(line.split())
        nn = int(matrix[0][0])
        aa = np.array(matrix[1:nn+1], dtype=float)
        bb = np.array(matrix[-1], dtype=float)
        return aa, bb


aa, bb = input_matrix("lindqd.in")
xx_gauss = solvers.gaussian_eliminate(aa, bb)
print(xx_gauss)
if xx_gauss is None:
    f = open('linsolve.out', 'w')
    f.write('ERROR: LINDEP')
    f.close()
else:
    np.savetxt('linsolve.out', xx_gauss)
