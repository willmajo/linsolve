# -*- coding: utf-8 -*-

"""Routines for solving a linear system of equations."""
import numpy as np


def gaussian_eliminate(aa, bb):
    """Solves a linear system of equations (Ax = b) by Gauss-elimination

    Args:
        aa: Matrix with the coefficients. Shape: (n, n).
        bb: Right hand side of the equation. Shape: (n,)

    Returns:
        Vector xx with the solution of the linear equation or None
        if the equations are linearly dependent.
    """
    nn = aa.shape[0]
    xx = np.zeros((nn, ), dtype=float)

    # do gauss algorithm and swap rows if necessary
    for i in range(nn):
        for j in range(i+1, nn):
            if abs(aa[i, i]) < abs(aa[j, i]):
                aa[[i, j]] = aa[[j, i]]
                bb[[i, j]] = bb[[j, i]]
            faktor = aa[j, i]/aa[i, i]
            bb[j] -= bb[i]*faktor
            for k in range(i, nn):
                aa[j, k] -= aa[i, k]*faktor

    # check for linear dependency
    if aa[-1, -1] == 0:
        xx = None
    else:
        # do back substition to get result xx
        for i in range(nn-1, -1, -1):
            zaehler = bb[i]
            for j in range(i+1, nn):
                zaehler -= aa[i, j] * xx[j]
            xx[i] = zaehler / aa[i, i]
    return xx
