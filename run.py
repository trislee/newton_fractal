# -*- coding: UTF-8 -*-

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import os

import numpy as np

from newton_fractal import NewtonFractal

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

N_POINTS = 1000
N_ITER = 100
BOUND = 2.0
METHOD = 'fortran'

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == '__main__':

  poly_params = np.array( [ 1, 1, 1, 1, 1, 1 ] )

  nf = NewtonFractal(
    poly_params = poly_params[ :: -1 ],
    N_points = N_POINTS,
    N_iter = N_ITER,
    bound = BOUND,
    method = METHOD )

  nf.visualize(
    colormap = 'tab10',
    filename = f'big_{METHOD}.png' )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#