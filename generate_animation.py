# -*- coding: UTF-8 -*-

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import os

import numpy as np

from newton_fractal import NewtonFractal

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

N_POINTS = 512
N_ITER = 100
BOUND = 2.0

N_T_VALS = 101

OUTPUT_DIR = '2_to_4'

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

if __name__ == '__main__':

  os.makedirs( OUTPUT_DIR, exist_ok = True )

  t_vals = np.linspace( 0, 1, N_T_VALS )

  for i, t in enumerate( t_vals ):

    poly_params = np.array( [ t, 1 - t, 0, 1 ] )

    nf = NewtonFractal(
      poly_params = poly_params,
      N_points = N_POINTS,
      N_iter = N_ITER,
      bound = BOUND )

    nf.visualize(
      colormap = 'viridis',
      filename = os.path.join( OUTPUT_DIR, f'{i:03d}.png' ) )

  for i, t in enumerate( t_vals ):

    poly_params = np.array( [ t, 1 - t, 0, 0, 1 ] )

    nf = NewtonFractal(
      poly_params = poly_params,
      N_points = N_POINTS,
      N_iter = N_ITER,
      bound = BOUND )

    nf.visualize(
      colormap = 'viridis',
      filename = os.path.join( OUTPUT_DIR, f'{i + N_T_VALS:03d}.png' ) )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#