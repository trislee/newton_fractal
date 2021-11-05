# -*- coding: UTF-8 -*-

"""Class for generating and visualizing Newton Fractals
"""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

import numpy as np
from matplotlib.cm import get_cmap
from skimage.io import imsave

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

class NewtonFractal:

  """Class for generating and visualizing Newton Fractals.

  Parameters
  ----------
  poly_params : numpy.ndarray<float>
    NumPy array corresponding to coefficients of the polynomial to generate
    Newton fractal with.
  N_points : int or None
    Number of points to compute fractal using, in both x and y dimensions.
    Defaults to ``512``.
  N_iter : int or None
    Number of Newton-Raphson iterations to perform for each point.
    Defaults to ``100``.
  bound : float or None
    Lower and upper bounds of grid.
    Defaults to ``2.0``.

  """

  #---------------------------------------------------------------------------#

  def __init__( self,
    poly_params,
    N_points = None,
    N_iter = None,
    bound = None,
    method = None ):

    # Convert input arguments to standard form
    #.........................................................................#

    self._poly_params = np.asarray( poly_params, dtype = float )

    if N_points is None:
      self._N_points = 512
    else:
      self._N_points = N_points + N_points % 2

    if N_iter is None:
      self._N_iter = 10
    elif N_iter <= 0:
      raise ValueError( '`N_iter` must be greater than or equal to 1' )
    else:
      self._N_iter = int( N_iter )

    if bound is None:
      self._bound = 2
    elif bound <= 0:
      raise ValueError( '`bound` must be greater than 0' )
    else:
      self._bound = float( bound )

    if method is None:
      self._method = 'fortran'
    else:
      if method in ( 'fortran', 'numpy' ):
        self._method = str( method )
      else:
        msg = "`method` argument must be either `'fortran'` or `'numpy'`."
        raise ValueError( msg )

    # Initialize auxilliary variables
    #.........................................................................#

    self._poly_der_params = np.polyder( p = self._poly_params, m = 1 )

    self._roots = np.roots( p = self._poly_params )

    _x = np.linspace( -self._bound, self._bound, self._N_points )
    X, Y = np.meshgrid( _x, _x )
    Z_0 = X + 1j * Y

    self.Z_history = np.zeros( ( self._N_iter, self._N_points, self._N_points ), dtype = complex )
    self.Z_history[ 0 ] = Z_0

    # Generate fractal
    #.........................................................................#

    self.execute( )

  #---------------------------------------------------------------------------#

  def _execute_numpy( self ):

    # Iterate over all iterations using Newton-Raphson method
    #.........................................................................#

    for i in range( 1, self._N_iter ):

      Z_i = self.Z_history[ i - 1 ]

      step = np.polyval( p = self._poly_params, x = Z_i ) \
        / np.polyval( p = self._poly_der_params, x = Z_i )

      self.Z_history[ i ] = Z_i - step

    # Compute index of closest root to each point
    #.........................................................................#

    Z = self.Z_history[ -1 ]

    distance = np.absolute(
      Z[ ..., np.newaxis ] - self._roots[ np.newaxis, np.newaxis, : ] )

    self._root_indices = np.argmin( distance, axis = -1 ).reshape( self._N_points, self._N_points )

    # Compute shading factor based on convergence speed
    #.........................................................................#

    all_distances = np.abs(
      self.Z_history[ ..., np.newaxis ] - self._roots[ np.newaxis, np.newaxis, np.newaxis, : ] )

    closest_distances = np.take_along_axis(
      arr = all_distances,
      indices = self._root_indices[ np.newaxis, :, :, np.newaxis ],
      axis = -1 ).squeeze( )

    self._convergence = np.argmin( ( closest_distances - 5e-16 ) ** 2, axis = 0 )

  #---------------------------------------------------------------------------#

  def _execute_fortran( self ):

    _root_indices = np.zeros( ( self._N_points, self._N_points ), dtype = np.int32 )
    root_indices = np.asfortranarray( _root_indices )

    _convergence = np.zeros( ( self._N_points, self._N_points ), dtype = np.int32 )
    convergence = np.asfortranarray( _convergence )

    self._fortran_function( self._poly_params, self._roots, self._N_iter, root_indices, convergence )

    self._root_indices = root_indices.T - 1
    self._convergence = convergence.T

  #---------------------------------------------------------------------------#

  def execute( self ):

    if self._method == 'fortran' :
      try:
        import newton_fractal_f2py
        self._fortran_function = newton_fractal_f2py.newton.main
      except ImportError:
        msg = '`newton_fractal f2oy extension module has not been created. Run the `make` command to create it.`'
        raise ImportError( msg )

      self._execute_fortran( )

    else:

      self._execute_numpy( )

  #---------------------------------------------------------------------------#

  def visualize( self, colormap = None, filename = None ):

    """Generate and save bitmap corresponding to Newton fractal.

    Parameters
    ----------
    index : int or None
      Iteration index of Newton-Raphson method to generate bitmap from.
      Defaults to ``-1`` (last iteration).
    colormap : str or None
      Matplotlib colormap to use to color different regions of fractal.
      Defaults to ``'tab10'``.
    filename : str or None
      File path to save bitmap of Newton fractal to.
      Defaults to ``'newton_fractal.png'``

    """

    # Convert input arguments to standard form
    #.........................................................................#

    if colormap is None:
      _colormap = 'tab10'
    else:
      _colormap = str( colormap )

    if filename is None:
      _filename = 'newton_fractal.png'
    else:
      _filename = str( filename )

    cmap = get_cmap( _colormap )

    _shades = self._convergence - np.min( self._convergence )
    _shades = _shades / np.max( _shades )
    shades = 1 - _shades
    shades **= 4

    # If the colormap is not qualitative, scale index values to range [ 0, 1].
    #.........................................................................#

    if cmap.N < 256:
      root_indices = self._root_indices
    else:
      root_indices = self._root_indices / np.max( self._root_indices )

    # Generate RGB bitmap from closest root indices
    #.........................................................................#

    _img_arr = cmap( root_indices )
    _img_arr = _img_arr[ :, :, :3 ] * shades[ ..., np.newaxis ]
    img_arr = np.asarray( _img_arr * 255, dtype = np.uint8 )

    imsave(
      fname = _filename,
      arr = img_arr )

    #.........................................................................#

  #---------------------------------------------------------------------------#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#