import math, pprint

def dot( a, b ):
  return a[ 0 ] * b[ 0 ] + a[ 1 ] * b[ 1 ] + a[ 2 ] * b[ 2 ] 
# end def

def norm_cross( a, b ):
  x = a[1] * b[2] - a[2] * b[1]
  y = a[2] * b[0] - a[0] * b[2]
  z = a[0] * b[1] - a[1] * b[0]
  n = (x*x + y*y + z*z)**0.5
  if n > 0:
    return [ x / n, y / n, z / n]
  else:
    return [ 0, 0, 0 ]
# end def

def x( u, v, d = 0, dir = 'u' ):
  if d == 0:
    return math.sin( u ) * math.cos( v )
  elif d == 1:
    if dir == 'u':
      return math.cos( u ) * math.cos( v )
    else:
      return -( math.sin( u ) * math.sin( v ) )
  else:
    if dir == 'u':
      return -( math.sin( u ) * math.cos( v ) )
    elif dir == 'v':
      return -( math.sin( u ) * math.cos( v ) )
    elif dir == 'uv':
      return -math.cos( u ) * math.sin( v )
    else:
      return -( math.cos( u ) * math.sin( v ) )

# end def

def y( u, v, d = 0, dir = 'u' ):
  if d == 0:
    return math.sin( u ) * math.sin( v )
  elif d == 1:
    if dir == 'u':
      return math.cos( u ) * math.sin( v )
    else:
      return math.sin( u ) * math.cos( v )
  else:
    if dir == 'u':
      return -math.sin( u ) * math.sin( v )
    elif dir == 'v':
      return -math.sin( u ) * math.sin( v )
    elif dir == 'uv':
      return math.cos( u ) * math.cos( v )
    else:
      return math.cos( u ) * math.cos( v )

# end def

def z( u, v, d = 0, dir = 'u' ):
  if d == 0:
    return math.cos( u )
  elif d == 1:
    if dir == 'u':
      return -math.sin( u )
    else:
      return 0
  else:
    if dir == 'u':
      return -math.cos( u )
    elif dir == 'v':
      return 0
    elif dir == 'uv':
      return 0
    else:
      return 0
      
# end def

def r( u, v, d = 0, dir = 'u' ):
  return [ x( u, v, d ), y( u, v, d ), z( u, v, d ) ]
# end def

def I1( u, v ):
  du = r( u, v, 1, 'u' )
  dv = r( u, v, 1, 'v' )
  E = dot( du, du )
  F = dot( du, dv )
  G = dot( dv, dv )
  return [ [ E, F ] , [ F, G ] ]
# end def

def I2( u, v ):
  duu = r( u, v, 2, 'u' )
  dvv = r( u, v, 2, 'v' )
  duv = r( u, v, 2, 'uv' )
  n = norm_cross( r( u, v, 1, 'u' ), r( u, v,1, 'v' ) )
  L = dot( n, duu )
  M = -dot( n, duv )
  N = dot( n, dvv )
  return [ [ L, M ] ,  [ M, N ] ]
# end def

pprint.pprint( I1( 0.1, 0.3 ) )
pprint.pprint( I2( 0.1, 0.3 ) )



