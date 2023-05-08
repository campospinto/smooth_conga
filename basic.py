import numpy as np

def create_open_knot_vector(p, n):
    """
    Create an open knot vector in [0,1] for splines with degree p and n control points.

    open knots = interior points and (p+1)-repated boundary points
    (by default the points are regular)
    
    example, for p = 1 and n = 5:

    X       X
    X X X X X
    0       1

    """
    m = n + p + 1
    npoints = n - p + 1
    ratio = 1./(npoints-1)
    knot_vector = np.zeros(m)
    # for i in range(p, n+1):
        # knot_vector[i] = (i - p + 1)*ratio
    knot_vector[p:p+npoints] = np.array(range(npoints))*ratio
    knot_vector[n+1:m] = 1.0
    return knot_vector

def construct_grid_from_knots(p, n, knots):
    """
    return the grid made by interior knots

    p            spline degree 
    n            number of control points
    knots        knot vector 
    """
    return knots[p:n+1].copy()
    