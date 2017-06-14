# integration routine, simpson's rule. integrates a function over a pre-set 1-D grid.
def integrate(integrand, x_grid, **kwargs):
        simpson = 0.0
        ii = 0
        for ii in range(0, len(x_grid)-1):
           start_pt = x_grid[ii]
           end_pt =  x_grid[ii+1]
           simpson = simpson + (end_pt - start_pt)/6.0*(integrand(start_pt, **kwargs) + 4.0*integrand(0.5*\
                                       (start_pt+end_pt), **kwargs) + integrand(end_pt, **kwargs))
           ii = ii+1

        return simpson