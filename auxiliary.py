# integration routine, simpson's rule. integrates a function over a pre-set 1-D grid.
#add option to run it with x_min, x_max npts as well.
def integrate(integrand, x_grid=[], a=0, b=0, npts=0, **kwargs):
        simpson = 0.0

        #if the user supplies the limits then use those
        #otherwise use the grid supplied by the user.
        if (a!=b and a!=npts):
            x_grid = [float(a)+ (float(b)-float(a))*k/npts for k in range(0, npts+1)]

        if x_grid==[]:
            print("ERROR: Integration grid in integrate() function not set up! Returning 0.")
            return simpson

        for ii in range(0, len(x_grid)-1):
           start_pt = x_grid[ii]
           end_pt =  x_grid[ii+1]
           simpson = simpson + (end_pt - start_pt)/6.0*(integrand(start_pt, **kwargs) + 4.0*integrand(0.5*\
                                       (start_pt+end_pt), **kwargs) + integrand(end_pt, **kwargs))


        return simpson

def write_data_to_file(x_data, y_data, filename=''):

    if len(x_data)!= len(y_data):
        print("ERROR: x and y columns don't have the same dimensions. Will not write the data.")
        return False

    try:
        with open(filename, 'w+') as f:
            for i, x in enumerate(x_data):
                f.write('%.5e     %.5e' %(x, y_data[i]))

    except OSError:
        print('ERROR: write_data_to_file can not open the file!')
        return
