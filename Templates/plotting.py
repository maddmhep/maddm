import logging
import sys
logging.basicConfig(filename='plotting.log',level=logging.DEBUG, filemode='w')

try:
    import matplotlib.pyplot as pl
except Exception, error:
    logging.debug(error)
    print "ERROR: MadDM plotting routines require matplotlib, which is not installed"
    sys.exit(1)

try:
    import numpy as np
except Exception, error:
    logging.debug(error)
    print "ERROR: MadDM plotting routines require numpy, which is not installed"
    sys.exit(1)

from matplotlib.mlab import griddata
import matplotlib.colors as colors

#----------------------------------------------------------------------------
# NOTE: All plotting functions will take either a list of numpy arrays or a list of file names.
#----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# This function plots a line out of data points,
# inputfiles - a list of input files (including a path if they are in a different
# folder), or a list of numpy arrays.
# Data files must be formatted as two column!
# colors - list of colors for histograms
# styles - list of styles ('solid', 'dashed', dotted' available)
# xlab - label for x-axis (can be latex formatted)
# ylab - label for y-axis (can be latex formatted)
# line_width - thickness of the lines
# linelabels - labels for individual histograms (can be latex formatted)
# title - plot title
# columns - which columns to use as x and y respectively
# logplot - True for log y axis, False for linear y axis
#-----------------------------------------------------------------------------

def LinePlot(inputfiles, columns=[0,1],  colors=['darkblue'], styles =['solid'],\
						line_width=2, labels=['data set 1'], xlab='', ylab='', title='',\
						logplot=False):

	if len(colors)!=len(labels)!=len(styles)!=len(inputfiles):
		print "ERROR: You have to specify colors, styles and labels for your histograms"
		exit

	i =0
	for input_file in inputfiles:
		if type(input_file)==type(str()):
			try:
				xdata, ydata= np.loadtxt(input_file, unpack=True, usecols=columns)
			except:
				print "Could not open input files!"
				exit
		elif type(input_file)==np.ndarray:
			xdata = input_file[0]
			ydata = input_file[1]
		else:
			print "ERROR: You tried to pass data to the LinePlot() function in a format it cannot read"
			exit

    		pl.plot(xdata, ydata,label=labels[i], alpha=1, linewidth=line_width,\
    	 		linestyle=styles[i],color=colors[i])

		i=i+1

	pl.grid(True)
	pl.legend()
	if logplot:
		pl.yscale('log')
	pl.title(title, fontsize=17)
	pl.tick_params(axis='x', labelsize=20)
	pl.tick_params(axis='y', labelsize=20)
	pl.xlabel(xlab, fontsize=20)
	pl.ylabel(ylab, fontsize=20)
	pl.show()

#-----------------------------------------------------------------------------
# This function creates a contour plot out of a file specified by the
# datafile variable (requires three column input format). You can also feed in a list
# of numpy data arrays as well.
# xlabel and ylabel are strings which will be placed
# on the x and y axes respectively. It is allowed to use latex labels.
# contours is either a number of contours one wishes to plot or an array of
# specific contours which should be plotted. colormap is the choice of coloring
# scheme for the contour plot. See matplotlib contour() and contourf() function
# documentation for more info.
# columns - which columns in the input file to use as x,y,z respectively
#-----------------------------------------------------------------------------
def ContourPlot(datafile,columns=[0,1,2], nXvals=100, nYvals=100, title='', xlab='', ylab='', colormap='rainbow', contours=10):

	if type(datafile)==type(str()):
		try:
			x, y, z = np.loadtxt(datafile, dtype='float', unpack=True, usecols=columns)
		except:
			print "Can't open the input file!"
			exit
	elif type(datafile)==np.ndarray:
			x = datafile[0]
			y = datafile[1]
			z = datafile[2]
	else:
			print "ERROR: You tried to pass data to the ContourPlot() function in a format it cannot read"
			exit

	print type(x)

	xi = np.linspace(np.amin(x), np.amax(x), nXvals)
	yi = np.linspace(np.amin(y), np.amax(y), nYvals)
	zi = griddata(x, y, z, xi, yi)
	norm = colors.Normalize(vmin = np.min(z), vmax = np.max(z), clip = False)
	pl.figure()
	pl.contourf(xi, yi, zi, 30, cmap = pl.get_cmap(colormap), norm =norm)
	CS = pl.contour(xi, yi, zi, contours, colors = 'k',lw = 3)
	pl.clabel(CS, inline=1, fontsize=10, fmt='%.1e')
	pl.tick_params(axis='x', labelsize=20)
	pl.tick_params(axis='y', labelsize=20)
	pl.title(title, fontsize=17)
	pl.xlabel(xlab, fontsize=20)
	pl.ylabel(ylab, fontsize=20)
	pl.show()


#-----------------------------------------------------------------------------
# This function plots a histogram using an arbitrary number of input files,
# inputfiles - a list of input files (including a path if they are in a different
# folder), or a list of numbers to be plotted (in numpy array format).
# Data files must be formatted as one column!
# nbins - number of histogram bins
# colors - list of colors for histograms
# styles - list of styles ('solid', 'dashed', dotted' available)
# xlab - label for x-axis (can be latex formatted)
# ylab - label for y-axis (can be latex formatted)
# line_width - thickness of the histogram lines
# linelabels - labels for individual histograms (can be latex formatted)
# title - plot title
# logplot - True for log y axis, False for linear y axis
# norm - 1 for histograms normalized to unit area, 0 for normalized to total number
# 				of data points.
#
#-----------------------------------------------------------------------------

def Histogram(inputfiles, nbins = 30, colors=['darkblue'], styles =['solid'],\
						line_width=2, labels=['data set 1'], xlab='', ylab='', title='',\
						logplot=False, norm=1):

	if len(colors)!=len(labels)!=len(styles)!=leb(inputfiles):
		print "ERROR: You have to specify colors, styles and labels for your histograms"
		exit

	i =0
	for input_file in inputfiles:
		#If the user fed in a filename, open the file name and read in the data.
		if type(input_file)==type(str()):
			try:
				observable= np.loadtxt(input_file, unpack=True, usecols=[0])
			except:
				print "Could not open input files!"
				exit
		#Else treat the input as a numpy array
		elif type(input_file)==np.ndarray:
			observable = input_file
		else:
			print "ERROR: You tried to pass data to the Histogram() function in a format it cannot read"
			exit


		pl.hist(observable, bins=nbins, histtype='step', normed=norm,\
				label=labels[i], alpha=1, linewidth=line_width, log=logplot, linestyle=styles[i],\
				color=colors[i])

		i=i+1

	pl.grid(True)
	pl.legend()

	pl.title(title, fontsize=17)
	pl.tick_params(axis='x', labelsize=20)
	pl.tick_params(axis='y', labelsize=20)
	pl.xlabel(xlab, fontsize=20)
	pl.ylabel(ylab, fontsize=20)
	pl.show()