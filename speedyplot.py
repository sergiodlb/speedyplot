#!/usr/bin/python3
# 2016-12-15 Sergio de la Barrera
# - forked from plotmacdata.py
# - 2017-01-19 added custom delimiter
# - 2017-01-28 fixed column count to ignore whitespace columns created by extra delimiters
# - 2017-01-31 added 'chain' option to allow plotting a series of files chained together in time
#              also added vertical option to force subplots to stack vertically
# - 2017-02-02 added horizontal option, like vertical
#              also added 'arrows' option to indicate direction of sweeping data
# - 2017-02-03 added glob function to allow expansion of filename wildcards
#              modified arrows to take optional argument of arrows length in data points
# - 2017-03-24 added ymult option to scale y-values
# - 2017-03-25 added active ymult to axis label
#              added logx option
# - 2017-03-30 added number of rows to stdout
#              added interpolation feature for better symmetrization (and plotting limits)
#              added boxcar average feature
# - 2017-04-04 added title for (anti)symmetrized panels
#              added sharey toggle to use same y-axis range
# - 2017-04-05 added divcol, column to divide data by
# - 2017-04-06 added yshift to allow constant shift of data before other computations
# - 2017-04-08 added option to fold negative x-data into positive axis
# - 2017-04-13 added option to superimpose all data in the same plot
#              enabled file/column specific ymult values: f1c1,f1c2;f2c1,f2c2;...
# - 2017-06-01 added manual axis label options
# - 2017-06-08 made legend draggable and moved legend options to rcParams
# - 2017-06-12 allowed for files with differing number of columns (usecols must mutually exist)
#              added xlim option for manual axis zooming
#              integrated speedy2d colorplotting into speedyplot given two x-columns
#              added semicolon specification for different ycolumns from different files: f1y1,f1y2;f2y1f2y2;...
# - 2017-06-13 added semicolon specification for different x-columns from different files: f1y1,f1y2;f2y1f2y2;...
#              added save command option to batch file (by pressing 'b' key)
# - 2017-06-26 changed default delimiter value to None, which should still use whitespace unless specified otherwise
# - 2017-06-29 corrected previous change to default delimiter value; now checks list of common delimiters until
#              totcols > 1 unless delimiter is explicitly specified; also seeks back to zero for 0 line file headers
# - 2017-11-07 incorporated tight_layout() and font.size, mathtext.default rcParams
# - 2018-01-12 added comma in list of trial delimiters and replaced row number replacement for first column of files with string columns to use np.genfromtxt(),
#              replacing the first column but leaving others as np.nan values
#              prints number of columns for each file
#              removed stripping of digits and underscores from column labels; now only left-strips '#'
#              created --list option to print column labels/headers and quit (useful for picking columns to plot)
# - 2018-01-16 added -r option ("combinecols") which will perform RMS on two columns to produce magnitude of Re and Im parts of lock-in signal
#              also added a basic "output" function (mapped to "o" key) to save the plotted/manipulated data from the last column of last file
# - 2018-01-26 added ability to invert x- or y-axes (e.g. for estimate of conductance)
#              extended yshifts to accept single argument and use that value to separate all curves as n*yshifts
#              remapped "output" function to "e" key (for "export") to prevent clash with "o" zoom toggle in view window
# - 2018-01-30 fixed ylabel in combinecols mode to print sqrt(x^2 + y^2)
# - 2018-02-07 added waterfall plot mode which shifts ydata by value of second x-axis with optional multiplier
# - 2018-03-22 changed order of divcol and invert so columns are inverted last (allows conversion of V-->R-->conductance)
#              added optional argument for --marker which allows specification of marker size for both line and colorplots
#              added --square marker option
#              created --derivative option (specify columns or perform on all y-columns)
#              added --integrate option which uses cumulative trapezoidal summation
#              added basic --datetime capability specifically for reading special-measure output (ISO 8601 timestamps)
# - 2018-03-23 implemented file number ranges using bash syntax e.g. {001..021} (zero padding conveys expected width of file numbers)
# - 2018-03-27 extended xhift --> xshifts, allowing semicolon-separated listof shifts for multiple columns/files
#              same for xmult
# - 2018-04-04 enabled individual colorplot colorbars with separate automatic scaling for each colorplot
#              added colorbar labels to colorplot axes
#              improved mode switching (added colorplot_mode, lineplot_mode with to change behavior based on mode throughout)
#              fixed labeling system to start with defaults (col 1, col 2, ...), replace with header labels if available, then x, y, c labels for each if provided (using defaults or header labels otherwise)
#              added interpolation mode for colorplots
# - 2018-04-05 added ylims option for manual y-axis zooming
#              fixed custom label truncation issue by using python list for labels until after manual label replacement
# - 2018-04-09 extended crange to allow semicolon-separated values for different panels
#              fixed refresh handling of file_or_glob by moving call to within refresh loop (no longer used as pre-parser for command-line argument input)
# - 2018-04-12 moved zeroy until after applying ymult
# - 2018-04-16 moved raw 2D scatter plot to after main file loop, similar to interpolated plot, in order to obtain better colorscale limits
# - 2018-04-17 added continue statement to skip empty data files; required empty elements to be added to data list and logic for handling those in concatenation
# - 2018-04-26 enabled two-column datetime input, specified by "-d c1,c2" where c1 is date column, c2 time (in BlueFors format)
# 2018-07-12    * added --deleterows option which masks selected rows (comma list) with NaN in plotted data
#               * added --lw (linewidth) option as well
# 2018-07-13 created --monotonic option which delets elements for which the sweep variable is non-monotonic (e.g. spikes in temperature)
# 2018-08-02 modified order of fold+xnorm to follow xmult+xshifts (used to be reverse; present order allows a shift in fold origin"
# 2018-08-03 created colorbyorder option to color curves according to colormap
# 2018-08-13    * introduced engineering notation for axis tick labels
#               * fixed bug which did not display y-label when using datetime mode
#               * changed x-label to 'points' when plotting versus row number
#               * added timedelta (-t) plotting option, which operates similarly to datetime (-d) but versus elapsed time instead of absolute time
#               * made color iterable
#               * added nolegend option to prevent display of legend
# 2018-08-21    * added xboxcar and x2boxcar to allow averaging of x-data (on both axes in colorplots)
#				* also moved boxcar averaging to occur before interpolation
# 2018-08-29 changed colorplot xdata to a mutable reference instead of a copy; this allows xdata manipulation to propogate to colorplot section
# 2018-09-05 gave colorbyorder an optional argument (cmap reference)
# 2018-09-05 gave line an optional argument (linewidth)
# 2018-09-16 added trimspikes to remove data spikes (2% deviation from rolling average, excluding self)
# 2018-09-16 added traces mode which plots specified linecuts from first x-axis versus the second x-axis
# 2018-09-16 added xdivcol (same functionality as divcol but for x-data)
# 2018-09-16 added ifinlabels option: checks file header for presence of a specific string label and skips file if not present
# 2018-12-11 edited trimspikes to allow multiple columns (comma-separated list)
# 2019-01-06 added even/odd flags to select only even/odd numbers in file ranges specified with {..} notation
# 2019-02-06 added 'orezy' option, which essentially works like zeroy, but uses the max value
# 2019-02-27 changed csints to allow range specification with colon-slices (e.g. colummns 1:9, or rows 100:1001)
# 2019-03-21 added csigma option which sets colorplot crange to a multiple of the z-data standard devation around the mean
#
# * add fft
# * add better refresh, autoupdate
# * static panel aspect ratio
# * quit using np array for labels; use list

import numpy as np
import scipy.interpolate as sci
import scipy.integrate as scg
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LinearSegmentedColormap
import itertools as it
import argparse, os, glob, re
import sys, datetime

def csints(string):
    if ':' in string:
        l = []
        for substring in string.split(','):
            if ':' in substring:
                r1, r2 = map(int, substring.split(':'))
                l.extend(range(r1, r2+1))
            else:
                l.append(int(substring))
    else:
        l = list(map(int, string.split(',')))
    return l

def csfloats(string):
    return list(map(float, string.split(',')))

def ssints(string):
    return list(map(csints, string.split(';')))

def ssfloats(string):
    return list(map(csfloats, string.split(';')))

def ssstrs(string):
    return list(map(str, string.split(';')))

def itstrs(string):
    return it.cycle(map(str, string.split(',')))

def file_or_glob(string):
    if os.path.isfile(string):
        fname = string
    else:
        pattern = r'{(\d+)\.\.(\d+)}'
        match = re.search(pattern, string)
        if match: # use file number range
            mgroup = match.group(1), match.group(2)
            if int(mgroup[0]) < int(mgroup[1])+1:
                fnum_range = range(int(mgroup[0]), int(mgroup[1])+1)
            else:
                fnum_range = range(int(mgroup[0]), int(mgroup[1])-1, -1)
            fnum_width = max(len(mgroup[0]), len(mgroup[1]))
            repl = '{{:0{}}}'.format(fnum_width)
            if args.even:
                strings = [re.sub(pattern, repl.format(fnum), string, count=1) for fnum in fnum_range if (fnum+1) % 2]
            elif args.odd:
                strings = [re.sub(pattern, repl.format(fnum), string, count=1) for fnum in fnum_range if fnum % 2]
            else:
                strings = [re.sub(pattern, repl.format(fnum), string, count=1) for fnum in fnum_range]
            fname = []
            for s in strings:
                fname.extend(glob.glob(s))
        else: # regular glob syntax
            fname = glob.glob(string)
    return fname

def on_key_press(event):
    if type(event.key) is str:
        if event.key == 'b':
            command_args = sys.argv.copy()
            command_args[0] = os.path.basename(command_args[0])
            batch_fname = os.path.splitext(command_args[0])[0] + '-' + datetime.datetime.now().strftime('%Y-%m-%d-%Hh%Mm%Ss') + '.bat'
            with open(batch_fname, 'w') as batch_f:
                print(' '.join(command_args), file=batch_f)
                print('command saved to --> {}'.format(os.path.join(os.getcwd(), batch_fname)))
        elif event.key == 'e':
            # saves a tab-delimited ascii file of plotted data
            # only makes sense for a single file one column
            command_args = sys.argv.copy()
            command_args[0] = os.path.basename(command_args[0])
            output_fname = os.path.join(os.path.dirname(fname), os.path.splitext(command_args[0])[0] + '-' + datetime.datetime.now().strftime('%Y-%m-%d-%Hh%Mm%Ss') + '-' + os.path.basename(fname))

            output_header = '{:>19}\t{:>20}'.format(ax.get_xlabel(), ax.get_ylabel())
            output = np.column_stack((xdata, col))
            np.savetxt(output_fname, output, fmt='%20.12e', delimiter='\t', newline='\n', header=output_header, comments='#')
            print('saved data to -->', output_fname)
            if len(files) > 1 or len(args.usecols[0]) > 1:
                print('(only saves final column from last file)')
        elif event.key == 'r':
            plt.close()
            main()
            #main(event.canvas.figure)

# Define input arguments
parser = argparse.ArgumentParser(description='plot somewhat arbitrary data stored in ascii-format column data')
#parser.add_argument('files', nargs='+', type=file_or_glob, help='input files')
parser.add_argument('files', nargs='+', type=str, help='input files')
parser.add_argument('--even', action='store_true', help='select only even number files in range, specified by {n1..n2} notation')
parser.add_argument('--odd', action='store_true', help='select only odd number files in range, specified by {n1..n2} notation')
parser.add_argument('--hysteresis', action='store_true', help='subtract every other curve')
parser.add_argument('--header', type=int, default=1, help='number of header lines to skip', metavar='N')
parser.add_argument('--listcols', action='store_true', help='print list of column headers or labels, if available (first file only)')
parser.add_argument('-x', '--plotvs', type=ssints, default=[[0,]], help='column index or indices to plot against (zero-indexed)', metavar='COLS')
parser.add_argument('-y', '--usecols', type=ssints, default=None, help='columns to plot', metavar='COLS')
parser.add_argument('-r', '--combinecols', type=csints, default=None, help='x,y columns --> R = (x**2 + y**2)**0.5', metavar='X,Y')
parser.add_argument('-d', '--datetime', nargs='?', default=False, const=(0,), type=csints, help='column index to interpret as datetime (replaces usual x-axis; combine two cols using comma list)', metavar='COL')
parser.add_argument('-t', '--timedelta', nargs='?', default=False, const=(0,), type=csints, help='column index to interpret as datetime (replaces usual x-axis)', metavar='COL')
parser.add_argument('--ifinlabels', type=str, default=None, help='filter file list based on check for a specific column header', metavar='STR')
parser.add_argument('--crange', type=ssfloats, default=None, help='range of 2D data values (z-axis; can be semicolon list)')
parser.add_argument('--csigma', type=float, default=None, help='multiple of z-axis stddev to use to set color range')
parser.add_argument('--cmap', default=mcm.inferno, help='colormap for 2d plots')
parser.add_argument('--cornplot', action='store_true', help='break up total field into parallel and perp components')
parser.add_argument('--yshift', type=float, default=None, help='constant shift to y-data before other manipulations')
parser.add_argument('--interpolate', type=ssfloats, default=False, metavar='X1,X2,Nx;Y1,Y2,Ny',
        help='limits and number of elements to interpolate with respect to plotting axis; single value uses data limits')
parser.add_argument('--boxcar', type=int, default=False, help='perform boxcar average over N neighbors', metavar='N')
parser.add_argument('--xboxcar', type=int, default=False, help='perform boxcar average on x-data over N neighbors', metavar='N')
parser.add_argument('--x2boxcar', type=int, default=False, help='perform boxcar average on second x-axis of colorplot over N neighbors', metavar='N')
parser.add_argument('--lorentz', type=int, default=False, help='perform lorentzian blur over N neighbors', metavar='N')
parser.add_argument('--symmetrize', type=csints, default=[], metavar='COLS',
        help='columns to symmetrize; requires rows to be symmetric about zero and equally spaced')
parser.add_argument('--antisymmetrize', type=csints, default=[], metavar='COLS',
        help='columns to antisymmetrize; requires rows to be symmetric about zero and equally spaced')
parser.add_argument('--chain', action='store_true', help='chain time series plots together from multiple files')
parser.add_argument('--line', nargs='?', type=float, default=False, const=True, help='plot solid lines (may also give linewidth)')
#parser.add_argument('--marker', action='store_true', help='plot using point markers')
parser.add_argument('--marker', nargs='?', default=False, const=True, type=int, help='plot using point markers (can take number argument which changes colorplot marker size; default s=200)')
parser.add_argument('--square', action='store_true', help='plot using square markers')
#parser.add_argument('--arrows', action='store_true', help='add arrows to indicate sweep direction')
parser.add_argument('--arrows', nargs='?', default=False, const=1, type=int,
                    help='add arrows to indicate sweep direction')
parser.add_argument('--dashes', action='store_true', help='use dashes to indicate negative sweep direction')
parser.add_argument('--fold', action='store_true', help='fold negative x-values into positive axis')
parser.add_argument('--logx', action='store_true', help='use log scale on x-axis')
parser.add_argument('--logy', action='store_true', help='use log scale on y-axis')
parser.add_argument('--xmult', type=ssfloats, default=None, help='multiplier to scale x-axes (can be list)')
parser.add_argument('--xshifts', type=ssfloats, default=None, help='offset to shift x-axes (can be list)')
parser.add_argument('--xinvert', action='store_true', help='invert x-axis')
parser.add_argument('--xlabel', type=str, default=None, help='manually specify x-axis label (string)')
parser.add_argument('--ylabels', type=ssstrs, default=None, help='manually specify y-axis labels (string; semicolons for multiple axes)')
parser.add_argument('--clabels', type=ssstrs, default=None, help='manually specify c-axis labels for colorplots (string; semicolons for multiple axes)')
parser.add_argument('--divcol', type=int, default=None, help='divide all columns by the data in this column', metavar='COL')
parser.add_argument('--xdivcol', type=int, default=None, help='divide x-column by the data in this column', metavar='COL')
parser.add_argument('--normy', action='store_true', help='normalize y-data to range between zero and one')
parser.add_argument('--xnorm', action='store_true', help='normalize x-data to range between zero and one')
parser.add_argument('--ymult', type=ssfloats, default=None, help='multiplier to scale y-axes (can be list)')
parser.add_argument('--zeroy', action='store_true', help='shift y-data down to zero by minimum value')
parser.add_argument('--orezy', action='store_true', help='shift y-data to zero by MAXimum value')
parser.add_argument('--meansubtract', action='store_true', help='shift y-data down to zero by average value')
parser.add_argument('--yshifts', type=csfloats, default=None, help='constant shifts to y-data applied after other manipulations; list or single value (creates crude waterfall plot ordered by file sequence)')
parser.add_argument('--waterfall', nargs='?', type=float, default=False, const=1, metavar='MULT', help='toggle waterfall plot; requires two x-columns; optional value specifies multiplier for y-separation')
parser.add_argument('--traces', type=csfloats, default=None, metavar='VALS', help='trace/linecut mode; comma list of values (from x-axis 1) along which to plot traces (vs x-axis 2)')
parser.add_argument('--derivative', nargs='?', type=csints, default=False, const=True, metavar='COLS', help='columns to take derivative/central differences (wrt x-axis; can be list)')
parser.add_argument('--integrate', nargs='?', type=csints, default=False, const=True, metavar='COLS', help='columns to integrate (cumulative trapezoidal summation wrt x-axis; can be list)')
parser.add_argument('--invert', nargs='?', type=csints, default=False, const=True, metavar='COLS', help='columns to invert (can be list)')
parser.add_argument('--sharey', action='store_true', help='use same y-range in all panels')
#parser.add_argument('--delimiter', type=str, default='\t', help='delimiter separating columns')
parser.add_argument('--delimiter', type=str, default=None, help='delimiter separating columns')
parser.add_argument('--superimpose', action='store_true', help='draw all columns together')
parser.add_argument('--vertical', action='store_true', help='draw subplots vertically')
parser.add_argument('--horizontal', action='store_true', help='draw subplots horizontally')
parser.add_argument('--xlim', type=csfloats, default=None, help='specified limits for x-axes')
parser.add_argument('--ylims', type=ssfloats, default=None, help='specified limits for y-axes (can specify panel with semicolons)')
parser.add_argument('--deleterows', type=csints, default=None, help='row numbers to delete (mask as NaN) in all loaded data (comma list)')
parser.add_argument('--trimspikes', type=csints, default=None, help='column numbers to detect and remove data spikes (mask as NaN)')
parser.add_argument('--monotonic', action='store_true', help='remove points with opposite sweep direction to the average (cuts data with non-monotonic direction in the x-variable)')
parser.add_argument('-c', '--color', type=itstrs, default=None, help='override default line colors (may be comma list)')
parser.add_argument('--colorbyorder', nargs='?', type=str, default=False, const=True, help='use colormap to define line colors based on order of files in list (may also give colormap name)')
parser.add_argument('--lw', type=float, default=None, help='override default linewidth')
parser.add_argument('--ls', type=itstrs, default=None, help='override default linestyles (may be comma list)')
parser.add_argument('--nolegend', action='store_true', help='toggle display of legend')
args = parser.parse_args()

# Plot settings
mpl.rcParams['keymap.quit'] = 'q'
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.handletextpad'] = 0.3
mpl.rcParams['legend.frameon'] = False

mpl.rcParams['lines.linewidth'] = 2
lineprops = dict(ls='none', marker='.')
if args.square:
    lineprops['marker'] = 's'
if args.line:
    lineprops['ls'] = '-'
    if not args.marker:
        lineprops['marker'] = ''
    if args.line is not True:
        lineprops['lw'] = args.line # overrides mpl.rcParams['lines.linewidth'], but can be overridden by args.lw
if args.marker:
    colorplot_marker_size = args.marker
    if args.marker is not True:
        lineprops['ms'] = args.marker
else:
    colorplot_marker_size = 200
if args.lw:
    lineprops['lw'] = args.lw # overrides mpl.rcParams['lines.linewidth']

# custom colormaps
zidict = {'red':    ((0.0, 0.3, 0.3),
                     (0.25,0.2, 0.2),
                     (0.5, 0.0, 0.0),
                     (0.75,0.9, 0.9),
                     (1.0, 1.0, 1.0)),
          'green':  ((0.0, 0.75, 0.75),
                     (0.25,0.3, 0.3),
                     (0.5, 0.0, 0.0),
                     (0.75,0.13, 0.13),
                     (1.0, 0.9, 0.9)),
          'blue':   ((0.0, 0.9, 0.9),
                     (0.25,0.6, 0.6),
                     (0.5, 0.0, 0.0),
                     (0.75,0.13, 0.13),
                     (1.0, 0.04, 0.04))}
zibrov   = LinearSegmentedColormap('zibrov', zidict) # similar to 'cool' (but ending in black) + 'hot' colormaps
zibrov_r = zibrov.reversed()
plt.register_cmap(cmap=zibrov)
plt.register_cmap(cmap=zibrov_r)

            
##### functions
def lorentzian(width, x0=0):
    # returns lorenzian with FWHM = width/2, over domain = width
    x = np.linspace(-width/2, width/2, width)
    hwhm = width/4
    f = 1/(1 + ((x-x0)/hwhm)**2)
    return f/np.sum(f)
#####

def main(fig=None):
    # Flatten list of files
    files = []
    #for sublist in args.files:
    globbed_files = [file_or_glob(arg) for arg in args.files]
    for sublist in globbed_files:
        if type(sublist) is list:
            files.extend([fn for fn in sublist])
        else:
            files.append(sublist)

    # Determine number of columns in first file
    with open(files[0], 'r') as first_file:
        # Seek to first data row
        for line in range(args.header+1):
            first_line = first_file.readline()
        #totcols = len(first_line.split(args.delimiter))
        if args.delimiter:
            if args.delimiter == 'None':
                args.delimiter = None
            totcols = len([col.strip() for col in first_line.split(args.delimiter) if col.strip()])
        else:
            for delim in ['\t', None, ',', ';']:
                totcols = len([col.strip() for col in first_line.split(delim) if col.strip()])
                if totcols > 1:
                    args.delimiter = delim
                    break
    if args.combinecols:
        args.usecols = [args.combinecols[0:1]]
        ncols = 1
    elif args.usecols:
        ncols = len(args.usecols[0])
    else:
        ncols = totcols - 1
        args.usecols = [list(range(totcols)),]
        for n, col in enumerate(args.plotvs[0]):
            args.usecols[0].pop(col-n)
    sq = np.sqrt(ncols)

    #if fig is None:
    # Create figure
    if args.superimpose:
        n, m = 1, 1
    elif args.vertical:
        n, m = ncols, 1
    elif args.horizontal:
        n, m = 1, ncols
    else:
        n, m = int(round(sq)), int(np.ceil(ncols/sq))
    fprops = dict(sharex=True)
    if args.sharey:
        fprops.update(dict(sharey=True))
    fig, axes = plt.subplots(n, m, **fprops)
    axes = np.atleast_1d(axes).ravel()

    # Create event callbacks
    fig.canvas.mpl_connect('key_press_event', on_key_press)
    #else:
        #fig.clf()

    # Load the data
    data = []
    first_file = True
    for n, (fname, plotvs, usecols) in enumerate(zip(files, it.cycle(args.plotvs), it.cycle(args.usecols))):
        with open(fname, 'r') as f:
            if args.header:
                for line in range(args.header):
                    header = f.readline()
            else:
                header = f.readline()
                f.seek(0)

            # Update number of columns in case files differ
            totcols = len([col.strip() for col in header.split(args.delimiter) if col.strip()])

            # Mode switching
            if len(plotvs) > 1 and not (args.waterfall or args.traces):
                lineplot_mode = False
                colorplot_mode = True
            else:
                lineplot_mode = True
                colorplot_mode = False
            if len(plotvs) == 1:
                if args.waterfall: # waterfall plot mode selected without second x-column
                    print('waterfall plot requires two x-axes') # display error message
                if args.traces: # traces mode selected without second x-column
                    print('linecuts require two x-axes (first is dimension from which to extract traces, second is plotting axis)') # display error message
            #if len(plotvs) == 1 and args.colorbyvalue:
            #    print('color-by-value plot requires two x-axes (avg value of 2nd x-column one determines color)')

            # Deal with axes labels
            labels = np.asarray(['col {}'.format(col) + 20*'\0' for col in range(totcols)], dtype=str)
            if args.header > 0 or args.listcols: # replace defaults with header labels if they exist
                # Header from last file in list used; only last line of header used
                #labels = np.asarray([element.strip().strip('#_0123456789') for element in header.split(args.delimiter)], dtype=str)
                #labels = np.asarray([element.strip().lstrip('#') for element in header.split(args.delimiter)], dtype=str)
                labels = [element.strip().lstrip('#') for element in header.split(args.delimiter)] # keep as list until AFTER label replacement
            if args.ifinlabels and args.ifinlabels not in labels:
                print(args.ifinlabels, 'not in header; skipping...')
                continue
            if args.listcols: # moved above label replacement to prevent overwrite before listing column headers
                print('list of columns:')
                for nl, label in enumerate(labels):
                    print('{}:\t{}'.format(nl, label))
                sys.exit(0)
            if args.ylabels: # replace header labels with custom labels if provided
                if lineplot_mode:
                    #labels[usecols] = args.ylabels
                    for col, ylabel in zip(usecols, args.ylabels):
                        labels[col] = ylabel
                elif colorplot_mode:
                    labels[plotvs[1]] = args.ylabels[0]
            if args.xlabel:
                labels[plotvs[0]] = args.xlabel
            if colorplot_mode and args.clabels:
                #labels[usecols] = args.clabels
                for col, clabel in zip(usecols, args.clabels):
                    labels[col] = clabel
            # convert to numpy array of strings AFTER replacing labels (prevents truncation of long custom labels)
            labels = np.asarray(labels)

            # Load the data
            try:
                if args.timedelta is not False:
                    args.datetime = args.timedelta # handle both in the same way until plotting
                if args.datetime is not False:
                    date_fmt = '%Y-%m-%dT%H:%M:%S.%f'
                    str2date = lambda s: datetime.datetime.strptime(s.decode('utf-8').strip(), date_fmt).timestamp() # key missing ingredient was whitespace at the end of string
                    #str2date = lambda s: float(np.datetime64(s))
                    #print(str2date('2018-03-22T11:57:03.456'))
                    if len(args.datetime) > 1: # date and time in separate columns
                        # read and combine datetime columns
                        date_fmt = '%d-%m-%y %H:%M:%S' # BlueFors format; note that this is not the same as our special measure format
                        date_columns = np.loadtxt(fname, usecols=args.datetime, delimiter=args.delimiter, dtype=str)
                        datetimes = np.array([datetime.datetime.strptime(' '.join((d, t)).strip(), date_fmt).timestamp() for d, t in date_columns])

                        # read data columns and fill plotvs[0] with datetimes
                        data_columns = np.genfromtxt(fname, delimiter=args.delimiter, skip_header=args.header)
                        data_columns[:, plotvs[0]] = datetimes[:data_columns.shape[0]] # places datetimes in first plotvs columns (default 0); may want to use first datetime column in future; try modifying plotvs[0] itself to store column of args.datetime[0]; this will preserve labeling and such
                    else:
                        data_columns = np.loadtxt(fname, delimiter=args.delimiter, converters={args.datetime[0]: str2date})
                    #data_columns = np.genfromtxt(fname, delimiter=args.delimiter, skip_header=args.header, converters={args.datetime: str2date})
                    #print(data_columns[0, :])
                else:
                    data_columns = np.loadtxt(f, delimiter=args.delimiter)
            except ValueError:
                # prepare for chaining
                #if n > 0:
                #    nlast += data_columns.shape[0]
                #else:
                #    nlast = 0
                if first_file: # added in haste on 11-15-2018 to allow file filtering based on presence of specific column header
                    nlast = 0
                    first_file = False
                else:
                    nlast += data_columns.shape[0]
                # method 1: replace first column with row number
                #data_columns = np.loadtxt(f, usecols=range(1, totcols), delimiter=args.delimiter)
                #data_columns = np.column_stack((np.arange(data_columns.shape[0])+nlast*args.chain, data_columns))

                # method 2: looping try block to skip columns from the left until a valid dataset is loaded (a hack to skip two-column date formats)
                #skipcols = 1
                #while True:
                    #try:
                        #data_columns = np.loadtxt(f, usecols=range(skipcols, totcols), delimiter=args.delimiter)
                        #data_columns = np.column_stack(skipcols*[np.arange(data_columns.shape[0])+nlast*args.chain] + [data_columns])
                        #break
                    #except ValueError:
                        #skipcols += 1
                        #continue

                # method 3: skips string values but can do much more with customization; only fill first column with row number
                data_columns = np.genfromtxt(fname, delimiter=args.delimiter, skip_header=args.header)
                data_columns[:, 0] = np.arange(data_columns.shape[0])+nlast*args.chain
                if not (args.xlabel and plotvs[0]==0):
                    labels[0] = 'points'

                # method 4: currently doesn't work due to matplotlib problem dealing with iterables
                #from matplotlib.dates import datestr2num
                #data_columns = np.genfromtxt(fname, delimiter=args.delimiter, converters={0: datestr2num})
            rows = data_columns.shape[0]
            stdoutstr = '{}\t{} rows\t{} columns'
            if rows == 0:
                stdoutstr += '\tskipping empty file'
                print(stdoutstr.format(fname, rows, totcols))
                data.append(None) # need empty element for correct file counting
                continue

            if args.deleterows:
                print('deleting specified rows.. ', end='')
                for row in args.deleterows:
                    data_columns[row, usecols] = np.nan
                print('{} rows deleted'.format(len(args.deleterows)))
            if args.trimspikes:
                for trimcol in args.trimspikes:
                    print('trimming spikes in column {}.. '.format(trimcol), end='')
                    def nan_get(y): # returns boolean array and function yielding indices of NaNs
                        return np.isnan(y), lambda z: z.nonzero()[0]
                    col = data_columns[:, trimcol]
                    nans, idx = nan_get(col) # find the NaNs
                    col[nans] = np.interp(idx(nans), idx(~nans), col[~nans]) # replace NaNs with interpolated values
                    avg = np.convolve(col, np.array([1,1,0,1,1])/4, 'same')
                    # remove data which differs from the rolling average
                    spikes = np.fabs((col-avg)/avg) > 0.02 # by 2%
                    #plt.plot(col)
                    #plt.plot(avg, '.')
                    #plt.plot((col-avg)/avg, '-.')
                    data_columns[spikes, trimcol] = np.nan
                    print('{} elements removed'.format(np.sum(spikes)))
            if args.monotonic:
                print('clipping data with non-monotonic sweep direction.. ', end='')
                xdata = data_columns[:, plotvs[0]]
                direction = np.sign(np.mean(np.diff(xdata)))
                if direction > 0:
                    cut_rows = xdata != np.maximum.accumulate(xdata)
                else:
                    cut_rows = xdata != np.minimum.accumulate(xdata)
                for col in usecols:
                    data_columns[cut_rows, col] = np.nan
                print('{} elements deleted'.format(cut_rows.size))
            if args.combinecols:
                data_columns[:, usecols] = np.sqrt(np.sum(data_columns[:, args.combinecols]**2, axis=1, keepdims=True))
                if not args.ylabels:
                    newlabels = list(labels) # convert to list for mutability
                    newlabels[args.combinecols[0]] = '$\sqrt{{[{}]^2+[{}]^2}}$'.format(*labels[args.combinecols]) # overwrite new label
                    labels = np.asarray(newlabels, dtype=str) # convert back to array (cannot otherwise change size of array or string length)
            if args.yshift:
                data_columns[:, usecols] += args.yshift
            if args.boxcar:
                for col in usecols:
                    data_columns[:, col] = np.convolve(data_columns[:, col], np.ones((args.boxcar,))/args.boxcar, 'same')
                    data_columns[:args.boxcar, col] = np.nan # mask elements outside valid boxcar range
                    data_columns[-args.boxcar:, col] = np.nan # mask elements outside valid boxcar range
            if args.xboxcar:
                data_columns[:, plotvs[0]] = np.convolve(data_columns[:, plotvs[0]], np.ones((args.xboxcar,))/args.xboxcar, 'same')
                data_columns[:args.xboxcar, plotvs[0]] = np.nan # mask elements outside valid boxcar range
                data_columns[-args.xboxcar:, plotvs[0]] = np.nan # mask elements outside valid boxcar range
            if args.x2boxcar:
                data_columns[:, plotvs[1]] = np.convolve(data_columns[:, plotvs[1]], np.ones((args.x2boxcar,))/args.x2boxcar, 'same')
                data_columns[:args.x2boxcar, plotvs[1]] = np.nan # mask elements outside valid boxcar range
                data_columns[-args.x2boxcar:, plotvs[1]] = np.nan # mask elements outside valid boxcar range
            if args.lorentz: # lorentzian blur of data
                for col in usecols:
                    data_columns[:, col] = np.convolve(data_columns[:, col], lorentzian(args.lorentz), 'same')
                    data_columns[:args.lorentz, col] = np.nan # mask elements outside valid lorentz range
                    data_columns[-args.lorentz:, col] = np.nan # mask elements outside valid lorentz range
            if args.interpolate and (lineplot_mode or len(args.interpolate) == 1): # colorplot interpolation must occur after all file data has been collected
                interpX = args.interpolate[0] # only use X1, X2, Nx group (ignore beyond semicolon in argument)
                if len(interpX) < 3: # use data limits for interpolation limits
                    X1, X2, Nx = np.nanmin(data_columns[:, plotvs[0]]), np.nanmax(data_columns[:, plotvs[0]]), int(interpX[0])
                else: # use user-provided limits
                    X1, X2, Nx = interpX
                    Nx = int(Nx)
                data_fn = sci.interp1d(data_columns[:, plotvs[0]], data_columns[:, usecols], axis=0)
                        #axis=0, fill_value='extrapolate')
                raw_data_columns = data_columns
                data_columns = np.empty((Nx, totcols))
                data_columns[:, plotvs[0]] = np.linspace(X1, X2, Nx)
                data_columns[:, usecols] = data_fn(data_columns[:, plotvs[0]])
                stdoutstr += '\tinterpolated and resampled to {} rows'.format(Nx)
            for col in args.symmetrize:
                data_columns[:, col] = 0.5*(data_columns[:, col] + data_columns[::-1, col])
                axes[usecols.index(col)].set_title('symmetrized')
            for col in args.antisymmetrize:
                data_columns[:, col] = 0.5*(data_columns[:, col] - data_columns[::-1, col])
                axes[usecols.index(col)].set_title('antisymmetrized')
            if args.divcol:
                data_columns[:, usecols] *= 1/data_columns[:, args.divcol][:, None]
            if args.xdivcol:
                data_columns[:, plotvs[0]] *= 1/data_columns[:, args.xdivcol]
            if args.derivative:
                if type(args.derivative) is list:
                    derivcols = args.derivative
                elif args.derivative is True:
                    derivcols = usecols.copy()
                data_columns[:, derivcols] = np.gradient(data_columns[:, derivcols], np.squeeze(data_columns[:, plotvs[0]]), axis=0)
                if not args.ylabels:
                    for derivcol in [derivcol for derivcol in derivcols if derivcol in usecols]:
                        labels[derivcol] = 'derivative of {}'.format(labels[derivcol])
#            if args.meansubtract:
#                data_columns[:, usecols] -= np.nanmean(data_columns[:, usecols], axis=0)
            if args.integrate:
                if type(args.integrate) is list:
                    intcols = args.integrate
                elif args.integrate is True:
                    intcols = usecols.copy()
                data_columns[:, intcols] = scg.cumtrapz(data_columns[:, intcols], data_columns[:, plotvs[0]], axis=0, initial=0)
                if not args.ylabels:
                    for intcol in [intcol for intcol in intcols if intcol in usecols]:
                        labels[intcol] = 'integral of {}'.format(labels[intcol])
            if args.invert is not False or args.xinvert:
                if type(args.invert) is list:
                    invcols = args.invert
                elif args.invert is True:
                    invcols = usecols.copy()
                else:
                    invcols = []
                if args.xinvert:
                    invcols.extend(plotvs)
                    invcols = list(set(invcols))
                    if not args.xlabel:
                        for invcol in [invcol for invcol in invcols if invcol in plotvs]:
                            labels[invcol] = '1 / {}'.format(labels[invcol])
                data_columns[:, invcols] = 1/data_columns[:, invcols]
                if not args.ylabels:
                    for invcol in [invcol for invcol in invcols if invcol in usecols]:
                        labels[invcol] = '1 / {}'.format(labels[invcol])
            data.append(data_columns)
        print(stdoutstr.format(fname, rows, totcols))

        for p, (ax, col, label) in enumerate(zip(it.cycle(axes),
                                  #data[n][:, args.usecols[:data[n].shape[1]-3]].T,
                                  data[n][:, usecols].T,
                                  labels[usecols])):
            formatter = EngFormatter(sep=u'\N{THIN SPACE}') # format axes using engineering notation
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)

            xdata = data[n][:, plotvs[0]].copy()
            if colorplot_mode or args.waterfall:
                ydata = data[n][:, plotvs[1]].copy()
                xdata = data[n][:, plotvs[0]] # don't use a copy; this allows xdata manipulations through to colorplot
                col   = data[n][:, usecols[p]] # don't copy; allows direct manipulation to propogate to colorplot
            if args.traces:
                tr_data = data[n][:, plotvs[0]]
                xdata = data[n][:, plotvs[1]] #don't use a copy; this allows xdata manipulations through to colorplot
            if args.xmult:
                if len(args.xmult) > 1 and len(args.xmult[0]) > 1:
                    xmult = args.xmult[n][p]
                elif len(args.xmult) > 1:
                    xmult = args.xmult[n][0]
                elif len(args.xmult[0]) > 1:
                    xmult = args.xmult[0][p]
                else:
                    xmult = args.xmult[0][0]
                xdata *= xmult
            if args.xshifts:
                if len(args.xshifts) > 1 and len(args.xshifts[0]) > 1:
                    xshifts = args.xshifts[n][p]
                elif len(args.xshifts) > 1:
                    xshifts = args.xshifts[n][0]
                elif len(args.xshifts[0]) > 1:
                    xshifts = args.xshifts[0][p]
                else:
                    xshifts = args.xshifts[0][0]
                xdata += xshifts
            if args.fold:
                xdata = np.fabs(xdata)
            if args.xnorm:
                xmin, xmax = np.nanmin(xdata), np.nanmax(xdata)
                xdata -= xmin
                xdata *= 1/(xmax-xmin)
            if not args.ylabels and args.divcol:
                label += ' / {}'.format(labels[args.divcol])
            if args.normy:
                ymin, ymax = np.nanmin(col), np.nanmax(col)
                col -= ymin
                col *= 1/(ymax-ymin)
                if not args.ylabels:
                    label = 'normalized ' + label
            if args.ymult:
                if len(args.ymult) > 1 and len(args.ymult[0]) > 1:
                    ymult = args.ymult[n][p]
                elif len(args.ymult) > 1:
                    ymult = args.ymult[n][0]
                elif len(args.ymult[0]) > 1:
                    ymult = args.ymult[0][p]
                else:
                    ymult = args.ymult[0][0]
                col *= ymult
                if not args.ylabels:
                    label += r'$\times {:.3g}$'.format(ymult)
            if args.zeroy:
                col -= np.nanmin(col)
                if not args.ylabels:
                    label += ' - ymin'
            if args.orezy:
                col -= np.nanmax(col)
                if not args.ylabels:
                    label += ' - ymax'
            if args.meansubtract:
                col -= np.nanmean(col)
                if not args.ylabels:
                    label += ' - ymean'
            if args.yshifts:
                if len(args.yshifts)>1:
                    col += args.yshifts[n]
                else: # Waterfall plot ordered by file sequence
                    if colorplot_mode:
                        data[n][:, plotvs[1]] += n*args.yshifts[0]
                    else:
                        col += n*args.yshifts[0]
            if args.hysteresis:
                if n % 2: # odd, and therefore second file of each pair
                    if colorplot_mode:
                        data[n][:, usecols]   = data[n][:, usecols] - np.flip(data[n-1][:, usecols])
                        data[n-1][:, usecols] = np.nan
                    else:
                        col = col - lastcol
                else:
                    lastcol = np.flip(col)
                    continue
            if args.waterfall: # Waterfall plot mode
                col += args.waterfall*ydata # separate lines by second x-value times multiplier (default 1)
            if args.traces: # trace plot mode
                cidx = np.asarray([np.nanargmin(np.fabs(tr_value-tr_data)) for tr_value in args.traces])
                traces = col[cidx]
                col = np.full_like(col, np.nan)
                col[cidx] = traces
            if colorplot_mode: # Colorplot mode
                xlabel = labels[plotvs[0]]
                ylabel = labels[plotvs[1]]
                clabel = label
                if args.cornplot: # Special colorplot for data versus rotating field
                    xdata, ydata = ydata*np.sin(xdata*np.pi/180), ydata*np.cos(xdata*np.pi/180)
                    ylabel = 'parallel magnetic field (T)'
                    xlabel = 'perpendicular magnetic field (T)' 
                if (n == 0 or (args.hysteresis and n==1)) and p == 0:
                    crange = len(axes)*[()]
                    gcmin, gcmax = np.full(len(axes), np.nan), np.full(len(axes), np.nan)
                if not args.crange:
                    cmin, cmax = np.nanmin(col), np.nanmax(col)
                    gcmin[p], gcmax[p] = np.nanmin([cmin, gcmin[p]]), np.nanmax([cmax, gcmax[p]]) # set global min/max for this panel
                    crange[p] = gcmin[p], gcmax[p] # use global upper and lower limits considering all files
                    print('{} data ranges from approx {:.6g} to {:.6g}; global limits {:.6g} to {:.6g}'.format(clabel, cmin, cmax, *crange[p]))
                elif n == 0:
                    if len(args.crange) > 1:
                        crange[p] = args.crange[p]
                    else:
                        crange[p] = args.crange[0]
                # Ends up labeling based on last file in list
                ax.set_ylabel(ylabel)
                ax.set_xlabel(xlabel)
                
                #if not args.interpolate or len(args.interpolate) == 1: # colorplot based on raw data (scattered points); interpolation handled after collecting all files
                    #cplot_props = dict(s=colorplot_marker_size, edgecolors='none')
                    #if args.square:
                        #cplot_props['marker'] = 's'
                    #cplot = ax.scatter(xdata, ydata, c=col, cmap=args.cmap, vmin=crange[p][0], vmax=crange[p][1], **cplot_props)
                    ##ax.autoscale(enable=True, axis='both', tight=True)
                    #if n == 0: # add colorbar first time around
                        #cbar = fig.colorbar(cplot, ax=ax)
                        #cbar.set_label(clabel)
                    #cbar.set_clim(*crange[p]) # update colorbar limits (only scales colormap)
                    #cbar.draw_all() # will not redraw ticks/limits without this
            if lineplot_mode:
                if args.colorbyorder:
                    if type(args.colorbyorder) is str:
                        cm = mcm.get_cmap(args.colorbyorder)
                    else:
                        cm = mcm.get_cmap(args.cmap)
                    lineprops['c'] = cm(n/(len(files)-1)) 
                elif args.color:
                    lineprops['c'] = next(args.color)
                elif args.ls:
                    lineprops['ls'] = next(args.ls)
                if args.datetime is not False: # date_plot mode
                    datetimes = [datetime.datetime.fromtimestamp(ts) for ts in xdata]
                    ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%b %d, %H:%M')) # month, day hours:minutes
                    if args.timedelta is not False:
                        d0 = datetimes[0]
                        dt0 = datetime.timedelta(days=d0.day, hours=d0.hour, minutes=d0.minute, seconds=d0.second, microseconds=d0.microsecond)
                        datetimes = [dt - dt0 for dt in datetimes]
                        ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M:%S')) # ignore dates and just plot versus elapsed time
                    datetimes = mpl.dates.date2num(datetimes)
                    line = ax.plot_date(datetimes, col, label=fname, **lineprops)
                    for l in ax.get_xticklabels():
                        l.set_rotation(45)
                else:
                    line = ax.plot(xdata, col, label=fname, **lineprops)
                                  
                # Ends up labeling based on last file in list
                ax.set_ylabel(label)
                if args.traces:
                    ax.set_xlabel(labels[plotvs[1]])
                else:
                    ax.set_xlabel(labels[plotvs[0]])
            if args.logx:
                ax.set_xscale('log', nonposx='clip')
            if args.logy:
                ax.set_yscale('log', nonposy='clip')
            if args.arrows:
                ind = col.shape[0]//2
                ax.annotate('', xy=(xdata[ind+args.arrows], col[ind+args.arrows]),
                            xytext=(xdata[ind], col[ind]),
                            arrowprops=dict(arrowstyle='->', ec=line[-1].get_color(), linewidth=4), size=24)
            if args.dashes and np.mean(np.diff(xdata)) < 0:
                line[-1].set_ls('--')
            ax.autoscale(enable=True, axis='x', tight=True)
            if args.xlim:
                ax.set_xlim(args.xlim)
            if args.ylims:
                if len(args.ylims) > 1:
                    ylim = args.ylims[p]
                else:
                    ylim = args.ylims[0]
                ax.set_ylim(ylim)

    if colorplot_mode:
        # use data from all files for 2D interpolation
        raw_data = np.concatenate([d for d in data if d is not None]) # skip empty entries (from emtpy files)

        if args.interpolate and len(args.interpolate) > 1:
            # set interpolation limits
            interpX, interpY = args.interpolate # ((X1, X2, Nx), (Y1, Y2, Ny))
            if len(interpX) < 3: # use data limits for interpolation limits
                X1, X2, Nx = np.nanmin(raw_data[:, plotvs[0]]), np.nanmax(raw_data[:, plotvs[0]]), int(interpX[0])
            else: # use user-provided limits
                X1, X2, Nx = interpX
                Nx = int(Nx)
            if len(interpY) < 3: # use data limits for interpolation limits
                Y1, Y2, Ny = np.nanmin(raw_data[:, plotvs[1]]), np.nanmax(raw_data[:, plotvs[1]]), int(interpY[0])
            else: # use user-provided limits
                Y1, Y2, Ny = interpY
                Ny = int(Ny)
            gridY, gridX = np.mgrid[Y1:Y2:Ny*1j, X1:X2:Nx*1j]

            # interpolate and plot
            cplot_props = dict(cmap=args.cmap, origin='lower', aspect='auto', interpolation='none')
            for p, (ax, usecol, clabel) in enumerate(zip(it.cycle(axes), usecols, labels[usecols])):
                interp_data = sci.griddata(raw_data[:, plotvs], raw_data[:, usecol], (gridX, gridY))
                if args.crange:
                    cplot_props.update(vmin=crange[p][0], vmax=crange[p][1])
                if args.csigma:
                    cmean, cstd = interp_data.mean(), interp_data.std()
                    cplot_props.update(vmin=cmean-args.csigma*cstd, vmax=cmean+args.csigma*cstd)
                    print('column \'{:s}\':\tmean = {:.5g}\tstddev = {:.5g}\tcrange = {:.5g} to {:.5g}'.format(clabel, cmean, cstd, cplot_props['vmin'], cplot_props['vmax']))
                cplot = ax.imshow(interp_data, extent=(X1, X2, Y1, Y2), **cplot_props)
                cbar = fig.colorbar(cplot, ax=ax)
                cbar.set_label(clabel)
        else:
            cplot_props = dict(s=colorplot_marker_size, edgecolors='none')
            if args.square:
                cplot_props['marker'] = 's'
            for p, (ax, usecol, clabel) in enumerate(zip(it.cycle(axes), usecols, labels[usecols])):
                if args.csigma:
                    cmean, cstd = raw_data[:, usecol].mean(), raw_data[:, usecol].std()
                    vmin, vmax = cmean-args.csigma*cstd, cmean+args.csigma*cstd
                    print('column \'{:s}\':\tmean = {:.5g}\tstddev = {:.5g}\tcrange = {:.5g} to {:.5g}'.format(clabel, cmean, cstd, vmin, vmax))
                else:
                    vmin, vmax = crange[p][0], crange[p][1]
                cplot = ax.scatter(raw_data[:, plotvs[0]], raw_data[:, plotvs[1]], c=raw_data[:, usecol], cmap=args.cmap, vmin=vmin, vmax=vmax, **cplot_props)
                cbar = fig.colorbar(cplot, ax=ax)
                cbar.set_label(clabel)
                #ax.autoscale(enable=True, axis='both', tight=True)

    fig.tight_layout()
    if lineplot_mode and not args.nolegend:
        leg = ax.legend(loc=0)
        #leg = fig.legend(loc='lower center')
        try:
            leg.set_draggable(True)
        except:
            leg.draggable()
    plt.show()

# call first instance of main()
main()
