#
# Plots pressure distributions from pressure data .txt as image montages for
# each 3 dimensions.
#

IMG_WSPACE = 0.07       # space between montage columns
IMG_HSPACE = 0.20       # space between montage rows
IMG_WIDTH = 1.0         # width of individual images
MONTAGE_DPI = 180       # montage DPI

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys


nargs = len(sys.argv)
if (nargs != 2 and nargs != 4):
    print("Usage: python [script] [pressure data] [max pressure (optional)] [min pressure (optional)]")
    exit()
fileIn = sys.argv[1]
highP = np.nan
lowP = np.nan
if (nargs == 4):
    highP = float(sys.argv[2])
    lowP = float(sys.argv[3])

data = np.loadtxt(fname=fileIn, delimiter='\t', skiprows=1, dtype=np.float32)
X = data[:,0]
Y = data[:,1]
Z = data[:,2]
P = data[:,3]
nx = len(np.unique(X))
ny = len(np.unique(Y))
nz = len(np.unique(Z))


#
# Returns figure object of pressure distribution montage for given axis tuple D,
# node coordinates C and pressure P. 
# D is one of [nx,ny,nz], [ny,nz,nx], [nz,ny,nx], where n* is the number of nodes
# along given axis.
#
def getImages(D, C, P):
    # Compute proper number of rows and columns for the montage
    rows = np.floor(np.sqrt(D[0]))
    cols = D[0] / rows
    while True:
        if (cols.is_integer()):
            rows = int(rows)
            cols = int(cols)
            break
        else:
            rows = rows - 1
            cols = D[0] / rows
    print('  Montage size: %d x %d' %(rows, cols))
    
    # Fill images array
    images = np.zeros(D)
    Cu = np.unique(C)
    for i in range(0, len(Cu)):
        I = np.where(C == Cu[i])[0]
        images[i] = np.reshape(P[I], [D[1],D[2]], order='F')
    
    # Calculate full montage figure dimensions.
    aspectRatio = D[1] / D[2]
    width = cols*IMG_WIDTH + (cols-1)*IMG_WSPACE
    height = rows * aspectRatio * IMG_WIDTH + (rows-1)*IMG_HSPACE
    fig = plt.figure(figsize = (width, height))
    matplotlib.rc('axes', linewidth=0.0)    # hide image borders
    print("  Full figure dimensions %f x %f." %(width, height))
    
    maxP = highP if not np.isnan(highP) else np.max(P)
    minP = lowP if not np.isnan(lowP) else np.min(P)
    print("  Pressure max. %f (red), min. %f (blue)." %(maxP, minP))
    
    for i in range(0, len(images)):
        ax = fig.add_subplot(rows, cols, i+1)
        ax.imshow(images[i], cmap=plt.get_cmap('jet'), interpolation='bilinear', vmin=minP, vmax=maxP)
        # ax.imshow(images[i], cmap=plt.get_cmap('jet'), interpolation='bilinear')
        plt.tick_params(axis='both', which='both', left=False, right=False, bottom=False, labelbottom=False, labelleft=False)
        s = 'x = %f' %np.unique(X)[i]
        plt.xlabel(s, fontsize=(int(40/cols)))      # under image
        # plt.title(s, fontsize=(int(40/cols)))     # above image
    
    fig.subplots_adjust(left=0.0, right=1.0, bottom=2*IMG_WSPACE/height, top=1.0, wspace=IMG_WSPACE, hspace=IMG_HSPACE)
    return fig, maxP, minP


# X slices
print('- Creating montage for slices along x axis.')
fig, maxP, minP = getImages([nx, ny, nz], X, P)
fileOut = fileIn.split('.')[0] + '_x_maxminP_%.2f_%.2f.eps' %(maxP, minP)
fig.savefig(fileOut, dpi=MONTAGE_DPI)
print('')

# Y slices
print('- Creating montage for slices along y axis.')
fig, maxP, minP = getImages([ny, nz, nx], Y, P)
fileOut = fileIn.split('.')[0] + '_y_maxminP_%.2f_%.2f.eps' %(maxP, minP)
fig.savefig(fileOut, dpi=MONTAGE_DPI)
print('')

# Z slices
print('- Creating montage for slices along z axis.')
fig, maxP, minP = getImages([nz, ny, nx], Z, P)
fileOut = fileIn.split('.')[0] + '_z_maxminP_%.2f_%.2f.eps' %(maxP, minP)
fig.savefig(fileOut, dpi=MONTAGE_DPI)
print('')
