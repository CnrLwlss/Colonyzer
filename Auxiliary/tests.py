import os,numpy,itertools
from colonyzer2 import *
from scipy import optimize as op
from scipy import ndimage

#LATESTIMAGE=os.path.realpath("Data\\384\\DLR00012647-2009-07-04_09-35-20.jpg")
#LATESTIMAGE=os.path.realpath("F:\\Downloads\\scan_images_registered_cropped_partial\\p0000071_1N_Hap_0144.tif")
#LATESTIMAGE=os.path.realpath("F:\\PLOSGeneticsImages\\pdump_cdc13\\DLR00012257-2009-11-27_17-22-47.JPG")
#EARLIESTIMAGE=os.path.realpath("F:\\PLOSGeneticsImages\\pdump_cdc13\\DLR00012257-2009-11-23_18-33-06.JPG")
EARLIESTIMAGE=os.path.realpath("Data\\384\\DLR00012255-2009-11-23_18-42-06.JPG")
LATESTIMAGE=os.path.realpath("Data\\384\\DLR00012255-2009-11-27_17-26-18.JPG")

imN,arrN=openImage(LATESTIMAGE)

if "EARLIESTIMAGE" not in locals():
    im0,arr0=imN,arrN
    arr=arrN
else:
    im0,arr0=openImage(EARLIESTIMAGE)
    arr=numpy.maximum(0,arrN-arr0)

#nx,ny=12,8
nx,ny=24,16
windowFrac=0.25
smoothWindow=0.13
showPlt=True
pdf=None

acmedian=True
rattol=0.1
nsol=144
verbose=True

### 1: Estimate height and width of spots by examining inter-peak distances in autocorrelation function

# Generate windowed mean intensities, scanning along x and y axes
# Estimate spot diameter, assuming grid takes up most of the plate
diam=min(float(arr.shape[0])/ny,float(arr.shape[1])/nx)
window=int(round(diam*windowFrac))
sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-window]):numpy.min([arr.shape[1],dx+window])]) for dx in range(0,arr.shape[1])],dtype=numpy.float)
sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-window]):numpy.min([arr.shape[0],dy+window]),0:arr.shape[1]]) for dy in range(0,arr.shape[0])],dtype=numpy.float)
# Smooth intensities to help eliminate small local maxima
#sumx=ndimage.gaussian_filter1d(sumx,2.5)
#sumy=ndimage.gaussian_filter1d(sumy,2.5)

# Look at autocorrelation for first estimate of distance between spots
maximay=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumy))))==-2)[0]
maximax=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0]
if acmedian:
    # Note that median inter-peak distance is more robust here
    # Mean is thrown by outliers: gives poor initial guess for optimisation routine
    dy=int(round(numpy.median(numpy.diff(maximay))))
    dx=int(round(numpy.median(numpy.diff(maximax))))
else:
    dy=int(round(maximay[0]))
    dx=int(round(maximax[0]))

# Some plate inoculation patterns (e.g. alternate columns of fit and sick strains in miniQFA) break this ACF procedure (skip every second column)
# However, the row ACF estimate is still robust.  First estimate: choose the smallest value
dmin=min(dy,dx)
dy,dx=dmin,dmin

# Given shape of grid and size of tiles, search range for x0 and y0, assuming grid parallel to plate edges
ry=arr.shape[0]-ny*dy
rx=arr.shape[1]-nx*dx

### 2: Find x0,y0 that maximises contrast between grid points and intermediate points (e.g. aligns grid with peaks and troughs),
### in horizontal and vertical ACF.  Optimise with fixed dx, for a range of dx.

dd=int(round(0.01*dx))
xtest=list(range(dx-dd,min(int(round(arr.shape[1]/nx))-1,dx+dd)))
ytest=list(range(dy-dd,min(int(round(arr.shape[0]/ny))-1,dy+dd)))
xvals=[fitProjection(sumx,int(round(dxval)),nx,1.0,1.0,True) for dxval in xtest]
yvals=[fitProjection(sumy,int(round(dyval)),ny,1.0,1.0,True) for dyval in ytest]
xind=numpy.argmin([x[1] for x in xvals])
yind=numpy.argmin([y[1] for y in yvals])

xbest=xvals[xind][0]
dx=xtest[xind]
ybest=yvals[yind][0]
dy=ytest[yind]

grd,gp=makeGrid((ybest,xbest),ny,nx,dy=dy,dx=dx,theta=0,makeGaps=False)
candy,candx=list(zip(*grd))
plotAC(sumy,sumx,candy,candx,maximay,maximax,pdf=pdf,main="Projection Estimate")

### 3: Optimise positions, first just optimise x0,y0 then update, optimising x0,y0,dx and theta

smarr=ndimage.filters.gaussian_filter(arr,dmin/10.0)
showIm(smarr)

checkvecs=[list(range(ry)),list(range(rx))]
checkpos=list(itertools.product(*checkvecs))

# Assume we can see the edges of the plate in the image (bright enough to make a peak in the smoothed intensities
peaksy=numpy.where(numpy.diff(numpy.sign(numpy.diff(sumy)))==-2)[0]
peaksx=numpy.where(numpy.diff(numpy.sign(numpy.diff(sumx)))==-2)[0]
corner=[peaksy[0],peaksx[0]]

com=ndimage.measurements.center_of_mass(arr)
com=[int(round(x)) for x in com]

#bounds=[(peaksy[0]+dy,ry),(peaksx[0]+dx,rx),(0.8*min(dy,dx),1.2*max(dy,dx)),(-5,5)]
bounds=[(max(int(round(dy/2.0)),ybest-2.0*dy),min(int(round(arr.shape[0]-(ny-1)*dy)),ybest+2.0*dy)),
         (max(int(round(dx/2.0)),xbest-2.0*dx),min(int(round(arr.shape[1]-(nx-1)*dx)),xbest+2.0*dx)),
         (0.8*min(dy,dx),1.2*max(dy,dx)),
         (-5,5)]

def makeOptAll(arr,ny,nx,bounds,sampfrac=0.35):
    def optfun(xvs):
        xrs=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,xvs)]
        #res=-1*checkPos(arr,ny,nx,xrs[0:2],xrs[2],xrs[2],xrs[3],sampfrac=sampfrac)
        res=-1*checkPoints(arr,ny,nx,xrs[0:2],xrs[2],xrs[2],xrs[3],True)
        return (res)
    return optfun

def makeOptPos(arr,ny,nx,dx_norm,theta_norm,bounds,sampfrac=0.35):
    dx=bounds[2][0]+dx_norm*(bounds[2][1]-bounds[2][0])
    theta=bounds[3][0]+theta_norm*(bounds[3][1]-bounds[3][0])
    def optfun(xvs):
        xrs=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds[0:2],xvs)]
        #res=-1*checkPos(arr,ny,nx,xrs,dx,dx,theta,sampfrac=sampfrac)
        res=-1*checkPoints(arr,ny,nx,xrs,dx,dx,theta,True)
        return (res)
    return optfun

### 3a: Optimise x0,y0 assuming grid parallel to edges of image and tile dimensions known (ignore theta, dx and dy)
optmess=False

dx_norm=(min(dx,dy)-bounds[2][0])/(bounds[2][1]-bounds[2][0])
optpos=makeOptPos(smarr,ny,nx,dx_norm,0.5,bounds)

# For even sampling: nsol = Nsamps**Ndim   
x0s=[sobol.i4_sobol(2,i)[0] for i in range(nsol)]
firstsol=[op.minimize(optpos,x0=x,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds[0:2]],jac=False,options={'eps':0.005,'disp':optmess,'maxiter':5}) for x in x0s]
solpos=firstsol[numpy.argmin([sol.fun for sol in firstsol])]

solnpos=solpos.x
solnpos=numpy.append(solnpos,[dx_norm,0.5])

soln=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,solnpos)]
candy,candx=grid(soln,ny,nx)
plotAC(sumy,sumx,candy,candx,maximay,maximax,pdf=pdf,main="First pass")
ybest,xbest,dx,theta=soln

### 3b: Optimise all parameters, using above as initial guess
optall=makeOptAll(smarr,ny,nx,bounds)

sol=op.minimize(optall,x0=solnpos,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds],jac=False,options={'eps':0.005,'disp':optmess})
soln=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,sol.x)]

##    x0s.append(xguess[0:2])
##    x0s.append(sol1.x[0:2])
##    firstpass=[optpos(x) for x in x0s]
##    firstguess=x0s[numpy.argmin(firstpass)]
##    sol2=op.minimize(optpos,x0=firstguess,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds[0:2]],jac=False,options={'eps':0.005,'disp':optmess,'gtol':0.1})
##    soln2=sol2.x
##    soln2=numpy.append(soln2,sol1.x[2:])
##    sol=op.minimize(optall,x0=soln2,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds],jac=False,options={'eps':0.005,'disp':optmess,'gtol':0.1})

if verbose:
    print(("Optimisation: "+sol.message+"\n\n"))
    
candy,candx=grid(soln,ny,nx)

# Output some plots
if showPlt:
    plotAC(sumy,sumx,candy,candx,maximay,maximax,pdf=pdf,main="Solution")

xguess=list([0.5 for b in bounds])
init=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,xguess)]

locationsN=locateCultures([int(round(cx-dx/2.0)) for cx in candx],[int(round(cy-dy/2.0)) for cy in candy],dx,dy,arr,nx,ny,update=True)
imloc=threshPreview(arr,-99,locationsN)
imloc.show()
