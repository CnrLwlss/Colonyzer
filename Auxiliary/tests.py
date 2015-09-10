import os,numpy,itertools
from colonyzer2 import *
from scipy import optimize as op

LATESTIMAGE=os.path.realpath("F:\\Colonyzer\\Auxiliary\\Data\\endpoints\\384\\X000001_001_004_2014-06-19_09-28-46.JPG")
im,arr=openImage(LATESTIMAGE)

nx,ny=24,16
windowFrac=0.25
smoothWindow=0.13
showPlt=True
pdf=None

acmedian=True
rattol=0.1
nsol=36
verbose=True

'''Automatically search for best estimate for location of culture array (based on culture centres, not top-left corner).'''
# Generate windowed mean intensities, scanning along x and y axes
# Estimate spot diameter, assuming grid takes up most of the plate
diam=min(float(arr.shape[0])/ny,float(arr.shape[1])/nx)
window=int(round(diam*windowFrac))
sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-window]):numpy.min([arr.shape[1],dx+window])]) for dx in range(0,arr.shape[1])],dtype=numpy.float)
sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-window]):numpy.min([arr.shape[0],dy+window]),0:arr.shape[1]]) for dy in range(0,arr.shape[0])],dtype=numpy.float)
# Smooth intensities to help eliminate small local maxima
sumx=ndimage.gaussian_filter1d(sumx,2.5)
sumy=ndimage.gaussian_filter1d(sumy,2.5)

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
# However, the row ACF estimate is still robust.  In the case of disagreement, choose the smallest value
rat=float(dy)/float(dx)
if rat > (rattol+1.0) or rat < 1.0/(rattol+1.0):
    dmin=min(dy,dx)
    dy,dx=dmin,dmin
    
ry=arr.shape[0]-((ny)*dy)
rx=arr.shape[1]-((nx)*dx)

checkvecs=[range(ry),range(rx)]
checkpos=list(itertools.product(*checkvecs))

# Assume we can see the edges of the plate in the image (bright enough to make a peak in the smoothed intensities
peaksy=numpy.where(numpy.diff(numpy.sign(numpy.diff(sumy)))==-2)[0]
peaksx=numpy.where(numpy.diff(numpy.sign(numpy.diff(sumx)))==-2)[0]
corner=[peaksy[0],peaksx[0]]

com=ndimage.measurements.center_of_mass(arr)
com=[int(round(x)) for x in com]

bounds=[(0,ry),(0,rx),(0.8*min(dy,dx),1.2*max(dy,dx)),(-5,5)]

def makeOptAll(arr,ny,nx,bounds,sampfrac=0.35):
    def optfun(xvs):
        xrs=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,xvs)]
        res=-1*checkPos(arr,ny,nx,xrs[0:2],xrs[2],xrs[2],xrs[3],sampfrac=sampfrac)
        return (res)
    return optfun

def makeOptPos(arr,ny,nx,dx_norm,theta_norm,bounds,sampfrac=0.35):
    dx=bounds[2][0]+dx_norm*(bounds[2][1]-bounds[2][0])
    theta=bounds[3][0]+theta_norm*(bounds[3][1]-bounds[3][0])
    def optfun(xvs):
        xrs=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds[0:2],xvs)]
        res=-1*checkPos(arr,ny,nx,xrs,dx,dx,theta,sampfrac=sampfrac)
        return (res)
    return optfun

optmess=False

dx_norm=(min(dx,dy)-bounds[2][0])/(bounds[2][1]-bounds[2][0])
optpos=makeOptPos(arr,ny,nx,dx_norm,0.5,bounds)

# For even sampling: nsol = Nsamps**Ndim   
x0s=[sobol.i4_sobol(2,i)[0] for i in range(nsol)]
firstpass=[optpos(x) for x in x0s]
#firstguess=x0s[numpy.argmin(firstpass)]
firstsol=[op.minimize(optpos,x0=x,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds[0:2]],jac=False,options={'eps':0.01,'disp':optmess,'gtol':0.1,'maxiter':3}) for x in x0s]
solpos=firstsol[numpy.argmin([sol.fun for sol in firstsol])]

#solpos=op.minimize(optpos,x0=firstguess,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds[0:2]],jac=False,options={'eps':0.005,'disp':optmess,'gtol':0.1})
solnpos=solpos.x
solnpos=numpy.append(solnpos,[dx_norm,0.5])

soln=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,solnpos)]
candy,candx=grid(soln,ny,nx)
plotAC(sumy,sumx,candy,candx,maximay,maximax)

optall=makeOptAll(arr,ny,nx,bounds)
xguess=list([0.5 for b in bounds])

sol=op.minimize(optall,x0=solnpos,method="L-BFGS-B",bounds=[(0.0,1.0) for b in bounds],jac=False,options={'eps':0.001,'disp':optmess,'gtol':0.1})
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
    print("Optimisation: "+sol.message)
    print(sol)
    
candy,candx=grid(soln,ny,nx)

# Output some plots
if showPlt:
    plotAC(sumy,sumx,candy,candx,maximay,maximax,pdf)

init=[b[0]+xv*(b[1]-b[0]) for b,xv in zip(bounds,xguess)]
