import os,numpy,itertools
from colonyzer2 import *
from scipy import optimize as op

LATESTIMAGE=os.path.realpath("F:\Colonyzer\Auxiliary\Data\DLR00012647-2009-07-04_09-35-20.jpg")
im,arr=openImage(LATESTIMAGE)

nx,ny=24,16
windowFrac=0.25
smoothWindow=0.13
showPlt=True
pdf=None
diam=20

# Generate windowed mean intensities, scanning along x and y axes
window=int(round(diam*windowFrac))
sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-window]):numpy.min([arr.shape[1],dx+window])]) for dx in xrange(0,arr.shape[1])],dtype=numpy.float)
sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-window]):numpy.min([arr.shape[0],dy+window]),0:arr.shape[1]]) for dy in xrange(0,arr.shape[0])],dtype=numpy.float)
# Smooth intensities to help eliminate small local maxima
sumx=ndimage.gaussian_filter1d(sumx,diam*smoothWindow)
sumy=ndimage.gaussian_filter1d(sumy,diam*smoothWindow)
# First peak in autocorrelation function is best estimate of distance between spots
dx=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0][0]
dy=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumy))))==-2)[0][0]

def sampleArr(arr,pos,svals):
    '''Sum pixel intensities from arr in a rectangle centred on y,x'''
    y,x=pos
    sx,sy=svals
    minx,miny=max(0,x-sx),max(0,y-sy)
    maxx,maxy=min(arr.shape[1]-1,x+sx), min(arr.shape[0]-1,y+sy)
    val=scipy.stats.nanmedian(arr[miny:maxy,minx:maxx],axis=None)
    return(val)

def makeGrid(pos0,ny,nx,dy,dx):
    '''Generate grid coordinates from top left position, grid dimension and gap sizes.'''
    y0,x0=pos0
    vecs=[range(ny),range(nx)]
    gpos=list(itertools.product(*vecs))
    pos=[(int(round(y0+gp[0]*dy)),int(round(x0+gp[1]*dx))) for gp in gpos]
    return(pos)
    
def checkPos(arr,ny,nx,pos0,dy,dx,theta=0,sampfrac=0.1):
    '''Return sum of pixel intensities in arr around grid points'''
    sx=int(round(dx*sampfrac))
    sy=int(round(dy*sampfrac))
    pos=makeGrid(pos0,ny,nx,dy,dx)
    vals=[sampleArr(arr,p,(sx,sy)) for p in pos]
    print(vals)
    return(sum(vals))

# Find all maxima
maxx=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(sumx)))==-2)[0]
maxy=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(sumy)))==-2)[0]
# Find the nspots maxima whose mean intermaximum distance is most internally consistent
varx,vary=[],[]
for i in xrange(0,len(maxx)-nx+1):
    varpos=numpy.var(numpy.diff(maxx[i:(i+nx)]))
    # Small penalty for deviations from centre of image
    symmpen=10*abs(maxx[i]-(arr.shape[1]-maxx[i+nx-1]))/dx
    varx.append(varpos+symmpen)
for i in xrange(0,len(maxy)-ny+1):
    # Small penalty for deviations from centre of image
    varpos=numpy.var(numpy.diff(maxy[i:(i+ny)]))
    symmpen=10*abs(maxy[i]-(arr.shape[0]-maxy[i+ny-1]))/dy
    vary.append(varpos+symmpen)

candx=maxx[numpy.argmin(varx):(numpy.argmin(varx)+nx)]
candy=maxy[numpy.argmin(vary):(numpy.argmin(vary)+ny)]

maximay=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumy))))==-2)[0]
maximax=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0]

# Update spot sizes based on grid location estimates
#dx=int(numpy.round(numpy.mean(numpy.diff(candx))))
#dy=int(numpy.round(numpy.mean(numpy.diff(candy))))

dy=int(round(numpy.mean(numpy.diff(maximay))))
dx=int(round(numpy.mean(numpy.diff(maximax))))

ry=arr.shape[0]-((ny-1)*dy)
rx=arr.shape[1]-((nx-1)*dx)

checkvecs=[range(ry),range(rx)]
checkpos=list(itertools.product(*checkvecs))

#checkvals=[checkPos(arr,nx,ny,p,dx,dy) for p in checkpos[0:1000]]

def optfun(x):
    if(x[0]<0 or x[0]>rx or x[1]<0 or x[1]>ry):
        res=300
    else:
        res=-1*checkPos(arr,ny,nx,x,dy,dx,sampfrac=0.55)
    return(res)

def minfun(x):
    return(op.minimize(optfun,x0=x,method="L-BFGS-B",options={'eps':1.0}))

#sol=op.basinhopping(optfun,niter=10,x0=(ry/2,rx/2),T=10000,minimizer_kwargs={'method':'L-BFGS-B','options':{'eps':1.0}})
#sol=op.minimize(optfun,x0=(ry/2,rx/2),method="L-BFGS-B",options={'eps':5.0})

ds=3

yguess,xguess=numpy.linspace(0,ry,ds),numpy.linspace(0,rx,ds)
guesses=itertools.product(yguess,xguess)
gsol=[minfun(guess) for guess in guesses]
fgsol=[sol.fun for sol in gsol]

indbest=numpy.argmin(fgsol)

print(optfun((100,100)))
print(optfun((100,1000)))
print(optfun((ry/2,rx/2)))
print(optfun(gsol[indbest].x))

grid=makeGrid(gsol[indbest].x,ny,nx,dy,dx)
lgrid=zip(*grid)
candy=list(set(lgrid[0]))
candx=list(set(lgrid[1]))
candy.sort()
candx.sort()

# Output some plots
if showPlt:
    fig,ax=plt.subplots(2,2,figsize=(20,10))
    
    ax[0,0].plot(sumx)
    for cand in candx:
        ax[0,0].axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    ax[0,0].set_xlabel('x coordinate (px)')
    ax[0,0].set_ylabel('Mean Intensity')

    ax[0,1].plot(autocor(sumx))
    for cand in maximax:
        ax[0,1].axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    ax[0,1].set_xlabel('Offset dx (px)')
    ax[0,1].set_ylabel('Autocorrelation')
        
    ax[1,0].plot(sumy)
    for cand in candy:
        ax[1,0].axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    ax[1,0].set_xlabel('y coordinate (px)')
    ax[1,0].set_ylabel('Mean Intensity')

    ax[1,1].plot(autocor(sumy))
    for cand in maximay:
        ax[1,1].axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    ax[1,1].set_xlabel('Offset dy (px)')
    ax[1,1].set_ylabel('Autocorrelation')
    plt.show()

