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

acmean=True

# Generate windowed mean intensities, scanning along x and y axes
# Estimate spot diameter, assuming grid takes up most of the plate
diam=min(float(arr.shape[0])/ny,float(arr.shape[1])/nx)
window=int(round(diam*windowFrac))
sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-window]):numpy.min([arr.shape[1],dx+window])]) for dx in range(0,arr.shape[1])],dtype=numpy.float)
sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-window]):numpy.min([arr.shape[0],dy+window]),0:arr.shape[1]]) for dy in range(0,arr.shape[0])],dtype=numpy.float)
# Smooth intensities to help eliminate small local maxima
sumx=ndimage.gaussian_filter1d(sumx,2.5)
sumy=ndimage.gaussian_filter1d(sumy,2.5)

# First peak in autocorrelation function is best estimate of distance between spots
maximay=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumy))))==-2)[0]
maximax=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0]
if acmean:
    dy=int(round(numpy.mean(numpy.diff(maximay))))
    dx=int(round(numpy.mean(numpy.diff(maximax))))
else:
    dy=int(round(maximay[0]))
    dx=int(round(maximax[0]))
    
ry=arr.shape[0]-((ny-1)*dy)
rx=arr.shape[1]-((nx-1)*dx)

checkvecs=[range(ry),range(rx)]
checkpos=list(itertools.product(*checkvecs))

def optfun(x):
    if(x[0]<0 or x[0]>rx or x[1]<0 or x[1]>ry or x[2]<=0):
        res=300
    else:
        try:
            res=-1*checkPos(arr,ny,nx,x[0:2],dx=x[2],dy=x[2],theta=x[3],sampfrac=0.35)
        except:
            res=300
    return(res)

sol=op.minimize(optfun,x0=(ry/2,rx/2,(dx+dy)/2,0),method="L-BFGS-B",options={'eps':[1.0,1.0,1.0,0.05]})

if len(sol.x)==4:
    theta=sol.x[3]
else:
    theta=0
grid=makeGrid(sol.x[0:2],ny,nx,dy=sol.x[2],dx=sol.x[2],theta=theta)

##lgrid=list(zip(*grid))
##candy=list(set(lgrid[0]))
##candx=list(set(lgrid[1]))
##candy.sort()
##candx.sort()

candy,candx=list(zip(*grid))

# Output some plots
if showPlt:
    fig,ax=plt.subplots(2,2,figsize=(15,15))
    
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



