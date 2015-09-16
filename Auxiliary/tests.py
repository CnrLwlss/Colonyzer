import os,numpy,itertools
from colonyzer2 import *
from scipy import optimize as op

LATESTIMAGE=os.path.realpath("Data\\endpoints\\384\\DLR00012647-2009-07-04_09-35-20.jpg")
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
# However, the row ACF estimate is still robust.  In the case of disagreement, choose the smallest value
rat=float(dy)/float(dx)
#if rat > (rattol+1.0) or rat < 1.0/(rattol+1.0):
dmin=min(dy,dx)
dy,dx=dmin,dmin
    
ry=arr.shape[0]-ny*dy
rx=arr.shape[1]-nx*dx

def fitProjection(proj,delt,n,sp=0.0,st=0.0):
    '''Find grid position that best fits a 1D projection of intensity by brute force'''
    checkinds=range(int(round(delt/2.0)),int(round(len(proj)-delt*n)))
    def getObj(i,proj,sp,st):
        peaks=proj[i:int(round((i+delt*n))):int(round(delt))]
        troughs=proj[int(round((i-delt/2.0))):int(round((i+delt*(n+0.5)))):int(round(delt))]
        return(numpy.median(peaks)-numpy.median(troughs)-sp*numpy.std(peaks)-st*numpy.std(troughs))
    grds=[getObj(i,proj,sp,st) for i in checkinds]
    maxind=numpy.argmax(grds)
    return((checkinds[maxind],grds[maxind]))

#xvals=[fitProjection(sumx,int(round(dxval)),nx) for dxval in numpy.linspace(0.95*dx,min(int(round(arr.shape[1]/nx))-1,1.05*dx),100)]
#yvals=[fitProjection(sumy,int(round(dyval)),ny) for dyval in numpy.linspace(0.95*dy,min(int(round(arr.shape[0]/ny))-1,1.05*dy),100)]
dd=int(round(0.01*dx))
dd=1
xtest=range(dx-dd,min(int(round(arr.shape[1]/nx))-1,dx+dd))
ytest=range(dy-dd,min(int(round(arr.shape[0]/ny))-1,dy+dd))
#xtest=numpy.linspace(dx-dd,min(int(round(arr.shape[1]/nx))-1,dx+dd),20)
#ytest=numpy.linspace(dy-dd,min(int(round(arr.shape[0]/ny))-1,dy+dd),20)
xvals=[fitProjection(sumx,int(round(dxval)),nx,1,1) for dxval in xtest]
yvals=[fitProjection(sumy,int(round(dyval)),ny,1,1) for dyval in ytest]
xind=numpy.argmin([x[1] for x in xvals])
yind=numpy.argmin([y[1] for y in yvals])

xbest=xvals[xind][0]
dx=xtest[xind]
ybest=yvals[yind][0]
dy=ytest[yind]

#xbest=fitProjection(sumx,dx,nx,1,1)[0]
#ybest=fitProjection(sumy,dy,ny,1,1)[0]

grd=makeGrid((ybest,xbest),ny,nx,dy=dy,dx=dx,theta=0)
candy,candx=list(zip(*grd))
plotAC(sumy,sumx,candy,candx,maximay,maximax,pdf=pdf,main="Projection Estimate")
