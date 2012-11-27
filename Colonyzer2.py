import numpy,pandas,PIL,math,os,sys
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
from scipy import stats, optimize, ndimage, signal
import scipy

def contiguous_regions(condition):
    '''Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array'''
    # Find the indicies of changes in "condition"
    d = numpy.diff(condition)
    idx, = d.nonzero() 
    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = numpy.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = numpy.r_[idx, condition.size] # Edit
    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

def autocor(x):
    s = numpy.fft.fft(x)
    res=numpy.real(numpy.fft.ifft(s*numpy.conjugate(s)))/numpy.var(x)
    res=res[0:len(res)/2]
    return(res)

def getMaxima(intensity):
    '''Numerical method to find local maxima in a 1D list with plateaus'''
    npoints=len(intensity)
    diffs=numpy.diff(intensity)
    zeroregions=contiguous_regions(diffs==0)
    maxima=[]
    for z in zeroregions:
        if z[0]>0 and z[1]<npoints-2 and diffs[z[0]-1]>0 and diffs[z[1]]<0:
            maxima.append(numpy.mean(z)+1)
    
    return(maxima)

def plotSpectrum(y,Fs):
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    n = len(y) # length of the signal
    k = numpy.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range

    Y = scipy.fft(y)/n # fft computing and normalization
    Y = Y[range(n/2)]

    plt.plot(frq[1:],abs(Y)[1:]) # plotting the spectrum
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')
    plt.show()

def estimateOffsets(arr,diam=20,limFrac=1.075,showPlt=True,pdfPlt=False):
    '''Sum intensities along edge of array and search for the first maximum which is greater than limFrac*smoothed sum as offset estimate'''
    sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-diam/4]):numpy.min([arr.shape[1],dx+diam/4])]) for dx in xrange(0,arr.shape[1])],dtype=numpy.int)
    sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-diam/4]):numpy.min([arr.shape[0],dy+diam/4]),0:arr.shape[1]]) for dy in xrange(0,arr.shape[0])],dtype=numpy.int)
    smoothedx=signal.medfilt(sumx,2*diam+1)
    smoothedy=signal.medfilt(sumy,2*diam+1)
    limx,limy=limFrac*smoothedx,limFrac*smoothedy
    maxx,maxy=getMaxima(sumx),getMaxima(sumy)
##    candx=maxx[sumx[maxx]>limx[maxx]]
##    candy=maxy[sumy[maxy]>limy[maxy]]
    candx,candy=maxx,maxy
    candx=signal.find_peaks_cwt(sumx,numpy.arange(0,arr.shape[1]/24))
    candy=signal.find_peaks_cwt(sumy,numpy.arange(0,arr.shape[1]/16))
    if showPlt:
        plt.plot(sumx)
        plt.plot(smoothedx)
        plt.plot(limx)
        for cand in candx:
            plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
        plt.xlabel('x coordinate')
        plt.ylabel('Mean Intensity')
        if pdfPlt:
            pdf.savefig()
            plt.close()
        else:
            plt.show()
        plt.plot(sumy)
        plt.plot(smoothedy)
        plt.plot(limy)
        for cand in candy:
            plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
        plt.xlabel('y coordinate')
        plt.ylabel('Mean Intensity')
        if pdfPlt:
            pdf.savefig()
            plt.close()
        else:
            plt.show()
    resx,resy=0,0
    if len(candx)==0:
        resx=-1
    else:
        resx=candx[0]
    if len(candy)==0:
        resy=-1
    else:
        resy=candy[0]
    return((resx,resy))

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def SetUp(instructarr):
    '''Set up plate description (i.e. 1536 or 384 format, top-left and bottom-right coordinates etc.)'''
    # Check if the first element is an integer:
    if is_number(instructarr[0]):
        NoSpots=int(instructarr[0])
    else:
        NoSpots=instructarr[0]
    TLX,TLY,BRX,BRY=instructarr[1],instructarr[2],instructarr[3],instructarr[4]
    if NoSpots==384:
        nocols,norows = 24,16
    elif NoSpots==1536:
        nocols,norows = 48,32
    elif NoSpots==768:
        nocols,norows = 48,32
    elif NoSpots==96:
        nocols,norows = 12,8
    elif NoSpots==48:
        nocols,norows = 6,8
    elif not is_number(NoSpots) and 'x' in NoSpots and len(NoSpots.split("x"))==2:
        tmp=NoSpots.split("x")
        nocols,norows = int(tmp[0]),int(tmp[1])
    else:
        nocols,norows=0,0
        print "WARNING: Incorrect spot number specified!"
    tlx,tly=TLX,TLY
    brx,bry=BRX,BRY

    # Best estimates for tile dimensions
    xdimf=float(abs(brx-tlx))/float(nocols-1)
    ydimf=float(abs(bry-tly))/float(norows-1)
    xdim=int(round(xdimf))
    ydim=int(round(ydimf))

    # Parameter for specifying the search area (area is (NSearch*2+1)^2)
    NSearch = int(round(float(min(xdim,ydim))/2.0))

    # Best estimates for the starting coordinates
    xstart=max(0,int(round(float(tlx)-0.5*xdimf)))
    ystart=max(0,int(round(float(tly)-0.5*ydimf)))

    print "Instructions: ",xstart,ystart,xdim,ydim,NSearch
    return((nocols,norows, tlx,tly,brx,bry,xdim,ydim,xdimf,ydimf,xstart,ystart,NSearch))

# Current directory
syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

# Read in the Colonyzer input file
Instructions=open(os.path.join(fullpath,'Colonyzer.txt'),'r')
InsData={}
InsTemp=Instructions.readlines()

for x in xrange(0,len(InsTemp)):
    if InsTemp[x][0]!="#" and InsTemp[x][0]!="\n":
        tlist=InsTemp[x].split(',')
        if len(tlist)>1:
            InsData[tlist[0]]=[tlist[1],int(tlist[2]),int(tlist[3]),int(tlist[4]),int(tlist[5])]

if 'default' in InsData:
    (nocols,norows, tlx,tly,brx,bry,xdim,ydim,xdimf,ydimf,xstart,ystart,NSearch)=SetUp(InsData['default'])
else:
    print "ERROR: No default instructions"
    sys.exit()

# Create directories for storing output data and preview images
try:
    os.mkdir(os.path.join(fullpath,"Output_Images"))
    os.mkdir(os.path.join(fullpath,"Output_Data"))
except:
    print("Output directories already present")

outputimages=os.path.join(fullpath,"Output_Images")
outputdata=os.path.join(fullpath,"Output_Data")

allfiles=os.listdir(fullpath)
alldats=os.listdir(os.path.join(fullpath,"Output_Data"))
barcRange=(0,11) # Read this in from Colonyzer.txt?
barcsDone=list(numpy.unique([dat[barcRange[0]:barcRange[1]] for dat in alldats]))
barcdict={}
for filename in allfiles:
    barc=filename[barcRange[0]:barcRange[1]]
    if filename[-4:] in ('.jpg','.JPG') and barc not in barcsDone:
        if barc not in barcdict:
            barcdict[barc]=[filename]
        else:
            barcdict[barc].append(filename)
for b in barcdict:
    barcdict[b].sort(reverse=True)

BARCODE=barcdict.keys()[0]
LATESTIMAGE=barcdict[BARCODE][0]
im=Image.open(LATESTIMAGE)
img=im.convert("F")
arr=numpy.array(img,dtype=numpy.float)
#(dx,dy)=estimateOffsets(arr,2*NSearch,showPlt=True)
for FILENAME in barcdict[BARCODE]:
    print FILENAME
    

diam=2*NSearch
limFrac=1.075
showPlt=True
pdfPlt=False
sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-diam/4]):numpy.min([arr.shape[1],dx+diam/4])]) for dx in xrange(0,arr.shape[1])],dtype=numpy.int)
sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-diam/4]):numpy.min([arr.shape[0],dy+diam/4]),0:arr.shape[1]]) for dy in xrange(0,arr.shape[0])],dtype=numpy.int)
# First peak in autocorrelation function is best estimate of distance between spots
dx=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0][0]
dy=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0][0]
# Minimum filter to find best estimate of background
smoothedx=ndimage.filters.minimum_filter(sumx,2*dx+1)
smoothedy=ndimage.filters.minimum_filter(sumy,2*dy+1)
background=numpy.median(numpy.append(smoothedx,smoothedy))

limx,limy=limFrac*smoothedx,limFrac*smoothedy
maxx,maxy=getMaxima(sumx),getMaxima(sumy)
candx,candy=maxx,maxy
if showPlt:
    plt.plot(sumx)
    plt.plot(smoothedx)
    plt.plot(limx)
    for cand in candx:
        plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    plt.xlabel('x coordinate')
    plt.ylabel('Mean Intensity')
    if pdfPlt:
        pdf.savefig()
        plt.close()
    else:
        plt.show()
    plt.plot(autocor(sumx))
    maxima=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0]
    for cand in maxima:
        plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    plt.show()
    plt.plot(sumy)
    plt.plot(smoothedy)
    plt.plot(limy)
    for cand in candy:
        plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
    plt.xlabel('y coordinate')
    plt.ylabel('Mean Intensity')
    if pdfPlt:
        pdf.savefig()
        plt.close()
    else:
        plt.show()
resx,resy=0,0
if len(candx)==0:
    resx=-1
else:
    resx=candx[0]
if len(candy)==0:
    resy=-1
else:
    resy=candy[0]

		
