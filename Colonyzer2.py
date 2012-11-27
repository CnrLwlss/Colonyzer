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

def optimiseSpot(arr,x,y,rad,RAD):
    '''Search in a square of width 2*RAD for the centre of a smaller square (of width 2*rad) which maximises the sum of signal inside the small square'''
    x_grid, y_grid = numpy.meshgrid(numpy.arange(x-(RAD-rad),x+(RAD-rad)),numpy.arange(y-(RAD-rad),y+(RAD-rad)))
    x_grid=x_grid.flatten()
    y_grid=y_grid.flatten()
    signal=numpy.array([numpy.sum(arr[(y_grid[i]-rad):(y_grid[i]+rad),(x_grid[i]-rad):(x_grid[i]+rad)]) for i in xrange(0,len(x_grid))])
    iopt=numpy.argmax(signal)
    return(x_grid[iopt],y_grid[iopt])

def autocor(x):
    s = numpy.fft.fft(x)
    res=numpy.real(numpy.fft.ifft(s*numpy.conjugate(s)))/numpy.var(x)
    res=res[0:len(res)/2]
    return(res)

def showIm(arr,returnIm=False):
    '''Quick 8-bit preview images from float arrays, useful for debugging'''
    imarr=numpy.array(arr,dtype=numpy.uint8)
    imnew=Image.fromarray(imarr,"L")
    if returnIm:
        return(imnew)
    else:
        imnew.show()

def estimateLocations(arr,diam=20,showPlt=True,pdfPlt=False):
    '''Automatically search for best estimate for location of culture array'''
    # Generate windowed mean intensities, scanning along x and y axes
    sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],numpy.max([0,dx-diam/4]):numpy.min([arr.shape[1],dx+diam/4])]) for dx in xrange(0,arr.shape[1])],dtype=numpy.float)
    sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-diam/4]):numpy.min([arr.shape[0],dy+diam/4]),0:arr.shape[1]]) for dy in xrange(0,arr.shape[0])],dtype=numpy.float)
    # First peak in autocorrelation function is best estimate of distance between spots
    dx=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0][0]
    dy=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumx))))==-2)[0][0]
    # Find all maxima
    maxx=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(sumx)))==-2)[0]
    maxy=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(sumy)))==-2)[0]
    # Find the nspots maxima whose mean intermaximum distance is most internally consistent
    varx,vary=[],[]
    for i in xrange(0,len(maxx)-nx-1):
        varx.append(numpy.var(numpy.diff(maxx[i:(i+nx)])))
    for i in xrange(0,len(maxy)-ny-1):
        vary.append(numpy.var(numpy.diff(maxy[i:(i+ny)])))
    candx=maxx[numpy.argmin(varx):(numpy.argmin(varx)+nx)]
    candy=maxy[numpy.argmin(vary):(numpy.argmin(vary)+ny)]
    # Output some plots
    if showPlt:
        plt.plot(sumx)
        for cand in candx:
            plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
        plt.xlabel('x coordinate (px)')
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
        plt.xlabel('Offset dx (px)')
        plt.ylabel('Autocorrelation')
        if pdfPlt:
            pdf.savefig()
            plt.close()
        else:
            plt.show()
            
        plt.plot(sumy)
        for cand in candy:
            plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
        plt.xlabel('y coordinate (px)')
        plt.ylabel('Mean Intensity')
        if pdfPlt:
            pdf.savefig()
            plt.close()
        else:
            plt.show()
        plt.plot(autocor(sumy))
        maxima=numpy.where(numpy.diff(numpy.sign(numpy.diff(autocor(sumy))))==-2)[0]
        for cand in maxima:
            plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
        plt.xlabel('Offset dy (px)')
        plt.ylabel('Autocorrelation')
        if pdfPlt:
            pdf.savefig()
            plt.close()
        else:
            plt.show()
    return((candx,candy,dx,dy))

def initialGuess(intensities,counts):
    '''Construct non-parametric guesses for distributions of two components and use these to estimate Gaussian parameters'''
    # Get all maxima
    maxima=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(counts)))==-2)[0]
    maxima=maxima[counts[maxima]>0.01*numpy.sum(counts)]
    # Use first maximum of distribution as estimate of mean of first component
    mu1=intensities[maxima[0]]
    # Mirror curve from 0...mu1 to estimate distribution of first component
    P1=numpy.concatenate((counts[0:mu1+1],counts[mu1-1::-1],numpy.zeros(len(counts)-2*mu1-1,dtype=numpy.int)))
    # Use last maximum of distribution as estimate of mean of second component
    mu2=intensities[maxima[-1]]
    # Mirror curve for second peak also
    halfpeak=counts[mu2:]
    peak=numpy.concatenate((halfpeak[-1::-1],halfpeak[1:]))
    P2=numpy.concatenate((numpy.zeros(len(counts)-len(peak),dtype=numpy.int),peak))
    bindat=pandas.DataFrame(intensities,columns=["intensities"])
    bindat["counts"]=counts
    bindat["P1"]=P1
    bindat["P2"]=P2
    # Calculate standard deviation of (binned) observations from first and second components
    sigma1=numpy.sqrt(numpy.sum(P1*(numpy.array(intensities-mu1,dtype=numpy.float)**2)/numpy.sum(P1)))/2
    sigma2=numpy.sqrt(numpy.sum(P2*(numpy.array(intensities-mu2,dtype=numpy.float)**2)/numpy.sum(P2)))/2
    # Estimate component weighting
    theta=float(numpy.sum(P1))/float(numpy.sum(P1)+numpy.sum(P2))
    # Discard empty bins
    bindat=bindat[bindat.counts>0]
    bindat["frac"]=numpy.array(numpy.cumsum(bindat.counts),dtype=numpy.float)/numpy.sum(bindat.counts)
    bindat["freq"]=numpy.array(bindat.counts,dtype=numpy.float)/numpy.sum(bindat.counts)
    return((bindat,[theta,mu1,mu2,sigma1,sigma2]))

def totFunc(x,p):
    '''Probability density function for a 2-component mixed Gaussian model'''
    [theta,mu1,mu2,sigma1,sigma2]=p
    if mu2-mu1<2:
        candidate=1e-100
    else:
        candidate=theta*stats.norm.pdf(x,mu1,sigma1)+(1.0-theta)*stats.norm.pdf(x,mu2,sigma2)
##    if candidate<1e-100:
##        candidate=1e-100
    return(candidate)

def makeObjective(ints,cnts,PDF):
    '''Returns a function for (log likelihood)*-1 (suitable for minimisation), given a set of binned observations and a PDF'''
    ints=numpy.array(ints,dtype=numpy.int)
    cnts=numpy.array(cnts,dtype=numpy.int)
    def logL(p):
        modeldens=numpy.array([PDF(x,p) for x in ints],dtype=numpy.float)
        lik=numpy.sum(cnts*numpy.log(modeldens))
        return(-1*lik)
    return(logL)

def getRoot(p,ints):
    '''Get the point at which two component Gaussians intersect.  Specifically looking for root with highest probability.'''
    [theta,mu1,mu2,sigma1,sigma2]=p
    ints=numpy.array(ints,dtype=numpy.int)
    def diffFunc(x):
        return(theta*stats.norm.pdf(x,mu1,sigma1)-(1.0-theta)*stats.norm.pdf(x,mu2,sigma2))
    # Find pairs of points in truncated, filtered intensity list which bracket any roots
    diffs=[numpy.sign(diffFunc(x)) for x in ints]
    switches=[]
    for i in xrange(1,len(diffs)):
        if abs((diffs[i]-diffs[i-1]))==2:
            switches.append((i,i-1))
    # Fine-tune root locations
    threshlist=[]
    for switch in switches:
        thresh=optimize.brentq(diffFunc,ints[switch[0]],ints[switch[1]])
        threshlist.append(thresh)
    # Get root which gives the highest probability for peak 1 (or peak 2, either is fine)
    p1=[stats.norm.pdf(thresh,mu1,sigma1) for thresh in threshlist]
    thresh=threshlist[numpy.argmax(p1)]
    return(thresh)

def thresholdImage(arr,thresh,saveim=True):
    '''Thresholding array representation of an image'''
    arr[arr<thresh]=0
    arr[arr>=thresh]=255
    arr=numpy.array(arr,dtype=numpy.uint8)
    arr=arr.reshape(im.size[::-1])
    imnew=Image.fromarray(arr, "L")
    if saveim:
        imnew.save(os.path.join(fullpath,fname+".png"))
    return(imnew)

def plotGuess(bindat,label="",pdf=None):
    '''Plot intensity frequency histogram and non-parametric estimates of component distributions'''
    plt.figure()
    plt.plot(bindat.intensities,bindat.counts,color="black")
    plt.plot(bindat.intensities,bindat.P1,color="red")
    plt.plot(bindat.intensities,bindat.P2,color="blue")
    plt.xlabel('Intensity')
    plt.ylabel('Frequency')
    plt.suptitle(label)
    if pdf==None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()

def plotModel(bindat,label="",pdf=None):
    '''Plot intensity density histogram, modelled distribution, component distributions and threshold estimate.'''
    plt.figure()
    plt.plot(bindat.intensities,bindat.freq,color="black")
    plt.plot(bindat.intensities,bindat.gauss1,color="red")
    plt.plot(bindat.intensities,bindat.gauss2,color="blue")
    plt.plot(bindat.intensities,bindat.mixed,color="green")
    plt.xlabel('Intensity')
    plt.ylabel('Density')
    plt.suptitle(label)
    plt.axvline(x=thresh,linestyle='--',linewidth=0.5,color="darkorchid")
    if pdf==None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()

# Current directory
syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

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
EARLIESTIMAGE=barcdict[BARCODE][-1]
im=Image.open(LATESTIMAGE)
img=im.convert("F")
arr=numpy.array(img,dtype=numpy.float)
    
nx,ny=24,16
diam=int(round(min(float(arr.shape[0])/ny,float(arr.shape[1])/nx)))
(candx,candy,dx,dy)=estimateLocations(arr,diam,showPlt=False)
xloc,yloc=numpy.meshgrid(candx,candy)
cols,rows=numpy.meshgrid(numpy.arange(1,nx+1),numpy.arange(1,ny+1))
d={"Row":rows.flatten(),"Column":cols.flatten(),"y":yloc.flatten(),"x":xloc.flatten()}
locations=pandas.DataFrame(d)
rad=int(round(float(min(dx,dy))/2.0))
RAD=1.2*rad
for i in xrange(0,len(locations.x)):
    (dx,dy)=optimiseSpot(arr,locations.x[i],locations.y[i],rad,RAD)
    locations.x[i]=dx
    locations.y[i]=dy

# Analyse first image to allow lighting correction
im0=Image.open(EARLIESTIMAGE)
img=im0.convert("F")
arr0=numpy.array(img,dtype=numpy.float)
smoothed_arr=ndimage.gaussian_filter(arr0,arr.shape[1]/500)
average_back=numpy.mean(smoothed_arr)

# Apply lighting correction to last image
corrected_arr=arr*average_back/smoothed_arr

# Threshold corrected image
# Initial guess for mixed model parameters for thresholding lighting corrected image
# Trim outer part of image to remove plate walls
trimmed_arr=corrected_arr[(min(locations.y)-dy):(max(locations.y)+dy),(min(locations.x)-dx):(max(locations.x)+dx)]
(counts,intensities)=numpy.histogram(trimmed_arr,bins=2**8,range=(0,2**8))
intensities=numpy.array(intensities[0:-1],dtype=numpy.int)
(bindat,[theta,mu1,mu2,sigma1,sigma2])=initialGuess(intensities,counts)

# Maximise likelihood of 2-component mixed Gaussian model parameters given binned observations by constrained optimisation
logL=makeObjective(bindat.intensities,bindat.counts,totFunc)
b=[(0.5,1.0),(float(mu1)/5.0,5*float(mu1)),(float(mu2)/5.0,5.0*float(mu2)),(float(sigma1)/5.0,5.0*float(sigma1)),(float(sigma2)/5.0,5.0*float(sigma2))]
opt=optimize.fmin_l_bfgs_b(logL,[theta,mu1,mu2,sigma1,sigma2],bounds=b,approx_grad=True)
[theta_opt,mu1_opt,mu2_opt,sigma1_opt,sigma2_opt]=opt[0]

# Best estimate for threshold is point of intersection of two fitted component Gaussians
thresh=getRoot(opt[0],intensities)

# Modelled densities
bindat["mixed"]=numpy.array([totFunc(x,opt[0]) for x in bindat.intensities],dtype=numpy.float)
bindat["gauss1"]=numpy.array([theta_opt*stats.norm.pdf(x,mu1_opt,sigma1_opt) for x in bindat.intensities],dtype=numpy.float)
bindat["gauss2"]=numpy.array([(1.0-theta_opt)*stats.norm.pdf(x,mu2_opt,sigma2_opt) for x in bindat.intensities],dtype=numpy.float)

# Threshold image and clean up noise
imthresh=thresholdImage(numpy.copy(corrected_arr),thresh,saveim=False)
    
imthresh.show()
plotGuess(bindat)
plotModel(bindat)

for FILENAME in barcdict[BARCODE]:
    print FILENAME



		
