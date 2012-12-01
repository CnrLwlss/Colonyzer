import numpy,pandas,PIL,math,os,sys, time
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
    '''Search from x-RAD to x+RAD for pixel range dx-rad to dx+rad with the greatest mean intensity'''
    xmin,xmax=max(0,x-RAD),min(arr.shape[1],x+RAD)
    ymin,ymax=max(0,y-RAD),min(arr.shape[0],y+RAD)
    # Generate windowed mean intensities, scanning along x and y axes
    sumx=numpy.array([numpy.mean(arr[ymin:ymax,numpy.max([0,dx-rad]):numpy.min([arr.shape[1],dx+rad])]) for dx in xrange(xmin,xmax)],dtype=numpy.float)
    sumy=numpy.array([numpy.mean(arr[numpy.max([0,dy-rad]):numpy.min([arr.shape[0],dy+rad]),xmin:xmax]) for dy in xrange(ymin,ymax)],dtype=numpy.float)
    # Find all maxima
    maxx=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(sumx)))==-2)[0]
    maxy=1+numpy.where(numpy.diff(numpy.sign(numpy.diff(sumy)))==-2)[0]
    # Get maxima with highest peak
    if len(maxx)>0:
        bestx=maxx[0]
        for dx in maxx:
            if sumx[dx]>sumx[bestx]:
                best=dx
        bestx=xmin+bestx
    else:
        bestx=x
    if len(maxy)>0:
        besty=maxy[0]
        for dy in maxy:
            if sumy[dy]>sumy[besty]:
                best=dy
        besty=ymin+besty
    else:
        besty=y    
    return(bestx,besty)

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
    sigma1=numpy.sqrt(numpy.sum(P1*(numpy.array(intensities-mu1,dtype=numpy.float)**2)/numpy.sum(P1)))
    sigma2=numpy.sqrt(numpy.sum(P2*(numpy.array(intensities-mu2,dtype=numpy.float)**2)/numpy.sum(P2)))
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

def thresholdArr(arrim,thresh):
    '''Thresholding array representation of an image'''
    arrim[arrim<thresh]=0
    arrim[arrim>=thresh]=255
    arrim=numpy.array(arrim,dtype=numpy.uint8)
    imnew=Image.fromarray(arrim, "L")
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

def plotModel(bindat,thresholds=(),label="",pdf=None):
    '''Plot intensity density histogram, modelled distribution, component distributions and threshold estimate.'''
    plt.figure()
    plt.plot(bindat.intensities,bindat.freq,color="black")
    plt.plot(bindat.intensities,bindat.gauss1,color="red")
    plt.plot(bindat.intensities,bindat.gauss2,color="blue")
    plt.plot(bindat.intensities,bindat.mixed,color="green")
    plt.xlabel('Intensity')
    plt.ylabel('Density')
    plt.suptitle(label)
    for thresh in thresholds:
        plt.axvline(x=thresh,linestyle='--',linewidth=0.5,color="darkorchid")
    if pdf==None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()

def getEdges(arr,cutoff=0.9975):
    '''Sobel edge detection for 2d array using scipy functions'''
    sx = ndimage.sobel(arr, axis=0)
    sy = ndimage.sobel(arr, axis=1)
    sob = numpy.hypot(sx, sy)
    sob[sob<stats.mstats.mquantiles(sob,cutoff)[0]]=0
    sob[sob>0]=1
    return(numpy.array(sob,dtype=numpy.bool))  

def sizeSpots(locations,arr,thresharr,edge,background=0):
    '''Adds intensity measures columns to locations dataFrame'''
    intMax=255.0
    # http://en.wikipedia.org/wiki/Shape_factor_(image_analysis_and_microscopy)#Circularity
    # Calculate area, intensity and trimmed intensity for each spot
    sumInt,sumArea,trim,fMed,bMed,circ,fVar,perim=[],[],[],[],[],[],[],[]
    for i in xrange(0,len(locations.x.values)):
        x,y,rad=locations.x.values[i],locations.y.values[i],int(math.ceil(max(locations.Diameter.values)/2.0))
        tile=arr[y-rad:(y+rad+1),x-rad:(x+rad+1)]-background
        threshtile=thresharr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        edgetile=edge[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        perimeter=numpy.sum(edgetile)
        area=numpy.sum(threshtile)
        if perimeter>0:
            circularity=4*math.pi*area/(perimeter)**2
        else:
            circularity=0
        featureMedian=numpy.median(tile[threshtile])/intMax
        backgroundMedian=numpy.median(tile[numpy.logical_not(threshtile)])/intMax
        sumInt.append(float(numpy.sum(tile))/(float(tile.size)*intMax))
        sumArea.append(float(area)/float(tile.size))
        trim.append(float(numpy.sum(tile[threshtile]))/(float(tile.size)*intMax))
        fMed.append(featureMedian/intMax)
        bMed.append(backgroundMedian/intMax)
        circ.append(circularity)
        fVar.append(numpy.var(tile[threshtile]/intMax))
        perim.append(float(perimeter)/float(tile.size))
    locations["Intensity"]=sumInt
    locations["Area"]=sumArea
    locations["Trimmed"]=trim
    locations["FeatureMedian"]=fMed
    locations["FeatureVariance"]=fVar
    locations["BackgroundMedian"]=bMed
    locations["Circularity"]=circ
    locations["Perimeter"]=perim
    return(locations)

def getColours(im,locations,thresharr):
    '''Extract feature and background mean and median Red Green and Blue channel values for a given 24 bit image'''
    (red,green,blue)=im.split()
    redarr,greenarr,bluearr=numpy.array(red,dtype=numpy.uint8),numpy.array(green,dtype=numpy.uint8),numpy.array(blue,dtype=numpy.uint8)
    r,g,b,rB,gB,bB,rm,gm,bm,rmB,gmB,bmB=[],[],[],[],[],[],[],[],[],[],[],[]
    store=numpy.zeros((len(locations.x.values),12),numpy.float)
    for i in xrange(0,len(locations.x.values)):
        x,y,rad=locations.x.values[i],locations.y.values[i],int(math.ceil(max(locations.Diameter.values)/2.0))
        redtile=redarr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        greentile=greenarr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        bluetile=bluearr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        threshtile=thresharr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        rMean,gMean,bMean=numpy.mean(redtile[threshtile]),numpy.mean(greentile[threshtile]),numpy.mean(bluetile[threshtile])
        rMed,gMed,bMed=numpy.median(redtile[threshtile]),numpy.median(greentile[threshtile]),numpy.median(bluetile[threshtile])
        rMeanBk,gMeanBk,bMeanBk=numpy.mean(redtile[numpy.logical_not(threshtile)]),numpy.mean(greentile[numpy.logical_not(threshtile)]),numpy.mean(bluetile[numpy.logical_not(threshtile)])
        rMedBk,gMedBk,bMedBk=numpy.median(redtile[numpy.logical_not(threshtile)]),numpy.median(greentile[numpy.logical_not(threshtile)]),numpy.median(bluetile[numpy.logical_not(threshtile)])
        store[i]=[rMean,gMean,bMean,rMeanBk,gMeanBk,bMeanBk,rMed,gMed,bMed,rMedBk,gMedBk,bMedBk]
    locations["redMean"]=store[:,0]
    locations["greenMean"]=store[:,1]
    locations["blueMean"]=store[:,2]
    locations["redMeanBack"]=store[:,3]
    locations["greenMeanBack"]=store[:,4]
    locations["blueMeanBack"]=store[:,5]
    locations["redMedian"]=store[:,6]
    locations["greenMedian"]=store[:,7]
    locations["blueMedian"]=store[:,8]
    locations["redMedianBack"]=store[:,9]
    locations["greenMedianBack"]=store[:,10]
    locations["blueMedianBack"]=store[:,11]
    return(locations)      

def saveColonyzer(filename,locs,thresh,dx,dy):
    '''Save output data in original Colonyzer format'''
    # FILENAME ROW COLUMN TOPLEFTX TOPLEFTY WHITEAREA(px) TRIMMED THRESHOLD INTENSITY EDGEPIXELS COLR COLG COLB BKR BKG BKB EDGELEN XDIM YDIM
    df={}
    df["FILENAME"]=locs["Filename"].values
    df["ROW"]=locs["Row"].values
    df["COLUMN"]=locs["Column"].values
    df["TOPLEFTX"]=locs["x"].values-locations["Diameter"].values/2.0
    df["TOPLEFTY"]=locs["y"].values-locations["Diameter"].values/2.0
    df["WHITEAREA"]=locs["Area"].values
    df["TRIMMED"]=locs["Trimmed"].values
    df["THRESHOLD"]=thresh
    df["INTENSITY"]=locs["Intensity"].values
    df["EDGEPIXELS"]=locs["FeatureMedian"].values ### NOTE LABEL INCORRECT!
    df["COLR"]=locs["redMedian"].values
    df["COLG"]=locs["greenMedian"].values
    df["COLB"]=locs["blueMedian"].values
    df["BKR"]=locs["redMedianBack"].values
    df["BKG"]=locs["greenMedianBack"].values
    df["BKB"]=locs["blueMedianBack"].values
    df["EDGELEN"]=locs["Perimeter"].values
    df["XDIM"]=dx
    df["YDIM"]=dy
    colorder=("FILENAME","ROW","COLUMN","TOPLEFTX","TOPLEFTY","WHITEAREA","TRIMMED","THRESHOLD","INTENSITY","EDGEPIXELS","COLR","COLG","COLB","BKR","BKG","BKB","EDGELEN","XDIM","YDIM")
    dataf=pandas.DataFrame(df)
    dataf.to_csv(filename,"\t",index=False,header=False)
    return(dataf)

start=time.time()

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
#barcsDone=list(numpy.unique([dat[barcRange[0]:barcRange[1]] for dat in alldats]))
barcsDone=[]
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
arrN=numpy.array(img,dtype=numpy.float)
    
nx,ny=24,16
diam=int(round(min(float(arrN.shape[0])/ny,float(arrN.shape[1])/nx)))
(candx,candy,dx,dy)=estimateLocations(arrN,diam,showPlt=False)
xloc,yloc=numpy.meshgrid(candx,candy)
cols,rows=numpy.meshgrid(numpy.arange(1,nx+1),numpy.arange(1,ny+1))
d={"Row":rows.flatten(),"Column":cols.flatten(),"y":yloc.flatten(),"x":xloc.flatten()}
locations=pandas.DataFrame(d)
rad=int(round(float(min(dx,dy))/2.0))
RAD=int(round(1.2*rad))
RAD=rad/6
for i in xrange(0,len(locations.x)):
    (x,y)=optimiseSpot(arrN,locations.x[i],locations.y[i],rad,RAD)
    locations.x[i]=x
    locations.y[i]=y
locations["Diameter"]=min(dx,dy)

print("Cultures located")

# Analyse first image to allow lighting correction
im0=Image.open(EARLIESTIMAGE)
img=im0.convert("F")
arr0=numpy.array(img,dtype=numpy.float)
smoothed_arr=ndimage.gaussian_filter(arr0,arr0.shape[1]/500)
average_back=numpy.mean(smoothed_arr[numpy.min(locations.y):numpy.max(locations.y),numpy.min(locations.x):numpy.max(locations.x)])
correction_map=average_back/smoothed_arr

print("Lighting correction map constructed")

# Apply lighting correction to last image
corrected_arr=arrN*correction_map

# Threshold corrected image
# Initial guess for mixed model parameters for thresholding lighting corrected image
# Trim outer part of image to remove plate walls
trimmed_arr=corrected_arr[(min(locations.y)-dy):(max(locations.y)+dy),(min(locations.x)-dx):(max(locations.x)+dx)]
(counts,intensities)=numpy.histogram(trimmed_arr,bins=2**8,range=(0,2**8))
intensities=numpy.array(intensities[0:-1],dtype=numpy.int)
(bindat,[theta,mu1,mu2,sigma1,sigma2])=initialGuess(intensities,counts)

# Maximise likelihood of 2-component mixed Gaussian model parameters given binned observations by constrained optimisation
logL=makeObjective(bindat.intensities,bindat.counts,totFunc)
b=[(0.0,1.0),(float(mu1)/5.0,5*float(mu1)),(float(mu2)/5.0,5.0*float(mu2)),(float(sigma1)/5.0,5.0*float(sigma1)),(float(sigma2)/5.0,5.0*float(sigma2))]
opt=optimize.fmin_l_bfgs_b(logL,[theta,mu1,mu2,sigma1,sigma2],bounds=b,approx_grad=True)
[theta_opt,mu1_opt,mu2_opt,sigma1_opt,sigma2_opt]=opt[0]

# Best estimate for threshold is point of intersection of two fitted component Gaussians
thresh=getRoot(opt[0],intensities)
thresh0=int(round(thresh))
thresh1=thresh0

# Make threshold as low as possible, for maximum sensitivity
smoothcounts=ndimage.gaussian_filter1d(counts,1)
while smoothcounts[thresh1]>=smoothcounts[thresh1-1]:
    thresh1-=1

# Save final mask for cutting out all cell signal from earlier images
finalMask=numpy.ones(arrN.shape,dtype=numpy.bool)
finalMask[arrN<thresh1]=False

# Modelled densities
bindat["mixed"]=numpy.array([totFunc(x,opt[0]) for x in bindat.intensities],dtype=numpy.float)
bindat["gauss1"]=numpy.array([theta_opt*stats.norm.pdf(x,mu1_opt,sigma1_opt) for x in bindat.intensities],dtype=numpy.float)
bindat["gauss2"]=numpy.array([(1.0-theta_opt)*stats.norm.pdf(x,mu2_opt,sigma2_opt) for x in bindat.intensities],dtype=numpy.float)
#plotGuess(bindat)
#plotModel(bindat,(thresh0,thresh1))

# Find culture area as a function of threshold value
#cellarea=[arr[arr>i].size for i in xrange(0,256)]

print("Threshold located")
barcdict[BARCODE].sort()
for FILENAME in barcdict[BARCODE]:
    print FILENAME
    im=Image.open(FILENAME)
    img=im.convert("F")
    arr=numpy.array(img,dtype=numpy.float)
    # Correct spatial gradient
    arr=arr*correction_map
    arrsm=arr[numpy.min(locations.y):numpy.max(locations.y),numpy.min(locations.x):numpy.max(locations.x)]
    masksm=finalMask[numpy.min(locations.y):numpy.max(locations.y),numpy.min(locations.x):numpy.max(locations.x)]
    meanPx=numpy.mean(arrsm[numpy.logical_not(masksm)])
    # Correct lighting differences
    arr=arr+(average_back-meanPx)
    mask=numpy.ones(arr.shape,dtype=numpy.bool)
    mask[arr<thresh1]=False
    imthresh=thresholdArr(numpy.copy(arr),thresh1)
    edge=getEdges(arr,0.925)
    locations=sizeSpots(locations,arr,mask,edge,average_back)
    locations=getColours(im,locations,mask)
    locations["Barcode"]=BARCODE
    locations["Filename"]=FILENAME[0:-4]
    locations.to_csv(FILENAME[0:-4]+".out","\t",index=False)
    dataf=saveColonyzer(FILENAME[0:-4]+".dat",locations,thresh1,dx,dy)
    draw=ImageDraw.Draw(im)
    for i in xrange(0,len(locations.x)):
        x,y,r=int(round(locations.x[i])),int(round(locations.y[i])),int(round(float(locations.Diameter[i])/2.0))
        draw.rectangle((x-r,y-r,x+r,y+r),outline=(255,255,0))
    im.save(FILENAME[0:-4]+".png")
    #im.show()
    #imthresh.show()

print("Finished: "+str(time.time()-start)+" s")

		
