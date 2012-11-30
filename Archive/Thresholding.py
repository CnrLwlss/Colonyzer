import numpy,pandas,PIL,math,os,sys
from PIL import Image, ImageDraw, ImageFont
from scipy import stats, optimize, ndimage, signal
from GridLocations import *

import matplotlib.pyplot as plt
# http://matplotlib.sourceforge.net/examples/pylab_examples/multipage_pdf.html
from matplotlib.backends.backend_pdf import PdfPages

def totFunc(x,p):
    '''Probability density function for a 2-component mixed Gaussian model'''
    [theta,mu1,mu2,sigma1,sigma2]=p
    if mu2-mu1<70:
        candidate=1e-100
    else:
        candidate=theta*stats.norm.pdf(x,mu1,sigma1)+(1.0-theta)*stats.norm.pdf(x,mu2,sigma2)
    if candidate<1e-100:
        candidate=1e-100
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

def initialGuess(intensities,counts):
    '''Construct non-parametric guesses for distributions of two components and use these to estimate Gaussian parameters'''
    # Use mode of distribution as estimate of mean of first component
    mu1=intensities[numpy.argmax(counts)]
    # Mirror curve from 0...mu1 to estimate distribution of first component
    P1=numpy.concatenate((counts[0:mu1+1],counts[mu1-1::-1],numpy.zeros(len(counts)-2*mu1-1,dtype=numpy.int)))
    # Subtract first component from full distribution to estimate distribution of second component
    P2=numpy.maximum(0,counts-P1)
    # Calculate mean of observations from second component
    mu2=numpy.sum(numpy.array(P2,dtype=numpy.float)*numpy.array(intensities,dtype=numpy.float))/numpy.sum(P2)
    # Second component in particular has a heavy right tail, discard tail
    P2[2*int(mu2):]=0
    # Calculate standard deviation of (binned) observations from first and second components
    sigma1=numpy.sqrt(numpy.sum(P1*(numpy.array(intensities-mu1,dtype=numpy.float)**2)/numpy.sum(P1)))
    sigma2=numpy.sqrt(numpy.sum(P2*(numpy.array(intensities-mu2,dtype=numpy.float)**2)/numpy.sum(P2)))
    # Estimate component weighting
    theta=float(sum(P1))/float(sum(P1)+sum(P2))

    bindat=pandas.DataFrame(intensities,columns=["intensities"])
    bindat["counts"]=counts
    bindat["P1"]=P1
    bindat["P2"]=P2
    
    # Discard empty bins
    bindat=bindat[bindat.counts>0]
    # Discard observations beyond the 99th quantile
    bindat["frac"]=numpy.array(numpy.cumsum(bindat.counts),dtype=numpy.float)/numpy.sum(bindat.counts)
    bindat=bindat[bindat.frac<=0.99]
    bindat["freq"]=numpy.array(bindat.counts,dtype=numpy.float)/numpy.sum(bindat.counts)
    return((bindat,[theta,mu1,mu2,sigma1,sigma2]))

def plotGuess(bindat,pdf=None):
    '''Plot intensity frequency histogram and non-parametric estimates of component distributions'''
    plt.figure()
    plt.plot(bindat.intensities,bindat.counts,color="black")
    plt.plot(bindat.intensities,bindat.P1,color="red")
    plt.plot(bindat.intensities,bindat.P2,color="blue")
    plt.xlabel('Intensity')
    plt.ylabel('Frequency')
    plt.suptitle(fname)
    if pdf==None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()

def plotModel(bindat,pdf=None):
    '''Plot intensity density histogram, modelled distribution, component distributions and threshold estimate.'''
    plt.figure()
    plt.plot(bindat.intensities,bindat.freq,color="black")
    plt.plot(bindat.intensities,bindat.gauss1,color="red")
    plt.plot(bindat.intensities,bindat.gauss2,color="blue")
    plt.plot(bindat.intensities,bindat.mixed,color="green")
    plt.xlabel('Intensity')
    plt.ylabel('Density')
    plt.suptitle(fname)
    plt.axvline(x=thresh,linestyle='--',linewidth=0.5,color="darkorchid")
    if pdf==None:
        plt.show()
    else:
        pdf.savefig()
        plt.close()

def thresholdImage(arr,thresh,saveim=True):
    '''Thresholding array representation of a 16-bit image'''
    arr[arr<thresh]=0
    arr[arr>=thresh]=255
    arr=numpy.array(arr,dtype=numpy.uint8)
    arr=arr.reshape(im.size[::-1])
    imnew=Image.fromarray(arr, "L")
    if saveim:
        imnew.save(os.path.join(fullpath,fname+".png"))
    return(imnew)

def sizeSpots(locations,arr,thresharr):
    '''Adds intensity measures columns to locations dataFrame'''
    # Calculate area, intensity and trimmed intensity for each spot
    sumInt,sumArea,trim=[],[],[]
    for i in xrange(0,len(locations.x.values)):
        x,y,rad=locations.x.values[i],locations.y.values[i],int(math.ceil(max(locations.Diameter.values)))
        rad=1.5*rad # Scale this up for the moment, values from .gal file are underestimates
        tile=arr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        threshtile=thresharr[y-rad:(y+rad+1),x-rad:(x+rad+1)]
        sumInt.append(float(numpy.sum(tile))/float(tile.size))
        sumArea.append(float(numpy.sum(threshtile))/float(tile.size))
        trim.append(float(numpy.sum(tile[threshtile]))/float(tile.size))
    locations["Intensity"]=sumInt
    locations["Area"]=sumArea
    locations["Trimmed"]=trim
    return(locations)

def scanLocations(locations,arr,searchx,searchy):
    '''Brute force search for best location of grid specified in locations data frame which maximises the signal from arr at each gridpoint'''
    coords=numpy.array([(x,y) for x in xrange(searchx[0],searchx[1]) for y in xrange(searchy[0],searchy[1])],dtype=numpy.int)
    up=numpy.array([(x,max(searchy[0],y-1)) for (x,y) in coords],dtype=numpy.int)
    down=numpy.array([(x,min(searchy[1],y+1)) for (x,y) in coords],dtype=numpy.int)
    left=numpy.array([(max(searchx[0],x-1),y) for (x,y) in coords],dtype=numpy.int)
    right=numpy.array([(min(searchx[1],x+1),y) for (x,y) in coords],dtype=numpy.int)
    opt=numpy.array([numpy.sum(arr[numpy.array(locations.y+dy,dtype=numpy.int),numpy.array(locations.x+dx,dtype=numpy.int)]) for [dx,dy] in coords],dtype=numpy.int)
    opt=opt+numpy.array([numpy.sum(arr[numpy.array(locations.y+dy,dtype=numpy.int),numpy.array(locations.x+dx,dtype=numpy.int)]) for [dx,dy] in up],dtype=numpy.int)
    opt=opt+numpy.array([numpy.sum(arr[numpy.array(locations.y+dy,dtype=numpy.int),numpy.array(locations.x+dx,dtype=numpy.int)]) for [dx,dy] in down],dtype=numpy.int)
    opt=opt+numpy.array([numpy.sum(arr[numpy.array(locations.y+dy,dtype=numpy.int),numpy.array(locations.x+dx,dtype=numpy.int)]) for [dx,dy] in left],dtype=numpy.int)
    opt=opt+numpy.array([numpy.sum(arr[numpy.array(locations.y+dy,dtype=numpy.int),numpy.array(locations.x+dx,dtype=numpy.int)]) for [dx,dy] in right],dtype=numpy.int)
    soln=coords[numpy.argmax(opt)]
    return(tuple(soln))

def scanLocations(locations,arr,searchx,searchy,delta=0):
    '''Brute force search for best location of grid specified in locations data frame which maximises the signal from arr at each gridpoint'''
    x_grid, y_grid = numpy.meshgrid(numpy.arange(searchx[0],searchx[1]),numpy.arange(searchy[0],searchy[1]))
    x_grid, y_grid = x_grid.flatten(), y_grid.flatten()
    locs=[(y,x) for x in xrange(int(locations.x.values[i])-delta,int(locations.x.values[i])+delta+1) for y in xrange(int(locations.y.values[i])-delta,int(locations.y.values[i])+delta+1) for i in xrange(0,len(locations.x.values))]
    opt=numpy.array([numpy.sum(arr[(locations.y.values+y_grid[i]-delta):(locations.y.values+y_grid[i]-delta),(locations.x.values+x_grid[i]-delta):(locations.x.values+x_grid[i]+delta)]) for i in xrange(len(x_grid))],dtype=numpy.int)
    soln=coords[numpy.argmax(opt)]
    return(tuple(soln))

def getMaxima(intensity):
    '''Numerical method to find local maxima in a 1D list'''
    diffs=numpy.insert(numpy.diff(intensity),0,0)
    signs=numpy.sign(diffs)
    signdiffs=numpy.diff(signs)
    signdiffs=numpy.insert(signdiffs,0,0)
    maxima=numpy.where(signdiffs<0)[0]
    return(maxima)
    
def estimateOffsets(arr,diam=20,limFrac=1.075,showPlt=True,pdfPlt=False):
    '''Sum intensities along edge of array and search for the first maximum which is greater than limFrac*smoothed sum as offset estimate'''
    sumx=numpy.array([numpy.mean(arr[0:arr.shape[0],dx]) for dx in xrange(0,arr.shape[1])],dtype=numpy.int)
    sumy=numpy.array([numpy.mean(arr[dy,0:arr.shape[1]]) for dy in xrange(0,arr.shape[0])],dtype=numpy.int)
    smoothedx=signal.medfilt(sumx,2*diam+1)
    smoothedy=signal.medfilt(sumy,2*diam+1)
    limx,limy=limFrac*smoothedx,limFrac*smoothedy
    maxx,maxy=getMaxima(sumx),getMaxima(sumy)
    candx=maxx[sumx[maxx]>limx[maxx]]
    candy=maxy[sumy[maxy]>limy[maxy]]
    if showPlt:
        plt.plot(sumx)
        plt.plot(smoothedx)
        plt.plot(limx)
        for cand in candx:
            plt.axvline(x=cand,linestyle='--',linewidth=0.5,color="black")
        plt.suptitle(fname+"dx")
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
        plt.suptitle(fname+"dy")
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
    
def showIm(arr,sigMax=2000.0,returnIm=False):
    '''Quick 8-bit preview images from float arrays, useful for debugging'''
    imarr=numpy.minimum(255.0*arr/numpy.float(sigMax),255.0)
    imarr=numpy.array(imarr,dtype=numpy.uint8)
    imnew=Image.fromarray(imarr,"L")
    if returnIm:
        return(imnew)
    else:
        imnew.show()

def optimiseSpot(arr,x,y,rad,RAD):
    '''Search in a square of width 2*RAD for the centre of a smaller square (of width 2*rad) which maximises the sum of signal inside the small square'''
    x_grid, y_grid = numpy.meshgrid(numpy.arange(x-(RAD-rad),x+(RAD-rad)),numpy.arange(y-(RAD-rad),y+(RAD-rad)))
    x_grid=x_grid.flatten()
    y_grid=y_grid.flatten()
    signal=numpy.array([numpy.sum(arr[(y_grid[i]-rad):(y_grid[i]+rad),(x_grid[i]-rad):(x_grid[i]+rad)]) for i in xrange(0,len(x_grid))])
    iopt=numpy.argmax(signal)
    return(x_grid[iopt],y_grid[iopt])

# Find what directory we're in
syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

# Get all histogram files in current directory
tiflist,gallist=[],[]
allfiles=os.listdir(fullpath)
for filename in allfiles:
    if filename[-4:] in ('.tif','.TIF','.tiff','.TIFF'):
        tiflist.append(filename[0:-4])
    if filename[-4:] in ('.gal','.GAL'):
        gallist.append(filename)

if len(gallist)>1:
    print "Warning!  Multiple .gal files found.  Only first one will be used..."

pdf=PdfPages(os.path.join(fullpath,'HistogramReport.pdf'))

for fname in tiflist:
    print fname
    # Build pixel intensity histogram from image
    im=Image.open(os.path.join(fullpath,fname+".tif"))
    # Convert to array, forcing to row-column format (i.e y-coordinate, then x-coordinate)
    arr=numpy.array(im.getdata(),dtype=numpy.uint16).reshape(im.size[::-1])

    # Truncate at 99th percentile
    sigMax=1.5*stats.mstats.mquantiles(arr,0.99)
    arr_trunc=arr[:,:]
    arr_trunc[arr>sigMax]=sigMax
    
    # Read in spot locations from existing file
    if os.path.exists(fname+".dat"):
        locations=pandas.read_csv(os.path.join(fullpath,fname+".dat"),sep="\t")
        search=int(numpy.mean(locations.Diameter))
    else:    
        galDict=parseGal(gallist[0])
        locations=zeroLocations(predictLocations(galDict))
        spots=galDict["spotIDs"]
        locations=locations.merge(spots,on=("Block","Row","Column"))
        search=int(numpy.mean(locations.Diameter))
        (dx,dy)=estimateOffsets(arr_trunc,diam=search,showPlt=True,pdfPlt=True)
        
        if dx<0 or dy<0:
            raise Exception("Automatic spot location has failed.  Is there any signal in image "+fname+".tif ?")
            break
        locations.x+=dx
        locations.y+=dy

        # Update each block in the grid individually, over a smaller search space
        for bno in numpy.unique(locations.Block):
            scaleBlock=1
            block=locations[locations.Block==bno]
            minx,maxx=max(0,min(block.x)-int(round(scaleBlock*search))),min(im.size[0],max(block.x)+int(round(scaleBlock*search)))
            miny,maxy=max(0,min(block.y)-int(round(scaleBlock*search))),min(im.size[1],max(block.y)+int(round(scaleBlock*search)))
            block_arr=arr_trunc[miny:maxy,minx:maxx]
            (dx,dy)=estimateOffsets(block_arr,diam=search,limFrac=1.0725,showPlt=True,pdfPlt=True)
            # Reject updates if optimum has not been found 
            if dx<0:
                dx,dy=0,0
            else:
                dx,dy=(dx+minx)-min(block.x),(dy+miny)-min(block.y)
            # Reject updates if they are too large
            if dx>int(round(scaleBlock*search)):
                dx=0
            if dy>int(round(scaleBlock*search)):
                dy=0
            locations.x[locations.Block==bno]+=dx
            locations.y[locations.Block==bno]+=dy
        
        # Update each spot individually
        rad=int(round(float(search/2.0)))
        RAD=search    
        for i in xrange(0,len(locations.x)):
            (dx,dy)=optimiseSpot(arr_trunc,locations.x[i],locations.y[i],rad,RAD)
            locations.x[i]=dx
            locations.y[i]=dy

    # New arrays for lighting correction
    cutout_arr=numpy.array(arr,dtype=numpy.float)
    mask_arr=numpy.zeros(arr.shape,dtype=numpy.bool8)

    # Increase estimate of spot radius by 10% (it seems to be an understimate) before cutting out spot locations
    # Important to cut out all signal, otherwise holes have bright edges, and will be filled by bright pixels instead of backround
    bigsearch=1.5*search
    
    for i in xrange(0,len(locations.x)):
        x_grid, y_grid = numpy.meshgrid(numpy.arange(locations.x.values[i]-float(bigsearch)/2,locations.x.values[i]+float(bigsearch)/2), numpy.arange(locations.y.values[i]-float(bigsearch)/2,locations.y.values[i]+float(bigsearch)/2))
        disk = ((x_grid-locations.x.values[i])**2 + (y_grid-locations.y.values[i])**2) <= (float(bigsearch)/2)**2
        x_res,y_res=numpy.array(x_grid[disk],dtype=numpy.int),numpy.array(y_grid[disk],dtype=numpy.int)
        mask_arr[y_res,x_res]=True

    (y_list,x_list)=numpy.where(mask_arr)
    cutout_arr[mask_arr]=numpy.nan
    hole_arr=numpy.copy(cutout_arr)
    old=numpy.zeros(cutout_arr.shape,dtype=float)
    
    # Tolerance for average pixel intensity difference between iterations to declare convergence of Markov update
    tol=5
    while numpy.sum(numpy.abs(old-cutout_arr))/numpy.sum(mask_arr)>tol or numpy.isnan(numpy.sum(numpy.abs(old-cutout_arr))/numpy.sum(mask_arr)):
        old=numpy.copy(cutout_arr)
        # Markov field update
        for i in xrange(0,len(x_list)):
            plist=[cutout_arr[y_list[i],x_list[i]+1],cutout_arr[y_list[i]+1,x_list[i]],cutout_arr[y_list[i],x_list[i]-1],cutout_arr[y_list[i]-1,x_list[i]]]
            cutout_arr[y_list[i],x_list[i]]=stats.nanmean(plist)
    
    # Flatten intensities to correct for lighting gradient
    smoothed_arr=ndimage.gaussian_filter(cutout_arr,search/2)
    average_back=numpy.mean(smoothed_arr)
    # Assume that differences between mean plate background intensities are not interesting
    average_back=100.0
    corrected_arr=arr*average_back/smoothed_arr

    # Initial guess for mixed model parameters for thresholding lighting corrected image
    (counts,intensities)=numpy.histogram(corrected_arr,bins=2**16,range=(0,2**16))
    intensities=numpy.array(intensities[0:-1],dtype=numpy.int)
    (bindat,[theta,mu1,mu2,sigma1,sigma2])=initialGuess(intensities,counts)

    # Maximise likelihood of 2-component mixed Gaussian model parameters given binned observations by constrained optimisation
    logL=makeObjective(bindat.intensities,bindat.counts,totFunc)
    b=[(0.5,0.85),(float(mu1)/5.0,5*float(mu1)),(float(mu2)/5.0,5.0*float(mu2)),(float(sigma1)/5.0,5.0*float(sigma1)),(float(sigma2)/5.0,5.0*float(sigma2))]
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
    thresharr=numpy.array(numpy.array(imthresh,dtype=numpy.int)/255,dtype=numpy.bool8)
    erodearr=ndimage.binary_erosion(thresharr,iterations=2)

    # Estimate and store spot sizes, areas, IODs together with previously estimated coordinates, ids and locations
    locations=sizeSpots(locations,corrected_arr,erodearr)
    locations.to_csv(fname+".out","\t",index=False)

    # Reporting results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Write plots to a .pdf report
    plotGuess(bindat,pdf)
    plotModel(bindat,pdf)

    # Stretch histogram of original image and create a preview showing spot locations
    stretch=im.point(lambda i:i*(255.0/sigMax)).convert('RGB')
    draw=ImageDraw.Draw(stretch)
    for i in xrange(0,len(locations.x)):
        x,y,r=locations.x[i],locations.y[i],float(locations.Diameter[i])/2.0
        draw.ellipse((x-r,y-r,x+r,y+r),outline=(255,0,0))
    stretch.save(fname+"_LOCATE.jpg")

    # Write image with holes cutout to file
    imcut=showIm(hole_arr,sigMax,returnIm=True)
    imcut.save(fname+"_CUT.jpg")

    # Write corrected image preview to file
    imnew=showIm(corrected_arr,sigMax,returnIm=True)
    imnew.save(fname+"_CORR.jpg")

    # Write thresholded image to file
    reconim=Image.fromarray(numpy.array(erodearr*255,dtype=numpy.uint8),mode="L")
    reconim.save(fname+"_THRESH.png")
    
    # Label spots and intensities on corrected image
    r=int(math.ceil(max(locations.Diameter.values)))
    imnew=imnew.convert("RGB")
    draw=ImageDraw.Draw(imnew)
    for i in xrange(0,len(locations.x)):
        lab=str(locations.Name.values[i])
        if lab=="nan":
            lab=""
        intense=str(round(locations.Intensity.values[i],2))
        draw.rectangle((locations.x.values[i]-r,locations.y.values[i]-r,locations.x.values[i]+r,locations.y.values[i]+r),outline="red")
        draw.text((locations.x.values[i]+20,locations.y.values[i]-5),lab,fill="green")
        draw.text((locations.x.values[i]+20,locations.y.values[i]+5),intense,fill="yellow")
    imnew.save(fname+"_LAB.png")
    
    print fname,thresh,theta_opt,mu1_opt,mu2_opt,sigma1_opt,sigma2_opt,average_back

pdf.close()
