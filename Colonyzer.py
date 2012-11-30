# Set this to a positive value (between 0 and 255) to turn off automatic thresholding of image
# For example, 125 is a good value for the demo image
fixThreshold=-99
# Set this to False to disable background lighting correction
lightCorrect=True
# Set this to False to disable search for optimum colony location (i.e. just accept regular grid guess)
colonySearch=True
# Set this to False to shrink images before analysis (saves memory, but can be less precise)
originalSize=True

try:
    import psyco
    psyco.full()
except ImportError:
    pass

import numpy,os,sys,time,math,gc
from PIL import Image, ImageDraw, ImageFilter, ImageOps

#### Let's import all the functions in R
### If you have multiple copies of R installed, you may need to set RHOME like this:
from rpy_options import set_options
set_options(RHOME="C:\\Program Files\\R\\R-2.12.1\\")
#These lines are just included to test your RPy environment
#If you're happy it's working, you can comment out the next 4 lines
from rpy import r
r.library("stats")
r.library("rgenoud")
del r
gc.collect()

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def SetUp(instructarr):
    # Check if the first element is an integer:
    if is_number(instructarr[0]):
        NoSpots=int(instructarr[0])
    else:
        NoSpots=instructarr[0]
    TLX,TLY,BRX,BRY=instructarr[1],instructarr[2],instructarr[3],instructarr[4]
    global nocols,norows, tlx,tly,brx,bry,xdim,ydim,xdimf,ydimf,xstart,ystart,NSearch
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

    # Best estimates for the starting coordinates
    xstart=max(0,int(round(float(tlx)-0.5*xdimf)))
    ystart=max(0,int(round(float(tly)-0.5*ydimf)))

    # Parameter for specifying the search area (area is (NSearch*2+1)^2)
    NSearch = int(round(3.0*float(min(xdim,ydim))/8.0))
    print "Instructions: ",xstart,ystart,xdim,ydim,NSearch

# Let's find what directory we're in
syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)

# Let's try to read in the Colonyzer input file
Instructions=open(os.path.join(fullpath,'Colonyzer.txt'),'r')
InsData={}
InsTemp=Instructions.readlines()

for x in xrange(0,len(InsTemp)):
    if InsTemp[x][0]!="#" and InsTemp[x][0]!="\n":
        tlist=InsTemp[x].split(',')
        if len(tlist)>1:
            InsData[tlist[0]]=[tlist[1],int(tlist[2]),int(tlist[3]),int(tlist[4]),int(tlist[5])]

if 'default' in InsData:
    SetUp(InsData['default'])
else:
    print "ERROR: No default instructions"
    sys.exit()

# If you call the script via :
# python ImageAnalysis1.6.py ImageFile.jpg
# it will process the single file ImageFile.jpg
# Otherwise it will process all .jpg files in the current directory

singlefile=False
# Switch for "verbose" graphics output - currently fails for very large image input
graphicsout=False
params=sys.argv[1:]

if len(params)>0 and params[0]!='ALL':
    singlefile=True
    imagetoscan=params[0]

# Width of strip to take off edge to correct for lighting problems
LightStrip=0
MaximiseOption=1

# Make this a global parameter for use in output
global medintensity
# Give filename for output file
if singlefile:
    outfile=imagetoscan.split('.')[0]+"OUT.dat"
else:
    outfile="Output.dat"

# Create directories for storing output data and preview images
try:
    os.mkdir(os.path.join(fullpath,"Output_Images"))
    os.mkdir(os.path.join(fullpath,"Output_Data"))
except:
    print("Output directories already present")

# We could set outputimages (directory for storing images) to be different to current directory
if graphicsout:
    outputimages=fullpath
else:
    # Print where python thinks it lives 
    print "Working Directory: ",os.getcwd()
    outputimages=os.path.join(fullpath,"Output_Images")
outputdata=os.path.join(fullpath,"Output_Data")

def smooth(x,window_len=10,window='hanning'):
    # This is a general smoothing function
    # From http://www.scipy.org/Cookbook/SignalSmooth
    """smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    see also: 
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    x=numpy.array(x)
    #if x.ndim != 1:
    #    raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def getMin(list):
    # Return the position of the minimum in a list
    # Ok, in Python 2.5 this works:
    #best = min(xrange(len(list)),key=list.__getitem__)
    #return best
    # For Python 2.3 (on cisbclust) we will have to do this manually
    minplace=len(list)+1
    minimum=min(list)
    for element in xrange(len(list)):
        if list[element]==minimum:
            minplace=element
            return minplace

def getMax(list):
    # Return the position of the minimum in a list
    # Ok, in Python 2.5 this works:
    #best = max(xrange(len(list)),key=list.__getitem__)
    #return best
    # For Python 2.3 (on cisbclust) we will have to do this manually
    maxplace=len(list)+1
    maximum=max(list)
    for element in xrange(len(list)):
        if list[element]==maximum:
            maxplace=element
            return maxplace

# http://bitecode.co.uk/2008/07/edge-detection-in-python/
def sobel(img):
    starting=time.time()
    if img.mode != "RGB":
        img = img.convert("RGB")
    out_image = Image.new(img.mode, img.size, None)
    imgdata = img.load()
    outdata = out_image.load()
    gx = [[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]
    gy = [[-1, -2, -1], [0, 0, 0], [1, 2, 1]]
    for row in xrange(1, img.size[0]-1):
        for col in xrange(1, img.size[1]-1):
            pixel_gx = pixel_gy = 0
            pxval = sum(imgdata[row,col])/3
            for i in range(-1, 2):
                for j in range(-1, 2):
                    val = sum(imgdata[row+i,col+j])/3
                    pixel_gx += gx[i+1][j+1] * val
                    pixel_gy += gy[i+1][j+1] * val
            newpixel = math.sqrt(pixel_gx * pixel_gx + pixel_gy * pixel_gy)
            newpixel = 255 - int(newpixel)
            outdata[row, col] = (newpixel, newpixel, newpixel)
    finishing=time.time()
    print "Sobel Edge Detection Took: ", finishing-starting, " seconds."
    del imgdata
    del outdata
    gc.collect()
    return out_image

# http://alwaysmovefast.com/2007/12/05/basic-edge-detection-in-python/
# Uses hashes of tuples to simulate 2-d arrays for the masks.  
def get_prewitt_masks():
    xmask = {}
    ymask = {}  
    xmask[(0,0)] = -1  
    xmask[(0,1)] = 0  
    xmask[(0,2)] = 1  
    xmask[(1,0)] = -1  
    xmask[(1,1)] = 0  
    xmask[(1,2)] = 1  
    xmask[(2,0)] = -1  
    xmask[(2,1)] = 0  
    xmask[(2,2)] = 1  
    
    ymask[(0,0)] = 1  
    ymask[(0,1)] = 1  
    ymask[(0,2)] = 1  
    ymask[(1,0)] = 0  
    ymask[(1,1)] = 0  
    ymask[(1,2)] = 0  
    ymask[(2,0)] = -1  
    ymask[(2,1)] = -1  
    ymask[(2,2)] = -1  

    return (xmask, ymask)

def prewitt(img):
    starting=time.time()
    if img.mode != 'L': img = img.convert('L')
    pixels = list(img.getdata())
    width,height=int(im.size[0]),int(im.size[1])
    xmask,ymask = get_prewitt_masks()  

    # create a new greyscale image for the output  
    outimg = Image.new('L', (width, height))  
    outpixels = list(outimg.getdata())  

    for y in xrange(height):  
        for x in xrange(width):  
            sumX, sumY, magnitude = 0, 0, 0  

            if y == 0 or y == height-1: magnitude = 0  
            elif x == 0 or x == width-1: magnitude = 0  
            else:  
                for i in xrange(-1, 2):  
                    for j in xrange(-1, 2):  
                        # convolve the image pixels with the Prewitt mask, approximating dI / dx  
                        sumX += (pixels[x+i+(y+j)*width]) * xmask[i+1, j+1]
                for i in xrange(-1, 2):  
                    for j in xrange(-1, 2):  
                        # convolve the image pixels with the Prewitt mask, approximating dI / dy  
                        sumY += (pixels[x+i+(y+j)*width]) * ymask[i+1, j+1]  
                # approximate the magnitude of the gradient  
                magnitude = abs(sumX) + abs(sumY)
            if magnitude > 255: magnitude = 255
            if magnitude < 0: magnitude = 0  
            outpixels[x+y*width] = 255 - magnitude
    outimg.putdata(outpixels)
    finishing=time.time()
    print "Prewitt Edge Detection Took: ", finishing-starting, " seconds."
    return outimg  

def automaticSetThreshold(histim):
    # Skip automatic thresholding...
    if fixThreshold >=0:
        return(fixThreshold)
    
    # Let's just try to import r
    from rpy import r
    r.library("stats")
    r.library("rgenoud")
    # Cut out the inner tiles out of the image
    x0=xstart+LightStrip*xdim
    y0=ystart+LightStrip*ydim
    xf=xstart+xdim*(nocols-LightStrip)
    yf=ystart+ydim*(norows-LightStrip)
    histpix=histim.load()
    whist,hhist=xf-x0,yf-y0
    # Set these to be negative so we can test for sensible values later
    darkest,brightest=-99,-99
    onepeak=False
    # Convert image to grayscale
    newim=histim.convert('L')
    box = (x0, y0, xf, yf)
    region = newim.crop(box)
    histlist=region.histogram()
    # Smooth the histogram, and find the darkest and brightest maxima in the cropped image
    listnew=list(smooth(histlist))
    # Convert to frequency to allow for differently-sized image files
    sumnew=float(sum(listnew))
    scalednew=[]
    for elem in listnew:
        scalednew.append(float(elem)/sumnew)
    # Moving from dark to white, get a list of maxima
    global maxima
    maxima=[]
    # Check for truncated peaks (maxima at edges)
    if listnew[0]>listnew[1]:
        maxima.append((1,listnew[1]))
    for pixintensity in xrange(1,254):
        # The last condition is to help prevent minor wiggles being identified as extrema
        if listnew[pixintensity]>=listnew[pixintensity-1] and listnew[pixintensity+1]<listnew[pixintensity] and scalednew[pixintensity]>0.000141285:
            maxima.append((pixintensity,listnew[pixintensity]))
    # Check for truncated peaks (maxima at edges)
    if listnew[255]>listnew[254]:
        maxima.append((255,listnew[255]))
    # Now we have our list of maxima, let's find the two biggest ones, and get ready to threshold
    # somewhere between those
    highest=0
    highind=-99
    secondhighest=0
    secondhighind=-99
    for maximum in maxima:
        if maximum[1]>highest:
            highest=maximum[1]
            highind=maximum[0]
    for maximum in maxima:
        if maximum[1]>secondhighest and maximum[1]!=highest:
            secondhighest=maximum[1]
            secondhighind=maximum[0]
    # If we've identified more than one peak
    if len(maxima)>1:
        # Our estimate of mu1 should be the smallest of these two
        mu1est=min(highind,secondhighind)
        mu2est=max(highind,secondhighind)
    else:
        # We better flag this as it may affect our optimisation choices below
        onepeak=True
        # Assume dark background (highest peak is dark peak)
        mu1est=maxima[0][0]
        # Guess next peak is 35 units along
        mu2est=min(mu1est+35,255)
    thetabest=float(listnew[int(round(mu2est))])/float(listnew[int(round(mu1est))])
    # Let's try fitting two Gaussian curves by least squares in R
    # First build a list of x,y frequency data
    freq=[]
    intense=[]
    for i in xrange(0,256):
        intense.append(i)
        freq.append(float(histlist[i])/float(sum(histlist)))
    # Pass the original histogram to R
    r.assign('histlist',histlist)
    # Pass the smoothed one used for peak estimation
    r.assign('smoothlist',listnew)
    # Make sure expected values are away from the edges (for saturated image, although
    # the maximum value is on the edge, the expected value is somewhat in from the edge)
    E1=max(mu1est,10.0)
    E2=min(mu2est,250.0)
    print thetabest, mu1est, mu2est,
    # Increment these by one to account for the converting 0:255 -> 1:256
    E1=E1+1
    E2=E2+1
    r.assign('E1',E1)
    r.assign('E2',E2)
    r.assign('thetabest',thetabest)
    r.assign('histname',os.path.join(fullpath,'RHistogram'+FILENAME +'.PNG'))
    # This dictates whether or not the mixed model parameter estimate optimisation procedure
    # is global (rgenoud)
    if MaximiseOption==0:
        r.assign('globopt',0)
    if MaximiseOption==1:
        r.assign('globopt',1)
    if MaximiseOption==2:
        if onepeak==False:
            r.assign('globopt',0)
        else:
            r.assign('globopt',1)
    # Execute the fitting function (sum of two Gaussians, second weighted by alpha
    # to give a small second peak).  The fit seems fairly robust to guesses for mu1 and mu2
    try:
        r("""
freq<-c(histlist[1:256]/sum(histlist))
intlist<-c(1:256)
# Define the normal and lognormal distributions
Norm<-function(x,sigma,mu){
	(1.0/(sigma*sqrt(2*pi)))*exp((-(x-mu)^2)/(2*sigma^2))
}
# Define the bimodal function which will approximate the histogram
totFunc<-function(x,p){
     theta<-p[1]
     mu1<-p[2]
     mu2<-p[3]
     sigma1<-p[4]
     sigma2<-p[5]
     # We need to penalise conditions where the two peaks are lined up
     # on top of each other (effectively one peak!)
     if (abs(mu2-mu1)<7){0}else{
     Norm(x,sigma1,mu1)/(1+theta)+theta*Norm(x,sigma2,mu2)/(1+theta)}
}
diffFunc<-function(x,p){
     theta<-p[1]
     mu1<-p[2]
     mu2<-p[3]
     sigma1<-p[4]
     sigma2<-p[5]
     Norm(x,sigma1,mu1)/(1+theta)-theta*Norm(x,sigma2,mu2)/(1+theta)
}
LogL<-function(p){
   # I've added 0.00000001 here to avoid problems with taking log of zero (zero to machine precision)
    sum(-1*histlist*log(0.00000001+
    totFunc(intlist,p)
    ))
}
# We need to estimate the variances...
# Roughly, if we just split the histogram into high and low (at the average of E1 and E2)
# we can calculate means and variances...
lowlist=histlist[1:E1]
highlist=histlist[E2:256]

# Get "half-peaks" to estimate variances
lowvar=sum(histlist[1:E1]*(E1-((1:E1)-1))^2)/sum(histlist[1:E1])
highvar=sum(histlist[E2:256]*(E2-((E2:256)-1))^2)/sum(histlist[E2:256])

theta.start<-thetabest

# Estimate LogNormal parameters from guesses for Expected values and Variances above
mu1.start<-E1
mu2.start<-E2
sigma1.start<-sqrt(lowvar)
sigma2.start<-sqrt(highvar)
# For the normal distribution, we can improve the guess for theta by incorporating our
# guesses for sigma1 and sigma2.  Technically, the estimate which has been passed to RPy
# is just the ratio:
# function at x=mu2 / function at x=mu1
# If we substitute those two values into the total function, we find that we should
# multiply this ratio by sigma2/sigma1 to have an estimate of theta
# In other words, using this ratio directly assumes that sigma1=sigma2
theta.start=theta.start*sigma2.start/sigma1.start
bestinit<-c(theta.start,mu1.start,mu2.start,sigma1.start,sigma2.start)

dommat<-matrix(c(
0.001,
1.0,
20,
0.1,
0.1, 
50.0,
235.0,
255.0,
7.0,
7.0
),nrow=5,ncol=2,byrow=FALSE)
if (globopt==0){    
# Note that we restrict mu1 to be the dark peak and mu2 to be the bright one
# Otherwise we lose the non-negative advantage of the two lognormals
out<-optim(
    c(theta.start,mu1.start,mu2.start,sigma1.start,sigma2.start),
    LogL,
    method="L-BFGS-B",
    lower=list(0.01,0.0,mu1.start,0.001,0.001),
    upper=list(200.0,mu2.start,255.0,75.0,75.0)
)
reslist<-out$par
}else{
pergen<-1500
cov<-0.5
init<-matrix(c(
rnorm(pergen/2,mean=theta.start,sd=cov*theta.start),
rnorm(pergen/2,mean=mu1.start,sd=cov*mu1.start),
rnorm(pergen/2,mean=mu2.start,sd=cov*mu2.start),
rnorm(pergen/2,mean=sigma1.start,sd=cov*sigma1.start),
rnorm(pergen/2,mean=sigma2.start,sd=cov*sigma2.start)
),nrow=5,ncol=pergen,byrow=FALSE)
globout<-genoud(
LogL,
nvars=5,
starting.values=init,
Domains=dommat,
hard.generation.limit=TRUE,
boundary.enforcement=2,
pop.size=pergen,
solution.tolerance=1,
print.level=0)
reslist<-globout$par
}

# Calculate final expected values (in case initial guesses were bad)
# Instead of the expected values, lets look at the medians!
EFinal1=reslist[2]
EFinal2=reslist[3]
#Reset root variables
r1=-99
r2=-99
r3=-99
root=-99

# Note that sometimes these roots don't exist so we need to check for that with try
# Search for the root before the darkest peak
try((r1=uniroot(diffFunc,c(1,EFinal1),p=reslist)$root),silent=TRUE)

# Let's make sure these are in the order we would expect!
EMin=min(EFinal1,EFinal2)
EMax=max(EFinal1,EFinal2)

# Search for the root between the two peaks
try((r2=uniroot(diffFunc,c(EMin,EMax),p=reslist)$root),silent=TRUE)

# Search for the root after the brightest peak
try((r3=uniroot(diffFunc,c(EFinal2,256),p=reslist)$root),silent=TRUE)

# Now find which one of these gives the largest value of totFunc
roots<-c(r1,r2,r3)
if(r1>0) roots<-append(roots,r1)
if(r2>0) roots<-append(roots,r2)
if(r3>0) roots<-append(roots,r3)
roots<-roots[4:length(roots)]
roots
values<-c(totFunc(roots,reslist))
values

# Finally find the index of the maximum of the list of values, this is our best root
root<-roots[rev(order(values))[1]]
# Remember that this is in 1:256 scale, may need to subtract 1 to change to 0:255
root

totpts<-totFunc(1:256,reslist)

# If searching by segment has failed just do one last search over the whole range
if(root<0) try((root=uniroot(diffFunc,c(1,256),p=reslist)$root),silent=TRUE)
##png(filename = histname, width = 700, height = 700, units = "px",
##pointsize = 12, bg = "white", res = NA, restoreConsole = TRUE)
##plot(totFunc(1:256,reslist),type="l",col="red",xlab="pixel intensity", ylab="frequency",xlim=c(1,256),ylim=c(0,max(freq)))
##points(freq)
##points(Norm(1:256,reslist[4],reslist[2])/(1+reslist[1]),type="l",col="blue")
##points(reslist[1]*Norm(1:256,reslist[5],reslist[3])/(1+reslist[1]),type="l",col="blue")
##points(c(root,root),c(0,max(freq)),type="l",col="green")
### Draw in the initial guesses for the peaks
##points(c(E1,E2),c(freq[E1],freq[E2]),col="red")
##dev.off()
""")
    except:
        print "RPy Exception"
    # This is to work around an rpy bug!
    r.histlist[0]+=1
    r.histlist[0]-=1
    del r.histlist

    SIGMA1=r.reslist[3]
    SIGMA2=r.reslist[4]
    MU1=r.reslist[1]
    MU2=r.reslist[2]
    THETA=r.reslist[0]
    del r.reslist
    print "Theta, Mu1, Mu2, Sigma1, Sigma2"
    print THETA, MU1,MU2,SIGMA1,SIGMA2
    print "Roots: ",r.roots
    del r.roots

    # Convert back to 
    optthresh=r.root-1
    print "Threshold:", optthresh
    print "Using threshold from double Gaussian fit: "+str(optthresh)
    if graphicsout:
        # Draw the histogram for trouble shooting
        plotHistograms(histlist,r.totpts,optthresh)
    #print histlist
    #print r.totpts
    try:
        del r.totpts, r.root
    except:
        print("RPy not deleting r.totpts and r.root properly again...")
    try:
        # Destroy the r object
        del r,newim
    except:
        print("RPy not deleting r object properly again...")
    gc.collect()
    return optthresh

def plotHistograms(histlist,funcpts,cut):
    # Let's make a quick plot of the image histograms which we can scan by eye
    # to make sure nothing weird is happening
    plotHeight=500
    widthFactor=4
    plotWidth=widthFactor*len(histlist)
    # Scale the frequencies to fit in the image
    freqmax=float(max(max(histlist),max(funcpts)))
    histcoords=[]
    funccoords=[]
    for i in xrange(len(histlist)):
        histcoords.append((i*widthFactor,int(round(plotHeight*float(histlist[i])/freqmax))))
        funccoords.append((i*widthFactor,int(round(plotHeight*float(funcpts[i])))))
    plot=Image.new("RGB",(plotWidth,plotHeight),'white')
    draw = ImageDraw.Draw(plot)

    # Draw the curves
    draw.line(histcoords,fill="blue")
    draw.line(funccoords,fill="red")
    # Draw the selected threshold
    draw.line(((widthFactor*cut,0),(widthFactor*cut,plotHeight-1)),fill="green")
    del draw

    # Flip the image about the centre horizontal axis to account for the coordinate system
    plot=plot.transpose(Image.FLIP_TOP_BOTTOM)
    # Save the image
    plot.save(os.path.join(outputimages,FILENAME +'Histogram.png'))
    del plot
    gc.collect()

def BackgroundCorrection(im,threshim):
    # Let's do the thresholding with the PIL point function
    # Strip the colour out of the image
    greyim=im.convert(mode="L")
    greypix=numpy.array(greyim,numpy.uint8)
    # Strongly smooth it before doing the initial thresholding
    for count in xrange(5):
        greyim=greyim.filter(ImageFilter.SMOOTH_MORE)
    # SAVE IMAGE
    if graphicsout:
        greyim.save(os.path.join(outputimages,'Smoothed.PNG'))
    threshpix=numpy.array(threshim,numpy.uint8)
    # Now we want to calculate the median intensity of pixels which are marked as background within the
    # relevant frame of tiles (Should this be the mode intensity?)
    intensities=[]
    for x in xrange(xstart+xdim*LightStrip,min(w-1,xstart+xdim*(nocols-LightStrip))):
        for y in xrange(ystart+ydim*LightStrip,min(h-1,ystart+ydim*(norows-LightStrip))):
            if threshpix[y,x]<128:
                # This is background
                intensities.append(greypix[y,x])
    global medintensity

    medintensity=numpy.median(intensities)

    # Skip background correction
    if lightCorrect==False:
        return(im)
    print "Median Background Intensity: ",medintensity
    # Disable the correction!!!!
    #return im
    # Now we copy the greyim, and smooth it to have less noisy interpolation across the holes
    greycopy=greyim.copy()
    # A clean copy of the image as an array, which will not be updated
    cleanpix=numpy.array(greycopy,numpy.float)
    # Convert image to numpy array
    numpixH=numpy.array(greycopy,numpy.float)
    # Convert image to numpy array
    numpixV=numpy.array(greycopy,numpy.float)
    # Ok, let's scan through this image, and everywhere there is a zero, let's
    # replace it with the median of its neighbours that are not zero
    # This parameter gives the size of the neighbourhood:(scandim*2+1)x(scandim*2+1)
    print "Horizontal scan"
    x=0
    y=0
    while y<h:
        if y%100==0: print y," out of ",h
        while x<w:
            # Check if we are at the beginning of a hole
            if threshpix[y,x]>0:
                if x==0:
                    P1=0
                    # This is the left edge of the image
                    # March along to the first background pixel
                    while threshpix[y,x]>0 and x<w-1:
                        x=x+1
                    # Now we've found the first background pixel
                    P2=x
                else:
                    # We are away from the left edge
                    # x0 is the last background pixel
                    P1=x-1
                    while threshpix[y,x]>0 and x<w-1:
                        x=x+1
                    P2=x
                # Now we have identified the two edges of the hole
                # where P1 and P2 are background pixels at the edge
                # we can search away from these and look for minimum backgrounds
                # within the range xdim/2 (assume no light gradient at spot scale)
                if P1<=0 and P2>=w-1:
                    # We haven't found any edges, just put in the median background
                    numpixH[y,0:w]=medintensity
                    bestbackground1=medintensity
                    bestbackground2=medintensity
                elif P1<=0:
                    # We only have one edge, copy best value for that edge across gap
                    rangemax=min(w-1,P2+xdim/2)
                    bestbackground2=min(cleanpix[y,P2:rangemax])
                    bestbackground1=bestbackground2
                elif P2>=w-1:
                    # We have only one edge, copy best value for that edge across gap
                    rangemin=max(0,P1-xdim/2)
                    bestbackground1=min(cleanpix[y,rangemin:P1])
                    bestbackground2=bestbackground1
                else:
                    # We are somewhere in the middle
                    rangemax=min(w-1,P2+xdim/2)
                    rangemin=max(0,P1-xdim/2)
                    bestbackground1=min(cleanpix[y,rangemin:P1])
                    bestbackground2=min(cleanpix[y,P2:rangemax])
                # Now we have our two points and our best background values for those points,
                # we can linearly interpolate between them
                m=(bestbackground2-bestbackground1)/(P2-P1)
                for P in xrange(P1+1,P2):
                    numpixH[y,P]=m*(P-P1)+bestbackground1
                # Now we can set x to be the next point after P2
                x=P2+1
            x+=1              
        y+=1
        x=0
    if graphicsout:
        numpix1=numpy.array(numpixH,numpy.uint8)
        pim=Image.fromarray(numpix1,'L')
        pim.save(os.path.join(outputimages,'Horizontal.PNG'))
        del pim,numpix1
        gc.collect()
    # Ok, let's scan through this image, and everywhere there is a zero, let's
    # replace it with the median of its neighbours that are not zero
    # This parameter gives the size of the neighbourhood:(scandim*2+1)x(scandim*2+1)
    print "Vertical scan"
    x=0
    y=0
    while x<w:
        if x%100==0: print x," out of ",w
        while y<h:
            # Check if we are at the beginning of a hole
            if threshpix[y,x]>0:
                if y==0:
                    P1=0
                    # This is the left edge of the image
                    # March along to the first background pixel
                    while threshpix[y,x]>0 and y<h-1:
                        y=y+1
                    # Now we've found the first background pixel
                    P2=y
                else:
                    # We are away from the left edge
                    # x0 is the last background pixel
                    P1=y-1
                    while threshpix[y,x]>0 and y<h-1:
                        y=y+1
                    P2=y
                # Now we have identified the two edges of the hole
                # where P1 and P2 are background pixels at the edge
                # we can search away from these and look for minimum backgrounds
                # within the range ydim (assume no light gradient at spot scale)
                if P1<=0 and P2>=h-1:
                    # We haven't found any edges, just put in the median background
                    numpixV[0:h,x]=medintensity
                    bestbackground1=medintensity
                    bestbackground2=medintensity
                elif P1<=0:
                    # We only have one edge, copy best value for that edge across gap
                    rangemax=min(h-1,P2+ydim/2)
                    bestbackground2=min(cleanpix[P2:rangemax,x])
                    bestbackground1=bestbackground2
                elif P2>=h-1:
                    # We have only one edge, copy best value for that edge across gap
                    rangemin=max(0,P1-ydim/2)
                    bestbackground1=min(cleanpix[rangemin:P1,x])
                    bestbackground2=bestbackground1
                else:
                    # We are somewhere in the middle
                    rangemax=min(h-1,P2+ydim/2)
                    rangemin=max(0,P1-ydim/2)
                    bestbackground1=min(cleanpix[rangemin:P1,x])
                    bestbackground2=min(cleanpix[P2:rangemax,x])
                # Now we have our two points and our best background values for those points,
                # we can linearly interpolate between them
                m=(bestbackground2-bestbackground1)/(P2-P1)
                for P in xrange(P1+1,P2):
                    numpixV[P,x]=m*(P-P1)+bestbackground1
                # Now we can set x to be the next point after P2
                y=P2+1
            y+=1              
        x+=1
        y=0
    if graphicsout:    
        numpix2=numpy.array(numpixV,numpy.uint8)
        pim=Image.fromarray(numpix2,'L')
        pim.save(os.path.join(outputimages,'Vertical.PNG'))
        del pim,numpix2
        gc.collect()
    print("Vertical & Horizontal scanning complete")

    # Now we have the vertically and horizontally filled images, we can average them, and then smooth the result
    # Include numpy.maximum term to artificially eliminate zeros (causes floating point error later on)
    #avpix=numpy.array(numpy.maximum(0.01,(numpixV+numpixH)/2.0),numpy.float)
    avpix=numpy.array(numpy.maximum(0.01,numpy.minimum(numpixV,numpixH)),numpy.float)
    del numpixV,numpixH
    
    #avpix=numpy.array(greycopy,numpy.float) # TEMPORARY LINE
    
    if graphicsout:
        # Let's save that and have a look at it
        avpix2=numpy.array(avpix,numpy.uint8)
        avim=Image.fromarray(avpix2,'L')
        avim.save(os.path.join(outputimages,'MinBack.PNG'))
        del avpix2,avim
        gc.collect()
    # Now, let's create an array of the adjustments which need to be done to each pixel
    adjpix=numpy.array(avpix/medintensity,numpy.float)
    # Make the adjustments, first open the original image
    finalim=im.convert(mode="RGB")
    # Split into channels to reduce probability of memory errors when opening large images
    chans = finalim.split()
    chanlist=[]
    for chan in chans:
        print("Correcting RGB channel...")
        # Correct intensity for each channel
        print("Creating float array version of channel")
        finalpix=numpy.array(chan,numpy.float)
        print("Adjusting for lighting problems")
        #finalpix=finalpix/adjpix
        #finalpix/=adjpix
        finalpix=numpy.divide(finalpix,adjpix)
        print("Converting adjusted float array to unint8")
        finalpix=numpy.array(finalpix,numpy.uint8)
        print("Converting array to greyscale image")
        chan=Image.fromarray(finalpix,'L')
        print("Storing corrected channel")
        chanlist.append(chan)
    # Merge three channels back together
    print("Merging three corrected channels")
    fim=Image.merge("RGB",chanlist)
    del avpix,adjpix,finalpix
    gc.collect()
    # After all that, this is what we want:
    return fim

# Scan along cell edges based on top left hand corner summing the white in the threshold array
# Remember, for these images, (0,0) is the top left hand corner (normally the origin is the bottom left)
def scanedge(x,y):
    bot=sum(thresharr[y,x:x+xdim-1,0])/255
    top=sum(thresharr[y+ydim-1,x:x+xdim-1,0])/255
    left=sum(thresharr[y+1:y+ydim-2,x,0])/255
    right=sum(thresharr[y+1:y+ydim-2,x+xdim-1,0])/255
    return top+bot+left+right

# Try to minimize the number of white pixels on the cell edge by just trying a NxN square of possibilities for the
# top left hand corner of the cell, and picking the coordinates which give the lowest number (or stopping if we find one with zero)
# Should try to make N less than the shortest distance from the edge of the grid to the edge of the image e.g. N<xstart

def spiralsearch(x0,y0,N):
    Direct=1
    sx,sy=x0,y0
    xnew,ynew=x0,y0
    outlist=[[x0,y0]]
    for rad in xrange(1,2*N+1):
        Direct=Direct*-1
        # Move along y axis
        for xnew in xrange(sx+Direct,sx+rad*Direct+Direct,Direct):
            outlist.append([xnew,ynew])
        sx=xnew
        # Move along x axis
        for ynew in xrange(sy+Direct,sy+rad*Direct+Direct,Direct):
            outlist.append([xnew,ynew])
        sy=ynew
    return outlist

def bruteforcesearch(x,y,N):
    coordlist=[]
    reslist=[]
    squarelen=2*xdim+2*ydim
    # Get a list of square spiral points
    sqlist=spiralsearch(x,y,N)
    # Take every xth element to make about 1000 points
    dropsq=max(1,int(round(float(len(sqlist))/1000.0)))
    sqlist=sqlist[0::dropsq]
    # Calculate the white values on the edge for each corner in the search area
    for pt in sqlist:
        [xnew,ynew]=pt
        # Make sure that the full tile is on the image
        xadj=min(max(0,xnew),w-xdim-1)
        yadj=min(max(0,ynew),h-ydim-1)
        whiteedge=scanedge(xadj,yadj)
        # Zero is the best result possible.  If we find it, just go home and be happy
        fracedge=float(whiteedge)/float(squarelen)
        if fracedge==0.0:
            return (xadj,yadj)
        coordlist.append((xadj,yadj))
        reslist.append(whiteedge)
    # Now we get the index of the best result
    best=getMin(reslist)
    # Now we can return the best coordinates
    return coordlist[best]

def tileBorder(bigArray,TLx,TLy,colourno):
    # Pick a colour
    colours=([255,0,0],[0,255,0],[0,0,255],[255,255,0],[0,255,255],[255,0,255])
    colour=colours[colourno]
    # Colour the lines:
    bigArray[TLy,TLx:TLx+xdim]=colour#Top
    bigArray[TLy+ydim,TLx:TLx+xdim]=colour#Bottom
    bigArray[TLy:TLy+ydim,TLx]=colour#Left
    bigArray[TLy:TLy+ydim,TLx+xdim]=colour#Right

#--------------------------------------------------
# This is where the main part of the program starts
#--------------------------------------------------

# Open the output file for writing.  Anything in this already will be destroyed at this stage.
if not singlefile:
    outputfile=open(outfile,'w')

# Now, we can just run through the list of .jpg files that we have:
#for filename in filelist:
while 1:
    # Now, we'll search in this directory for all .jpg files, and put their filenames in filelist
    filelist=[]
    if singlefile==False:
        allfiles=os.listdir(fullpath)
        alldats=os.listdir(os.path.join(fullpath,"Output_Data"))
        for filename in allfiles:
            if filename[-4:] in ('.jpg','.JPG') and filename[0:-4]+"OUT.dat" not in alldats:
                filelist.append(filename)
    else:
        if imagetoscan[-4:] in ('.jpg','.JPG'):
            filelist.append(imagetoscan)

    # If we haven't found any jpg files, we'd better raise an error
    if len(filelist)==0:
        raise Exception("There aren't any .jpg files in this directory which have not been analysed!")
        break
    filename=filelist[0]
    # If we have special instructions for the filename, use those, otherwise use default
    if filename in InsData:
        SetUp(InsData[filename])
    else:
        SetUp(InsData['default'])
    singlef=os.path.join(fullpath,"Output_Data",filename.split('.')[0]+"OUT.dat")
    singleout=open(singlef,'w')

    # Make the filename and the barcode name
    FILENAME=filename[0:-4]
    if len(FILENAME)>11:
        BARCODE=filename[0:11]
    else:
        BARCODE=FILENAME
    # Get a list of all files and directories in the outputimages directory
    #allfolders=os.listdir(outputimages)
    # This is the name of the directory we expect to be putting this image
    outpath=os.path.join(outputimages,BARCODE)
    # If there's no directory called BARCODE, then we need to make it
    #if BARCODE not in allfolders:
    #    os.mkdir(outpath)        
    # Some timing variables
    start=time.time()
    totalstart=start
    # Open the image
    im=Image.open(filename)
    # Get the file dimensions
    w,h=int(im.size[0]),int(im.size[1])
    # TEMPORARY FIX MEMORY ERRORS FOR IMAGES > 10 Megapixels
    if(w*h>10000000) or originalSize==False:
        # Half size of image
        im=im.resize((w/2,h/2),Image.ANTIALIAS)
        w,h=w/2,h/2
        tlx,tly,brx,bry=tlx/2,tly/2,brx/2,bry/2
        xdim,ydim=xdim/2,ydim/2
        xdimf,ydimf=xdimf/2.0,ydimf/2.0
        xstart,ystart=xstart/2,ystart/2
        NSearch=NSearch/2
    
    print filename+' loaded...  '
    # Make up a mask for doing the background correction
    # Do edge detection
    edgeim=sobel(im)
    # SAVE IMAGE
    if graphicsout:
        edgeim.save(os.path.join(outputimages,'Gradient.PNG'))
    invedge=ImageOps.invert(edgeim)
    if graphicsout:
        invedge.save(os.path.join(outputimages,'InvertedGradient.PNG'))

    # Get the gradient histogram, but only for the inner tiles
    x0=max(0,xstart+LightStrip*xdim)
    y0=max(0,ystart+LightStrip*ydim)
    xf=min(w-1,xstart+xdim*(nocols-LightStrip))
    yf=min(h-1,ystart+ydim*(norows-LightStrip))
    whist,hhist=xf-x0,yf-y0

    # Convert image to grayscale
    newim=edgeim.convert('L')
    box = (x0, y0, xf, yf)

    # Crop inner tiles
    region = newim.crop(box)

    # Get histogram
    histlist=region.histogram()
    print "Histogram Ready!"

    # Get the top qq*100% highest gradients
    qq=0.05
    cutind=0
    while sum(histlist[0:cutind])<qq*whist*hhist:
        cutind+=1
    # Interpolate to get best threshold
    cutbest=int(round(cutind-1+(qq*whist*hhist-sum(histlist[0:cutind-1]))/(sum(histlist[0:cutind])-sum(histlist[0:cutind-1]))))
    # Now define a map of the edges
    edgethresh=newim.point(lambda p:(p>cutbest) and 255)
    edgethresh=ImageOps.invert(edgethresh)
    edgearr=numpy.array(edgethresh,numpy.uint8)

    if graphicsout:
        # SAVE IMAGE
        tempim=Image.fromarray(edgearr,'L')
        tempim.save(os.path.join(outputimages,'ThresholdedGradient.PNG'))
        del tempim

    # Now, for every tile, we want to find the darkest edge pixel (where edge pixel is defined by the map above)
    # and threshold so that pixel is classified as foreground
    # First we paste the original image onto a white background using the edge mask
    blank=Image.new("L",(w,h),color=255)
    imcp=im.convert(mode="L")
    blank.paste(imcp,(0,0),edgethresh)
    # Now, let's slowly build the initial threshold map
    # First, let's do a rough thresholding, to have something in place for the parts of the image outside the tiles
    original=im.convert(mode="L")
    newmap=original.point(lambda p:(p>123) and 255)
    edgepix=numpy.array(blank)
    original=im.convert(mode="L")

    if graphicsout:
        # SAVE IMAGE
        blank.save(os.path.join(outputimages,'ThresholdedGradientMask.PNG'))

    # Now let's run through the tiles and threshold each one individually
    for ROW in xrange(norows):
        for COL in xrange(nocols):
            xnow,ynow = xstart+int(round(float(COL)*xdimf)), ystart+int(round(float(ROW)*ydimf))
            # Get the darkest pixel in the tile
            #darkpix=edgepix[ynow:ynow+ydim,xnow:xnow+xdim].min()-1
            darkim=blank.crop((xnow,ynow,xnow+xdim,ynow+ydim))
            # Get histogram for the tile
            histtile=darkim.histogram()
            # Chop off the brightest pixel (mask background)
            histtile=histtile[0:255]
            # Sum up the number of masked pixels
            totpix=sum(histtile)
            if totpix==0:
                cbest=123
            else:
                # Move from darkest pixel towards lightest, until we have counted
                # one third of the pixels
                frac=0.33
                dpix=0
                while dpix<256 and sum(histtile[0:dpix])<frac*totpix:
                    dpix+=1
                # Interpolate to get best threshold
                cbest=int(round(dpix-1+(frac*totpix-sum(histtile[0:dpix-1]))/(sum(histtile[0:dpix])-sum(histtile[0:dpix-1]))))
            # Cut out the tile from the original image
            tileim=original.crop((xnow,ynow,xnow+xdim,ynow+ydim))
            # Threshold it
            tileim=tileim.point(lambda p:(p>cbest) and 255)
            # Now paste it into the map
            newmap.paste(tileim,(xnow,ynow))
    finalim=BackgroundCorrection(im,newmap)

    if graphicsout:
        # SAVE IMAGE
        newmap.save(os.path.join(outputimages,'FirstPassThreshold.PNG'))

    cutoff=automaticSetThreshold(finalim)
    finish=time.time()
    print 'Thresholding complete '+str(int(round(finish-start)))+'sec  Building arrays...'
    start=time.time()

    pixarr=numpy.array(finalim,numpy.uint8)
    origarr=numpy.array(im,numpy.uint8)
    greyim=finalim.convert(mode="L")

    # These are kept as 8-bit so that we can have colour later
    threshim1chan=greyim.point(lambda p:(p>cutoff) and 255)
    # We also need to check if the thresholding has done something stupid
    # occasionally it does since the genetic algorithm is stochastic
    # Let's just find the fraction of white pixels in the image
    tarr=numpy.array(threshim1chan,numpy.uint8)
    sumtotpix=0
    for x in xrange(w):
        sumtotpix+=sum(tarr[0:,x])/255
    fracpix=float(sumtotpix)/float(w*h)
    # If this fraction is ridiculously small (black) or ridiculously large (white)
    # then let's try the thresholding again!
    mutcount=0
    while (fracpix<0.05 or fracpix>0.80) and mutcount<10:
        print "Analyzing again....",mutcount
        cutoff=automaticSetThreshold(finalim)
        pixarr=numpy.array(finalim,numpy.uint8)
        origarr=numpy.array(im,numpy.uint8)
        greyim=finalim.convert(mode="L")
        threshim1chan=greyim.point(lambda p:(p>cutoff) and 255)
        tarr=numpy.array(threshim1chan,numpy.uint8)
        sumtotpix=0
        for x in xrange(w):
                sumtotpix+=sum(tarr[0:,x])/255
        fracpix=float(sumtotpix)/float(w*h)    
        mutcount+=1
    thresharr1chan=numpy.array(threshim1chan,numpy.uint8)
    thresharr=numpy.zeros((h,w,3), numpy.uint8)
    thresharr[0:,0:,0]=thresharr1chan
    thresharr[0:,0:,1]=thresharr1chan
    thresharr[0:,0:,2]=thresharr1chan    

    greyarr1chan=numpy.array(greyim,numpy.uint8)
    greyarr=numpy.zeros((h,w,3), numpy.uint8)
    greyarr[0:,0:,0]=greyarr1chan
    greyarr[0:,0:,1]=greyarr1chan
    greyarr[0:,0:,2]=greyarr1chan

    tilearr=numpy.zeros((ydim,xdim,3),numpy.uint8)
    finish=time.time()            
    print 'Arrays built...'+str(int(round(finish-start)))+'sec  Performing Analysis and Saving Tiles...'
    start=time.time()

    if graphicsout:
        # Testing to see if I can explain Darren's observation of different levels of growth
        Untrimmed=numpy.zeros((norows,nocols),numpy.float)
        Trimmed=numpy.zeros((norows,nocols),numpy.float)
        # Expect the untrimmed area to be the most reliable
        FracDiff=numpy.zeros((norows,nocols),numpy.float)
    # We want to calculate the area of the culture in each cell
    # Looking at the output.dat files that are in the robot folders I guess that they need to be in this format:
    # FILENAME ROW COLUMN TOPLEFTX TOPLEFTY WHITEAREA(px) TRIMMED THRESHOLD INTENSITY EDGEPIXELS COLR COLG COLB BKR BKG BKB EDGELEN XDIM YDIM
    tilecoords=[]

    for ROW in xrange(norows):
        for COL in xrange(nocols):
	    # Initialise counters for this cell
            # Declaring the types (int) to stop these counters getting typecast as uint later (problems with -1=255)
            WHITEAREA=int(0)
            TOTALINTENSITY=int(0)
            COLR,COLG,COLB=int(0),int(0),int(0)
            BKR,BKG,BKB=int(0),int(0),int(0)
            EDGE=int(0)
            IPWA=int(0)
	    # Make a first guess at the coordinates of the cell's top left hand corner
            xnow,ynow = min(w-xdim-1,xstart+int(round(float(COL)*xdimf))), min(h-ydim-1,ystart+int(round(float(ROW)*ydimf)))
            # If there are white pixels on the edge of the cell, then search for a better pair
            # But only do this if we are away from the edges of the plate (since these are controls, and are blown out by lighting issues)
            if ROW>LightStrip-1 and COL>LightStrip-1 and ROW<norows-LightStrip and COL<nocols-LightStrip and scanedge(xnow,ynow)!=0:
                # Here we are scanning over a square, looking for an improvement
                if colonySearch:
                    (xnow,ynow)=bruteforcesearch(xnow,ynow,NSearch)
            # Now that we have a square that contains the colony, minimizing the number of white pixels on the edge,
            # it might be nice to try to centre the image on the colony (without increasing the no. of edge pixels).
            edgetest=scanedge(xnow,ynow)
            xmid,ymid,npoints=0,0,0
            # Add up the x and y coordinates (and count the number of white pixels inside the cell)
            for x in xrange (xnow,xnow+xdim):
                for y in xrange (ynow,ynow+ydim):
                    if thresharr[y,x][0]>128:
                        xmid+=x
                        ymid+=y
                        npoints+=1
            # Calculate "centre of mass" if there are white pixels in the cell (if we're away from the edge of the plate)
            if npoints>0 and ROW>LightStrip-1 and COL>LightStrip-1 and ROW<norows-LightStrip and COL<nocols-LightStrip:
                xmid=int(round(float(xmid)/float(npoints)))-xdim/2
                ymid=int(round(float(ymid)/float(npoints)))-ydim/2
                # If the edge test is unaffected (increased by less than 2% of the max edge length), let's use the centre of mass instead
                if colonySearch and xmid>=0 and xmid<w-xdim and ymid>=0 and ymid<h-ydim and scanedge(xmid,ymid)-edgetest<=int(round(2*(xdim+ydim)*0.02)):
                    xnow,ynow=xmid,ymid                        

            for x in xrange(xnow,xnow+xdim):
                for y in xrange(ynow,ynow+ydim):
                    # Fill the tile array for outputting later
                    tilearr[y-ynow,x-xnow]=origarr[y,x]
                    # Add total intensities
                    TOTALINTENSITY+=int(round(float(greyarr[y,x][0])-medintensity))
                    if edgearr[y,x]>128:
                        EDGE+=1
                    if thresharr[y,x][0]>128:
                        # Add white pixels
                        WHITEAREA+=1
                        # For every white pixel, add the pixel intensity from the greyscale image
                        # Subtract the background value to account for small (thin) colonies
                        # which are not completely white:
                        IPWA+=int(round(float(greyarr[y,x][0])-medintensity))
                        # Store the colours from the original image
                        COLR+=origarr[y,x][0]
                        COLG+=origarr[y,x][1]
                        COLB+=origarr[y,x][2]
                    else:
                        # Store the colours from the backround (original image)
                        BKR+=origarr[y,x][0]
                        BKG+=origarr[y,x][1]
                        BKB+=origarr[y,x][2]
            # Average the colours for storing
            if WHITEAREA>0:
                COLR=int(round(float(COLR)/float(WHITEAREA)))
                COLG=int(round(float(COLG)/float(WHITEAREA)))
                COLB=int(round(float(COLB)/float(WHITEAREA)))
            else:
                COLR=""
                COLG=""
                COLB=""
            if xdim*ydim-WHITEAREA>0:
                BKR=int(round(float(BKR)/float(xdim*ydim-WHITEAREA)))
                BKG=int(round(float(BKG)/float(xdim*ydim-WHITEAREA)))
                BKB=int(round(float(BKB)/float(xdim*ydim-WHITEAREA)))
            else:
                BKR=""
                BKG=""
                BKB=""          
            # Now we can write this stuff to the output file
            # Remember: ROW goes from top to bottom (because of the coordinate system)
            outstr=FILENAME+"\t"+str(ROW+1)+"\t"+str(COL+1)+"\t"+str(xnow)+"\t"+str(ynow)+"\t"+str(WHITEAREA)+"\t"+str(IPWA)+"\t"
            outstr+=str(cutoff)+"\t"+str(TOTALINTENSITY)+"\t"+str(scanedge(xnow,ynow))+"\t"+str(COLR)+"\t"+str(COLG)+"\t"
            outstr+=str(COLB)+"\t"+str(BKR)+"\t"+str(BKG)+"\t"+str(BKB)+"\t"+str(EDGE)+"\t"+str(xdim)+"\t"+str(ydim)+"\n"
	    #outstr+=str(int(round(medintensity)))+"\n"    
            if not singlefile:
                outputfile.write(outstr)
                singleout.write(outstr)
            else:
                singleout.write(outstr)
            # Save the current tile
            #tilefile=os.path.join(outpath,FILENAME+'R%02dC%02d.png'%(ROW+1,COL+1))
            #pilImage=Image.fromarray(tilearr,'RGB')
            #pilImage.save(tilefile)
            # Record the top left hand corner for the preview
            tilecoords.append((xnow,ynow))
            if graphicsout:
                # Save some data to arrays for testing Darren's growth gradient
                Trimmed[ROW,COL]=IPWA
                Untrimmed[ROW,COL]=TOTALINTENSITY
        # Just to let us know how we're getting on
        print "Row: "+str(ROW+1)
    finish=time.time()
    print 'Analysis Finished and Tiles Saved...'+str(int(round(finish-start)))+'sec  Building Previews...'
    start=time.time()

    if graphicsout:
        pixdim=50
        FracMax=0.1
        GrowthImArr=numpy.zeros((norows*pixdim,nocols*pixdim,3), numpy.uint8)
        # For Darren's growth gradient, let's build a difference array
        for ROW in xrange(norows):
            for COL in xrange(nocols):
                Frac=(Untrimmed[ROW,COL]-Trimmed[ROW,COL])/Untrimmed[ROW,COL]
                LookupColour=int(round((Frac/FracMax)*255.0))
                FracDiff[ROW,COL]=Frac
                GrowthImArr[ROW*pixdim:(ROW+1)*pixdim,COL*pixdim:(COL+1)*pixdim]=[255,0,LookupColour]
        # Save the output file for this analysis
        GrowthIm=Image.fromarray(GrowthImArr,'RGB')
        GrowthIm.save(os.path.join(outputimages,FILENAME+'DarrenGrowth.png'))
        del GrowthIm,Trimmed,Untrimmed,FracDiff

        # How to print midslices of intensity from the image
        sl1=open(os.path.join(outputimages,'SLICE1.dat'),'w')
        bstring=""
        for x in greyarr[1152,0:,0]:
            bstring+=str(x)+"\n"
        sl1.write(bstring[0:-1])
        sl1.close()
        sl2=open(os.path.join(outputimages,'SLICE2.dat'),'w')
        bstring=""
        for x in greyarr[0:,1536,0]:
            bstring+=str(x)+"\n"
        sl2.write(bstring[0:-1])
        sl2.close()
        del GrowthIm,GrowthImArr
        gc.collect()

    if graphicsout:
        # Finally write preview images in lossless .png format
        pilCorr=Image.fromarray(outputimages,'RGB')
        pilCorr.save(os.path.join(outputimages,FILENAME +'OriginalCorrected.png'))
        im.save(os.path.join(outputimages,FILENAME +'Original.png'))

    # Let's draw the tile borders on the images
    tileno=0
    for coord in tilecoords:
        tileBorder(pixarr,coord[0],coord[1],tileno%5)
        tileBorder(greyarr,coord[0],coord[1],tileno%5)
        tileBorder(thresharr,coord[0],coord[1],tileno%5)
        tileno+=1

    if graphicsout:
        pilGridCorr=Image.fromarray(pixarr,'RGB')
        pilGridCorr.save(os.path.join(outputimages,FILENAME +'GriddedCorrected.png'))

    pilThresh=Image.fromarray(thresharr,'RGB')
    pilThresh.save(os.path.join(outputimages,FILENAME +'GriddedThresh.png'))

##    # Now let's make some small previews
##    smw=200
##    smh=int(round(float(h)*float(smw)/float(w)))
##    smallCorr=pilGridCorr.resize((smw,smh),Image.ANTIALIAS)
##    smallOriginal=im.resize((smw,smh),Image.ANTIALIAS)
##    smallThresh=pilThresh.resize((smw,smh),Image.ANTIALIAS)
##    smallCorr.save(os.path.join(outputimages,FILENAME +'CorrectedPreview.jpeg'),quality=100)
##    smallOriginal.save(os.path.join(outputimages,FILENAME +'OriginalPreview.jpeg'),quality=100)
##    smallThresh.save(os.path.join(outputimages,FILENAME +'ThreshPreview.jpeg'),quality=100)

    finish=time.time()
    print 'Previews created...'+str(int(round(finish-start)))+'sec  Total time: '+str(int(round(finish-totalstart)))+' sec'
    singleout.close()
    # Let's clear some memory for the next round
    del im,edgeim,invedge,newim,region,edgethresh,edgearr,blank,imcp,original
    del newmap,edgepix,darkim,tileim,finalim,pixarr,origarr,greyim,threshim1chan
    del thresharr,greyarr1chan,greyarr,tilearr,pilThresh
    #del smallCorr,smallOriginal,smallThresh,pilImage,pilCorr,pilGridCorr
    gc.collect()

# Now we can close the output files
if not singlefile:
    outputfile.close()  
