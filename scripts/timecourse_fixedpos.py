import colonyzer2 as c2
import json
import argparse
import shutil
import string
import os
import time
import numpy
from matplotlib.backends.backend_pdf import PdfPages

def parseArgs(dummyRun=False,dummyInp='-d testdir -n -t'):
    '''Define console script behaviour, hints and documentation for setting off Colonyzer analysis.'''
    parser=argparse.ArgumentParser(description="Analyse timeseries of QFA images: locate cultures on plate, segment image into agar and cells, apply lighting correction, write report including cell density estimates for each location in each image.")
    parser.add_argument("-d","--dir", type=str, help="Directory in which to search for image files that have not been analysed (current directory by default).")
    parser.add_argument("-n","--nolc", help="Disable lighting correction?", action="store_true")
    parser.add_argument("-f","--fixthresh", type=float, help="Value to fix image segmentation threshold at (default is automatic thresholding).")
    parser.add_argument("-t","--threshplots", help="Plot pixel intensity distributions and segmentation thresholds?", action="store_true")
    parser.add_argument("-u","--usedict", type=str, help="Load .json file specifying images to analyse from HTS-style directory structure.  See C2Find.py in HTSauto package.")
    if dummyRun:
        args = parser.parse_args(dummyInp.split())
    else:
        args = parser.parse_args()
    return(args)

def checkImages(fdir,fdict=None,barcRange=(0,-24),verbose=False):
    '''Discover barcodes in current working directory (or in fdir or in those specified in fdict) for which analysis has not started.'''
    if fdict!=None:
        with open(fdict, 'rb') as fp:
            barcdict = json.load(fp)
            # Drop any barcodes that are currently being analysed/already analysed
            barcdict={x:barcdict[x] for x in barcdict.keys() if not c2.checkAnalysisStarted(barcdict[x][-1])}
    else:
        # Find image files which have yet to be analysed
        # Lydall lab file naming convention (barcRange)
        # First 15 characters in filename identify unique plates
        # Remaining charaters can be used to store date, time etc.
        barcdict=c2.getBarcodes(fdir,barcRange,verbose=verbose)
    return(barcdict)

def buildVars(dummyRun=False,dummyInp='-d testdir -n -t',verbose=False):
    '''Read user input, set up flags for analysis, report on options chosen and find files to be analysed.'''
    inp=parseArgs(dummyRun,dummyInp)

    if inp.dir==None:
        fdir=os.getcwd()
    else:
        fdir=inp.dir
    
    correction = not inp.nolc
    if verbose and not correction: print "Lighting correction will be disabled."
    
    if inp.fixthresh!=None:
        fixedThresh=inp.fixthresh
        if verbose: print "Segmentation threshold will be fixed at "+str(fixedThresh)
    else:
        fixedThresh=0

    threshplots=inp.threshplots
    if verbose and threshplots: print "Pixel intensity distributions will be plotted."
    
    expt=inp.usedict
    if expt==None:
        fdict=None
    else:
        exptType=expt[0:-4]
        fdict=os.path.join(rootDir,exptType+"_EXPERIMENTS",expt,"AUXILIARY",expt+"_C2.json")
        if verbose: print "Preparing to load barcodes for experiment "+expt+" from "+fdict+", assuming HTS directory structure."
    return((correction,fixedThresh,threshplots,expt,fdict,fdir))

def prepareTimecourse(barcdict,verbose=False):
    '''In timecourse mode, prepares "next" batch of images for analysis from dictionary of image names (unique image barcodes are dictionary keys).'''
    BARCs=barcdict.keys()
    BARCs.sort()
    BARCODE=BARCs[0]
    imdir=os.path.dirname(barcdict[BARCODE][0])
    InsData=c2.readInstructions(imdir)
    IMs=barcdict[BARCODE]
    LATESTIMAGE=IMs[0]
    EARLIESTIMAGE=IMs[-1]
    imRoot=EARLIESTIMAGE.split(".")[0]
    if verbose:
        print("Analysing images labelled with the barcode "+BARCODE+" in "+imdir)
        print("Earliest image: "+EARLIESTIMAGE)
        print("Latest image: "+LATESTIMAGE)
    return((BARCODE,imdir,InsData,LATESTIMAGE,EARLIESTIMAGE,imRoot))

def locationGuesses(IMAGE,InsData):
    '''Set up initial guesses for location of spots on image by parsing data from Colonyzer.txt
    TODO: What about situtation when there is no Colonyzer.txt?  i.e. InsData=None?'''
    # If we have ColonyzerParametryzer output for this filename, use it for initial culture location estimates
    if os.path.basename(IMAGE) in InsData:
        (candx,candy,dx,dy)=c2.SetUp(InsData[os.path.basename(IMAGE)])
    # If there are multiple calibrations available, choose the best one based on date of image capture
    elif any(isinstance(el, list) for el in InsData['default']):
        imname=os.path.basename(IMAGE).split(".")[0]
        imdate=imname[-19:-9]
        (candx,candy,dx,dy)=c2.SetUp(InsData['default'],imdate)
    else:
        (candx,candy,dx,dy)=c2.SetUp(InsData['default'])
    return((candx,candy,dx,dy))

def cutEdgesFromMask(mask,locations):
    '''Mask for identifying culture areas (edge detection). Set all pixels outside culture grid to background, to aid binary filling later.'''
    mask[0:min(locations.y-dy/2),:]=False
    mask[max(locations.y+dy/2):mask.shape[0],:]=False
    mask[:,0:min(locations.x-dx/2)]=False
    mask[:,max(locations.x+dx/2):mask.shape[1]]=False
    return(mask)

def edgeFill(arr,locations,cutoff=0.8):
    edgeN=getEdges(arr,cutoff=0.8)
    dilN=ndimage.morphology.binary_dilation(edgeN,iterations=2)
    erodeN=ndimage.morphology.binary_erosion(dilN,iterations=1)
    dil2N=ndimage.morphology.binary_dilation(dilN,iterations=3)
    
    fillN=ndimage.morphology.binary_fill_holes(cutEdgesFromMask(dil2N,locations))
    maskN=ndimage.morphology.binary_erosion(fillN,iterations=7)
    return(maskN)


def main(verbose=False):
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Note that this script requires a Colonyzer.txt file (as generated by ColonyzerParametryzer) describing initial guess for culture array")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    correction,fixedThresh,threshplots,expt,fdict,fdir=buildVars(verbose=verbose,dummyRun=True,dummyInp='-d ../Auxiliary/Data -n -t')
    barcdict=checkImages(fdir,fdict,verbose=verbose)
    rept=c2.setupDirectories(barcdict,verbose=verbose)

    start=time.time()

    while len(barcdict)>0:
        BARCODE,imdir,InsData,LATESTIMAGE,EARLIESTIMAGE,imRoot=prepareTimecourse(barcdict,verbose=True)  
        
        # Create empty file to indicate that barcode is currently being analysed, to allow parallel analysis (lock files)
        tmp=open(os.path.join(os.path.dirname(EARLIESTIMAGE),"Output_Data",os.path.basename(EARLIESTIMAGE).split(".")[0]+".out"),"w").close()

        # Get latest image for thresholding and detecting culture locations
        imN,arrN=c2.openImage(LATESTIMAGE)
        # Get earliest image for lighting gradient correction
        im0,arr0=c2.openImage(EARLIESTIMAGE)

        (candx,candy,dx,dy)=locationGuesses(LATESTIMAGE,InsData)

        # Update guesses and initialise locations data frame
        locationsN=c2.locateCultures([int(round(cx-dx/2.0)) for cx in candx],[int(round(cy-dy/2.0)) for cy in candy],dx,dy,arrN)

        cutFromFirst=False
        cythonFill=False
        
        if cutFromFirst:
            mask=edgeFill(arr0,locationsN,0.8)
            startFill=time.time()
            if cythonFill:
                pseudoempty=maskAndFillCython(arr0,maskN,0.005)
                print("Inpainting using Cython & NumPy: "+str(time.time()-startFill)+" s")
            else:
                pseudoempty=maskAndFill(arr0,mask,0.005)
                print("Inpainting using NumPy: "+str(time.time()-start)+" s")
        else:
            pseudoempty=arr0
            
        # Smooth (pseudo-)empty image 
        (correction_map,average_back)=c2.makeCorrectionMap(pseudoempty,locationsN,printMess=correction)

        # Correct spatial gradient in final image
        corrected_arrN=arrN*correction_map

        # Trim outer part of image to remove plate walls
        trimmed_arrN=arrN[max(0,int(round(min(locationsN.y)-dy/2.0))):min(arrN.shape[0],int(round((max(locationsN.y)+dy/2.0)))),max(0,int(round(min(locationsN.x)-dx/2.0))):min(arrN.shape[1],int(round((max(locationsN.x)+dx/2.0))))]
        
        if fixedThresh!=0:
            thresh=fixedThresh
        else:
            if threshplots:
                pdf=PdfPages(BARCODE+'_HistogramReport.pdf')
                (thresh,bindat)=c2.automaticThreshold(trimmed_arrN,BARCODE,pdf)
                c2.plotModel(bindat,label=BARCODE,pdf=pdf)
                pdf.close()
            else:
                (thresh,bindat)=c2.automaticThreshold(trimmed_arrN)

        # Mask for identifying culture areas
        maskN=numpy.ones(arrN.shape,dtype=numpy.bool)
        maskN[arrN<thresh]=False

        for FILENAME in barcdict[BARCODE]:
            im,arr=c2.openImage(FILENAME)
            if correction:
                arr=arr*correction_map
            
            # Correct for lighting differences between plates
            arrsm=arr[max(0,int(round(min(locationsN.y)-dy/2.0))):min(arrN.shape[0],int(round((max(locationsN.y)+dy/2.0)))),max(0,int(round(min(locationsN.x)-dx/2.0))):min(arrN.shape[1],int(round((max(locationsN.x)+dx/2.0))))]
            masksm=maskN[max(0,int(round(min(locationsN.y)-dy/2.0))):min(arrN.shape[0],int(round((max(locationsN.y)+dy/2.0)))),max(0,int(round(min(locationsN.x)-dx/2.0))):min(arrN.shape[1],int(round((max(locationsN.x)+dx/2.0))))]
            meanPx=numpy.mean(arrsm[numpy.logical_not(masksm)])

            #arr=arr+(average_back-meanPx)
            #threshadj=thresh+(average_back-meanPx)
            threshadj=thresh

            mask=numpy.ones(arr.shape,dtype=numpy.bool)
            mask[arrN<threshadj]=False

            # Measure culture phenotypes
            locations=c2.measureSizeAndColour(locationsN,arr,im,mask,average_back,BARCODE,FILENAME[0:-4])

            # Write results to file
            locations.to_csv(os.path.join(os.path.dirname(FILENAME),"Output_Data",os.path.basename(FILENAME).split(".")[0]+".out"),"\t",index=False,engine='python')
            dataf=c2.saveColonyzer(os.path.join(os.path.dirname(FILENAME),"Output_Data",os.path.basename(FILENAME).split(".")[0]+".dat"),locations,threshadj,dx,dy)

            # Visual check of culture locations
            imthresh=c2.threshPreview(arr,threshadj,locations)
            imthresh.save(os.path.join(os.path.dirname(FILENAME),"Output_Images",os.path.basename(FILENAME).split(".")[0]+".png"))

        # Get ready for next image
        print("Finished: "+FILENAME+" "+str(time.time()-start)+" s")

        barcdict={x:barcdict[x] for x in barcdict.keys() if not c2.checkAnalysisStarted(barcdict[x][-1])}

    print("No more barcodes to analyse... I'm done.")

if __name__ == '__main__':
    main(True)
