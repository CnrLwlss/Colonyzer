from colonyzer2.functions import *
import time, sys

def main(fmt="384"):
    # Plate format names and dimension definitions
    formats=["48","96","117","384","768","1536"]
    dims=[(6,8),(12,8),(13,9),(24,16),(48,32),(48,32)]
    
    # Lydall lab file naming convention
    # First 15 characters in filename identify unique plates
    # Remaining charaters can be used to store date, time etc.
    barcRange=(0,15)

    # Parse arguments... (use library for this...)
    correction=True
    if len(sys.argv)>1:
        # Disabling lighting correction
        if '--nolc' in sys.argv:
            print "No lighting correction..."
            correction=False
        # Non-default plate formats
        for fmat in formats:
            if fmat in sys.argv:
                fmt=fmat    

    nx,ny=dims[formats.index(fmt)] 

    start=time.time()

    # Find image files which have yet to be analysed
    barcdict=getBarcodes(os.getcwd(),barcRange)
    # Setup output directories if not already present
    rept=setupDirectories(barcdict)
    if len(rept)>0:
        print ("Newly created directories:")
        for line in rept:
            print rept

    while len(barcdict)>0:
        BARCODE=barcdict.keys()[0]
        print(BARCODE)
        LATESTIMAGE=barcdict[BARCODE][0]
        EARLIESTIMAGE=barcdict[BARCODE][-1]
        imRoot=os.path.basename(EARLIESTIMAGE).split(".")[0]
        
        # Indicate that barcode is currently being analysed, to allow parallel analysis
        tmp=open(os.path.join(os.path.dirname(EARLIESTIMAGE),"Output_Data",os.path.basename(EARLIESTIMAGE).split(".")[0]+".out"),"w").close()

        # Get latest image for thresholding and detecting culture locations
        imN,arrN=openImage(LATESTIMAGE)
        # Get earliest image for lighting gradient correction
        im0,arr0=openImage(EARLIESTIMAGE)

        # Automatically generate guesses for gridded array locations
        diam=int(1.05*round(min(float(arrN.shape[0])/(ny+1),float(arrN.shape[1])/(nx+1))))
        (candx,candy,dx,dy)=estimateLocations(arrN,nx,ny,diam,showPlt=False)

        # Update guesses and initialise locations data frame
        locationsN=locateCultures([int(round(cx-dx/2.0)) for cx in candx],[int(round(cy-dy/2.0)) for cy in candy],dx,dy,arrN)

        # Smooth (pseudo-)empty image 
        (correction_map,average_back)=makeCorrectionMap(arr0,locationsN,printMess=correction)

        # Correct spatial gradient in final image
        corrected_arrN=arrN*correction_map

        # Trim outer part of image to remove plate walls
        trimmed_arr=arrN[max(0,int(round(min(locationsN.y)-dy/2.0))):min(arrN.shape[0],int(round((max(locationsN.y)+dy/2.0)))),max(0,int(round(min(locationsN.x)-dx/2.0))):min(arrN.shape[1],int(round((max(locationsN.x)+dx/2.0))))]
        (thresh,bindat)=automaticThreshold(trimmed_arr)

        # Mask for identifying culture areas
        maskN=numpy.ones(arrN.shape,dtype=numpy.bool)
        maskN[arrN<thresh]=False

        for FILENAME in barcdict[BARCODE]:
            im,arr=openImage(FILENAME)
            if correction:
                arr=arr*correction_map
            
            # Correct for lighting differences between plates
            arrsm=arr[max(0,int(round(min(locationsN.y)-dy/2.0))):min(arrN.shape[0],int(round((max(locationsN.y)+dy/2.0)))),max(0,int(round(min(locationsN.x)-dx/2.0))):min(arrN.shape[1],int(round((max(locationsN.x)+dx/2.0))))]
            masksm=maskN[max(0,int(round(min(locationsN.y)-dy/2.0))):min(arrN.shape[0],int(round((max(locationsN.y)+dy/2.0)))),max(0,int(round(min(locationsN.x)-dx/2.0))):min(arrN.shape[1],int(round((max(locationsN.x)+dx/2.0))))]
            meanPx=numpy.mean(arrsm[numpy.logical_not(masksm)])

            arr=arr+(average_back-meanPx)
            threshadj=thresh+(average_back-meanPx)

            mask=numpy.ones(arr.shape,dtype=numpy.bool)
            mask[arrN<threshadj]=False

            # Measure culture phenotypes
            locations=measureSizeAndColour(locationsN,arr,im,mask,average_back,BARCODE,FILENAME[0:-4])

            # Write results to file
            locations.to_csv(os.path.join(os.path.dirname(FILENAME),"Output_Data",os.path.basename(FILENAME).split(".")[0]+".out"),"\t",index=False,engine='python')
            dataf=saveColonyzer(os.path.join(os.path.dirname(FILENAME),"Output_Data",os.path.basename(FILENAME).split(".")[0]+".dat"),locations,threshadj,dx,dy)

            # Visual check of culture locations
            imthresh=threshPreview(arr,threshadj,locations)
            imthresh.save(os.path.join(os.path.dirname(FILENAME),"Output_Images",os.path.basename(FILENAME).split(".")[0]+".png"))

        # Get ready for next image
        print("Finished: "+FILENAME+" "+str(time.time()-start)+" s")
        barcdict=getBarcodes(os.getcwd(),barcRange)
        # Setup output directories if not already present
        rept=setupDirectories(barcdict)
    print("No more barcodes to analyse... I'm done.")

if __name__ == '__main__':
    main()
