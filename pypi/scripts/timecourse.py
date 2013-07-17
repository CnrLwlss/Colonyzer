from colonyzer2 import *
import time, sys

def main(fmt="384"):
    # Lydall lab file naming convention
    # First 15 characters in filename identify unique plates
    # Remaining charaters can be used to store date, time etc.
    barcRange=(0,15)
    if len(sys.argv)>1:
        fmt=sys.argv[1]
    
    # Format names and dimension definitions
    formats=["48","96","117","384","768","1536"]
    dims=[(6,8),(12,8),(13,9),(24,16),(48,32),(48,32)]
    nx,ny=dims[formats.index(fmt)] 

    start=time.time()

    # Find image files which have yet to be analysed
    (fullpath,outputimages,outputdata)=setupDirectories()
    barcdict=getBarcodes(outputimages,outputdata,fullpath,barcRange)

    while len(barcdict)>0:
        BARCODE=barcdict.keys()[0]
        print(BARCODE)
        LATESTIMAGE=barcdict[BARCODE][0]
        EARLIESTIMAGE=barcdict[BARCODE][-1]
        imRoot=EARLIESTIMAGE.split(".")[0]
        
        # Indicate that barcode is currently being analysed, to allow parallel analysis
        tmp=open(os.path.join(outputdata,imRoot+".dat"),"w").close()

        # Get latest image for thresholding and detecting culture locations
        imN,arrN=openImage(LATESTIMAGE)
        # Get earliest image for lighting gradient correction
        im0,arr0=openImage(EARLIESTIMAGE)

        # Automatically generate guesses for gridded array locations
        diam=int(1.05*round(min(float(arrN.shape[0])/(ny+1),float(arrN.shape[1])/(nx+1))))
        (candx,candy,dx,dy)=estimateLocations(arrN,nx,ny,diam,showPlt=False)

        # Update guesses and initialise locations data frame
        locationsN=locateCultures(candx,candy,dx,dy,arrN)

        # Smooth (pseudo-)empty image 
        (correction_map,average_back)=makeCorrectionMap(arr0,locationsN)

        # Correct spatial gradient in final image
        corrected_arrN=arrN*correction_map

        # Trim outer part of image to remove plate walls
        trimmed_arr=arrN[max(0,min(locationsN.y)-dy):min(arr0.shape[0],(max(locationsN.y)+dy)),max(0,(min(locationsN.x)-dx)):min(arr0.shape[1],(max(locationsN.x)+dx))]
        (thresh,bindat)=automaticThreshold(trimmed_arr)

        # Mask for identifying culture areas
        maskN=numpy.ones(arrN.shape,dtype=numpy.bool)
        maskN[arrN<thresh]=False

        for FILENAME in barcdict[BARCODE]:
            im,arr=openImage(FILENAME)
            arr=arr*correction_map
            
            # Correct for lighting differences between plates
            arrsm=arr[numpy.min(locationsN.y):numpy.max(locationsN.y),numpy.min(locationsN.x):numpy.max(locationsN.x)]
            masksm=maskN[numpy.min(locationsN.y):numpy.max(locationsN.y),numpy.min(locationsN.x):numpy.max(locationsN.x)]
            meanPx=numpy.mean(arrsm[numpy.logical_not(masksm)])

            arr=arr+(average_back-meanPx)
            threshadj=thresh+(average_back-meanPx)

            mask=numpy.ones(arr.shape,dtype=numpy.bool)
            mask[arrN<threshadj]=False

            # Measure culture phenotypes
            locations=measureSizeAndColour(locationsN,arr,im,mask,average_back,BARCODE,FILENAME[0:-4])

            # Write results to file
            locations.to_csv(os.path.join(outputdata,FILENAME[0:-4]+".out"),"\t",index=False,engine='python')
            dataf=saveColonyzer(os.path.join(outputdata,FILENAME[0:-4]+".dat"),locations,threshadj,dx,dy)

            # Visual check of culture locations
            imthresh=threshPreview(arr,threshadj,locations)
            imthresh.save(os.path.join(outputimages,FILENAME[0:-4]+".png"))

        # Get ready for next image
        print("Finished: "+FILENAME+" "+str(time.time()-start)+" s")
        barcdict=getBarcodes(outputimages,outputdata,fullpath,barcRange)
    print("No more barcodes to analyse... I'm done.")

if __name__ == '__main__':
    main()