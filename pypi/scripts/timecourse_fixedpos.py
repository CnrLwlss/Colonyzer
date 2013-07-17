from colonyzer2 import *
import time, sys, os, numpy, PIL

def main():
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Note that this script requires a Colonyzer.txt file (as generated by ColonyzerParametryzer) describing initial guess for culture array")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    # Lydall lab file naming convention
    # First 15 characters in filename identify unique plates
    # Remaining charaters can be used to store date, time etc.
    barcRange=(0,15)
    correction=True
    if len(sys.argv)>1:
        if '--nolc' in sys.argv:
            print "No lighting correction..."
            correction=False

    start=time.time()

    # Find image files which have yet to be analysed
    (fullpath,outputimages,outputdata)=setupDirectories()
    barcdict=getBarcodes(outputimages,outputdata,fullpath,barcRange)
    InsData=readInstructions(fullpath)

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

        # If we have ColonyzerParametryzer output for this filename, use it for initial culture location estimates
        if LATESTIMAGE in InsData:
            (candx,candy,dx,dy)=SetUp(InsData[LATESTIMAGE])
        else:
            (candx,candy,dx,dy)=SetUp(InsData['default'])

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
            if correction:
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
            locations.to_csv(os.path.join(outputdata,FILENAME[0:-4]+".out"),"\t",index=False)
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
