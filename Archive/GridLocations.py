import sys, os, pandas, StringIO, numpy, math, pygame
from PIL import Image, ImageDraw, ImageFilter, ImageOps
from pygame.locals import *
#import matplotlib.pyplot as plt

def getSpot(block,row,col,dat):
    '''Pulls data for an individual feature from the pandas version of .gal file data frame (dat)'''
    return(dat[(dat.Block==block)&(dat.Row==row)&(dat.Column==col)])

def parseGal(filepath):
    '''Parses a .gal microarray description file, returning a dictionary containing file contents, particularly
    the dataFrame element, which is the list of names of objects at each microarray location'''
    # Parsing .gal microarray files
    gal=open(filepath,'r')
    galLines=gal.readlines()
    galDict={}
    # WTF is ATF?
    tmp=galLines[0].rstrip().split("\t")
    galDict[tmp[0]]=float(tmp[1])
    # Skipping 2nd line, not sure what these integers mean
    for i in xrange(2,9):
        tmp=galLines[i].rstrip().split("=")
        galDict[tmp[0]]=tmp[1]
    galDict["BlockCount"]=int(galDict["BlockCount"])
    # Get block coordinates
    for i in xrange(9,9+galDict["BlockCount"]):
        tmp=galLines[i].replace('\"','').replace(' ','').rstrip().split("=")
        lst=[int(l) for l in tmp[1].split(",")]
        galDict[tmp[0]]=lst
    # Read in block, row, column data
    stream=StringIO.StringIO("".join(galLines[9+galDict["BlockCount"]:]))
    spots=pandas.read_table(stream)
    galDict["spotIDs"]=spots
    return(galDict)

def predictLocations(galDict,micronsPerPixel=5):
    '''Assuming 8 chips arranged 4 wide 2 high, numbered left to right then top to bottom.
    Returns pandas data frame with estimate for relative feature locations (px), based on .gal file instructions.'''
    dat=[]
    b=0
    for crow in xrange(0,2):
        for ccol in xrange(0,4):
            chipDat=galDict["Block"+str(b+1)]
            tlx,tly=float(chipDat[0])/float(micronsPerPixel),float(chipDat[1])/float(micronsPerPixel)
            diam=float(chipDat[2])/float(micronsPerPixel)
            ny,nx=chipDat[5],chipDat[3]
            dy=float(chipDat[6])/float(micronsPerPixel)
            dx=float(chipDat[4])/float(micronsPerPixel)
            for y in xrange(0,ny):
                for x in xrange(0,nx):
                    dat.append([b+1,y+1,x+1,crow+1,ccol+1,crow*ny+y+1,ccol*nx+x+1,diam,tlx+x*dx,tly+y*dy])
            b+=1
    dat=pandas.DataFrame(dat)
    dat.columns=["Block","Row","Column","ChipRow","ChipColumn","ImageRow","ImageColumn","Diameter","x","y"]
    return(dat)

def zeroLocations(locations):
    '''Set top left hand feature to origin.'''
    locations.x=locations.x-locations.x[(locations.Block==1)&(locations.Row==1)&(locations.Column==1)].values[0]
    locations.y=locations.y-locations.y[(locations.Block==1)&(locations.Row==1)&(locations.Column==1)].values[0]
    return(locations)

def pointsInRect(xCoords,yCoords,x1,y1,x2,y2):
    tlx,tly=min(x1,x2),min(y1,y2)
    brx,bry=max(x1,x2),max(y1,y2)
    return([tlx<=xCoords[i]<=brx and tly<=yCoords[i]<=bry for i in xrange(0,len(xCoords))])

if __name__ == "__main__":
    # Find what directory we're in
    syspath = os.path.dirname(sys.argv[0])
    fullpath = os.path.abspath(syspath)

    # Get current monitor resolution
    pygame.init()
    display=pygame.display.Info()
    scldown=0.9
    dispW,dispH=int(round(scldown*display.current_w)),int(round(scldown*display.current_h))

    # Now, we'll search in this directory for all .gal files, and put their filenames in filelist
    GALlist=[]
    TIFlist=[]
    allfiles=os.listdir(fullpath)
    for filename in allfiles:
        if filename[-4:] in ('.gal','.GAL'):
            GALlist.append(filename)
        if filename[-4:] in ('.tif','.TIF') or filename[-5:] in ('.tiff','.TIFF'):
            TIFlist.append(filename)
    if len(GALlist)>1:
        print "Warning!  Multiple .gal files found.  Only first one will be used..."

    # If we haven't found appropriate files, we'd better raise an error
    if len(GALlist)==0:
        raise Exception("There aren't any .gal files in this directory...  Stopping...")
    if len(TIFlist)==0:
        raise Exception("There aren't any .tif files in this directory...  Stopping...")

    GALfile=GALlist[0]
    TIFfile=TIFlist[0]

    # Check parsing .gal file
    galDict=parseGal(os.path.join(fullpath,GALfile))
    spots=galDict["spotIDs"]
    test=getSpot(1,2,3,spots)
    test.Name.values[0]

    # Get guesses for feature locations from .gal file
    locations=predictLocations(galDict)
    locations=locations.merge(spots,on=("Block","Row","Column"))

    for tif in TIFlist:
        locations=zeroLocations(locations)

        # Open .tif file
        im=Image.open(os.path.join(fullpath,tif))

        # Convert to 3 channel 8-bit image before viewing.
        preview=im.point(lambda i:i*(1./256)).convert('RGB')
        # Stretch contrast to allow us to see signal in darkest images
        sigMax=600.0
        stretch=im.point(lambda i:i*(255.0/sigMax)).convert('RGB')
        # Resize images for display on screen?
        if dispW*float(im.size[1])/float(im.size[0])<=dispH:
            dispH=int(round(im.size[1]*float(dispW)/float(im.size[0])))
        else:
            dispW=int(round(im.size[0]*float(dispH)/float(im.size[1])))
        smallStretch=stretch.resize((dispW,dispH),Image.ANTIALIAS)
        smallPreview=preview.resize((dispW,dispH),Image.ANTIALIAS)

        # What scaling did we actually do there?
        sclX=float(dispW)/float(im.size[0])
        sclY=float(dispH)/float(im.size[1])

        # Create and display a pygame image from PIL image
        surface = pygame.image.fromstring(smallStretch.tostring(),smallStretch.size,smallStretch.mode)
        pygame.display.set_icon(pygame.Surface((100,100)))
        screen = pygame.display.set_mode((dispW,dispH))
        pygame.display.set_caption("Avacta Microarray "+tif.split(".")[0])

        screen.blit(surface, (0, 0))
        pygame.display.flip()

        pointsPlaced=False
        pointsSelected=False
        firstRect=False

        xplus,yplus=0.0,0.0
        cxplus,cyplus=0.0,0.0

        running=True
        while running:
            for event in pygame.event.get():
                if event.type == QUIT or (event.type == KEYDOWN and event.key == K_q):
                    running=False
                    pygame.quit()
                    sys.exit(0)
                if pointsPlaced==False and event.type in (MOUSEBUTTONDOWN,MOUSEMOTION) and pygame.mouse.get_pressed()[0]:
                    # First, place grid estimate from .gal file over features (specfically, guide top left of grid to location of top left culture)
                    mx,my=event.pos
                    radii=[int(round(diam*min(sclX,sclY)/2.0)) for diam in locations.Diameter.values]
                    xcoords=[int(round(mx+x*sclX)) for x in locations.x.values]
                    ycoords=[int(round(my+y*sclY)) for y in locations.y.values]
                    drawgrid=pygame.Surface((dispW,dispH), pygame.SRCALPHA, 32)
                    drawgrid=drawgrid.convert_alpha()
                    for i in xrange(0,len(xcoords)):
                        pygame.draw.circle(drawgrid,(255,0,0,125),(xcoords[i],ycoords[i]),radii[i])
                    screen.blit(surface,(0,0))
                    screen.blit(drawgrid,(0,0))
                    pygame.display.flip()
                if pointsPlaced==False and event.type==KEYDOWN and event.key==K_SPACE:
                    # Once we're happy with the grid estimate, place it
                    pointsPlaced=True
                if pointsPlaced==True and pointsSelected==False and pygame.mouse.get_pressed()[0] and event.type == MOUSEBUTTONDOWN:
                    # Once initial, nudged grid is in place, we can now select sections and move them further by mouse
                    # Click top left to select
                    rect1x,rect1y=event.pos
                    firstRect=True
                if pointsPlaced==True and pointsSelected==False and firstRect==True and pygame.mouse.get_pressed()[0] and event.type == MOUSEMOTION:
                    # Drag bottom right to appropriate location to define selection area
                    rect2x,rect2y=event.pos
                    drawbox=pygame.Surface((dispW,dispH), pygame.SRCALPHA, 32)
                    drawbox=drawbox.convert_alpha()
                    pygame.draw.rect(drawbox,(0,0,255,125),(min(rect1x,rect2x),min(rect1y,rect2y),abs(rect1x-rect2x),abs(rect1y-rect2y)),3)
                    screen.blit(surface,(0,0))
                    screen.blit(drawgrid,(0,0))
                    screen.blit(drawbox,(0,0))
                    pygame.display.flip()
                if firstRect==True and event.type==MOUSEBUTTONUP:
                    # Release mouse button to finalise selection
                    pointsSelected=True
                    firstRect=False
                if pointsSelected==True:
                    if event.type == MOUSEBUTTONDOWN and pygame.mouse.get_pressed()[0]:
                        # Once selection area defined, find points inside selected rectangle
                        startx,starty=event.pos
                        inside=pointsInRect(xcoords,ycoords,rect1x,rect1y,rect2x,rect2y)
                        oldx=xcoords[:]
                        oldy=ycoords[:]
                    if event.type == MOUSEMOTION and pygame.mouse.get_pressed()[0]:
                        # Drag selected points about
                        currx,curry=event.pos
                        drawgrid=pygame.Surface((dispW,dispH), pygame.SRCALPHA, 32)
                        drawgrid=drawgrid.convert_alpha()
                        for i in xrange(0,len(xcoords)):
                            if inside[i]:
                                xcoords[i]=oldx[i]+currx-startx
                                ycoords[i]=oldy[i]+curry-starty
                            pygame.draw.circle(drawgrid,(255,0,0,125),(xcoords[i],ycoords[i]),radii[i])
                        screen.blit(surface,(0,0))
                        screen.blit(drawgrid,(0,0))
                        pygame.display.flip()
                    if event.type==KEYDOWN and event.key==K_SPACE:
                        # Place selected points to allow us to make another selection 
                        pointsSelected=False
                if pointsPlaced==True and event.type==KEYDOWN and event.key == K_f:
                    # Once we're happy with the new grid, write updated guesses (appropriately rescaled) to data frame (and to file)
                    xfinal=[int(round(x/sclX)) for x in xcoords]
                    yfinal=[int(round(y/sclY)) for y in ycoords]
                    locations.x=xfinal
                    locations.y=yfinal
                    locations.to_csv(tif.split(".")[0]+".dat","\t",index=False)
                    pygame.quit()
                    running=False

        # Draw feature location estimates onto stretch and save to file
        stretchpoints=stretch.copy()
        draw=ImageDraw.Draw(stretchpoints)
        for i in xrange(0,len(locations.x)):
            x,y,r=locations.x[i],locations.y[i],float(locations.Diameter[i])/2.0
            draw.ellipse((x-r,y-r,x+r,y+r),outline=(255,0,0))
        stretchpoints.save(tif.split(".")[0]+".jpg")
                
                
            
            
             








