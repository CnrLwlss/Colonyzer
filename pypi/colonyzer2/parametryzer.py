import pygame
from pygame.locals import *
import os,sys,time

# Let's find what directory we're in
#syspath = os.path.dirname(sys.argv[0])
#fullpath = os.path.abspath(syspath)
fullpath=os.getcwd()

useRoot=False

filelist=[]
rootdict={}
allfiles=os.listdir(fullpath)
for filename in allfiles:
    if len(filename)>4 and filename[-4:] in ['.jpg','.JPG','.jpeg','.JPEG','.tif','.TIF','.tiff','.TIFF','.PNG','.png']:
        filelist.append(filename)
        # Old style DLR barcodes
        if len(filename)==35:
            root=filename[0:11]
            if root not in rootdict:
                rootdict[root]=[filename]
            else:
                rootdict[root].append(filename)
        # New imager filenames
        # Lydall lab file naming convention
        # First 15 characters in filename identify unique plates
        # Remaining charaters can be used to store date, time etc.
        else:
            root=filename[0:15]
            if root not in rootdict:
                rootdict[root]=[filename]
            else:
                rootdict[root].append(filename)

# If we haven't found any image files, we'd better raise an error
if len(filelist)==0:
    raise Exception("There aren't any suitable image files in this directory!")

filedict={}
# Attempt to parse the filenames and put them in a sensible order (reverse time order for each plate)
### Normal imager output
##if len(filelist[0])==35:
##    for x in filelist:
##        filedict[x[0:31]]=x
### Steve type output
##elif len(filelist[0])==41:
##    for x in filelist:
##        filedict[x[4:35]]=x
##else:
##    for x in filelist:
##        filedict[x]=x

# Attempt to parse the filenames and put them in a sensible order (reverse time order)
# Normal imager output
for x in filelist:
    # Old Barcode style
    if len(filelist[0])==35:
        filedict[x[12:31]]=x
    # New SPImager style output
    elif len(filelist[0])==39:
        filedict[x[16:35]]=x
    else:
        filedict[x]=x

# Get the keys
keylist=filedict.keys()
keylist.sort()
keylist.reverse()
newfiles=[]
for x in keylist:
    newfiles.append(filedict[x])

if useRoot:
    newfiles=[]
    for x in rootdict:
        rootdict[x].sort()
        newfiles.append(rootdict[x][-1])

# Typical AR for latest SPImager images
#w,h,=1000,667
w,h=1200,800
size=w,h
counter=0
LAST=()
Results=[]
for x in xrange(len(newfiles)):
    Results.append([])

BkCol=(0,0,0)
# Initialise screen
pygame.init()
screen = pygame.display.set_mode(size)

def InitializeImage():
    global mx,my,linew,dotradius,TL,BR, finished
    finished=False
    pygame.display.set_caption("ColonyzerParametryzer -"+newfiles[counter]+"- Click on centre of top left spot, press r to restart, t to go back one image, b to accept current points.")
    mx,my=0,0
    linew=1
    dotradius=2
    TL=()
    BR=()
    # Blit everything to the screen
    screen.blit(PhotoSurface, (0, 0))
    if LAST!=():
        LTL=pygame.draw.circle(screen, (0,0,255), LAST[0], dotradius)
        LBR=pygame.draw.circle(screen, (0,0,255), LAST[1], dotradius)
    pygame.display.flip()

def SaveResults(Results,Format):
    # Scale the results file up to the original image coordinates
    NewRes=[]
    for zz in xrange(0,len(Results)):
        x=Results[zz]
        filenm=newfiles[zz]
        PhotoName = os.path.join(fullpath,filenm)
        PhotoSurface = pygame.image.load(PhotoName)
        ow,oh=PhotoSurface.get_size()
        NewRes.append(((int(round(float(ow)*float(x[0][0])/float(w))),int(round(float(oh)*float(x[0][1])/float(h)))),(int(round(float(ow)*float(x[1][0])/float(w))),int(round(float(oh)*float(x[1][1])/float(h))))))
    #print NewRes

    defstring="default,%s,%s,%s,%s,%s,%s"%(Format,str(NewRes[0][0][0]),str(NewRes[0][0][1]),str(NewRes[0][1][0]),str(NewRes[0][1][1]),time.strftime("%Y-%m-%d"))

    specstr=""
    for x in xrange(0,len(Results)):
        if useRoot:
            filenm=rootdict.keys()[x]
        else:
            filenm=newfiles[x]
        specstr+="%s,%s,%s,%s,%s,%s\n"%(filenm,Format,str(NewRes[x][0][0]),str(NewRes[x][0][1]),str(NewRes[x][1][0]),str(NewRes[x][1][1]))

    OutString="""####################################
# Colonyzer information file
####################################
# Please input comma-separated default values for:
# No. of spots (e.g. 96,384 or 1536), top left spot x, top left spot y, bottom right spot x, bottom right spot y
# Where the x and y coordinates refer to the centres of the spots
# with y-axis increasing moving from top of image towards bottom
# with x-axis increasing moving from left of image towards right.
# This is the default reading you get hovering mouse pointer over
# an image in ImageJ (http://rsbweb.nih.gov/ij/) for example.
# Finally, add the calibration date in ISO format: YYYY-MM-DD to allow time-dependent calibrations
# Note that the comma-separated list must be preceded by:
# default
####################################

"""+defstring+"""

####################################
# Please input any image-specific differences here
# same format as before, but include the filename 
# of the image followed by spot no. and coordinates
# as above.  Calibration date should be omitted in this case.
####################################
"""+specstr
    #print OutString
    final=open(os.path.join(fullpath,"Colonyzer.txt"),"w")
    final.write(OutString)
    final.close()


PhotoName = os.path.join(fullpath,newfiles[counter])
PhotoSurface = pygame.image.load(PhotoName)
ow,oh=PhotoSurface.get_size()
PhotoSurface = pygame.transform.scale(PhotoSurface,(w,h))
background = pygame.Surface(size)
InitializeImage()

# Blit everything to the screen
screen.blit(PhotoSurface, (0, 0))
pygame.display.flip()

while 1:
        for event in pygame.event.get():
                if event.type == QUIT or (event.type ==KEYDOWN and event.key==K_q):
                        pygame.quit()
                        sys.exit(0)
                else:
                    #print event
                    if event.type == MOUSEMOTION:
                        mx,my=event.pos
                        screen.blit(PhotoSurface, (0, 0))
                        if TL==() and BR==():
                            if LAST!=():
                                LTL=pygame.draw.circle(screen, (0,0,255), LAST[0], dotradius)
                                LBR=pygame.draw.circle(screen, (0,0,255), LAST[1], dotradius)
                            L1=pygame.draw.line(screen, (255,0,0), (0,my), (w,my), linew)
                            L2=pygame.draw.line(screen, (255,0,0), (mx,0), (mx,h), linew)

                        elif BR==():
                            if LAST!=():
                                LTL=pygame.draw.circle(screen, (0,0,255), LAST[0], dotradius)
                                LBR=pygame.draw.circle(screen, (0,0,255), LAST[1], dotradius)
                            PTL=pygame.draw.circle(screen, (255,0,0), TL, dotradius)
                            L1=pygame.draw.line(screen, (0,255,0), (0,my), (w,my), linew)
                            L2=pygame.draw.line(screen, (0,255,0), (mx,0), (mx,h), linew)
                        else:
                            if LAST!=():
                                LTL=pygame.draw.circle(screen, (0,0,255), LAST[0], dotradius)
                                LBR=pygame.draw.circle(screen, (0,0,255), LAST[1], dotradius)
                            PTL=pygame.draw.circle(screen, (255,0,0), TL, dotradius)
                            PBR=pygame.draw.circle(screen, (0,255,0), BR, dotradius)
                            
                        pygame.display.flip()
                    if event.type == MOUSEBUTTONUP and TL==() and BR==():
                        TL=event.pos
                        pygame.display.set_caption("ColonyzerParametryzer -"+newfiles[counter]+"- Now click on centre of bottom right spot, press r to restart, t to go back one image.")
                    elif event.type == MOUSEBUTTONUP and BR==():
                        BR=event.pos
                        pygame.display.set_caption("ColonyzerParametryzer -"+newfiles[counter]+"- Now press spacebar to save, press r to restart.")
                    # Press r to reset current image
                    if event.type == KEYDOWN and event.key==K_r:
                        InitializeImage()
                    # Press t to go back one image
                    if event.type == KEYDOWN and event.key==K_t:
                        counter=max(0,counter-1)
                        PhotoName = os.path.join(fullpath,newfiles[counter])
                        PhotoSurface = pygame.image.load(PhotoName)
                        ow,oh=PhotoSurface.get_size()
                        PhotoSurface = pygame.transform.scale(PhotoSurface,(w,h))
                        background = pygame.Surface(size)
                        InitializeImage()
                    # Press b to accept the last points as valid for current image
                    if event.type == KEYDOWN and event.key==K_b and LAST!=():
                        TL=LAST[0]
                        BR=LAST[1]
                        PTL=pygame.draw.circle(screen, (255,0,0), TL, dotradius)
                        PBR=pygame.draw.circle(screen, (0,255,0), BR, dotradius)
                        pygame.display.flip()
                        pygame.display.set_caption("ColonyzerParametryzer -"+newfiles[counter]+"- Now press spacebar to save, press r to restart.")
                    if event.type == KEYDOWN and event.key==K_SPACE and TL!= () and BR!=():
                        # If not on last image yet...
                        if counter < len(newfiles)-1:
                            LAST=(TL,BR)
                            Results[counter]=[TL,BR]

                            counter+=1
                            PhotoName = os.path.join(fullpath,newfiles[counter])
                            PhotoSurface = pygame.image.load(PhotoName)
                            ow,oh=PhotoSurface.get_size()
                            PhotoSurface = pygame.transform.scale(PhotoSurface,(w,h))
                            background = pygame.Surface(size)
                            InitializeImage()
                        else:
                            LAST=(TL,BR)
                            Results[counter]=[TL,BR]
                            
                            pygame.display.set_caption("ColonyzerParametryzer - Press key to select format & save: 1536 (d), 384 (f), 117 (p), 96 (g), 48 (h)")
                            finished=True
                    if finished==True and event.type == KEYDOWN and (event.key==K_d or event.key==K_f or event.key==K_p or event.key==K_g or event.key==K_h):
                        if event.key==K_d:
                            Format='1536'
                        elif event.key==K_f:
                            Format='384'
                        elif event.key==K_p:
                            Format='117'
                        elif event.key==K_g:
                            Format='96'
                        elif event.key==K_h:
                            Format='48'
                        # write file...
                        SaveResults(Results,Format)
                        #print Results
                        pygame.quit()
                        sys.exit(0)
                            
       
                        
