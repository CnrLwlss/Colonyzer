# This set of functions archives Lydall style images in
# YEAR - MONTH directory structure. Useful for archiving
# genome-wide timeseries images
import os,sys, shutil, time

syspath = os.path.dirname(sys.argv[0])
fullpath = os.path.abspath(syspath)
today=time.gmtime()
fltoday=float(today.tm_year)+(float(today.tm_mon)-1.0)/12

def getfiles(fullpath):
    return os.listdir(fullpath)

def getfolders(fullpath):
    allfiles=getfiles(fullpath)
    allfolders=[]
    for f in allfiles:
        if "." not in f:
            allfolders.append(f)
    return allfolders

def getJPG(fullpath):
    allfiles=getfiles(fullpath)
    allJPG=[]
    for f in allfiles:
        if f[-4:]==".JPG" or f[-4:]==".jpg":
            allJPG.append(f[0:-4])
    return allJPG

def getDAT(fullpath):
    allfiles=getfiles(fullpath)
    allDAT=[]
    for f in allfiles:
        if f[-7:]=="OUT.dat":
            allDAT.append(f[0:-7])
    return allDAT

def getDate(fname):
    dt=fname.split("_")[-2]
    year=dt.split("-")[0]
    month=dt.split("-")[1]
    return([year,month])


# Get name of all unarchived JPGs in main directory
JPGs=getJPG(fullpath)

# Get all the .dat names in the Output_Data directory
datdirect=os.path.join(fullpath,"Output_Data")
predirect=os.path.join(fullpath,"Output_Images")
if (os.path.exists(datdirect)):
    DATs=getDAT(datdirect)
else:
    print("Output_Data directory does not exist!")
    print("Please analyse images before archiving")
    sys.exit(0)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Find years and months of images from filenames
# Must be in XXXXXXXXXXXXXX_YYYY-MM-DD_HH-MM-SS.JPG format
# where XXXXXXXXXXXXXX can be used to identify specific plates

YearMonths={}
for JPG in JPGs:
    ym=getDate(JPG)
    # Check if date is valid (if filename is parseable)
    if len(ym[0])==4 and len(ym[1])==2 and is_number(ym[0]) and is_number(ym[1]):
        if ym[0] not in YearMonths:
            YearMonths[ym[0]]=[ym[1]]
        if ym[1] not in YearMonths[ym[0]]:
            YearMonths[ym[0]].append(ym[1])

print("Building directories...")
# For every year month in list above, generate directory if missing
for year in YearMonths:
    if not os.path.exists(os.path.join(fullpath,"Year_"+year)):
        os.mkdir(os.path.join(fullpath,"Year_"+year))
    for month in YearMonths[year]:
        if not os.path.exists(os.path.join(fullpath,"Year_"+year,"Month_"+month)):
            os.mkdir(os.path.join(fullpath,"Year_"+year,"Month_"+month))
        if not os.path.exists(os.path.join(fullpath,"Year_"+year,"Month_"+month,"Output_Data")):
            os.mkdir(os.path.join(fullpath,"Year_"+year,"Month_"+month,"Output_Data"))
        if not os.path.exists(os.path.join(fullpath,"Year_"+year,"Month_"+month,"Output_Images")):
            os.mkdir(os.path.join(fullpath,"Year_"+year,"Month_"+month,"Output_Images"))
        #shutil.copyfile(os.path.join(fullpath,"Colonyzer.py"),os.path.join(fullpath,"Year_"+year,"Month_"+month,"Colonyzer.py"))
        shutil.copyfile(os.path.join(fullpath,"Colonyzer.txt"),os.path.join(fullpath,"Year_"+year,"Month_"+month,"Colonyzer.txt"))

print("Moving images and data files...")
# Copy JPGS, DATS and previews across
toDelete=[]
for JPG in JPGs:
    try:
        # Check if it is too soon to archive this file
        ym=getDate(JPG)
        # Check if date is valid (if filename is parseable)
        if len(ym[0])==4 and len(ym[1])==2 and is_number(ym[0]) and is_number(ym[1]):
            fldate=float(ym[0])+(float(ym[1])-1.0)/12.0
            if fltoday-fldate>3.0/12.0:
                fname=os.path.join(datdirect,JPG+".dat")
                iname=os.path.join(fullpath,JPG+".JPG")
                pname=os.path.join(predirect,JPG+".png")
                ftarg=os.path.join(fullpath,"Year_"+ym[0],"Month_"+ym[1],"Output_Data",JPG+".dat")
                itarg=os.path.join(fullpath,"Year_"+ym[0],"Month_"+ym[1],JPG+".JPG")
                ptarg=os.path.join(fullpath,"Year_"+ym[0],"Month_"+ym[1],"Output_Images",JPG+".png")
                xfname=os.path.exists(fname)
                xiname=os.path.exists(iname)
                xpname=os.path.exists(pname)
                xftarg=os.path.exists(ftarg)
                xitarg=os.path.exists(itarg)
                xptarg=os.path.exists(ptarg)
                # If files have already been backed up, delete them from folder
                #if xitarg and xiname:
                #    os.remove(iname)
                #if xftarg and xfname:
                #    os.remove(fname)
                #if xptarg and xpname:
                #    os.remove(pname)
                # If files have not already been backed up and analysis complete:
                if xfname and xiname and xpname:
                    shutil.copyfile(fname,ftarg)
                    shutil.copyfile(iname,itarg)
                    shutil.copyfile(pname,ptarg)
                    toDelete.append(JPG)
                else:
                    print(JPG, "Missing files.")
            else:
                print(JPG, " too recent")
    except:
        print(JPG, " failed")

# If that's been successful, let's delete the remaining files
for JPG in toDelete:
    fname=os.path.join(datdirect,JPG+".dat")
    iname=os.path.join(fullpath,JPG+".JPG")
    pname=os.path.join(predirect,JPG+".png")
    os.remove(fname)
    os.remove(iname)
    os.remove(pname)
