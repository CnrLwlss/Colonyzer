import numpy as np
import cv2
import os
import time
import shutil
from PIL import Image
from matplotlib import pyplot as plt
import matplotlib

def getMatches(fname,featlist,frects,orb,draw=False):
    '''Match list of ORB-detected features (featlist) against cropped sections (cropped according to frects) of image file (fname)'''
    if type(fname) in (str, unicode, np.string_):
        fullplate = cv2.imread(fname,0) # trainImage
    else: # assume PIL image object
        fname.load()
        fullplate = fname.convert('L')
        fullplate = np.array(fullplate)
    h,w=fullplate.shape

    bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)

    fres=[]
    for i,f in enumerate(featlist):
        fx0,fy0,fx1,fy1=frects[i]
        plate=fullplate[int(round(fy0*h)):int(round(fy1*h)),int(round(fx0*w)):int(round(fx1*w))]

        kp_plate, des_plate = orb.detectAndCompute(plate,None)
        kp_feature=f["kp"]
        des_feature=f["des"]
        # Match descriptors.
        matches = bf.match(des_feature,des_plate)
        # Sort them in the order of their distance.
        matches = sorted(matches, key = lambda x:x.distance)
        res={"matches":matches}
        if draw:
            res["plate"]=plate
            res["kp_plate"]=kp_plate
        fres.append(res)
    return(fres)

def testMatches(posmatches,negmatches,pmatchdist=20,nmatchdist=20):
    '''Compare lists of positive (lidded) and negative (unlidded) matches to decide whether image is more likely to be lidded or unlidded'''
    pmatches=[[m for m in match["matches"] if m.distance<pmatchdist] for match in posmatches]
    nmatches=[[m for m in match["matches"] if m.distance<nmatchdist] for match in negmatches]
    posHits=[len(match) for match in pmatches]
    negHits=[len(match) for match in nmatches]
    hit=sum(posHits)>sum(negHits)
    return({"hit":hit,"posLengths":posHits,"negLengths":negHits,"posMatches":pmatches,"negMatches":nmatches})

def checkMatches(fname,posfeats,negfeats,frects,orb,pmatchdist=20,nmatchdist=20,draw=False):
    '''Check whether image in fname is more similar to features in posfeats or features in negfeats'''
    posmatches=getMatches(fname,posfeats,frects,orb,draw)
    negmatches=getMatches(fname,negfeats,frects,orb,draw)
    testmatches=testMatches(posmatches,negmatches,pmatchdist,nmatchdist)
    if draw:
        for i,gm,feat in zip(range(0,len(posfeats)),posmatches,posfeats):
            plate=gm["plate"]
            kp_plate=gm["kp_plate"]
            feature=feat["plate"]
            kp_feature=feat["kp"]
            fmatches=testmatches["posMatches"][i]
            print("Number of matches: "+str(len(fmatches)))
            plt.figure(figsize=(20,20))
            img3 = cv2.drawMatches(feature,kp_feature,plate,kp_plate,fmatches, flags=2,outImg=None)
            plt.imshow(img3),plt.savefig(os.path.basename(fname)[0:-4]+"_PosMatches_{:03}.png".format(i),bbox_inches='tight', pad_inches=0)
        for i,gm,feat in zip(range(0,len(negfeats)),negmatches,negfeats):
            plate=gm["plate"]
            kp_plate=gm["kp_plate"]
            feature=feat["plate"]
            kp_feature=feat["kp"]
            fmatches=testmatches["negMatches"][i]
            print("Number of matches: "+str(len(fmatches)))
            plt.figure(figsize=(20,20))
            img3 = cv2.drawMatches(feature,kp_feature,plate,kp_plate,fmatches, flags=2,outImg=None)
            plt.imshow(img3),plt.savefig(os.path.basename(fname)[0:-4]+"_NegMatches_{:03}.png".format(i),bbox_inches='tight', pad_inches=0)
    return(testmatches["hit"])

def calcHits(matchDicts,actuallyLidded,actuallyUnlidded,pmatchdist=20,nmatchdist=20,report=True):
    '''Check automated classification against a set of manually classified images'''
    fnames=set([m[0:-4] for m in matchDicts.keys()])
    lidded={f:testMatches(matchDicts[f+"_POS"],matchDicts[f+"_NEG"],pmatchdist,nmatchdist)["hit"] for f in fnames}
    probs=[l for l in lidded.keys() if lidded[l]]
    noprobs=[l for l in lidded.keys() if not lidded[l]]
    
    truehits=[f for f in probs if f in actuallyLidded]
    falsehits=[f for f in probs if f in actuallyUnlidded]
    truemisses=[f for f in noprobs if f in actuallyUnlidded]
    falsemisses=[f for f in noprobs if f in actuallyLidded]

    fracTrueHits=float(len(truehits))/len(actuallyLidded)
    fracFalseHits=float(len(falsehits))/len(actuallyUnlidded)
    fracTrueMisses=float(len(truemisses))/len(actuallyUnlidded)
    fracFalseMisses=float(len(falsemisses))/len(actuallyLidded)

    if report:
        print("{} lidded plates were detected".format(len(probs)))
        print("{0:.2f}% true hits".format(100.0*fracTrueHits))
        print("{0:.2f}% false hits".format(100.0*fracFalseHits))
        print("{0:.2f}% true misses".format(100.0*fracTrueMisses))
        print("{0:.2f}% false misses".format(100.0*fracFalseMisses))
    return([fracTrueHits,fracFalseHits,fracTrueMisses,fracFalseMisses])

# find the keypoints and descriptors with SIFT
def getFeatures(imlist,orb):
    '''Get ORB features in each OpenCV image in imlist'''
    flist=[]
    for im in imlist:
       kp_feature, des_feature = orb.detectAndCompute(im,None)
       fdict={"kp":kp_feature,"des":des_feature,"plate":im}
       flist.append(fdict)
    return(flist)

def makeLidTest(posfiles,negfiles,frects,pmatchdist=50,nmatchdist=50,draw=False):
    '''Function closure returning function to check whether fname has lid'''
    posims = [cv2.imread(fname,0) for fname in posfiles]
    negims = [cv2.imread(fname,0) for fname in negfiles]
    # Initiate detector
    orb = cv2.ORB_create()
    posfeats=getFeatures(posims,orb)
    negfeats=getFeatures(negims,orb)
    def lidTest(fname):
        '''Does image in fname have a lid?'''
        return(checkMatches(fname,posfeats,negfeats,frects,orb,pmatchdist,nmatchdist,draw))
    return(lidTest)



  



