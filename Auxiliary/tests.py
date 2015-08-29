import os,numpy,itertools
from colonyzer2 import *
from scipy import optimize as op

LATESTIMAGE=os.path.realpath("F:\Colonyzer\Auxiliary\Data\DLR00012647-2009-07-04_09-35-20.jpg")
im,arr=openImage(LATESTIMAGE)

nx,ny=24,16
windowFrac=0.25
smoothWindow=0.13
showPlt=True
pdf=None
diam=20



