#cython.wraparound=False
#cython.boundscheck=False

import numpy as np
cimport numpy as np

DTYPEf = np.float
ctypedef np.float_t DTYPEf_t

DTYPEb = np.uint8
ctypedef np.uint8_t DTYPEb_t

def maskAndFill(arrNArg,finalMaskArg,np.float_t tol=5.0):
	'''Cut out masked pixels from image and re-fill using a Markov field update.'''
	cdef int h = arrNArg.shape[0]
	cdef int w = arrNArg.shape[1]
	
	cdef np.ndarray[np.float_t, ndim=2, mode="c"] arrN
	cdef np.ndarray[np.uint8_t, ndim=2, mode="c"] finalMask
	arrN = np.ascontiguousarray(arrNArg,dtype=np.float)
	cutout_arr = np.ascontiguousarray(arrNArg,dtype=np.float)
	finalMask = np.ascontiguousarray(finalMaskArg,dtype=np.uint8)
	
	# Unmask edges to allow Markov field update
	finalMask[0,:] = 0
	finalMask[-1,:] = 0
	finalMask[:,0] = 0
	finalMask[:,-1] = 0
	cdef float diff = 100*tol
	cdef float cumsum
	cdef float count
	cdef int x, y
	cdef np.ndarray old = np.zeros((h,w), dtype = DTYPEf)
	cdef np.ndarray filled = np.zeros((h,w), dtype = DTYPEb)
	for y in xrange(1,h-1):
		for x in xrange(1,w-1):
			filled[y,x]=1-finalMask[y,x]
	
	while diff>tol or np.isnan(diff):
		# Invert filling order at every pass to minimise bias towards a particular direction
		print diff>tol, np.isnan(diff)
		for y in xrange(1,h-1):
			for x in xrange(1,w-1):
				# Markov field update
				if finalMask[y,x]==1:
					# Store last image to allow assessment of convergence
					old[y,x]=cutout_arr[y,x]
					# Initialise counters for calculation of mean
					count=0.0
					cumsum=0.0
					if filled[y-1,x]>0:
						count=count+1.0
						cumsum=cumsum+float(cutout_arr[y-1,x])
					if filled[y+1,x]>0:
						count=count+1.0
						cumsum=cumsum+float(cutout_arr[y+1,x])
					if filled[y,x-1]>0:
						count=count+1.0
						cumsum=cumsum+float(cutout_arr[y,x-1])
					if filled[y,x+1]>0:
						count=count+1.0
						cumsum=cumsum+float(cutout_arr[y,x+1])
					cutout_arr[y,x]=cumsum/count
					filled[y,x]=1
		diff=float(np.sum(np.abs(old*finalMask-cutout_arr*finalMask)))/float(np.sum(finalMask))
	return(cutout_arr)