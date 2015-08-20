Colonyzer
=========

Image analysis software for estimating cell density in arrayed microbial cultures growing on solid agar.

Website describing current, stable release:
http://research.ncl.ac.uk/colonyzer/

Open access manuscript describing Colonyzer algorithms:
http://dx.doi.org/10.1186/1471-2105-11-287

Open access video and manuscript demonstrating the use of Colonyzer within a Quantitative Fitness Analysis workflow
http://www.jove.com/video/4018

Building this package generates a lot of intermediate files and directories.  Analysing the images in the data directory also generates output files that are not under version control.
To carry out a dry run of clearing out all unversioned files and directories from your working copy (remove '-n' argument to irreversibly delete files):

    git clean -f -d -x -n

To register project with PYPI:

	python setup.py register 
	
Once project has been registered, submit updates by:

    python setup.py sdist upload