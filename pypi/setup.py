from setuptools import setup, find_packages
setup(name='Colonyzer2',
      version='1.0.80',
      packages=['colonyzer2'],
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description=open('README.txt').read(),
      entry_points={"console_scripts":["independent = colonyzer2.independent:main",
                                       "timecourse = colonyzer2.timecourse:main",
                                       "fixedpos = colonyzer2.timecourse_fixedpos:main",
                                       "parametryzer = colonyzer2.parameteryzer_script:main"]},
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ]      
      )
