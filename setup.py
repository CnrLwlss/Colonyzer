from setuptools import setup, find_packages
import numpy
import os
import subprocess

# http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
# Note that this package requires 'numpy>=1.9.0','scipy>=0.14.1','pandas','matplotlib','pillow','sobol'
# These have not been added to the install_requires list in setup as that leads to very inefficient
# re-compilation of dependencies, even when compatible versions are already available

def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

version='1.1.20'
VERSION=version+"."+git_version()
f=open('colonyzer2/version.py',"w")
f.write("__version__='{}'".format(VERSION))
f.close()

# Nice tutorial on packaging python code:
# http://www.scotttorborg.com/python-packaging/

try:
    from Cython.Distutils import build_ext
except ImportError:
    ext_modules=[]
else:
    from Cython.Build import cythonize
    ext_modules=cythonize("./colonyzer2/maskfill.pyx")

setup(name='Colonyzer2',
      version=version,
      packages=['colonyzer2','scripts'],
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description=open('README.txt').read(),
      entry_points={"console_scripts":["colonyzer = scripts.parseAndRun3:main",
                                       "parametryzer = scripts.parameteryzer_script:main",
                                       "merge = scripts.findAndMerge:main"]},
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ],
      ext_modules=ext_modules,
      include_dirs=[numpy.get_include()],
      data_files=[('data',["./data/BottomRightLid.png","./data/BottomRightNoLid.png","./data/CornerLid.png","./data/CornerNoLid.png","./data/GreenLabLid.png","./data/GreenLabNoLid.png"])]
      )
