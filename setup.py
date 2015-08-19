from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy
import os
import subprocess

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

VERSION='1.0.88'+"."+git_version()[:7]
f=open('colonyzer2/version.py',"w")
f.write("__version__='{}'".format(VERSION))
f.close()


# Nice tutorial on packaging python code:
# http://www.scotttorborg.com/python-packaging/

setup(name='Colonyzer2',
      version=VERSION,
      packages=['colonyzer2','scripts'],
      description='Image analysis for microbial cultures growing on solid agar surfaces',
      long_description=open('README.txt').read(),
      entry_points={"console_scripts":["independent = scripts.independent:main",
                                       "timecourse = scripts.timecourse:main",
                                       "fixedpos = scripts.timecourse_fixedpos:main",
                                       "parametryzer = scripts.parameteryzer_script:main"]},
      author='Conor Lawless',
      author_email='conor.lawless@ncl.ac.uk',
      url='http://research.ncl.ac.uk/colonyzer/',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Intended Audience :: Science/Research'
        ],
      install_requires=['numpy','scipy','pandas','matplotlib','pillow'],
      ext_modules=cythonize("colonyzer2/maskfill.pyx"),
      include_dirs=[numpy.get_include()]
      )
