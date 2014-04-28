#!/usr/bin/env python

from distutils.core import setup
from setuptools.command.install import install
from setuptools.command.develop import develop

def installPeano3():
  p3Path = 'src/peano3'
  p3RepositorySubpath = '.'
  p3Revision = 'HEAD'
  p3Build = 'release'
  try:
    import peanoConfiguration
    p3Path = peanoConfiguration.getPeano3Path()
    p3RepositorySubpath = peanoConfiguration.getPeano3RepositorySubpath()
    p3Revision = peanoConfiguration.getPeano3Revision()
    p3Build = peanoConfiguration.getPeano3Build()
    p3Dimension = peanoConfiguration.getPeano3Dimension()
  except ImportError:
    pass
  
  from os.path import join
  from os.path import exists
  from subprocess import call
  if exists(join(p3Path, '.svn')):
    pass
    #print("Updating Peano3 Repository")
    #call("svn update -r" + p3Revision + " " + p3Path, shell=True)
  else:
    print("Checking out Peano3 Repository")
    call("svn checkout -r" + p3Revision + join(" svn://svn.code.sf.net/p/peano/code/trunk", p3RepositorySubpath) + " " + p3Path, shell=True)
  print("Building PeanoClaw")
  returnValue = call("scons build=" + str(p3Build) + " parallel=yes dim=2" + " -j2", shell=True)
  
  if returnValue != 0:
    raise Exception("PeanoClaw: Build failed.")


class Peano3Install(install):
  def run(self):
    installPeano3()
    install.run(self)
    
class Peano3Develop(develop):
  def run(self):
    installPeano3()
    develop.run(self)
    
#Read command line parameters
import sys

setup(name='PeanoClaw',
      version='0.1',
      description='PeanoClaw - AMR Extension for PyClaw/Clawpack',
      author='Kristof Unterweger',
      author_email='unterweg@in.tum.de',
      url='http://github.com/unterweg/peanoclaw',
      packages=['peanoclaw', 'peanoclaw.callbacks'],
      package_dir={'': 'src/python'},
      package_data={'clawpack.peanoclaw': ['libpeano-claw-*']},
      cmdclass={'install': Peano3Install, 'develop': Peano3Develop}
     )