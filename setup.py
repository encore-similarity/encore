#! /usr/bin/python
"""Setuptools-based setup script for ENCORE.

ENCORE requires:

	* A working python2.6 distribution or higher (2.X)
	* the MDAnalysis package ( http://https://code.google.com/p/mdanalysis/ )
	* a command-line C compiler (as gcc or clang)
	* a working Cython distribution

For a basic installation just type the command::

  python setup.py install


  http://groups.google.com/group/mdnalysis-discussion


By default we use setuptools <http://pypi.python.org/pypi/setuptools>.  The
details of such an "EasyInstall" installation procedure are shown on

  http://peak.telecommunity.com/DevCenter/EasyInstall

By changing the code below you can also switch to a standard distutils
installation.
"""

#------------------------------------------------------------
# selection of the installation system
#------------------------------------------------------------
#
# Standard distutils-based installation:
#
##from distutils.core import setup, Extension

# setuptools ("EasyInstall") installation:
#
# If you want EasyInstall features then enable the next three lines and comment
# out the preceding line 'from distutils.core import ...'
#

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension
from distutils.ccompiler import new_compiler
from distutils.version import LooseVersion
#
#------------------------------------------------------------

import sys, os
import shutil
import tempfile

REQUIRED_MDANALYSIS = '0.8.1'
REQUIRED_MDANALYSIS = LooseVersion(REQUIRED_MDANALYSIS)

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 6):
    print "ENCORE requires Python 2.6 or better. Python %d.%d detected" % sys.version_info[:2]
    print "Please upgrade your version of Python."
    sys.exit(-1)

try:
    # Find out if MDAnalysis is installed
    import MDAnalysis
except ImportError:
    print "*** package 'MDAnalysis' not found ***"
    sys.exit(-1)
    
    # Check MDAnalysis version
if LooseVersion(MDAnalysis.__version__) < REQUIRED_MDANALYSIS:
    print "*** package 'MDAnalysis' was found, but supported version is >= 0.8.1 ***"
    sys.exit(-1)

import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dirs = [numpy_include]

# Handle cython modules
try:
    from Cython.Distutils import build_ext
    use_cython = True
    cmdclass = {'build_ext': build_ext}
except ImportError:
    use_cython = False
    cmdclass = {}

def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            f = open(fname, 'w')
            if include is not None:
                f.write('#include %s\n' % include)
            f.write('int main(void) {\n')
            f.write('    %s;\n' % funcname)
            f.write('}\n')
            f.close()
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception as e:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)


if __name__ == '__main__':
    RELEASE = "0.1b"     # NOTE: keep in sync with MDAnalysis.version in __init__.py

    if 'DEBUG_CFLAGS' in os.environ:
        extra_compile_args = '\
            -std=c99 -pedantic -Wall -Wcast-align -Wcast-qual -Wpointer-arith \
            -Wchar-subscripts -Winline -Wnested-externs -Wbad-function-cast \
            -Wunreachable-code -Werror'
        define_macros = [('DEBUG', '1')]
    else:
        extra_compile_args = ''
        define_macros = []

    with open("SUMMARY.txt") as summary:
        LONG_DESCRIPTION = summary.read()
    
    extensions = [Extension('cutils',['src/cutils/cutils.%s' % ("pyx" if use_cython else "c")],
                            include_dirs = include_dirs,
                            extra_compile_args = ["-O3", "-ffast-math"]),
                  Extension('clustering.affinityprop', ['src/clustering/affinityprop.%s' % ("pyx" if use_cython else "c"), "src/clustering/ap.c"],
                            include_dirs = include_dirs.extend(['src/clustering']),
                            libraries=["m"],
                            extra_compile_args=["-O3", "-ffast-math","-std=c99"]),
                  Extension('dimensionality_reduction.stochasticproxembed', 
                            ['src/dimensionality_reduction/stochasticproxembed.%s' % ("pyx" if use_cython else "c"), "src/dimensionality_reduction/spe.c"],
                            include_dirs = include_dirs.extend(['src/dimensional_reduction']),
                            libraries=["m"],
                            extra_compile_args=["-O3", "-ffast-math","-std=c99"])
		]
                  

    setup(name              = 'ENCORE',
          version           = RELEASE,
          description       = 'A module built ',
          author            = '*',
          author_email      = '',
          url               = 'http://somewhere.over.the.rain.bow/',
          requires          = ['MDAnalysis'],
          provides          = ['encore'],
          license           = 'GPL 2',
          packages          = [ 'encore',
				'encore.clustering',
				'encore.dimensionality_reduction' ],
          package_dir       = {'encore': 'encore'},
          ext_package       = 'encore',
          ext_modules       = extensions,
          long_description  = LONG_DESCRIPTION,
          cmdclass          = cmdclass,
	  include_dirs	    = include_dirs,
          zip_safe          = False
          )
