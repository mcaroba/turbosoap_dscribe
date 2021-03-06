import sys
from setuptools import find_packages
from distutils.core import setup

try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError:
    raise RuntimeError('Numpy Needs to be installed '
                  'for extensions')

# Check python version
if sys.version_info[:2] < (3, 0):
    raise RuntimeError("Python version >= 3.0 required.")


import numpy.distutils.system_info as sysinfo
lapack_opt = sysinfo.get_info('lapack_opt')
if len(lapack_opt) == 0:
    raise RuntimeError('No Lapack found, please install Lapack libraries to install turbosoap_dscribe')
libraries = lapack_opt['libraries']
library_dirs = lapack_opt['library_dirs']

turbosoap_ext = Extension(name='turbosoap_ext',
        sources = ['turbosoap_ext.pyf', 
            'turbogap/src/angular.f90',
            'turbogap/src/functions.f90',
            'turbogap/src/neighbors.f90',
            'turbogap/src/radial.f90',
            'turbogap/src/soap.f90' ,
        ],
        libraries=libraries,
        library_dirs=library_dirs
    )


if __name__ == "__main__":
    setup(name="turbosoap_dscribe",
    version="0.1.2",
    url="",
    description="An optional dependency of DScribe which provides TurboSOAP machine learning descriptor",
    long_description="An optional dependency of DScribe which provides TurboSOAP machine learning descriptor",
    packages=find_packages(),
    ext_modules=[turbosoap_ext],
    setup_requires=["numpy"],
    install_requires=["dscribe"]
    )
