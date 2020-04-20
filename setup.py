import sys

try:
    from numpy.distutils.core import Extension
    from numpy.distutils.core import setup
except ImportError:
    raise RuntimeError('Numpy Needs to be installed '
                  'for extensions')

# Check python version
if sys.version_info[:2] < (3, 0):
    raise RuntimeError("Python version >= 3.0 required.")

from setuptools import find_packages

#Get lapack for turbosoap (numpy includes one)
#from numpy.distutils.system_info import get_info
#lapack = get_info("lapack_opt")
#if len(lapack) == 0:
#    raise RuntimeError("Numpy lapack not found. Necessary for turbosoap")
#lapack_link = lapack['libraries'] + ["-L" + d for d in lapack["library_dirs"]]

turbosoap_ext = Extension(name='turbosoap_ext',
        sources = ['turbosoap_ext.pyf', 
            'turbogap/src/angular.f90',
            'turbogap/src/functions.f90',
            'turbogap/src/neighbors.f90',
            'turbogap/src/radial.f90',
            'turbogap/src/soap.f90' ,
        ],
        libraries=['lapack']
        #extra_compile_args=["--link-lapack"],
        #extra_link_args = lapack['libraries'] + ["-L" + ]
    )


if __name__ == "__main__":
    setup(name="turbosoap_dscribe",
    version="0.1.0",
    #url="",
    description="A Python package for providing turbosoap machine learning descriptor for DScribe",
    long_description="A Python package for providing turbosoap machine learning descriptor for DScribe",
    packages=find_packages(),
    ext_modules=[turbosoap_ext],
    setup_requires=["numpy"],
    install_requires=["dscribe"]
    )
