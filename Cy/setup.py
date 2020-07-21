# all .pyx files in a folder
from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Acceleration',
  ext_modules = cythonize(["*.pyx"]),
)