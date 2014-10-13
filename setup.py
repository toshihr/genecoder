# -*- coding: utf-8 -*-
try:
    from Cython.Build import cythonize
    EXIST_CYTHON = True
except ImportError:
    EXIST_CYTHON = False
# http://wiki.python.org/moin/PortingPythonToPy3k
try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    # 2.x
    from distutils.command.build_py import build_py
from setuptools import find_packages
from setuptools.extension import Extension
try:
    from cx_Freeze import setup, Executable
    Executable  # for omitting the flake8 warning
    EXIST_CX_FREEZE = True
except ImportError:
    from setuptools import setup
    EXIST_CX_FREEZE = False
from setuptools.command.test import test as TestCommand
from genecoder.version import VERSION
import os.path
import glob
from itertools import chain

version = VERSION


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        pytest.main(self.test_args)


def _make_extensions(ext_name_with_wildcard):
    ext = '.pyx' if EXIST_CYTHON else '.c'
    filenames = glob.glob(ext_name_with_wildcard.replace('.', os.path.sep) + ext)

    return [Extension(
        name=filename.replace(os.path.sep, '.')[:-len(ext)],
        sources=[filename],
        extra_compile_args=["-Wno-unneeded-internal-declaration", "-Wno-unused-function"],
    ) for filename in filenames]


def _load_requires_from_file(filepath):
    return [pkg_name.rstrip('\r\n') for pkg_name in open(filepath).readlines()]


def _install_requires():
    return _load_requires_from_file('requirements.txt')


def _test_requires():
    return _load_requires_from_file('test-requirements.txt')


packages = find_packages(exclude=['tests'])

# gather package/*.[pyx|c]
extensions = list(chain.from_iterable(map(lambda s: _make_extensions(s + '.*'), packages)))

if EXIST_CYTHON:
    print('running cythonize')
    extensions = cythonize(extensions)

setup(
    name='genecoder',
    version=version,
    description='code analysis for genes',
    url='https://github.com/kerug/genecoder',
    long_description=open('README.rst').read() + '\n' + open('CHANGELOG.txt').read(),
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
    keywords='analyze gene',
    author='kerug',
    author_email='keru.work@gmail.com',
    license='MIT',
    packages=packages,
    package_data={'genecoder': ['genecoder.ini']},
    install_requires=_install_requires(),
    tests_require=_test_requires(),
    test_suite='nose.collector',
    zip_safe=False,
    ext_modules=extensions,
    scripts=[],
    entry_points={'console_scripts': ['genecoder = genecoder.main:main']},
    cmdclass={'build_py': build_py, 'test': PyTest},
)
