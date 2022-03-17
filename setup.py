from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.install import install as _install
from setuptools import setup
import setuptools
from setuptools.extension import Extension
from subprocess import run
import os
import glob
import shutil

extras = ["-Wno-sign-compare", "-Wno-unused-function", "-Wno-unused-result", '-Wno-ignored-qualifiers',
          "-Wno-deprecated-declarations"]

ext_modules = list()

root = os.path.abspath(os.path.dirname(__file__))
wfa = os.path.join(root, "pywfa/WFA2_lib")

libraries = [f"{wfa}/lib"]
library_dirs = [f"{wfa}/lib"]
include_dirs = [".", root, wfa, f"{wfa}/lib", f"{wfa}/utils", f"{wfa}/wavefront", f"{wfa}/bindings/cpp", f"{wfa}/system",
                f"{root}/pywfa"]

print("Libs", libraries)
print("Library dirs", library_dirs)
print("Include dirs", include_dirs)

ext_modules.append(Extension("pywfa.align",
                             ["pywfa/align.pyx"],
                             libraries=libraries,
                             library_dirs=library_dirs,
                             include_dirs=include_dirs,
                             extra_compile_args=extras,
                             language="c",
                             extra_objects=[f"{wfa}/lib/libwfa.a"]
                             ))


class MyBuild(_build_ext):
    # https://stackoverflow.com/questions/1754966/how-can-i-run-a-makefile-in-setup-py
    def build_extension(self, ext):
        # This is a hack to automate python setup.py build_ext --inplace during install
        # Note sure how robust this is
        if self.debug:
            ext.extra_compile_args.append("-O0")

        shared = glob.glob(f"{root}/build/lib*/pywfa/*.so") + glob.glob(f"{root}/build/lib*/pywfa/*.dll")
        _build_ext.build_extension(self, ext)
        [shutil.copy(i, f"{root}/pywfa") for i in shared]


class Build_ext_first(_install):
    def run(self):
        # Build WFA2-lib first
        # run("cd pywfa/WFA2_lib; make BUILD_CPP=0 clean all; cd ../../", shell=True)
        run("cd pywfa/WFA2_lib; make clean all; cd ../../", shell=True)
        return setuptools.command.install.install.run(self)

setup(
    name="pywfa",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    url="https://github.com/kcleal/pywfa",
    description="Align sequences using WFA2-lib",
    license="MIT",
    version='0.2.5',
    python_requires='>=3.7',
    install_requires=[  # runtime requires
            'cython',
        ],
    setup_requires=[
            'cython',
        ],
    tests_require=[
            'pysam', 'nose'
        ],
    packages=["pywfa", "pywfa/tests"],
    ext_modules=cythonize(ext_modules),
    cmdclass={'install': Build_ext_first, 'build_ext': MyBuild},
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
)
