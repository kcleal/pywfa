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
from distutils import ccompiler
import distutils.sysconfig


debug = False

cy_options = {
    'annotate': False,
    'compiler_directives': {
        'profile': debug,
        'linetrace': debug,
        'boundscheck': debug,
        'wraparound': debug,
        'nonecheck': debug,
        'initializedcheck': debug,
        'language_level': 3
    }
}

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.c') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def get_extra_args(flags):
    compiler = ccompiler.new_compiler()
    extra_compile_args = []
    for f in flags:
        if has_flag(compiler, f):
            extra_compile_args.append(f)
    return extra_compile_args


extras = ["-Wno-unused-function", "-Wno-unused-result", '-Wno-ignored-qualifiers', "-Wno-deprecated-declarations"]
extras_pywfa = get_extra_args(extras)

root = os.path.abspath(os.path.dirname(__file__))
wfa = os.path.join(root, "pywfa/WFA2_lib")
libraries = [f"{wfa}/lib"]
library_dirs = [f"{wfa}/lib"]
include_dirs = [".", root, wfa, f"{wfa}/lib", f"{wfa}/utils", f"{wfa}/wavefront", f"{wfa}/bindings/cpp", f"{wfa}/system",
                f"{wfa}/alignment", f"{root}/pywfa"]

print("Libs", libraries)
print("Library dirs", library_dirs)
print("Include dirs", include_dirs)

ext_modules = list()
ext_modules.append(Extension("pywfa.align",
                             ["pywfa/align.pyx"],
                             libraries=libraries,
                             library_dirs=library_dirs,
                             include_dirs=include_dirs,
                             extra_compile_args=extras_pywfa,
                             language="c",
                             extra_objects=[f"{wfa}/lib/libwfa.a"]
                             ))

openmp_support = False


class MyBuild(_build_ext):
    # https://stackoverflow.com/questions/1754966/how-can-i-run-a-makefile-in-setup-py
    def build_extension(self, ext):
        # This is a hack to automate python setup.py build_ext --inplace during install
        # Note sure how robust this is
        global openmp_support
        compiler = ccompiler.new_compiler()
        openmp_support = has_flag(compiler, '-fopenmp')

        shared = glob.glob(f"{root}/build/lib*/pywfa/*.so") + glob.glob(f"{root}/build/lib*/pywfa/*.dll")
        _build_ext.build_extension(self, ext)
        [shutil.copy(i, f"{root}/pywfa") for i in shared]


class Build_ext_first(_install):
    def run(self):
        # Build WFA2-lib first
        # run("cd pywfa/WFA2_lib; make clean all; cd ../../", shell=True)
        global openmp_support
        if openmp_support:
            run("cd pywfa/WFA2_lib; make clean all BUILD_WFA_PARALLEL=1 BUILD_MINIMAL=1; cd ../../", shell=True)
        else:
            run("pwd", shell=True)
            run("cd pywfa/WFA2_lib; make clean all BUILD_WFA_PARALLEL=0 BUILD_MINIMAL=1; cd ../../", shell=True)
        return setuptools.command.install.install.run(self)

setup(
    name="pywfa",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    url="https://github.com/kcleal/pywfa",
    description="Align sequences using WFA2-lib",
    license="MIT",
    version='0.4.0',
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
    ext_modules=cythonize(ext_modules, **cy_options),
    cmdclass={'install': Build_ext_first, 'build_ext': MyBuild},
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
)
