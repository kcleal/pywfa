import os
import glob
import shutil
import sysconfig
import setuptools
from subprocess import run
from distutils import ccompiler
from setuptools import Extension
from Cython.Build import cythonize
from setuptools.command.build_py import build_py as _build_py

###########
# Helpers #
###########
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

cfg_vars = sysconfig.get_config_vars()
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

##################
# WFA2_lib build #
##################
extras = ["-Wno-unused-function", "-Wno-unused-result",
          "-Wno-ignored-qualifiers", "-Wno-deprecated-declarations"]
extras_pywfa = get_extra_args(extras)

root = os.path.abspath(os.path.dirname(__file__))
wfa = os.path.join(root, "pywfa/WFA2_lib")
libraries = [f"{wfa}/lib"]
library_dirs = [f"{wfa}/lib"]
include_dirs = [".", root, wfa, f"{wfa}/lib", f"{wfa}/utils", f"{wfa}/wavefront",
                f"{wfa}/bindings/cpp", f"{wfa}/system", f"{wfa}/alignment", f"{root}/pywfa"]

print("Libs", libraries)
print("Library dirs", library_dirs)
print("Include dirs", include_dirs)

# Can't get openmp to behave. Whenever its on (1), the build fails.
# this has happened on multiple OS with/without `libomp-dev`
# compiler = ccompiler.new_compiler()
omp = 0 #1 if has_flag(compiler, "-fopenmp") else 0
ret = run(f"cd pywfa/WFA2_lib; make clean all BUILD_WFA_PARALLEL={omp} BUILD_TOOLS=0 BUILD_EXAMPLES=0",
          shell=True)
if ret.returncode != 0:
    print("Unable to build WFA2_lib")
    print(ret)
    exit(ret.returncode)

##################
# bindings build #
##################
m_ext_module = cythonize(Extension("pywfa.align",
                                ["pywfa/align.pyx"],
                                libraries=libraries,
                                library_dirs=library_dirs,
                                include_dirs=include_dirs,
                                extra_compile_args=extras_pywfa,
                                language="c",
                                extra_objects=[f"{wfa}/lib/libwfa.a"]
                                ),
                        **cy_options)

shared = glob.glob(f"{root}/build/lib*/pywfa/*.so") + glob.glob(f"{root}/build/lib*/pywfa/*.dll")
[shutil.copy(i, f"{root}/pywfa") for i in shared]

###################
# Basic build_ext #
###################
# Thanks to https://stackoverflow.com/questions/73800736/pyproject-toml-and-cython-extension-module
class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()

    def initialize_options(self):
        super().initialize_options()
        if self.distribution.ext_modules == None:
            self.distribution.ext_modules = []

        self.distribution.ext_modules.extend(
            m_ext_module
        )
