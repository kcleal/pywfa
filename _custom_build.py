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

root = os.path.abspath(os.path.dirname(__file__))
wfa = os.path.join(root, "pywfa/WFA2_lib")


def build_wfa2_lib():
    """Compile the vendored WFA2-lib (runs `make`, not needed for sdist/metadata builds)."""
    # Can't get openmp to behave. Whenever its on (1), the build fails.
    # this has happened on multiple OS with/without `libomp-dev`
    # compiler = ccompiler.new_compiler()
    omp = 0 #1 if has_flag(compiler, "-fopenmp") else 0
    # WFA2-lib's Makefile builds with -march=native, which ties the binary to the build CPU.
    # Set PYWFA_MARCH (e.g. x86-64-v3) to target a portable CPU level instead.
    march = os.environ.get("PYWFA_MARCH")
    cc_flags = f'"CC_FLAGS=-Wall -g -fPIC -O3 -march={march}"' if march else ""
    ret = run(f"cd pywfa/WFA2_lib; make clean all BUILD_WFA_PARALLEL={omp} BUILD_TOOLS=0 BUILD_EXAMPLES=0 {cc_flags}",
              shell=True)
    if ret.returncode != 0:
        raise RuntimeError(
            f"Unable to build WFA2_lib (make exited with code {ret.returncode}). "
            "A working C compiler and GNU make are required to build pywfa from source."
        )
    if not os.path.exists(f"{wfa}/lib/libwfa.a"):
        raise RuntimeError("WFA2_lib build did not produce lib/libwfa.a")


def make_ext_modules():
    """Cythonize the pywfa.align extension against the compiled WFA2-lib."""
    extras = ["-Wno-unused-function", "-Wno-unused-result",
              "-Wno-ignored-qualifiers", "-Wno-deprecated-declarations"]
    extras_pywfa = get_extra_args(extras)

    libraries = [f"{wfa}/lib"]
    library_dirs = [f"{wfa}/lib"]
    include_dirs = [".", root, wfa, f"{wfa}/lib", f"{wfa}/utils", f"{wfa}/wavefront",
                    f"{wfa}/bindings/cpp", f"{wfa}/system", f"{wfa}/alignment", f"{root}/pywfa"]

    return cythonize(Extension("pywfa.align",
                               ["pywfa/align.pyx"],
                               libraries=libraries,
                               library_dirs=library_dirs,
                               include_dirs=include_dirs,
                               extra_compile_args=extras_pywfa,
                               language="c",
                               extra_objects=[f"{wfa}/lib/libwfa.a"]
                               ),
                     **cy_options)


###################
# Basic build_ext #
###################
# Thanks to https://stackoverflow.com/questions/73800736/pyproject-toml-and-cython-extension-module
class build_py(_build_py):
    def run(self):
        build_wfa2_lib()
        if self.distribution.ext_modules is None:
            self.distribution.ext_modules = []
        self.distribution.ext_modules.extend(make_ext_modules())
        self.run_command("build_ext")
        shared = glob.glob(f"{root}/build/lib*/pywfa/*.so") + glob.glob(f"{root}/build/lib*/pywfa/*.dll")
        [shutil.copy(i, f"{root}/pywfa") for i in shared]
        return super().run()
