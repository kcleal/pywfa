from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import setup
from setuptools.extension import Extension
from subprocess import run
import os

# run("cd pywfa/WFA2-lib; make BUILD_TOOLS=0 BUILD_EXAMPLES=0 clean all; cd ../../", shell=True)

extras = ["-Wno-sign-compare", "-Wno-unused-function", "-Wno-unused-result", '-Wno-ignored-qualifiers',
          "-Wno-deprecated-declarations"]

ext_modules = list()

root = os.path.abspath(os.path.dirname(__file__))
wfa = os.path.join(root, "pywfa/WFA2-lib")

libraries = [f"{wfa}/lib"]
library_dirs = [f"{wfa}/lib"]
include_dirs = [root, wfa, f"{wfa}/lib", f"{wfa}/utils", f"{wfa}/wavefront", f"{wfa}/bindings/cpp", f"{wfa}/system",
                ]

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

setup(
    name="pywfa",
    author="Kez Cleal",
    author_email="clealk@cardiff.ac.uk",
    url="https://github.com/kcleal/pywfa",
    description="Align sequences using WFA2",
    license="MIT",
    version='0.2.0',
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
    packages=["pywfa", "tests"],
    ext_modules=cythonize(ext_modules),
    cmdclass={'build_ext': build_ext},
    include_package_data=True,
    zip_safe=False,
    test_suite='nose.collector',
)
