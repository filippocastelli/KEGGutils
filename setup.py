import os
from setuptools import setup
import re

# view filename
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()



PKG = "KEGGutils"
VERSIONFILE = os.path.join(PKG, "_version.py")
verstr = "unknown"
try:
    verstrline = open(VERSIONFILE, "rt").read()
except EnvironmentError:
    pass # Okay, there is no version file.
else:
    VSRE = r"^verstr = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        print("unable to find version in %s" % (VERSIONFILE,))
raise RuntimeError("if %s.py exists, it is required to be well-formed" % (VERSIONFILE,))



setup(
    name = "KEGGutils",
    version = "0.1.0",
    author = "Filippo Castelli",
    author_email = "filippocastelli42@gmail.com",
    description = ("Simple utils to work with KEGG data on NetworkX"),
    license = "Unlicense",
    keywords = "KEGG networkx utils",
    url = "https://github.com/filippocastelli/KEGGutils",
    packages=['KEGGutils'],
    install_requires=['networkx', 'requests', 'matplotlib', 'os'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: Public Domain",
    ],
)