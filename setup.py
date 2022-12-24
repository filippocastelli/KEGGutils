import os
from setuptools import setup
import re, pathlib

# view filename
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
    return result.group(1)

def get_requirements():
    reqpath = pathlib.Path("./requirements.txt")
    reqs = reqpath.read_text().splitlines()
    return reqs


project_name = "KEGGutils"

setup(
    name = project_name,
    version = get_property('__version__', project_name),
    python_reqyuires = '>3.8',
    author = "Filippo Castelli",
    author_email = "filippocastelli42@gmail.com",
    description = ("Simple utils to work with KEGG data on NetworkX"),
    license = "Unlicense",
    keywords = "KEGG networkx utils",
    url = "https://github.com/filippocastelli/KEGGutils",
    packages=['KEGGutils'],
    install_requires= get_requirements(),
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: Public Domain",
    ],
)