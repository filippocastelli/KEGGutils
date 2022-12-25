# KEGGutils: 

[![DOI](https://zenodo.org/badge/174000695.svg)](https://zenodo.org/badge/latestdoi/174000695)


Working with **KEGG** in **Python** and **NetworkX**
# <img src="/img/logo_cut.png" alt="Drawing" width = "630"></img>



**KEGGutils** is a toolkit designed for working with the *Kyoto Encyclopedia of Genes and Genome* database in *Python* with a quick and easy to use interface: in a single line you can download data from KEGG's REST API, organize in a graph-like format provided by *NetworkX* and immediately start exploring.

*KEGGutils* is much more than just an API interface: other than a series of tools to better interface yourself with the service, it provides expanded *NetworkX* classes and methods ( totally nx-compatible ) to handle different types of data and help you better exploit underlying structures. 
KEGGutils is easily expandable and can be immediately integrated anywhere you use *NetworkX* to process data.

## Current build status:

*master branch* : [![CircleCI](https://circleci.com/gh/filippocastelli/KEGGutils.svg?style=shield)](https://circleci.com/gh/filippocastelli/KEGGutils)

*devel branch* : [![CircleCI](https://circleci.com/gh/filippocastelli/KEGGutils/tree/devel.svg?style=shield)](https://circleci.com/gh/filippocastelli/KEGGutils/tree/devel)

## Installing KEGGutils
# <a href="https://pypi.org/"><img alt = PyPi src="https://pypi.org/static/images/logo-large.72ad8bf1.svg" height="100"></img></a>

You can just clone this repo and use it's content as any other python package, but if you just need a super easy drop-in solution *KEGGutils* is available as a *PyPi* package:

to install it you just need to run

`pip install KEGGutils`
 
 and that should be it!
 
## Dependencies
KEGGutils is tested and working against python `3.8`, `3.9`, `3.10` and `3.11`.

`3.6` and `3.7` are untested but should work.

To make sure KEGGutils works as it should, a few dependencies must be satisfied:
- `networkx`
- `matplotlib`
- `awesome-slugify`
- `requests`
- `Pillow`
- `scipy`

note: if you use pip to install KEGGutils, it should automatically get the needed dependencies for you!

note: you can create an environment with all the required dependencies using the incldued anaconda environment configurator, you just need to run

`conda env create -f keggutils_env.yml` 

to create a `keggutils_env` anaconda environment with all the required dependencies.

## Getting started

In this repo can find four dense yet easy to follow tutorials covering most of *KEGGutils*'s functionalities and more are coming in the next future.
To start just follow the links below: 
- [**Tutorial 0 - Basics and KEGG API**](https://github.com/filippocastelli/KEGGutils/blob/dev/tutorials/Tutorial%200%20-%20Basics%20and%20KEGG%20API.ipynb)
- [**Tutorial 1 - Enzymatic Graphs**](https://github.com/filippocastelli/KEGGutils/blob/dev/tutorials/Tutorial%201%20-%20EnzymeGraphs.ipynb)
- [**Tutorial 2 - KEGGgraphs, KGGlinkgraphs and KEGGchains**](https://github.com/filippocastelli/KEGGutils/blob/dev/tutorials/Tutorial%202%20-%20KEGGgraphs%2C%20KGGlinkgraphs%20and%20KEGGchains.ipynb)
- [**Tutorial 3 - KEGGpathways**](https://github.com/filippocastelli/KEGGutils/blob/dev/tutorials/Tutorial%203%20-%20KEGGpathways.ipynb)


## Contributing

The project is open to contributions and suggestions, just open an issue on the repo or contact me directly (contacts below).

## External links

Here are a few useful links
- [KEGG: Kyoto Encyclopedia of Genes and Genomes](https://www.kegg.jp/)
- [KEGG REST API reference page](https://www.kegg.jp/kegg/rest/keggapi.html)
- [KEGG KGML (KEGG Markup Language) reference page](https://www.kegg.jp/kegg/xml/)
- [Networkx Github IO](https://networkx.github.io/)

## How to Cite

KEGGutils DOI

[![DOI](https://zenodo.org/badge/174000695.svg)](https://zenodo.org/badge/latestdoi/174000695)

BibTeX-style citation:
```
@software{filippo_maria_castelli_2022_7482523,
  author       = {Filippo Maria Castelli},
  title        = {KEGGutils},
  month        = dec,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v0.4.1},
  doi          = {10.5281/zenodo.7482523},
  url          = {https://doi.org/10.5281/zenodo.7482523}
}
```

## Contacts

**Author:**

Filippo Maria Castelli  
castelli@lens.unifi.it  
LENS, European Laboratory for Non-linear Spectroscopy  
Via Nello Carrara 1  
50019 Sesto Fiorentino (FI), Italy


