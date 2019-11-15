# Overview
*pixsim* is a pure-Python package with a command line interface. The requirements are:

- numpy/matplotlib/
- BEM++
- click
- GMSH/pygmsh
- sqlalchemy/sqlite3

## Install prerequisites
```
$ sudo apt install gmsh python-sqlalchemy sqlite3 python-virtualenv \
                   python-matplotlib python-numpy
$ sudo apt install paraview paraview-python mayavi2
$ virtualenv --system-site-packages venv
$ source venv/bin/activate
$ pip install pygmsh meshio
```

## Install BEM++
See BEM++ install instructions *pixsim* needs the development version, so see “Building BEM++ from scratch”. Don’t forget to check out the development branch.
```
$ git clone https://github.com/bempp/bempp.git
$ cd bempp/
$ git checkout -b development origin/development
```
On Ubuntu, you will need, at least:
```
$ sudo apt install gmsh gmsh-doc libaec0 libann0 \
                   libfltk-gl1.3 libfltk-images1.3 \
                   libfltk1.3 libgl2ps0 libhdf5-openmpi-10 \
                   libmed1v5 liboce-foundation10 \
                   liboce-modeling10 libsz2 libtet1.5 \
                   python-decorator (python-scipy?)
```
This list is based on what is pulled in by the .deb package for the production version. Having these additional packages installed helps to reduce the amount of packages to compile:
```
$ sudo apt install libtbb-dev libtbb2 cython python-sphinx
```
Now, build:
```
$ mkdir bempp-build && cd bempp-build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/to/opt ../bempp
$ make -j8 && make install
```

## Install pixsim
It’s recommended to install in a virtualenv.
```
$ virtualenv --system-site-packages vevn
$ source venv/bin/activiate
$ git clone ...
$ cd larf
$ python setup.py install
$ pixsim --help
```