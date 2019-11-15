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
                   libfltk1.3 libgl2ps1.4 libhdf5-openmpi-100 \
                   libmed1v5 liboce-foundation11 \
                   liboce-modeling11 libsz2 libtet1.5 \
                   python-decorator python-scipy \
                   pck-config zlib1g-dev
```
This list is based on what is pulled in by the .deb package for the production version. Having these additional packages installed helps to reduce the amount of packages to compile:
```
$ sudo apt install libtbb-dev libtbb2 cython python-sphinx
```
Now, build:
```
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/to/bempp_install ../
$ make -j8 && make install
```
If you run into problems with the [install](http://roybijster.nl/2018/08/installing-bem-from-source/), you may need to install boost
```
$ sudo apt-get install libboost-all-dev
```
and Patchelf
```
$ cd ~
$ wget http://nixos.org/releases/patchelf/patchelf-0.8/patchelf-0.8.tar.bz2
$ tar xf patchelf-0.8.tar.bz2
$ cd patchelf-0.8/
$ ./configure --prefix="$HOME/.local"
$ make install
$ strip ~/.local/bin/patchelf
$ gzip -9 ~/.local/share/man/man1/patchelf.1
$ export PATH=$HOME/.local/bin:$PATH
```

## Install pixsim
It’s recommended to install in a virtualenv.
```
$ virtualenv --system-site-packages venv
$ source venv/bin/activate
$ pip install pygmsh meshio
$ git clone https://github.com/hcsullivan12/pixsim.git
$ cd pixsim
$ python setup.py install
$ export PYTHONPATH=/path/to/bempp-install/lib/python2.7/site-packages/
$ export LD_LIBRARY_PATH=/path/to/bempp/lib
$ pixsim --help
```