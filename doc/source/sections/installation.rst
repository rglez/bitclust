Installation
============

**MDTraj**
----------

Before installing **BitClust**, please be sure that ``mdtraj`` is available
in your system. It is recommended that you install ``mdtraj`` using conda. ::

  $ conda install -c conda-forge mdtraj

``conda`` is an open source package management system and environment management
system that runs on Windows, macOS and Linux. Conda quickly installs, runs and updates
packages and their dependencies.

You can install ``mdtraj`` with ``pip``, if you prefer. ::

  $ pip install mdtraj

``pip`` is the package installer for Python. You can use pip to install packages
from the Python Package Index and other indexes.

Any of them will install ``mdtraj`` along with all dependencies from a
pre-compiled binary. If you don't have Python or the ``conda`` package
manager, you can start with the `Anaconda Scientific Python
distribution <https://store.continuum.io/cshop/anaconda/>`_. ``pip`` homepage
is available `here <https://pip.pypa.io/en/stable/>`_


**BitClust**
------------

Via **pip**:


After succesful installation of ``mdtraj`` you can easily proceed to
install **BitClust** and the rest of its dependencies using pip. ::

  $ pip install bitclust

Then, you should be able to see **BitClust** help by typing in a console: ::

  $ bitclust -h


Via **GitHub**:


After succesful installation of ``mdtraj`` you can easily proceed to
install **BitClust** and the rest of its dependencies from GitHub. ::

  $ git clone https://github.com/rglez/bitclust
  $ cd bitclust
  $ python setup.py install

