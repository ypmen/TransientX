==================
Installation
==================

Using Docker is recommended, as it eliminates the need for manual installation. The Docker image for TransientX is available on `Docker Hub <https://hub.docker.com/r/ypmen/pulsarx/tags>`_. You can pull the image using the command: ``docker pull ypmen/pulsarx``.
If you prefer to install TransientX manually, please follow the instructions below, and the Dockerfile in the repository can serve as a reference `Dockerfile <https://github.com/ypmen/PulsarX/blob/main/Dockerfile>`_.

Dependencies
----------------

- boost > 1.56
- `PlotX <https://github.com/ypmen/PlotX>`_
- `XLibs <https://github.com/ypmen/XLibs>`_
- `sofa <https://www.iausofa.org/2020_0721_C/sofa_c-20200721.tar.gz>`_ or `erfa <https://github.com/liberfa/erfa>`_

Installation Steps
----------------------

- ``./bootstrap``
- ``./configure --prefix=[install_path] LDFLAGS="-L/path_to_sofa" CPPFLAGS="-I/path_to_sofa"``
- ``make``
- ``make install``

Bug Reports and Support
----------------------------

If you need support or report bugs, please open an issue on the `TransientX GitHub repository <https://github.com/ypmen/TransientX>`_ or contact me by email ypmen@mpifr-bonn.mpg.de.