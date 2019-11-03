Installation
============

There are two ways to install :code:`selectionfunctions`. If you are familiar with the installation of the :code:`dustmaps` module,
then this will be familiar.


1. Using :code:`pip`
--------------------

From the commandline, run

.. code-block :: bash

    pip install selectionfunctions

You may have to use :code:`sudo`.

Next, we'll configure the package and download the dust maps we'll want to use.
Start up a python interpreter and type:

.. code-block :: python
    
    from selectionfunctions.config import config
    config['data_dir'] = '/path/to/store/maps/in'
    
    import selectionfunctions.boubert_everall_2019
    selectionfunctions.boubert_everall_2019.fetch()

All the selection functions should now be in the path you gave to
:code:`config['data_dir']`. Note that these selection functions can be very large - some
are several Gigabytes! Only download those you think you'll need.


2. Using :code:`setup.py`
-------------------------

An alternative way to download :code:`selectionfunctions`, if you don't want to use
:code:`pip`, is to download or clone the respository from
https://https://github.com/DouglasBoubert/selectionfunctions.


In this case, you will have to manually make sure that the dependencies are
satisfied:

* :code:`numpy`
* :code:`scipy`
* :code:`astropy`
* :code:`h5py`
* :code:`healpy`
* :code:`requests`
* :code:`six`
* :code:`progressbar2`

These packages can typically be installed using the Python package manager,
:code:`pip`.

Once these dependencies are installed, run the following command from the root
directory of the :code:`selectionfunctions` package:

.. code-block :: bash
    
    python setup.py install --large-data-dir=/path/to/store/maps/in

Then, fetch the selection functions you'd like to use. Depending on which selection functions you choose
to download, this step can take up several Gigabytes of disk space. Be careful
to only download those you think you'll need:

.. code-block :: bash
    
    python setup.py fetch --map-name=boubert_everall_2019

That's it!
