# Modular Method SageMath Package
This package provides tools to work with Frey curves and Frey Q-curves
and apply the modular method to them. This includes:

- Frey (Q-)curves of which one can compute conductors

- (Frey) Q-curves of which one can compute various data such as
  splitting maps, levels and characters of the associated newforms, as
  well as most traces of Frobenius elements under associated Galois
  representations.

- Conditions that allow the parameters of a Frey (Q-)curve to be
  limited in the computations mentioned above, e.g. specifying that
  the parameters are coprime or are congruent to some value.

- Elimination strategies to eliminate newforms associated with one or
  multiple Frey (Q-)curves automatically.

An in depth explanation of the features provided in this package can
be found in Joey van Langen's PhD thesis (link will be added when
published by the VU university library).

## Requirements
This package requires [SageMath](http://ww.sagemath.org/). It was last
tested on version 9.1, but should work on version 9.0 or newer.

The package has the option to do newform computations in
[Magma](http://magma.maths.usyd.edu.au/magma) if it is installed. It
was last tested on version 2.25-7. Magma should be available through
[SageMath's Magma interface](http://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/magma.html)

## Installation

The easiest way to install the **modular_method** package is to use
the version of pip that comes with your distribution of SageMath. To
do this execute the following in a command line.

    sage -pip install --upgrade git+https://github.com/jmvlangen/modular-method-package.git@release

This will install the latest version of the package directly from
GitHub. Note that this does not download the files in the **examples**
directory for you.

### How to use

To use the package in a SageMath terminal you can simply import the
package after installation

    sage: import modular_method

Note that to use classes and functions of the package as global
variables you have to import them as global variables, e.g.

	sage: from modular_method.elliptic_curves.frey_curves import FreyCurve

To directly import the most commonly used classes and functions of the
package as global variables you can use

	sage: from modular_method import *

To see how most classes and functions can be used, have a look at the
examples in the **examples** directory.

### Source code

To get the raw code instead of the installation you can simply clone
the master branch of the GitHub repository.

    git clone https://github.com/jmvlangen/modular-method-package.git

If you want to turn this source code into a Python package, you can
use the compile scripts `compile.sh` and `compile-mac.sh`. The compile
script requires both *SageMath* and the *sed* command to be
installed. Use the compile script inside the top directory (the
directory this README is in), e.g.

    ./compile.sh

On Mac systems the standard compile script might give an error. Use

    ./compile-mac.sh

instead. Both scripts have not been tested on Windows based systems.

After using the script all source files are compiled into a directory
called *modular_method*. This is a Python package like the one
installed in [Installation](#Installation), which you can put either
in your Python path or current directory to use.

Whenever you update the source files you can update the Python package
by running the compile script again. Note that the *compile-mac.sh*
script does not keep track of timestamps and will update all files
regardless of whether they were changed.

## Directory layout

* **src** directory containing all SageMath source files

* **modular_method** directory containing the Python package. Should
  be compiled first from the source files, see [Setup](#Setup).

* **examples** A directory containing various examples for the
  framework provided by the modular method package. Most examples are
  written in a ReStructured Text, containing an explanation in text
  interjected with SageMath input and output. The output in each
  ReStructured Text file can be verified using SageMath's automated
  doctests, see [Testing](#Testing).

    * **examples/literature** A directory containing examples from the
      literature.

    * **examples/thesis** A directory containing examples from Joey van
      Langen's PhD thesis.

    * **examples/EDS** A directory containing the examples for the yet
      to be published article that is Chapter 5 of Joey van Langen's
      PhD thesis.

* **compile.sh** A script to compile the SageMath source code to a
  python package, see [Source code](#source-code). Should work on Linux based
  systems.

* **compile-mac.sh** A modified version of **compile.sh** that works
  on Mac based systems, see [Source code](#source-code).

* **test.sh** A script to easily use the SageMath automated doctest
  system with the correct settings for this package, see
  [Testing](#Testing) below.

* **load.sage** A SageMath file to easily load important functionality
  from the package

## Testing
Using the SageMath automated doctesting you can verify:

 - All examples and tests in the documentation of the code

 - The output of all the examples in the *examples* directory

Note that some of the examples take significantly more time and memory
than the SageMath doctesting framework allows by default. It is
therefore advised to run the doctests with the options *-T 0* and *-m
0* to disable the time and memory constraints
respectively. Furthermore most examples require the testing to be
executed with the top level directory (the directory this README file
is in) as the Python path for cross references to work correctly. For
convenience the script *test.sh* has been added that sets all these
options automatically when executed in the top level directory.

For example to test the entire Python package -- after it was compiled
as explained in [Source code](#source-code) -- you could type

    ./test.sh modular_method

You can also test a single file

    ./test.sh modular_method/elliptic_curves/Qcurves.py

It is also possible to test all the examples at once

	./test.sh examples/

or test a single one of them

	./test.sh examples/literature/Dieulefait-Freitas-2014.rst
