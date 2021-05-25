# Modular Method SageMath Package
This package provides tools to work with Frey curves and Frey Q-curves
and apply the modular method to them.

## Requirements
This package requires [SageMath](http://ww.sagemath.org/). It was last
tested on version 9.1, but should work on version 9.0 or newer.

The package has the option to do newform computations in
[MAGMA](http://magma.maths.usyd.edu.au/magma) if it is installed. It
was last tested on version 2.24-8. MAGMA should be available through
[SageMath's MAGMA interface](http://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/magma.html)

## Directory layout

* **src** directory containing all sage source files

* **modular_method** directory containing the python package. Should
  be compiled first from the source, see "How to use" : "Setup"

* **examples** A directory containing various examples for the
  framework provided by the modular method package. Each example is
  written in a restructured text format, which contains text
  explanation of what has to be computed and sage examples of how it
  should be computed. These files can be automatically tested using
  SageMath's automated doctests, see "How to use" : "Testing".

    * **examples/literature** 

    * **examples/thesis** A directory containing examples from Joey van
      Langen's PhD thesis.

* **compile.sh** A script to compile the sage source code to a python
  package, see [Setup](#Setup) below.
  
* **test.sh** A script to easily use the SageMath automated doctest
  system with the correct settings for this package, see
  [Testing](#Testing) below.

* **load.sage** A sage code file to easily load important
  functionality from the package
  
## How to use

### Setup
To be able to use the package, you should first compile all the sage
source files into a python package using the compile script inside the
top directory (the directory this README is in), e.g.

    > ./compile.sh

All source files are now compiled into a directory called
*modular_method*. You can move this directory wherever you like, but
it should be in your current working directory or python path in order
to use it as described below.

Whenever you update the source files you can update the python package
by running the compile script again.

### Testing
Using the SageMath automated doctesting you can verify all the code is
working correctly. This is completely optional. The testing script
allows you to easily test code and examples that are part of the
package.

*All testing should be done in a terminal in the top level directory,
i.e. the directory this README file is in*

For example to test the entire python package you could type

    >test.sh modular_method

You can also test a single file

    >test.sh modular_method/elliptic_curves/Qcurves.py
	
It is also possible to test the entire database of examples

	>test.sh database/
	
or test a single one of them

	>test.sh database/Dieulefait-Freitas-2014.rst
	
### Using in a sage terminal
To use the package in a sage terminal you can simply import the
package

    sage: import modular_method
	
Note that to use classes and functions of the package as global
variables you have to import them as global variables, e.g.

	sage: from modular_method.elliptic_curves.frey_curves import FreyCurve
	
To directly import the most commonly used classes and functions of the
package as global variables you can import * from modular_method

	sage: from modular_method import *

