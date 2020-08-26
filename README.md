# PhysicsUtils
A collection of tools for doing physics problems

## Installation
Clone the project into a folder named `utils`. Then, from any Python3 code, import the package `utils` as you would any other local package.

## Main Features

The whole package is sub-divided into several categories:
- **Number:** Overrides float and complex datatypes with versions that implement equality with tolerance to floating-point error, and print out in scientific notation by default.
- **Quantity:** Implements unit-aware calculations, tracking and converting units, including verifying dimensional compatibility, across calculations automatically. Includes definitions for all SI base units, and all the units I've used so far.
- **Constants:** Definitions for physics constants, appended to as I come across them. All are implemented as Quantities (i.e., have units attached).
- **Computation:** Routines for more complex computations, like numerically-stable Kahan summation and an adaptive integration routine. Compatible with purely numerical or Quantity inputs.
- **Graphics:** Visualization tools, like a utility for making animated line plots.
- **Evolution:** Special-purpose submodule for evolving a wavefunction according to the Schr&ouml;dinger equation (free or under a potential).
