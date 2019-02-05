# cu-astr2600
A collection of toys and tools for ASTR2600: Introduction to Scientific Computing at the University of Colorado Boulder.

## Installation
If you're working on the `scorpius` computers at SBO, you should be able to access these by default, or if necessary by running `source henrietta` from a UNIX prompt.

If you're working on your own computer, you should be able to install this by running
```
pip install cu-astr2600
```
from a UNIX prompt. As we add more code to this package, you may need to upgrade the version that is installed on your computer. To upgrade, run
```
pip install cu-astr2600 --upgrade
```
instead. That will make sure it always grabs the latest version of the code.

If you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```
git clone https://github.com/zkbt/cu-astr2600.git
cd cu-astr2600/
pip install -e .
```
This will link the installed version of the `astr2600` package to your local repository. Changes you make to the code in the repository should be reflected in the version Python sees when it tries to `import astr2600`.

## Contributors

This package was written by [Zach Berta-Thompson](https://github.com/zkbt).
