'''
This setup.py file sets up our package to be installable on any computer,
so that folks can `import henrietta` from within any directory.

Thanks to this file, you can...

...tell python to look for `henrietta` in the current directory (which you
can continue to edit), by typing *one* of the following commands:

`pip install -e .`
or
`python setup.py develop`

...move a copy of this code to your site-packages directory, where python will
be able to find it (but you won't be able to keep editing it), by typing *one*
of the following commands:

`pip install .`
or
`python setup.py install`

...upload the entire package to the Python Package Index, so that other folks
will be able to install your package via the simple `pip install henrietta`, by
running the following command:

`python setup.py release`


The template for this setup.py came was pieced together with help from
barentsen, christinahedges, timothydmorton, and dfm. Check them out on github
for more neat tricks!

[`python-packaging`](https://python-packaging.readthedocs.io/en/latest/index.html)
is a pretty useful resource too!
'''

# import our basic setup ingredients
from setuptools import setup, find_packages
import os,sys

# running `python setup.py release` from the command line will post to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    # uncomment the next line to test out on test.pypi.com/project/tess-zap
    #os.system("twine upload --repository-url https://test.pypi.org/legacy/ dist/*")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/cu-astr2600*")
    sys.exit()

# a little kludge to get the version number from __version__
exec(open('astr2600/version.py').read())

# run the setup function
setup(
    # people can type `import henrietta` to access this package
    name = "cu-astr2600",
    # this package will only be installed if the current version doesn't exist
    version = __version__,
    # what's a short description of the package?
    description = "A collection of toys and tools for ASTR2600 at CU Boulder.",
    # what's a more detailed description?
    long_description = open('README.md').read(),
    # who's the main author?
    author = "Zach Berta-Thompson",
    # what's the main author's email?
    author_email = "zach.bertathompson@colorado.edu",
    # what's the URL for the repository?
    url = "https://github.com/zkbt/cu-astr2600",
    # this figures out what subdirectories to include
    packages = find_packages(),
    # are the directories of data that should be accessible when installed?
    include_package_data=True,
    # where are those data directories?
    package_data = {'astr2600':['data/*.txt']},
    # any scripts will be copied into your $PATH, so that can run from the command line
    scripts = [],
    # some descriptions about this package (for searchability?)
    classifiers=[
      'Intended Audience :: Education',
      'Intended Audience :: Science/Research',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    # what other packages are required. these must be pip-installable
    install_requires=['numpy',
                      'scipy',
                      'ipython',
                      'matplotlib',
                      'jupyter',
                      'astropy'],

    # the packages in `key` will be installed if folks run `pip install henrietta[key]`
    extras_require={},
    # (I think just leave this set to False)
    zip_safe=False,
    # under what license is this code released?
    license='MIT')
