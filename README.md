# Laughlin-Metropolis
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1161969.svg)](https://doi.org/10.5281/zenodo.1161969)
[![Build Status](https://travis-ci.org/tcompa/Laughlin-Metropolis.svg?branch=master)](https://travis-ci.org/tcompa/Laughlin-Metropolis)

## What is this?
This is a simple python/cython code implementing the Metropolis Monte Carlo
algorithm to sample configurations from the Laughlin wave function with an
arbitrary number of quasiholes.

The program was developed by Tommaso Comparin and Elia Macaluso to produce some
of the results reported in the following paper: *Observing anyonic statistics
via time-of-flight measurements*  [[arXiv:1712.07940
cond-mat.quant-gas](https://arxiv.org/abs/1712.07940), by R. O. Umucalilar, E.
Macaluso, T. Comparin, and I. Carusotto].

If you use this code in a scientific project, please cite [the corresponding
Zenodo entry](https://zenodo.org/record/1161969):
```
@misc{tommaso_comparin_2018_1161969,
  author       = {Tommaso Comparin, Elia Macaluso},
  title        = {{tcompa/Laughlin-Metropolis: Laughlin-Metropolis v1.0}},
  year         = 2018,
  doi          = {10.5281/zenodo.1161969},
  url          = {https://doi.org/10.5281/zenodo.1161969}
}
```

## How to use it?
This code requires the [numpy](http://www.numpy.org/) and
[cython](http://cython.org/) libraries.  It works on python 2.7, 3.4, 3.5 and
3.6 (elementary tests are available in the `Tests` folder, and they are
performed at each commit - see the current status on
https://travis-ci.org/tcompa/Laughlin-Metropolis).
To run the examples, also the [future](https://pypi.python.org/pypi/future) and
[matplotlib](http://matplotlib.org/) libraries.

Before being imported in a python script, the module
`lib_laughlin_metropolis.pyx` has to be compiled through the command

    $ python setup_cython.py build_ext --inplace

After this step, it can be imported in ordinary python scripts.
Have a look at the two example files (in the `Examples` folder):
+ In `example_1_mean_square_radius.py`, the average square radius is computed
in several Monte Carlo runs.
+ In `example_2_density_profile.py`, the density profile is computed and shown,
for parameters corresponding to Fig. 2(b) of the reference above
[Umucalilar et al.].


## Note on physical units

Within our program, the unit of length is defined such that

<img src="https://latex.codecogs.com/gif.latex?l_B&space;\sqrt{2}&space;=&space;1" title="l_B \sqrt{2} = 1" />,

and the (non-normalized) probability distribution sampled with the Monte Carlo algorithm reads

<img src="https://latex.codecogs.com/gif.latex?p\left(\lbrace&space;\mathbf{z}&space;\rbrace,&space;\lbrace&space;\mathbf{\xi}&space;\rbrace\right)&space;=&space;\exp&space;\left\lbrace-&space;\left[&space;\sum_{i=1}^N&space;|z_i|^2&space;-&space;\sum_{i=1}^N&space;\sum_{j=1}^{N_\mathrm{qh}}&space;\log\left(&space;2\left|\mathbf{z}_i&space;-&space;\mathbf{\xi}_j&space;\right|^2&space;\right&space;)&space;-&space;m&space;\sum_{i=1}^N&space;\sum_{j<i}&space;\log\left(&space;2\left|\mathbf{z}_i&space;-&space;\mathbf{z}_j&space;\right|^2&space;\right&space;)&space;\right&space;]&space;\right&space;\rbrace" title="p\left(\lbrace \mathbf{z} \rbrace, \lbrace \mathbf{\xi} \rbrace\right) = \exp \left\lbrace- \left[ \sum_{i=1}^N |z_i|^2 - \sum_{i=1}^N \sum_{j=1}^{N_\mathrm{qh}} \log\left( 2\left|\mathbf{z}_i - \mathbf{\xi}_j \right|^2 \right ) - m \sum_{i=1}^N \sum_{j<i} \log\left( 2\left|\mathbf{z}_i - \mathbf{z}_j \right|^2 \right ) \right ] \right \rbrace" />.

## License
MIT License

Copyright (c) 2018 Tommaso Comparin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
