lamg_cpp
========

The following is a C++ porting of the Lean Algebraic Multigrid solver by O.E.Livne and A.Brandt, which was written in Matlab and available for download [here](https://code.google.com/p/lamg/), written using the new Smart Expression Templates library [Blaze-lib](https://code.google.com/p/blaze-lib/). The main purposes of this porting was to test this solution as an internal solver for the Interior Point algorithm's problem developed by "The Operation Research Group" of University of Pisa, exploiting the main ideas introduced by the authors in this [arxiv article](http://arxiv.org/abs/1108.1310).

__The code is still in a beta version and the results are very preliminary so don't expect it to work correctly__.

In particular you'll have to manually download all the required headers (MCFClass.h, IPClass.h, MCFLSSolver.h, etc.)  from the [ORGroup website](http://www.di.unipi.it/optimize/Software/MCF.html) in order to correctly compile the code.

The docs folder contains the source code of my thesis (written in Italian). To compile it all you need to do is install a basic LaTeX distribution and type
```
$ pdflatex Tesi.tex
```
in a command prompt.
