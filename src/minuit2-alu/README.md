## Taken from CERN Root 6.10

taken with: git clone http://root.cern.ch/git/root.git
on Oct 2 2017, commit fd317ca82612526ddcb6e9804cb88d3b029a50ec

made minimal modifications to:

- inc/Minuit2/FumiliMaximumLikelihoodFCN.h
  avoid linking to Root / mathcore by importing an inline function

- src/SimplexBuilder.cxx
  include one header file even if DEBUG is not set (otherwise compilation fails)

added make-related files from Minuit-5.24 release to allow standalone
compilation with make

## Cite Root authors
We are [![DOI](https://zenodo.org/badge/10994345.svg)](https://zenodo.org/badge/latestdoi/10994345)

Please cite us as

    Rene Brun and Fons Rademakers, ROOT - An Object Oriented Data Analysis Framework,
    Proceedings AIHENP'96 Workshop, Lausanne, Sep. 1996,
    Nucl. Inst. & Meth. in Phys. Res. A 389 (1997) 81-86.
    See also "ROOT" [software], Release vX.YY/ZZ, dd/mm/yyyy,
    https://doi.org/10.5281/zenodo.848819.

## Live Demo for CERN Users
[![](https://swanserver.web.cern.ch/swanserver/images/badge_swan_white_150.png)](http://cern.ch/swanserver/cgi-bin/go?projurl=https://github.com/cernphsft/rootbinder.git)
