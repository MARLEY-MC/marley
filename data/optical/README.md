Nuclear optical model potential inputs for MARLEY
=================================================

Updated 22 October 2026 by S. Gardiner

This directory contains inputs used to initialize the
`marley::KoningDelarocheOpticalModel` class in MARLEY. This class implements
the parameterization of the nuclear optical potential originally developed
by A. J. Koning & J. P. Delaroche in Nucl. Phys. A 713, 231-310 (2003),
https://doi.org/10.1016/S0375-9474(02)01321-0.

Starting with MARLEY v2.0.0, the default treatment of the optical model
parameters comes from the "KDUQ Federal" evaluation reported by
C. D. Pruitt, J. E. Escher, and R. Rahman in Phys. Rev. C 107, 014602 (2023),
https://doi.org/10.1103/PhysRevC.107.014602.

The recommended central values of the model parameters may be found in the file
`kduq_federal_cv.js`. An ensemble of 416 variations of the parameters (in files
named `kduq_federal_<X>.js`, where `<X>` is a three-digit number between `000`
and `415`) is also provided to quantify uncertainties. The central values of
the parameters from the original Koning-Delaroche fit are also tabulated in
`kd_global.js` and match those that were used (hard-coded) in MARLEY v1.
All of these files are slightly reformatted versions of originals taken from
the Supplemental Material included with
https://doi.org/10.1103/PhysRevC.107.014602.

If you use these input files in published calculations, please cite the papers
mentioned above to give proper credit to the data evaluators. For convenience,
BibTeX-format citations for these are provided in the file
[OMP-CITATION.bib](OMP-CITATION.bib).
