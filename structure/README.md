Recommended nuclear structure data for MARLEY
=============================================

Prepared 22 July 2016 by S. Gardiner

Introduction
------------

This directory contains the nuclear structure data recommended for use with the
MARLEY event generator. The files were originally produced by A. Koning, S.
Hilaire, and M. Duijvestijn for use with the nuclear code TALYS-1.6. They
appear in the TALYS tar.gz file under the directory
talys/structure/levels/exp/. If you use these data in published calculations,
please give proper attribution to the TALYS developers. For information about
TALYS, please see http://www.talys.eu.

License
-------

This subset of the discrete level data from the TALYS-1.6 nuclear code is
distributed under the terms of the GNU General Public License (see
http://tinyurl.com/zssovzw for details).

Even though this dataset (everything that appears in the structure/ folder of
the MARLEY git repository) is licensed under the GPL, the MARLEY source code
and other data distributed with MARLEY (e.g., reaction data files in the react/
folder) retain their BSD license (see the LICENSE file for details). Since
including these structure data with MARLEY merely aggregates a GPL-licensed
work with a separate, non-GPL-licensed work in the same distribution medium,
the GPL does not apply to the aggregate as a whole. For more information,
please see the text of the GNU General Public License, the "mere aggregation"
GPL-FAQ entry (http://www.gnu.org/licenses/gpl-faq.en.html#MereAggregation),
and a similar precedent involving GPL-licensed spell checking dictionaries for
Apache OpenOffice (http://tinyurl.com/j3gbr9m).

How to use these data
---------------------

Each data file is named zXXX, where XXX is a 3-digit number that corresponds to
the atomic number of all nuclides in the file. Many different mass numbers are
represented in these files.

To use these data in MARLEY, include the following lines in your run
configuration file:

    structure z019 talys all
    structure z018 talys all
    structure z017 talys all

If you prefer to use only a subset of these data, individual nuclides
may be selected on configuration file lines like this:

    structure z018 talys 37Ar 38Ar 39Ar 40Ar

Both ENSDF-style "nucids" (mass number + element symbol as a string) and
"ZAs" (integer codes equal to Z * 1000 + A) are allowed as nuclide identifiers.

If you run MARLEY from a directory other than the one where you're storing
these files, add the full data file path to the second entry in each line,
e.g., "structure /path/to/z019 talys ..."
