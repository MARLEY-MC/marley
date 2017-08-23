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

To use the data within MARLEY, you will need to include paths to each of these
files in the "structure" JSON array of your configuration file. For a working
example, please see the NUCLEAR STRUCTURE DATA section in examples/config.js.
