Recommended nuclear structure data for MARLEY
=============================================

Revised 18 November 2018 by S. Gardiner

Introduction
------------

This directory contains the nuclear structure data recommended for use with the
MARLEY event generator. The data for potassium-40 are based on the 2017 ENSDF
evaluation for A=40 (see https://bit.ly/2Y4nPdg for details). Data for all
other nuclides included in this collection are taken from files that were
originally produced by A. Koning, S. Hilaire, and M. Duijvestijn for use with
the nuclear code TALYS-1.6. The original versions (in a different format)
appear in the TALYS tar.gz file under the directory
talys/structure/levels/exp/. If you use these data in published calculations,
please give proper attribution to the TALYS developers. For information about
TALYS, please see http://www.talys.eu.

License
-------

These data are distributed under the terms of the GNU General Public License (see
http://tinyurl.com/zssovzw).

How to use these data
---------------------

Each data file is named zXXX, where XXX is a 3-digit number that corresponds to
the atomic number of all nuclides in the file. Many different mass numbers are
represented in these files.

To use the data within MARLEY, you will need to include paths to each of these
files in the "structure" JSON array of your configuration file. For a working
example, please see the NUCLEAR STRUCTURE DATA section in examples/config.js.
