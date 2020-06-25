.. Redirect to a local path trick taken from
   https://stackoverflow.com/a/37755644/4081973

   This gets around Sphinx's inability to handle relative links
   in the toctree (see https://github.com/sphinx-doc/sphinx/issues/701).
   This comes at the price of a manual redirect.

   An alternative you can consider is putting raw HTML in the toctree
   itself, see, e.g., https://stackoverflow.com/a/61506452/4081973

.. <meta http-equiv="refresh" content="0; url=./_static/doxygen/index.html" />

=======================
Developer documentation
=======================

MARLEY welcomes contributions of all kinds from new developers. This page
provides some basic resources for those who would like to get involved.

Starter guide
-------------

The MARLEY source code is written in `C++14
<https://en.wikipedia.org/wiki/C%2B%2B14>`__ with a small number of helper
scripts written for the `Bash <https://www.gnu.org/software/bash/>`__ shell.
`Git <https://git-scm.com/>`__ is used for version control, and the official
source code repository is hosted on GitHub. For new contributors who are
unfamiliar with Git, `this site
<https://hamwaves.com/collaboration/doc/rypress.com/index.html>`__ provides an
excellent tutorial.

A GitHub user account (available free of charge `here
<https://github.com/join>`__) is required to contribute changes to MARLEY via
the usual workflow. The instructions below assume that you already have an
account set up.

Ordinary MARLEY development is carried out using copies of the official
repository, called *forks*, which are managed by individual developers. After
code changes have been made and tested using a development branch in a fork,
they may be submitted for inclusion in the official repository (and thus in a
future MARLEY release) by means of a GitHub *pull request*. A brief tutorial
illustrating this approach to development on GitHub is available `here
<https://guides.github.com/activities/forking/>`__.

Clear communication is essential when contributing to an open source project
like MARLEY. This `guide
<https://opensource.guide/how-to-contribute/#how-to-submit-a-contribution>`__
presents some great tips on how to communicate most effectively.

Regular contributors to MARLEY may be invited to join the core development team
and be given write access to the official repository. Core developers are
responsible for reviewing and approving pull requests submitted from forks.

.. _code-checkout:

Checking out the source code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A development copy of the MARLEY source code may be obtained according to the
following steps. They need only be executed once.

1. Visit the `webpage <https://github.com/MARLEY-MC/marley>`__ for the official
   repository
2. Press the ``Fork`` button in the upper right corner. If prompted to choose a
   location for the fork, select your user account. A fork of MARLEY will be
   added to your GitHub user account at ``https://github.com/<username>/marley``
   (replace ``<username>`` with your GitHub username here and elsewhere).
3. Clone your GitHub fork to a folder on your local machine::

     git clone https://github.com/<username>/marley


4. Enter the folder for the cloned fork and check that your setup is correct:
   ::

     cd marley
     git remote -v

   The second command above should yield the following output:
   ::

     origin     https://github.com/<username>/marley (fetch)
     origin     https://github.com/<username>/marley (push)

If the steps above were followed successfully, then the MARLEY source code will
now be present in the ``marley`` folder.

MARLEY development guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Having checked out the source code in the manner described above, you may now
make changes, build, and test in the usual way. Frequent Git commits will make
your contributions easier to merge into the official ``master`` branch. Some
useful advice on best practices for making commits is available `here
<https://blog.hartleybrody.com/git-small-teams/>`__.

When preparing code for inclusion in MARLEY, please keep in mind the following
guidelines.

Dependencies
^^^^^^^^^^^^

External dependencies have deliberately been kept to a minimum throughout
MARLEY's history. Beyond the C++ Standard Library itself, the only required
external dependency is the `GNU Scientific Library
<https://www.gnu.org/software/gsl/>`__ (GSL). Even for GSL, only those portions
needed to compute the `Coulomb wavefunctions <https://dlmf.nist.gov/33.2>`__ are
actually used by MARLEY.

New features that impact core MARLEY functionality should avoid introducing new
external dependencies if at all possible. Use of `ROOT <https://root.cern.ch>`__
classes is acceptable only in analysis macros (e.g., ``examples/macros/``) and
in optionally-built portions of the code designed specifically to interface with
ROOT (e.g., ``src/RootEventFileReader.cc``).

Coding style
^^^^^^^^^^^^

MARLEY does not have a formal coding style guide. However, following the list of
suggestions below will help to ensure that new code contributions use a style
roughly consistent with previous work. Readability, clarity, and consistency are
more important than any particular style consideration. Above all, strive for
code that is "correct, beautiful, [and] fast (`in that order
<https://tinyurl.com/correct-beautiful-fast>`__)."

* When in doubt, imitate the conventions used in the existing MARLEY source
  code
* Member variable names end with an underscore character (``_``) 
* Class names and enumeration types are written in `UpperCamelCase
  <https://en.wikipedia.org/wiki/Camel_case>`__
* Variable and function names are written in `snake_case
  <https://en.wikipedia.org/wiki/Snake_case>`__
* Prefer `#pragma once <https://en.wikipedia.org/wiki/Pragma_once>`__ to
  `#include guards <https://en.wikipedia.org/wiki/Include_guard>`__ in header
  files
* Prefer an explicit namespace specifier (``std::``) to a
  `using-directive <https://tinyurl.com/cppref-using-directive>`__
  (``using namespace std``)
* Prefer ``constexpr`` variables to preprocessor
  `macros <https://en.cppreference.com/w/cpp/preprocessor/replace>`__.
  Use ``MACRO_CASE`` for the names of both of these entities.
* Prefer ``'\n'`` to ``std::endl`` when writing newlines to output
  streams (see `explanation <https://tinyurl.com/ccpcore-endl>`__ in the
  CppCoreGuidelines) 
* C++ source files have the filename extension ``.cc``. Header files have the
  filename extension ``.hh``. An exception to the latter rule occurs for
  header files original to another code base (e.g.,
  ``include/fftpack4/fftpack4.h``).
* Use whitespace to improve readability, e.g., ``std::cout << x << '\n';``
  rather than ``std::cout<<x<<'\n';``.
* Prefer a maximum line width of 80 columns. Tolerate modestly longer lines
  if doing so improves code readability.
* Use two spaces for each level of indentation. Never use tabs.
* Prefer multiple ``//`` comments to a single ``/* */`` comment block

As described in the :ref:`meta-doc` section, Doxygen comments should be used in
header files to document class methods and member variables. See `this section
<https://www.doxygen.nl/manual/docblocks.html>`__ of the Doxygen manual for a
description of the special comment blocks used for API documentation.

Pull requests
~~~~~~~~~~~~~

After your changes are finished, tested, committed, and pushed to a branch on
your GitHub fork, you are ready to submit a pull request to officially include
them in MARLEY. Please see the `instructions
<https://tinyurl.com/github-fork-pull-request>`__ for creating a GitHub pull
request if you are unfamiliar with the process. Under most circumstances, you
will want to use the official ``MARLEY-MC/marley`` repository as the "base
repository" and the ``master`` branch as the "base branch."

After your pull request has been created, the current MARLEY core developer(s)
will be notified on GitHub. After reviewing your changes, they may leave
comments on your pull request suggesting modifications. GitHub provides a
flexible `interface <https://tinyurl.com/github-pull-request-comments>`__ for
adding both overall and line-by-line comments to a submitted pull request.
If you do not receive any comments on your pull request within a day or
two, please `reach out <contact.html>`__ to the core developer(s) via email.

As the review process proceeds, you may respond to comments from the core
developer(s) by leaving comments of your own. New commits may be added to the
pull request by pushing them to the same branch on your GitHub fork of MARLEY.

When your changes have been approved, a core developer will merge your
development branch into the ``master`` branch of the official repository.
You have now completed a contribution which will be included in the next
MARLEY release. Congratulations and thanks for your hard work!

.. _meta-doc:

API and meta documentation
--------------------------

As a convenient reference for developers, a set of webpages that provide API
documentation for the MARLEY C++ classes, source files, etc. is available `here
<./_static/doxygen/index.html>`__. These webpages are generated automatically
from the source code using a tool called `Doxygen
<https://www.doxygen.nl/index.html>`__. Special `comment blocks
<https://www.doxygen.nl/manual/docblocks.html>`__ that are written in the MARLEY
header files are interpreted by Doxygen during the generation process.

For offline viewing, the API documentation may generated in any environment
in which Doxygen is installed. To create the HTML files, simply execute
::

   make docs

from within the ``build/`` folder. After Doxygen executes, open the file
``doxygen/html/index.html`` in a browser to view the local copy of the API
documentation website.

With the exception of the API webpages described above, all other content for
the official MARLEY website (https://marleygen.org) is produced from a set of
text files stored in the ``docs/`` folder of the source code tree. These text
files are written in the `reStructuredText
<https://docutils.sourceforge.io/rst.html>`__ (reST) markup language. The
`Sphinx <https://www.sphinx-doc.org/en/master/index.html>`__ documentation
generator is used with the `Guzzle theme
<https://github.com/guzzle/guzzle_sphinx_theme>`__ to produce HTML webpages
from the reST files. Two Sphinx extensions are required to fully build the
website. The `sphinxcontrib-bibtex
<https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`__ extension is used
to handle citations (see, e.g., the online `bibliography <pubs.html>`__). The
`sphinxcontrib-newsfeed <https://pypi.org/project/sphinxcontrib-newsfeed/>`__
is used to manage the posts on the `news webpage <news.html>`__.

Installation of the prerequisites needed to use Sphinx will vary somewhat
across different systems. Typically, however, the standard package manager
may be used to install Sphinx itself, and the remaining components may be
added using `pip <https://pypi.org/project/pip/>`__. For a computer
running macOS and Python 3, for instance, `Homebrew <https://brew.sh/>`__
may be used to install Sphinx and its extensions via the commands
::

  brew install sphinx-doc
  pip3 install guzzle-sphinx-theme sphinxcontrib-bibtex sphinxcontrib-newsfeed

After these components have been installed, one may build the MARLEY webpages
by navigating to the ``docs/`` folder and using the Makefile:
::

  cd docs
  make html

When the build completes, an offline copy of the MARLEY website may be viewed
by opening the file ``docs/_build/html/index.html`` in a browser.

Development wish list
---------------------

An informal list of possible new features that may be added to MARLEY in the
future is given below. `Feedback <contact.html>`__ from the community about the
contents of this list, including suggestions for new items, is welcome.

Physics
~~~~~~~

.. |nuebar| raw:: html

   &#x1d708;&#x304;<sub>e</sub>

* Additional reaction input files

  - New channels for :superscript:`40`\Ar: NC, |nuebar| CC

  - New nuclear targets: :superscript:`12`\C, :superscript:`16`\O,
    :superscript:`56`\Fe, :superscript:`63`\Cu, :superscript:`127`\I,
    :superscript:`208`\Pb, others?

* Implementation of an inclusive cross section model that includes forbidden
  nuclear transitions. A new class derived from ``marley::Reaction`` will
  likely be required.

* Handling of new job configuration file keys to vary the parameters used in the
  nuclear optical model, etc.

  - As an application of the new configuration options, event reweighting could
    be implemented to facilitate assessments of theoretical uncertainties on the
    MARLEY physics models. A prerequisite to the reweighting would be upgrades
    to the ``marley::Event`` class to allow storage of the full de-excitation
    history.

* Refinements of the nuclear de-excitation model

  - Pre-equilibrium particle emission

  - Internal conversion

  - Neutrino-induced fission 

  - Realistic angular distributions for evaporated particles

  - Finite particle emission times

* Non-neutrino projectiles (e.g., electrons, MeV-scale dark matter)

External Interfaces
~~~~~~~~~~~~~~~~~~~

* Interface to `NUISANCE <https://nuisance.hepforge.org/>`__ for comparisons
  to low-energy neutrino scattering data

* Interface to external flux and geometry `drivers
  <https://tinyurl.com/fnal-workshop-flux-geom>`__. This would enable
  simulations of non-uniform detector geometries with a spatially-varying
  neutrino flux.

Documentation and testing
~~~~~~~~~~~~~~~~~~~~~~~~~

* Full Doxygen documentation coverage

* A full suite of unit tests incorporated into the continuous integration
  system
