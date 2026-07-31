=======
Contact
=======

Correspondence addressed specifically to the current core developer(s) of MARLEY
may be sent to support@marleygen.org.

Mailing list
------------

.. |MARLEY-USERS| replace::
   `MARLEY-USERS <https://listserv.fnal.gov/scripts/wa.exe?A0=MARLEY-USERS>`__

.. Prevents Sphinx from automatically hyperlinking this email address
.. |FNAL-LISTSERV| raw:: html

   LISTSERV@FNAL.GOV

Release announcements and occasional MARLEY-related news items are sent to the
|MARLEY-USERS| mailing list hosted by Fermilab. The list also serves as a forum
for general discussion about MARLEY and low-energy neutrino physics more
broadly.

To `subscribe <mailto:LISTSERV@FNAL.GOV?body=SUBSCRIBE%20 MARLEY-USERS %20 your
%20 name>`__, send an email with a blank subject line to
|FNAL-LISTSERV|. The body of the email should contain only a single
line with the text

::

   SUBSCRIBE MARLEY-USERS your name

where ``your name`` should be replaced with your actual name. Subscribers are
allowed to send email messages to the list.

You may also `unsubscribe <mailto:LISTSERV@FNAL.GOV?body=SIGNOFF %20
MARLEY-USERS>`_ from the list by sending an email to |FNAL-LISTSERV|.
The subject line should once again be blank. In this case, the body of the
email should consist of the single line

::

  SIGNOFF MARLEY-USERS

More detailed instructions for users of Fermilab LISTSERV mailing lists are
available `here <https://listserv.fnal.gov/users.asp>`__.

Bug reports
-----------

The preferred method for reporting MARLEY bugs is via the official `issue
tracker <https://github.com/MARLEY-MC/marley/issues>`__ on GitHub. The ``New
issue`` button on the issue tracker webpage may be used by anyone with a (free)
GitHub `account <https://github.com/join>`__ to report a problem.

If you are unsure whether or not you have found a bug, please send an email to
support@marleygen.org and ask for help.

Specificity and detail in bug reports are greatly appreciated. Examples of
useful pieces of information to mention in a bug report include the following:

Bug description
  Describe the problem that you have encountered and the circumstances
  under which it arises. Are there situations in which the bug can be
  avoided? Include any warning or error messages issued by MARLEY.
  Screenshots may be helpful in some cases.

Expected behavior
  Describe how you would expect the code to behave if it were working
  correctly. Sometimes this is obvious (e.g., for bugs that cause MARLEY to
  crash), but physics-related bugs in particular may be subtle.

Minimal working example
  Provide instructions and any appropriate supplemental materials (e.g., the
  job configuration file used for the simulation) that will allow someone else
  to reproduce the problem. Supplemental files may be attached to bug reports
  posted on the GitHub issue tracker. A bug is typically easiest to fix when
  the instructions given to reproduce it are made as simple as possible.

Version information
  Mention the version of MARLEY (``marley --version``) used and the operating
  system (e.g., macOS Catalina 10.15.5) on which it was run. Version information
  for the compiler used to build MARLEY (e.g., Apple clang 11.0.3) is also
  frequently useful.

Stack trace
  For bugs that lead to a crash, including a stack trace in the report may
  facilitate diagnosis of the underlying cause. Users that feel comfortable
  building MARLEY with debugging symbols (``make debug``) and providing this
  information are encouraged to do so.
