New paper, DEAP-3600, and COHERENT CsI[Na]
==========================================

.. feed-entry::
   :date: 2021-11-04 17:30
   :author: Steven Gardiner

Since my last news post, the MARLEY v1.2.0 `implementation paper
<https://doi.org/10.1016/j.cpc.2021.108123>`__ was published by Computer
Physics Communications. The `recommended citations <../citing.html>`__ to use
for MARLEY have been updated accordingly.

.. cut::

New MARLEY-based simulations of |B8| solar neutrino events in the `DEAP-3600
<http://deap3600.ca/>`__ liquid argon dark matter experiment were reported in a
recent `poster
<https://indico.physics.ucsd.edu/event/1/contributions/73/attachments/36/55/LIDINE_poster2021_AE-2.pdf>`__
presented by Andrew Erlandson at `LIDINE 2021
<https://indico.physics.ucsd.edu/event/1/overview>`__. By combining MARLEY with
`RAT-PAC <https://rat-pac.readthedocs.io/en/latest/index.html>`__ to simulate
the detector response and apply a realistic event selection, the DEAP-3600
collaboration obtains an estimated rate of 7.34 ± 0.66 charged-current solar
|ve| events in a 7.20 tonne-year exposure assuming 100% acceptance. If
follow-up studies of backgrounds and acceptance conclude favorably, there could
be an exciting opportunity for DEAP-3600 to achieve a first detection of these
neutrinos in liquid argon.

The `COHERENT <https://sites.duke.edu/coherent/>`__ collaboration also reported
an interesting MARLEY calculation in a recent `preprint
<https://arxiv.org/abs/2110.07730>`__ describing a measurement of coherent
elastic neutrino-nucleus scattering (CEvNS) using a CsI[Na] detector. To
mitigate γ-ray backgrounds, the CsI[Na] crystal is shielded with low-activity
lead. Neutrino reactions occurring within the lead can eject neutrons,
potentially leading to spurious CEvNS-like signals in the detector. While a
data-driven estimate of the rate of these neutrino-induced neutrons (NINs) was
used in the CEvNS cross-section analysis, the NIN energy spectrum was
approximated using the results of a MARLEY simulation. No official reaction
input files have been released yet for |Pb208| (or any of the other Pb
isotopes), but user-defined ones may be provided (see appendix B of the MARLEY
v1.2.0 `paper <https://doi.org/10.1016/j.cpc.2021.108123>`__), as was done in
this analysis.

It continues to be fun to learn about new and exciting applications of MARLEY.
Please keep them coming!

.. |B8| raw:: html

   <sup>8</sup>B

.. |ve| raw:: html

   &nu;<sub>e</sub>

.. |Pb208| raw:: html

   <sup>208</sup>Pb
