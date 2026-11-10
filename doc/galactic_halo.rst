.. _galactic-halo:

Shared Galactic halo
====================

Set ``rhoDM`` in ``maddm_card.dat`` to choose the local dark matter
density, shared by direct detection and the Milky Way signals (default
:math:`0.3\,\mathrm{GeV\,cm^{-3}}`). At fixed halo shape, recoil rates
scale linearly with density and Galactic annihilation signals
quadratically. The original density used to publish a reference is
metadata; it does not override the card. The published direct-detection
upper limits, given at :math:`0.3\,\mathrm{GeV\,cm^{-3}}`, are rescaled by
``0.3/rhoDM``; dwarf J factors retain their original halo assumptions.

The default velocity distribution is the Standard Halo Model (SHM), a
Galactic Maxwellian truncated at ``vescape = 544`` km/s, with scale
``vMP = 238`` km/s. ``vMP`` is the mode of the untruncated Maxwellian
speed distribution. ``vMP`` and ``vescape`` are shared by the RAPIDD
recoil rates and the Milky Way annihilation. ``halo_velocity_profile``
selects the distribution for the Milky Way annihilation only: the direct
detection rates always use the SHM, and MadDM prints a warning when
another profile is selected together with direct detection.

.. list-table::
   :header-rows: 1

   * - Profile
     - Galactic velocity distribution
   * - ``SHM``
     - Native Maxwellian with a hard escape cutoff.
   * - ``SHM_numeric``
     - The same distribution, integrated numerically.
   * - ``SHM_beta``
     - Subtract a constant from the Maxwellian below the cutoff and
       renormalize. ``halo_beta = 1`` is the default; 0 recovers SHM.
   * - ``Lisanti``
     - Smoothly suppressed high-speed tail, with ``halo_lisanti_k = 1.5``
       by default.
   * - ``SHM_wLMC``
     - Add a second boosted, truncated Maxwellian, configured by
       ``halo_lmc_fraction``, ``halo_lmc_vb``, ``halo_lmc_sigb``,
       ``halo_lmc_vcut``, ``halo_lmc_cosb``, and ``halo_lmc_azimuth_deg``.

Below the escape speed, the unnormalized radial phase-space factors for
``SHM_beta`` and ``Lisanti`` are

.. math::

   f_\beta[v] &\propto e^{-v^2/v_0^2}-\beta e^{-v_{\rm esc}^2/v_0^2}, \\
   f_{\rm Lisanti}[v] &\propto
   \left[e^{-v^2/[k v_0^2]}-e^{-v_{\rm esc}^2/[k v_0^2]}\right]^k.

Both vanish outside the escape sphere. For the LMC option, the local
fraction defaults to 0.006, the stream speed to 570 km/s, the Maxwellian
scale to 100 km/s and its cutoff to 200 km/s in the stream rest frame.
The cosine (default −0.71) and azimuth define the stream direction
relative to the RAPIDD reference laboratory direction. Zero azimuth
points along the projection of the Galactic radial axis into the
transverse plane.

For annihilation, MadDM draws two independent velocities
from the same Galactic distribution and integrates the cross section
times their relative speed:

.. math::

   \langle\sigma u\rangle_f =
     \int\mathrm{d}u\,p_{\rm rel}[u],\sigma_f[u],u,
   \qquad u=|\mathbf v_1-\mathbf v_2|.

By default, ``mw_shm_npts = 16`` evaluates this integral as a weighted
sum at 16 relative speeds, with weights that sum to one. The two
incoming particles have speed :math:`u/2` in their centre-of-mass
frame. Increase the number of points and check convergence near
thresholds or narrow resonances. The LMC mixture has smooth--smooth,
smooth--LMC and LMC--LMC pair weights :math:`[1-f]^2`,
:math:`2f[1-f]` and :math:`f^2`. The bulk drift contributes to the
cross pairs but cancels in same-component pairs. The default
``mw_velocity_model = shm`` enables distribution averaging for all
the profiles above; ``fixed`` is the legacy single-speed option for
the native SHM. Indirect detection computes the velocity integral
without requiring an installed RAPIDD library. Dwarf analyses retain
their separate velocity and J-factor inputs.

Indirect detection runs record the choice in ``Halo_velocity.json``
and, when detailed outputs are retained, write the Galactic
relative-speed rule (``MW_relative_speed.dat``). A local velocity law is used for the
Galactic calculations; it does not model radial velocity evolution
or recompute cosmic-ray transport.
