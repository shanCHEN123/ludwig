
### Changes

version 0.9.2
- add back an operational surfactant free energy. There is no
  dynamics available yet.
- Add unit tests for the the same.
- Add a description of how to add a free energy to free_energy.h

version 0.9.1
- add option "fe_use_stress_relaxation" in the symmetric case
  to allow force via equilibrium stress.
- Symmetric free energy. Add option "fd_force_divergence 0" which
  allows computation of -phi grad mu in the case where solid is
  present with the approximation that grad mu at boundary is zero.



version 0.8.13
- add option for density output via "rho" commands in input

version 0.8.12
- allow force divergence method to see porous media

version 0.8.11
- updated util/length_from_sk.c to take double input to be consistent
  with extract output

version 0.8.10
- Add vtk format output for composition, velocity, in extract;
  automatically detect input format
- Replace dubious assertion in stats_symmetric_length()

version 0.8.09
- Avoid zero-sized allocations in wall.c

version 0.8.08
- Repaired util/colloid_init.c
- Separate vtk scalar order / director / biaxial order files

version 0.8.07
- Updated the extract program to include vtk headers
- Added some Bond number test examples regression/d3q19/serial-bond-c01.inp

Version 0.8.05
- Added the option to specify fixed liquid crystal anchoring orientation
  from the input.

Version 0.8.04
- Fixed bug in device halo swap.

Version 0.8.0

- The constraint that the number of MPI tasks divide exactly the system
  size in each direction has been relaxed. Logically rectangular, but
  uneven decompositions are computed automatically if required.

- Output format in parallel. Files for a given I/O group now appear with
  data in a format which is independent of parallel decomposition. This
  means an 'extract' step is no longer required. It also means files
  generated with previous versions will no longer work as input.

- Different collision relaxation time schemes are now available by using
  the 'lb_relaxation_scheme' input either ('bgk', 'trt', or 'm10'). The
  default is unchanged ('m10').

- An option for a 'split' treatment of symmetric and antisymmetric stress
  arising from the thermodynamic sector has been introduced. This is
  via the input key 'fe_use_stress_relaxation yes'. This introduces the
  symmetric part of the stress as a relaxation in the collision, while
  the anti-symmetric part is via the force (via divergence). The default
  is still to treat the whole stress via the divergence.

- A 'second active stress' is now available for active liquid crystals.
