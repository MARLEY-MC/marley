# Transmission-Coefficient Deduplication and Batch Refactor

This document describes the changes in this branch relative to the public
MARLEY code. The refactor removes repeated optical-model calculations from the
fragment-continuum decay path and prepares those calculations for a future
parallel or GPU backend without changing the underlying nuclear physics.

## Motivation

`FragmentContinuumExitChannel::differential_width()` previously called the
optical model inside the loop over the final nuclear spin \(J_f\). At a fixed
continuum excitation energy, the inputs to the transmission coefficient
\(T_{\ell j}\) do not change when only \(J_f\) changes. The same expensive
optical-model calculation was therefore repeated for every allowed \(J_f\).

The level density does depend on \(J_f\). This refactor reuses only
\(T_{\ell j}\); it continues to calculate the level density and width
contribution separately for every allowed final spin and parity.

## New Calculation Flow

For each continuum-energy point, MARLEY now:

1. builds one transmission-coefficient request for every required
   \((\ell,j)\) channel;
2. submits the ordered request vector through the optical-model batch
   interface;
3. associates each returned coefficient with the request at the same vector
   index; and
4. reuses that coefficient while evaluating the allowed \(J_f\)-dependent
   level densities and width contributions.

The calculation changes from the equivalent of

```cpp
for (const auto Jf : allowed_final_spins) {
  const double Tlj = solve_optical_model(same_l_and_j_inputs);
  width += Tlj * level_density(Jf);
}
```

to

```cpp
const double Tlj = solve_optical_model(l_and_j_inputs);
for (const auto Jf : allowed_final_spins) {
  width += Tlj * level_density(Jf);
}
```

The order of the \((\ell,j)\) channels, the \(J_f\) loop, the parity assigned
to each contribution, and the order in which width terms are accumulated are
preserved.

## Changes by File

### `include/marley/OpticalModel.hh`

- Adds `TransmissionCoefficientRequest`, a plain aggregate containing all
  inputs required for one transmission-coefficient calculation.
- Adds the virtual
  `OpticalModel::transmission_coefficients(const std::vector<...>&)` API.
- Defines the batch contract: the result at index `i` belongs to the request
  at index `i`.
- Provides a serial default implementation that calls the existing scalar
  `transmission_coefficient()` method for each request. Existing optical-model
  implementations therefore work without defining a custom batch backend.

### `src/ExitChannel.cc`

- Refactors `FragmentContinuumExitChannel::differential_width()` into request
  collection, ordered batch evaluation, and \(J_f\)-dependent accumulation.
- Calculates one \(T_{\ell j}\) for each required \((\ell,j)\) channel at the
  current excitation-energy point.
- Reuses that coefficient across the allowed \(J_f\) values.
- Derives final-state parity from the stored orbital angular momentum
  `request.l`, preserving the original parity selection.
- Checks that the optical model returns exactly one result for every request.
  A backend that violates the ordered batch contract now fails explicitly
  instead of producing incomplete or misaligned width contributions.

The discrete exit-channel path remains scalar because it does not contain the
same repeated-\(J_f\) optical solve.

### `include/marley/KoningDelarocheOpticalModel.hh`

- Introduces the private request-local `WorkingState` structure.
- Updates the Koning-Delaroche helper methods to receive that state explicitly
  instead of reading and writing temporary calculation members on the shared
  optical-model object.

### `src/KoningDelarocheOpticalModel.cc`

- Constructs a `WorkingState` for each scalar transmission-coefficient
  calculation.
- Stores kinematic quantities, target and fragment masses, potential
  parameters, and other temporary values in that state.
- Passes the state through the optical-potential, Schrödinger-equation, and
  S-matrix helper functions.
- Computes the charge-adjusted target mass without mutating the shared model
  object.

This state isolation makes the C++ calculation path suitable for a future
concurrent batch implementation. 

## Performance Scope

The current performance improvement comes from avoiding duplicate
Koning-Delaroche solves across final-spin values. The default batch
implementation is still serial: it organizes independent requests but does not
yet execute them concurrently or on a GPU.

The batch API supplies a stable boundary for a later optimized backend. Such a
backend must preserve request/result ordering and verify the thread safety of
all numerical-library calls used by the optical solver.


## Validation

The refactored code was checked with:

- a strict C++14 production build using GCC 14.2 and
  `-Wall -Wextra -Wpedantic -Werror`;
- focused scalar/batch agreement, batch-ordering, physical-bound, and
  request-state-isolation checks;
- an exact comparison of 60 Koning-Delaroche reference values;
- a fixed-seed comparison of 100 generated events, in which every event record
  was identical between the original and refactored calculations; and
- the combined native and focused test suite, which passed 30 assertions
  across 3 test cases.

## Compatibility

The public scalar optical-model API remains available, and the batch method has
a serial default implementation. However, this is an ABI-changing update:

- `OpticalModel` gains a virtual method; and
- `KoningDelarocheOpticalModel` has a different object layout after its
  temporary members move into `WorkingState`.


