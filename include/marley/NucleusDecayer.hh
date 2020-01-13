/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see \${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

namespace marley {

  // Forward-declare some needed classes
  class Event;
  class Generator;
  class Level;
  class MatrixElement;

  /// @brief EventProcessor that handles nuclear de-excitations
  class NucleusDecayer {

    public:

      inline NucleusDecayer() {}

      void deexcite_residue( marley::Event& event, const marley::Level* plevel,
        int twoJ, marley::Parity P, marley::Generator& gen );

      void deexcite_residue( marley::Event& event,
        const marley::MatrixElement& matrix_el, marley::Generator& gen );
  };

}
