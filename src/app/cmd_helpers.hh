#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
}

void for_each_event(
  const std::vector< std::string >& input_files,
  std::function< void( HepMC3::GenEvent&, bool,
    double, const std::shared_ptr< HepMC3::GenRunInfo >& ) > callback );
