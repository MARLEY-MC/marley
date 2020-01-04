#include "marley/Error.hh"

// Initialize the static member flag indicating whether error messages should
// be printed to the marley::Logger
bool marley::Error::log_them_ = true;
