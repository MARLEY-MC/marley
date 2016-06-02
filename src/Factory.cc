#include "Factory.hh"

std::unique_ptr<marley::Generator> marley::Factory::make_generator(
  const std::string& file_name)
{
  cf = std::make_unique<marley::ConfigFile>(file_name);
  return std::make_unique<marley::Generator>(*cf);
}
