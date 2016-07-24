// Ensure that ROOT (via rootcint) generates all of the dictionary entries
// needed to write MARLEY events and particles to a ROOT file
#ifdef __MAKECINT__
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ defined_in "marley/Event.hh";
#pragma link C++ defined_in "marley/Particle.hh";
#pragma link C++ class std::vector<marley::Particle*>+;
#pragma link C++ class std::vector<marley::Particle>::*+;
#pragma link C++ operators std::vector<marley::Particle>::iterator;
#pragma link C++ operators std::vector<marley::Particle>::const_iterator;
#pragma link C++ operators std::vector<marley::Particle>::reverse_iterator;
#endif
