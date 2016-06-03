// Ensure that ROOT (via rootcint) generates all of the dictionary entries
// needed to write MARLEY events and particles to a ROOT file
#ifdef __MAKECINT__
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ defined_in "./Event.hh";
#pragma link C++ defined_in "./Particle.hh";
#pragma link C++ class std::list<marley::Particle>+;
#pragma link C++ class std::list<marley::Particle>::*+;
#pragma link C++ operators std::list<marley::Particle>::iterator;
#pragma link C++ operators std::list<marley::Particle>::const_iterator;
#pragma link C++ operators std::list<marley::Particle>::reverse_iterator;
#endif
