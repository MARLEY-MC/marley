// Ensure that ROOT (via rootcint) generates all of the dictionary entries
// needed to write MARLEY events and particles to a ROOT file
#ifdef __MAKECINT__
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class marley::Particle+;
#pragma link C++ class std::vector<marley::Particle*>+;
#pragma link C++ class std::vector<marley::Particle>::*+;
#pragma link C++ operators std::vector<marley::Particle>::iterator;
#pragma link C++ operators std::vector<marley::Particle>::const_iterator;
#pragma link C++ operators std::vector<marley::Particle>::reverse_iterator;
#pragma link C++ class marley::Event+;
#endif
