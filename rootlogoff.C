void rootlogoff() {
  std::cout << std::endl << "Removing temporary MARLEY class dictionaries ..." << std::endl;
  gROOT->ProcessLine(".! rm -f TMarleyParticle_cc.* TMarleyROOTEvent_cc.*");
}
