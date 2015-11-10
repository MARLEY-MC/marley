// Add C++14 support (currently non-default) to ACLiC
void aclic_cpp14() {
  TString o;
  // Customize MakeSharedLib
  o = TString(gSystem->GetMakeSharedLib());
  o = o.ReplaceAll(" -c ", " -std=c++14 -c ");
  gSystem->SetMakeSharedLib(o.Data());
  // Customize MakeExe
  o = TString(gSystem->GetMakeExe());
  o = o.ReplaceAll(" -c ", " -std=c++14 -c ");
  gSystem->SetMakeExe(o.Data());
  // Customize interpreter's and compiler's flags
  // gSystem->AddIncludePath("-D__WHATEVER_YOU_NEED__");
}

void load_marley_dict() {
  // Instruct ACLiC to use C++14 compiler flags
  aclic_cpp14();

  // Load dictionaries for MARLEY classes that
  // can be written to ROOT files
  gROOT->ProcessLine(".L TMarleyParticle.cc+");
  gROOT->ProcessLine(".L TMarleyEvent.cc+");
}

void rootlogon() {
  std::cout << std::endl << "Creating temporary MARLEY class dictionaries ..." << std::endl;
  load_marley_dict();
}
