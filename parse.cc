#include <iostream>
#include <string>
#include <vector>
#include "TMarleyDecayScheme.hh"

int main() {

  std::vector<std::string> nuc_ids = { " 40AR",  " 40K ",
    " 40CA", " 40SC", " 40TI" };
  std::string filename = "ensdf.040";

  //// Create a decay scheme object to store data
  //// imported from the ENSDF file
  //TMarleyDecayScheme decay_scheme(nuc_id, filename);

  // Print a report describing the decay scheme
  //decay_scheme.print_report();
  //std::cout << std::endl;

  // Test the decay scheme object by simulating a sample gamma cascade
  //decay_scheme.do_cascade(6000);
  
  //decay_scheme.print_latex_table();

  //std::ofstream ofs;
  
  for (std::vector<std::string>::iterator i = nuc_ids.begin();
    i != nuc_ids.end(); ++i)
  {
    TMarleyDecayScheme ds = TMarleyDecayScheme((*i), filename); 
    //ofs.open(marley_utils::trim_copy((*i)) + ".txt", std::ofstream::trunc);
    ds.print_report();//(ofs);
    std::cout << std::endl;
    ds.do_cascade(6000);
    std::cout << std::endl;
    //ofs.close();
  }
}
