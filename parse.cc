#include <iostream>
#include <string>
#include <vector>
#include "DecayScheme.hh"
#include "marley_utils.hh"

int main() {

  std::vector<std::string> nuc_ids = { " 40AR",  " 40K ", " 40CA", " 40SC", " 40TI" };
  std::string filename = "ensdf.040";

  for (auto &id : nuc_ids) {
    marley::DecayScheme ds = marley::DecayScheme(id, filename);
    //ofs.open(marley_utils::trim_copy(id) + ".txt", std::ofstream::trunc);
    std::cout << ds; 
    //ds.print_report();//(ofs);
    //std::cout << std::endl;
    //ofs.close();
  }

}
