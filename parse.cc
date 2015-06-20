#include <iostream>
#include <string>
#include <vector>
#include "TMarleyDecayScheme.hh"
#include "marley_utils.hh"

int main() {

  std::vector<std::string> nuc_ids = { " 40AR",  " 40K ", " 40CA", " 40SC", " 40TI" };
  std::string filename = "ensdf.040";

  for (std::vector<std::string>::iterator i = nuc_ids.begin();
    i != nuc_ids.end(); ++i)
  {
    TMarleyDecayScheme ds = TMarleyDecayScheme((*i), filename); 
    //ofs.open(marley_utils::trim_copy((*i)) + ".txt", std::ofstream::trunc);
    ds.print_report();//(ofs);
    std::cout << std::endl;
    //ofs.close();
  }

}
