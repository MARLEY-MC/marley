#include <fstream>
#include <iostream>
#include <string>


class TEnsdfRecord {

  public:
    TEnsdfRecord(std::string rec_text="") {
      sRecordText = rec_text;
    }

  private:
    std::string sRecordText;
};


int main() {

  std::ifstream file_in("ensdf.040");

  std::string line; // String to store the current line
                    // of the ENSDF file during parsing

  std::string temp_line;

  bool no_advance = false; // Flag that prevents advancing through
                           // the ENSDF file when a continuation record
                           // is found

  while (!file_in.eof()) {
    // Get the next line of the file
    // unless we already did
    if (!no_advance) std::getline(file_in, line);
    no_advance = false;

    std::string nuc_id = line.substr(0,5);
    if (nuc_id == " 40K ") {

      // Characters that follow the NUCID entry
      // on each line. These indicate which type
      // of record is found on the line.
      std::string rec_type_id = line.substr(5,4);

      // Level Record
      if (rec_type_id.substr(1,3) == " L ") {
        temp_line = line;
        // Check for continuation records    
        while (!file_in.eof() && !no_advance) {
          // Get the next line of the file
          std::getline(file_in, line);
          // Check to see if the appropriate field
          // contains a continuation record character
          std::string cont_char = line.substr(5,1);
          
          if (cont_char != " " && cont_char != "1") {
            // If it does, add the next line to the current
            // record text.
            temp_line += "\n" + line;
          }

          else {
            // If it doesn't, set a flag so that the parsing
            // loop will restart without advancing to the next
            // line
            no_advance = true;
          }
          
        }

        std::cout << "level record:" << std::endl << temp_line << std::endl;
        
      }

      // Gamma Record
      else if (rec_type_id.substr(1,3) == " G ") {
        std::cout << "gamma record:" << std::endl << line << std::endl;
      }
     
      // Parent Record
      else if (rec_type_id.substr(0,3) == "  P") {
        std::cout << "parent record:" << std::endl << line << std::endl;
      }

      // Normalization Record
      else if (rec_type_id.substr(0,3) == "  N") {
        std::cout << "normalization record:" << std::endl << line << std::endl;
      }
      
      // Product Normalization Record
      else if (rec_type_id.substr(1,2) == "PN") {
        std::cout << "product normalization record:" << std::endl << line << std::endl;
      }

    }

  }

  file_in.close();

}
