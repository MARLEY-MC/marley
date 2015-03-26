#include <fstream>
#include <iostream>
#include <regex>
#include <string>


class TEnsdfRecord {

  public:
    TEnsdfRecord(std::string rec_text="") {
      sRecordText = rec_text;
    }

  private:
    std::string sRecordText;
};

// Regular expressions for identifying ensdf record types
std::string nuc_id = "";
std::string primary_record = "";
std::string continuation_record = "";
std::string record_data = "";
// g++ 4.8 doesn't yet support all C++11 features, including regex.
// When 4.9+ has been out longer, you 
//std::regex rx_primary_level_record(nuc_id + primary_record + " L " + record_data);
std::regex rx_primary_level_record("^[a-zA-Z0-9]\\{5\\}[ 1] L .\\{71\\}$");
std::regex rx_continuation_level_record(nuc_id + continuation_record + " L " + record_data);


int main() {

  std::cout << "blah!";

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
      if (std::regex_match(line, rx_primary_level_record)) {
        temp_line = line;
        // Check for continuation records    
        while (!file_in.eof() && !no_advance) {
          // Get the next line of the file
          std::getline(file_in, line);
         
          // Check to see if the next line is a continuation record 
          if (std::regex_match(line, rx_continuation_level_record)) {
            // If it is, add the next line to the current
            // record text.
            temp_line += "\n" + line;
          }

          else {
            // If it isn't, set a flag so that the parsing
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
