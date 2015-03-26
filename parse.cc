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
std::string nuc_id = "^[[:alnum:] ]{5}";
std::string primary_record = "[ 1]";
std::string continuation_record = "[^ 1]";
std::string record_data = ".{71}";

std::regex rx_primary_level_record(nuc_id + primary_record + " L " + record_data);
std::regex rx_continuation_level_record(nuc_id + continuation_record + " L " + record_data);
std::regex rx_primary_gamma_record(nuc_id + primary_record + " G " + record_data);
std::regex rx_continuation_gamma_record(nuc_id + continuation_record + " G " + record_data);

std::string process_continuation_records(std::ifstream &file_in,
  std::string &record, std::regex &rx_cont_record) {

  std::string line;
  while (!file_in.eof()) {
    // Get the next line of the file
    std::getline(file_in, line);
 
    // Check to see if the next line is a continuation record 
    if (std::regex_match(line, rx_cont_record)) {
      // If it is, add the next line to the current
      // record text.
      record += "\n" + line;
    }

    // If a non-continuation-record line is found,
    // stop reading from the file.
    else return line; 
  
  }

}



int main() {

  std::ifstream file_in("ensdf.040");

  std::string line; // String to store the current line
                    // of the ENSDF file during parsing

  std::string record;

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

      // Level Record
      if (std::regex_match(line, rx_primary_level_record)) {
        record = line;
        line = process_continuation_records(file_in, record, rx_continuation_level_record);
        no_advance = true;
          
        std::cout << "level record:" << std::endl << record << std::endl;
      }

      // Gamma Record
      else if (std::regex_match(line, rx_primary_gamma_record)) {
        record = line;
        line = process_continuation_records(file_in, record, rx_continuation_gamma_record);
        no_advance = true;

        std::cout << "gamma record:" << std::endl << record << std::endl;
      }
     
      //// Parent Record
      //else if (rec_type_id.substr(0,3) == "  P") {
      //  std::cout << "parent record:" << std::endl << line << std::endl;
      //}

      //// Normalization Record
      //else if (rec_type_id.substr(0,3) == "  N") {
      //  std::cout << "normalization record:" << std::endl << line << std::endl;
      //}
      //
      //// Product Normalization Record
      //else if (rec_type_id.substr(1,2) == "PN") {
      //  std::cout << "product normalization record:" << std::endl << line << std::endl;
      //}

    }

  }

  file_in.close();

}
