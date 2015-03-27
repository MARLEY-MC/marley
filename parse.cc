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

  std::string nuc_id = " 40K ";
  std::string filename = "ensdf.040";

  // Regular expressions for identifying ensdf record types
  //std::string generic_nuc_id = "^[[:alnum:] ]{5}";
  std::string primary_record = "[ 1]";
  std::string continuation_record = "[^ 1]";
  std::string record_data = ".{71}";
  
  std::regex rx_end_record("\\s*"); // Matches blank lines
  //std::regex rx_generic_primary_identification_record(nuc_id + primary_record + "   " + record_data);
  std::regex rx_primary_identification_record(nuc_id + primary_record
    + "   ADOPTED LEVELS, GAMMAS.{49}");
  std::regex rx_primary_level_record(nuc_id + primary_record + " L " + record_data);
  std::regex rx_continuation_level_record(nuc_id + continuation_record + " L " + record_data);
  std::regex rx_primary_gamma_record(nuc_id + primary_record + " G " + record_data);
  std::regex rx_continuation_gamma_record(nuc_id + continuation_record + " G " + record_data);

  // Open the ENSDF file for parsing
  std::ifstream file_in(filename);

  std::string line; // String to store the current line
                    // of the ENSDF file during parsing

  std::string record;

  bool found_cascade = false; // Flag that indicates whether or not
                              // gamma cascade data were found for the
                              // current nuc_id.
                              //
  while (!file_in.eof()) {
    std::getline(file_in, line);
    if (std::regex_match(line, rx_primary_identification_record)) {
      found_cascade = true;
      break;
    }
  }

  if (!found_cascade) {
    std::cout << "Gamma cascade data (adopted levels, gammas) for " + nuc_id << std::endl;
    std::cout << "could not be found in the file " + filename << std::endl;
  }
  else {
    std::cout << "Gamma cascade data for " + nuc_id << " found. Using ENSDF dataset" << std::endl;
    std::cout << line << std::endl;
  }

  bool no_advance = false; // Flag that prevents advancing through
                           // the ENSDF file when a continuation record
                           // is found

  while (!file_in.eof()) {
    // Get the next line of the file
    // unless we already did
    if (!no_advance) std::getline(file_in, line);
    no_advance = false;

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

      std::cout << "gamma record:" << std::endl;
      std::cout << record << std::endl;
      std::cout << "   energy: " << record.substr(9,9) << std::endl;
      std::cout << "   relative photon intensity: " << record.substr(21,7) << std::endl;
      std::cout << "   total conversion coefficient: " << record.substr(55,6) << std::endl;
      std::cout << "   relative total transition intensity: " << record.substr(64,9) << std::endl;
    }

    else if (std::regex_match(line, rx_end_record)) {
      std::cout << "Finished parsing gamma cascade data." << std::endl;
      break;
    }

  }

  file_in.close();

}
