#include <chrono>
#include <csignal>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "marley_utils.hh"

#ifdef USE_ROOT
#include "RootConfigFile.hh"
#else
#include "ConfigFile.hh"
#endif

#include "Generator.hh"
#include "Event.hh"
#include "Logger.hh"

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

// Flag that will alert the main event generation loop that
// the user has interrupted the program's execution (probably
// by pressing ctrl+c). This is a global variable (out of
// necessity, since it must be modified by our signal handler),
// but its scope is limited to this source file by the use
// of the static keyword in its declaration.
volatile static std::sig_atomic_t interrupted = false;

// Function that will be used to handle a SIGINT signal
// in our event generation loop. The method used here is
// taken from the second answer at
// http://stackoverflow.com/questions/19366503/in-c-when-interrupted-with-ctrl-c-call-a-function-with-arguments-other-than-s
// Note that, since we don't use the signal code, we leave the integer argument of the function nameless. This prevents
// "unused parameter" warnings from the compiler.
void signal_handler(int)
{
  interrupted = true;
}

// Global help message string. Its scope is limited to this
// source file via the static keyword
static const std::string help_message1 = "Usage: ";
static const std::string help_message2 = " [OPTION...] FILE\n"
"\n"
"  -h, --help     Print this help message\n"
"  -v, --version  Print version and exit\n";

inline void print_help(std::string executable_name) {
  std::cout << help_message1 + executable_name + help_message2;
  exit(0);
}

inline void print_version() {
  std::cout << "MARLEY (Model of Argon Reaction Low Energy Yields) v"
    << marley_utils::MARLEY_VERSION << std::endl;
  exit(0);
}

// Use this instead of std::put_time to allow this executable to be built
// with g++ 4.9 (issue fixed in 5.0). See discussion here:
// http://stackoverflow.com/a/14137287
inline std::string put_time(std::tm* time, const char* format)
{
  static constexpr size_t TIME_STR_SIZE = 100;
  std::string time_str(TIME_STR_SIZE, ' ');
  // Pass a pointer to the first character in the time string as the char*
  // argument of std::strftime(). This is well-defined behavior because the
  // storage of std::string is guaranteed by the C++11 standard to be
  // contiguous, and we have pre-allocated storage inside the string of the
  // necessary size using the constructor above.
  std::strftime(&time_str.front(), TIME_STR_SIZE, format, time);
  marley_utils::trim_right_inplace(time_str);
  return time_str;
}

inline bool check_if_file_exists(const std::string& filename) {
  static std::ifstream test_stream;
  test_stream.open(filename);
  bool file_exists = test_stream.good();
  test_stream.close();
  return file_exists;
}

int main(int argc, char* argv[]){

  // TODO: if you need command line parsing beyond the
  // trivial stuff used here, consider using a header-only
  // option parser library like cxxopts (https://github.com/jarro2783/cxxopts)
  std::string config_file_name;

  // If the user has not supplied any command-line
  // arguments, display the standard help message
  // and exit
  if (argc <= 1) {
    print_help(argv[0]);
  }
  // The first command-line argument does not begin with a
  // hyphen, so assume that it is the configuration file name.
  else if (std::string(argv[1]).substr(0,1) != "-") {
    config_file_name = argv[1];
  }
  // The first command-line argument begins with a hyphen,
  // so treat it as an option and react accordingly.
  // All of the current options cannot be combined, so
  // just parse the first one and ignore the others.
  else {
    std::string option = argv[1];
    if (option == "-h" || option == "--help") {
      print_help(argv[0]);
    }
    else if (option == "-v" || option == "--version") {
      print_version();
    }
    else {
     std::cout << argv[0] << ": unrecognized "
       << "command line option '" <<  option
       << "'" << std::endl;
     print_help(argv[0]);
    }
  }

  //std::ofstream log_file("marley.log", std::ofstream::out | std::ofstream::trunc);
  //marley::Logger::Instance().add_stream(log_file, marley::Logger::LogLevel::INFO);
  marley::Logger::Instance().add_stream(std::cout, marley::Logger::LogLevel::INFO);
  #ifdef USE_ROOT
  marley::RootConfigFile cf(config_file_name);
  #else
  marley::ConfigFile cf(config_file_name);
  #endif
  marley::Generator gen(cf);

  std::cout << std::endl << marley_utils::marley_logo << std::endl;
  std::cout << "\"Don't worry about a thing," << std::endl;
  std::cout << "'Cause every little thing gonna be all right.\"" << std::endl;
  std::cout << "-- Bob, \"Three Little Birds\"" << std::endl << std::endl;

  // Get the time that the program was started
  std::chrono::system_clock::time_point start_time_point
    = std::chrono::system_clock::now();

  std::time_t start_time = std::chrono::system_clock::to_time_t(start_time_point);


  std::cout << "MARLEY started on "
    << put_time(std::localtime(&start_time), "%c %Z")
    << std::endl;
  #ifndef USE_ROOT
  std::cout << "Seed for random number generator: "
    << gen.get_seed() << std::endl;
  #endif

  int num_old_events = 0;
  // Desired number of events to be generated in this run. This
  // will be read back from the ROOT file and overwritten
  // if we're doing a continuation run.
  int num_events = cf.get_num_events();

  // Whether to write the events to a HEPEvt file
  bool write_hepevt = cf.check_write_hepevt();
  std::ofstream hepevt_stream;
  // If writing a HEPEvt file is enabled, then prepare one
  if (write_hepevt) {
    std::string hepevt_file_name = cf.get_hepevt_filename();
    bool hepevt_file_exists = check_if_file_exists(hepevt_file_name);
    bool overwrite_check = cf.check_overwrite_hepevt();
    auto open_mode_flag = std::ofstream::ate; // add new events at the end by default
    if (hepevt_file_exists && overwrite_check) {
      std::string response;
      while (std::cout << "Overwrite HEPEvt file " << hepevt_file_name
        << " [y/n]? " && std::getline(std::cin, response)
        && !(response == "y" || response == "n" || response == "Y" || response == "N"));
      if (response == "y" || response == "Y") open_mode_flag = std::ofstream::trunc;
    }
    hepevt_stream.open(hepevt_file_name, std::ofstream::out
      | open_mode_flag);
    std::cout << "Events for this run will be ";
    if (hepevt_file_exists && open_mode_flag == std::ofstream::ate)
      std::cout << "appended";
    else std::cout << "written";
    std::cout << " to the HEPEvt format file "
      << hepevt_file_name << std::endl;
  }

  #ifdef USE_ROOT
  // Whether to write the events to a ROOT file
  bool write_root = cf.check_write_root();
  // Whether to check before overwriting the ROOT file
  bool check_overwrite_root = cf.check_overwrite_root();

  // Create a pointer to an event object. This will
  // be used to fill the event tree with the
  // events generated in the loop below
  marley::Event* p_event = nullptr;

  // Check if the ROOT file used to store the event tree already exists
  std::string tree_file_name = cf.get_root_filename();
  bool root_file_existed = check_if_file_exists(tree_file_name);

  // The amount of data written to the ROOT file
  // represented as a string (e.g., "10 MB")
  std::string data_written;

  std::unique_ptr<TFile> treeFile(nullptr);
  // This is a bare pointer, but ROOT will associate it with
  // treeFile, so we don't want to delete it ourselves or let
  // a smart pointer do it.
  TTree* event_tree = nullptr;

  if (write_root) {

    std::string tfile_open_mode("recreate");
    if (root_file_existed && check_overwrite_root) {
      std::string response;
      while (std::cout << "Overwrite ROOT file " << tree_file_name
        << " [y/n]? " && std::getline(std::cin, response)
        && !(response == "y" || response == "n" || response == "Y" || response == "N"));
      if (response == "n" || response == "N") tfile_open_mode = std::string("update");
    }

    // Create a ROOT file to store the event tree
    // if it does not already exist. Otherwise, open
    // it so that we can append to it.
    treeFile = std::make_unique<TFile>(tree_file_name.c_str(),
      tfile_open_mode.c_str());

    // Check if there was a problem opening the file
    // (e.g., pre-existing file that has the wrong
    // format, etc.). If so, complain and quit.
    if (treeFile->IsZombie()) {
      throw marley::Error(std::string("Invalid format or other error ")
        + "encountered while opening the ROOT file " + tree_file_name);
    }

    // Attempt to get the old ROOT event tree if the ROOT file already exists
    if (root_file_existed) {
      treeFile->GetObject("MARLEY Event Tree", event_tree);
    }

    // If the ROOT file didn't already exist, or if we had a problem
    // retrieving the old event tree, then make a new event tree
    if (!event_tree) {
      // Create a ROOT tree to store the events
      event_tree = new TTree("MARLEY Event Tree", "A tree of marley::Event objects");

      // Create a branch in this ROOT tree, and associate
      // it with the event pointer we made before
      event_tree->Branch("events", "marley::Event", &p_event);

      // Write the seed for the random number generator to the ROOT file
      // Since ROOT does not allow us to write a single integer to the file,
      // we will convert it to a std::string object first
      // TODO: consider writing both the RNG seed and the RNG state string
      // to the UserInfo array of TObjects (a member of the event TTree)
      std::string dummy_str = std::to_string(gen.get_seed());
      treeFile->WriteObject(&dummy_str, "MARLEY RNG Seed");

      // Write the desired number of events in this run to the ROOT file.
      // Use a string so that we can write a single int without much trouble.
      std::string str_num_events = std::to_string(num_events);
      treeFile->WriteObject(&str_num_events, "number of MARLEY events to generate");

      std::cout << "Events for this run will be written to the ROOT file "
        << tree_file_name << std::endl;

      // Notify the user of the seed that will be used for this run
      std::cout << "Seed for random number generator: "
      << gen.get_seed() << std::endl;
    }

    // If we were able to read in the event tree from the file, get
    // ready to add new events to it, and notify the user
    else {

      // Get previous number of events to generate from the ROOT file
      std::string* str_num_events;
      treeFile->GetObject("number of MARLEY events to generate", str_num_events);
      num_events = std::stoi(*str_num_events);

      //TODO: add check to handle cases where we can read an event
      //tree with the correct name from the ROOT file, but we can't find
      //a branch that matches the one we expect
      event_tree->SetBranchAddress("events", &p_event);
      num_old_events = event_tree->GetEntries();
      std::cout << "Continuing previous run from ROOT file " << tree_file_name
        << std::endl << "which contains " << num_old_events
        << " events." << std::endl;
      std::cout << "A total of " << num_events << " events will be generated."
        << std::endl;

      // Get previous random number generator state string from the ROOT file
      std::string* p_rng_state_string = nullptr;
      treeFile->GetObject("MARLEY RNG State String", p_rng_state_string);

      // Use the state string to reset the MARLEY RNG
      gen.seed_using_state_string(p_rng_state_string);

      // Get previous RNG seed from the ROOT file
      std::string* p_rng_seed_str = nullptr;
      treeFile->GetObject("MARLEY RNG Seed", p_rng_seed_str);
      // TODO: add error handling here (p_rng_seed_str may be nullptr)
      uint_fast64_t old_seed = static_cast<uint_fast64_t>(
        std::stoull(*p_rng_seed_str));

      // Notify the user of the seed that was used previously
      std::cout << "The previous run was initialized using the "
        << "random number generator seed " << old_seed
        << std::endl;
    }
  }
  #endif

  // TODO: debug numerical errors that arise when
  // Ea = E_threshold

  // Display all floating-point numbers without
  // using scientific notation and using
  // one decimal digit
  std::cout << std::fixed << std::setprecision(1);

  // Use the signal handler defined above to deal with
  // SIGINT signals (e.g., ctrl+c interruptions initiated
  // by the user). This will allow us to terminate the
  // loop gracefully, leaving a valid event tree file, etc.
  std::signal(SIGINT, signal_handler);

  // Generate all of the requested events. End the loop early
  // if the user interrupts execution (e.g., via ctrl+C)
  for (int i = 1 + num_old_events; i <= num_events && !interrupted; ++i) {

    // Create an event using the generator object
    marley::Event e = gen.create_event();

    #ifdef USE_ROOT
    // Get the address of this event object
    p_event = new marley::Event;
    *p_event = e;

    // Store this event in the ROOT tree
    if (write_root) event_tree->Fill();
    #endif

    if (write_hepevt && hepevt_stream.good()) {
      e.write_hepevt(i, hepevt_stream);
    }

    // Print status messages about simulation progress after every 100
    // events have been generated
    if ((i - num_old_events) % 100 == 1 || i == num_events) {
      // Print a status message showing the current number of events
      std::cout << "Event Count = " << i << "/" << num_events
        << " (" << i*100/static_cast<double>(num_events)
        << "% complete)" << std::endl;

      // Print timing information
      std::chrono::system_clock::time_point current_time_point
        = std::chrono::system_clock::now();
      std::cout << "Elapsed time: "
        << marley_utils::elapsed_time_string(start_time_point,
        current_time_point) << " (Estimated total run time: ";

      marley_utils::seconds<float> estimated_total_time =
        (current_time_point - start_time_point)*(static_cast<float>(num_events
        - num_old_events)/(i - num_old_events));

      std::cout << marley_utils::duration_to_string
        <marley_utils::seconds<float>>(estimated_total_time)
        << ")\033[K" << std::endl;

      #ifdef USE_ROOT
      if (write_root) {
        data_written = marley_utils::num_bytes_to_string(treeFile->GetBytesWritten(),2);
        std::cout << "Data written to ROOT file = " << data_written << "\033[K" << std::endl;
      }
      #endif

      if (write_hepevt) {
        hepevt_stream.flush();
        double num_bytes_to_hepevt_file = hepevt_stream.tellp();
        std::cout << "Data written to HEPEvt file = "
          << marley_utils::num_bytes_to_string(num_bytes_to_hepevt_file)
          << "\033[K" << std::endl;
      }

      std::time_t estimated_end_time = std::chrono::system_clock::to_time_t(
        start_time_point + std::chrono::duration_cast
        <std::chrono::system_clock::duration>(estimated_total_time));

      std::cout << "MARLEY is estimated to terminate on "
        << put_time(std::localtime(&estimated_end_time), "%c %Z")
        << std::endl;

      #ifdef USE_ROOT
      // Move up an extra line if we're using ROOT and
      // therefore displaying information about the
      // amount of data written to disk
      if (write_root) std::cout << "\033[F";
      #endif

      // Move up an extra line if we're writing data to a HEPEvt format file
      // and therefore displaying the amount of data written to it.
      if (write_hepevt) {
        std::cout << "\033[F";
      }

      // Move up three lines in std::cout
      std::cout << "\033[F\033[F\033[F";
    }
  }

  std::cout << "\033[E";

  #ifdef USE_ROOT
  // Write the event tree to the ROOT file, replacing the previous
  // version if one exists. Avoid data loss by not deleting the previous
  // version until the new version is completely written to disk.
  if (write_root) {
    event_tree->Write(event_tree->GetName(), TTree::kWriteDelete);

    // Write the internal state string of the random number generator to
    // disk. This will allow MARLEY to resume event generation from where
    // it left off with no loss of consistency. This trick is based on
    // http://stackoverflow.com/questions/18361050/saving-random-number-generator-state-in-c11
    std::string rng_state_string = gen.get_state_string();
    treeFile->WriteObject(&rng_state_string, "MARLEY RNG State String", "WriteDelete");

    // Close the ROOT file and print final information about
    // the amount of data written to disk
    treeFile->Close();
    data_written = marley_utils::num_bytes_to_string(treeFile->GetBytesWritten());
    std::cout << "Data written to ROOT file = " << data_written << "\033[K" << std::endl;
  }
  #endif

  if (write_hepevt) {
    hepevt_stream.flush();
    double num_bytes_to_hepevt_file = hepevt_stream.tellp();
    std::cout << "Data written to HEPEvt file = "
      << marley_utils::num_bytes_to_string(num_bytes_to_hepevt_file)
      << "\033[K" << std::endl;
  }

  // Display the time that the program terminated
  std::chrono::system_clock::time_point end_time_point
    = std::chrono::system_clock::now();
  std::time_t end_time = std::chrono::system_clock::to_time_t(end_time_point);

  if (!interrupted) {
    std::cout << "MARLEY terminated normally on ";
  }
  else {
    std::cout << "MARLEY was interrupted by the user on ";
  }
  std::cout << put_time(std::localtime(&end_time), "%c %Z")
    << "\033[K\033[E\033[K";

  return 0;
}
