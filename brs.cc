#include <chrono>
#include <csignal>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "marley_utils.hh"
#include "ConfigFile.hh"
#include "Generator.hh"
#include "Event.hh"
#include "NuclearPhysics.hh"

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

// Macro that we can use to suppress unused parameter warnings for our signal
// handler function. We don't need the signal code s, but we have to include it
// in the function declraration to keep std::signal happy. This trick was taken
// from the accepted answer at
// http://stackoverflow.com/questions/3599160/unused-parameter-warnings-in-c-code
#define UNUSED(x) (void)(x)

// Function that will be used to handle a SIGINT signal
// in our event generation loop. The method used here is
// taken from the second answer at
// http://stackoverflow.com/questions/19366503/in-c-when-interrupted-with-ctrl-c-call-a-function-with-arguments-other-than-s
void signal_handler(int s)
{
  UNUSED(s);
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
  std::cout << "react (MARLEY) v0.1" << std::endl;
  exit(0);
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

  marley::ConfigFile cf(config_file_name);
  marley::Generator gen(cf);

  std::cout << marley_utils::marley_logo << std::endl;
  std::cout << "\"Don't worry about a thing," << std::endl;
  std::cout << "'Cause every little thing gonna be all right.\"" << std::endl;
  std::cout << "-- Bob, \"Three Little Birds\"" << std::endl << std::endl;

  // Get the time that the program was started
  std::chrono::system_clock::time_point start_time_point
    = std::chrono::system_clock::now();

  std::time_t start_time = std::chrono::system_clock::to_time_t(start_time_point);


  std::cout << "MARLEY started on "
    << std::put_time(std::localtime(&start_time), "%c %Z")
    << std::endl;
  #ifndef USE_ROOT
  std::cout << "Seed for random number generator: "
    << gen.get_seed() << std::endl;
  #endif

  int num_old_events = 0;
  #ifdef USE_ROOT
  // Create a pointer to an event object. This will
  // be used to fill the event tree with the
  // events generated in the loop below
  marley::Event* p_event = nullptr;

  // Check if the ROOT file used to store the event tree already exists
  std::string tree_file_name("event_tree.root");
  std::ifstream test_stream(tree_file_name.c_str());
  bool root_file_existed = false;
  if (test_stream.good()) {
    root_file_existed = true;
  }
  test_stream.close();

  // Create a ROOT file to store the event tree
  // if it does not already exist. Otherwise, open
  // it so that we can append to it.
  TFile treeFile(tree_file_name.c_str(),"UPDATE");

  // Check if there was a problem opening the file
  // (e.g., pre-existing file that has the wrong
  // format, etc.). If so, complain and quit.
  if (treeFile.IsZombie()) {
    throw std::runtime_error(std::string("Invalid format or other error ")
      + "encountered while opening the ROOT file " + tree_file_name);
  }

  // Attempt to get the old ROOT event tree if the ROOT file already exists
  TTree* event_tree = nullptr;
  if (root_file_existed) {
    treeFile.GetObject("MARLEY Event Tree", event_tree);
  }

  // If the ROOT file didn't already exist, or if we had a problem
  // retrieving the old event tree, then make a new event tree
  if (!event_tree) {
    // Create a ROOT tree to store the events
    event_tree = new TTree("MARLEY Event Tree", "A tree of marley::Event objects");

    // Create a branch in this ROOT tree, and associate
    // it with the event pointer we made before
    event_tree->Branch("events", &p_event, 32000, 99);

    // Write the seed for the random number generator to the ROOT file
    // Since ROOT does not allow us to write a single integer to the file,
    // we will convert it to a std::string object first
    // TODO: consider writing both the RNG seed and the RNG state string
    // to the UserInfo array of TObjects (a member of the event TTree)
    std::string dummy_str = std::to_string(gen.get_seed());
    treeFile.WriteObject(&dummy_str, "MARLEY RNG Seed");

    // Notify the user of the seed that will be used for this run
    std::cout << "Seed for random number generator: "
    << gen.get_seed() << std::endl;
  }

  // If we were able to read in the event tree from the file, get
  // ready to add new events to it, and notify the user
  else {
    //TODO: add check to handle cases where we can read an event
    //tree with the correct name from the ROOT file, but we can't find
    //a branch that matches the one we expect
    event_tree->SetBranchAddress("events", &p_event);
    num_old_events = event_tree->GetEntries();
    std::cout << "Continuing previous run from ROOT file " << tree_file_name
      << std::endl << "which contains " << num_old_events
      << " events." << std::endl;

    // Get previous random number generator state string from the ROOT file
    std::string* p_rng_state_string = nullptr;
    treeFile.GetObject("MARLEY RNG State String", p_rng_state_string);

    // Use the state string to reset the MARLEY RNG
    gen.seed_using_state_string(p_rng_state_string);

    // Get previous RNG seed from the ROOT file
    std::string* p_rng_seed_str = nullptr;
    treeFile.GetObject("MARLEY RNG Seed", p_rng_seed_str);
    // TODO: add error handling here (p_rng_seed_str may be nullptr)
    uint_fast64_t old_seed = static_cast<uint_fast64_t>(
      std::stoull(*p_rng_seed_str));

    // Notify the user of the seed that was used previously
    std::cout << "The previous run was initialized using the "
      << "random number generator seed " << old_seed
      << std::endl;
  }

  // The amount of data written to the ROOT file
  // represented as a string (e.g., "10 MB")
  std::string data_written;
  #endif

  // TODO: debug numerical errors that arise when
  // Ea = E_threshold

  int n_events = 1e4;

  double Exi = 15;
  int Zi = 19;
  int Ai = 40;
  int twoJi = 2;
  marley::Parity Pi = 1;

  marley::StructureDatabase& db = gen.get_structure_db();

  marley::Particle initial, first, second;

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
  for (int i = 1 + num_old_events; i <= n_events && !interrupted; ++i) {

    bool good = true;
    double Ex = Exi;
    int Z = Zi;
    int A = Ai;
    int twoJ = twoJi;
    marley::Parity P = Pi;
    initial = marley::Particle(marley_utils::get_nucleus_pid(Z, A),
      marley::MassTable::get_atomic_mass(Z, A) + Ex);

    marley::Event e(Ex);
    e.add_initial_particle(initial, marley::Event::ParticleRole::pr_target);
    e.add_final_particle(initial, marley::Event::ParticleRole::pr_residue);

    //std::cout << "Decay #" << i << std::endl;
    //std::cout << "Beginning with " << initial.get_id() << " at Ex = " << Ex
    //  << " MeV" << std::endl;

    while (Ex > 0.01 && good) {
      good = marley::NuclearPhysics::hauser_feshbach_decay(Z, A, initial, first, second,
        Ex, twoJ, P, db, gen);
      //std::cout << "Decay to " << first.get_id() << " and " << second.get_id()
      //  << std::endl;
      //std::cout << second.get_id() << " is at Ex = " << Ex << " MeV."
      //  << std::endl;
      initial = second;
      Z = marley::MassTable::get_particle_Z(initial.get_id());
      A = marley::MassTable::get_particle_A(initial.get_id());

      e.add_final_particle(first);
      *(e.get_residue()) = second;
    }

    if (Ex > 0.01) {
      marley::DecayScheme* ds = db.get_decay_scheme(second.get_id());
      ds->do_cascade(ds->get_pointer_to_closest_level(Ex), &e, gen);
    }

    #ifdef USE_ROOT
    // Get the address of this event object
    p_event = new marley::Event;
    *p_event = e;

    // Store this event in the ROOT tree
    event_tree->Fill();
    #endif

    // Print status messages about simulation progress after every 10
    // events have been generated
    if ((i - num_old_events) % 10 == 1 || i == n_events) {
      // Print a status message showing the current number of events
      std::cout << "Event Count = " << i << "/" << n_events
        << " (" << i*100/static_cast<double>(n_events)
        << "% complete)" << std::endl;

      // Print timing information
      std::chrono::system_clock::time_point current_time_point
        = std::chrono::system_clock::now();
      std::cout << "Elapsed time: "
        << marley_utils::elapsed_time_string(start_time_point,
        current_time_point) << " (Estimated total run time: ";

      marley_utils::seconds<float> estimated_total_time =
        (current_time_point - start_time_point)*(static_cast<float>(n_events
        - num_old_events)/(i - num_old_events));

      std::cout << marley_utils::duration_to_string
        <marley_utils::seconds<float>>(estimated_total_time)
        << ")\033[K" << std::endl;

      #ifdef USE_ROOT
      data_written = marley_utils::num_bytes_to_string(treeFile.GetBytesWritten(),2);
      std::cout << "Data written = " << data_written << "\033[K" << std::endl;
      #endif

      std::time_t estimated_end_time = std::chrono::system_clock::to_time_t(
        start_time_point + std::chrono::duration_cast
        <std::chrono::system_clock::duration>(estimated_total_time));

      std::cout << "MARLEY is estimated to terminate on "
        << std::put_time(std::localtime(&estimated_end_time), "%c %Z")
        << std::endl;

      #ifdef USE_ROOT
      // Move up an extra line if we're using ROOT and
      // therefore displaying information about the
      // amount of data written to disk
      std::cout << "\033[F";
      #endif

      // Move up three lines in std::cout
      std::cout << "\033[F\033[F\033[F";
    }
  }

  std::cout << "\033[E";

  #ifdef USE_ROOT
  // Write the event tree to the ROOT file, replacing the previous
  // version if one exists. Avoid data loss by not deleting the previous
  // version until the new version is completely written to disk.
  event_tree->Write(event_tree->GetName(), TTree::kWriteDelete);

  // Write the internal state string of the random number generator to
  // disk. This will allow MARLEY to resume event generation from where
  // it left off with no loss of consistency. This trick is based on
  // http://stackoverflow.com/questions/18361050/saving-random-number-generator-state-in-c11
  std::string rng_state_string = gen.get_state_string();
  treeFile.WriteObject(&rng_state_string, "MARLEY RNG State String", "WriteDelete");

  // Close the ROOT file and print final information about
  // the amount of data written to disk
  treeFile.Close();
  data_written = marley_utils::num_bytes_to_string(treeFile.GetBytesWritten());
  std::cout << "Data written = " << data_written << "\033[K" << std::endl;
  #endif

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
  std::cout << std::put_time(std::localtime(&end_time), "%c %Z")
    << "\033[K\033[E\033[K";

  return 0;
}
