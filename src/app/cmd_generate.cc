// Standard library includes
#include <chrono>
#include <csignal>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/Logger.hh"
#include "marley/OutputFile.hh"
#include "marley/marley_utils.hh"

namespace {

  volatile static std::sig_atomic_t interrupted = false;

  void signal_handler( int )
  {
    interrupted = true;
  }

  constexpr int DEFAULT_STATUS_UPDATE_INTERVAL = 100;

  std::string format_number( double number ) {
    static std::stringstream temp_stream;
    static bool configured = false;
    if (!configured) {
      temp_stream << std::fixed << std::setprecision(1);
      configured = true;
    }

    temp_stream.str("");
    temp_stream.clear();
    temp_stream << number;
    return temp_stream.str();
  }

  std::string put_time( std::tm* time, const char* format )
  {
    constexpr size_t TIME_STR_SIZE = 100;
    std::string time_str( TIME_STR_SIZE, ' ' );
    std::strftime( &time_str.front(), TIME_STR_SIZE, format, time );
    marley_utils::trim_right_inplace( time_str );
    return time_str;
  }

  std::string makeStatusLines( long ev_count, long num_events,
    long num_old_events,
    std::chrono::system_clock::time_point start_time_point,
    const std::vector< std::shared_ptr<marley::OutputFile> >& output_files )
  {
    std::chrono::system_clock::time_point current_time_point
      = std::chrono::system_clock::now();

    auto elapsed_time = std::chrono::duration_cast<
      marley_utils::seconds<double> >( current_time_point - start_time_point );
    long events_since_start = ev_count - num_old_events;
    double elapsed_seconds = elapsed_time.count();
    double avg_event_rate = events_since_start / elapsed_seconds;

    std::ostringstream temp_oss;
    double percent_complete = static_cast<double>( ev_count )
      / num_events * 100.;
    temp_oss << "\nEvent Count = " << ev_count << "/" << num_events
      << " (" << format_number( percent_complete ) << "% complete, "
      << format_number( avg_event_rate ) << " events / s)\033[K\n";

    temp_oss << "Elapsed time: "
      << marley_utils::elapsed_time_string(start_time_point,
      current_time_point) << " (Estimated total run time: ";

    marley_utils::seconds< double > estimated_total_time =
      ( current_time_point - start_time_point )
      * ( static_cast<double>(num_events - num_old_events)
      / (ev_count - num_old_events) );

    temp_oss << marley_utils::duration_to_string
      < marley_utils::seconds<double> >( estimated_total_time )
      << ")\033[K\n";

    for ( const auto& file : output_files ) {
      temp_oss << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string( file->bytes_written(), 2 )
        << " (estimate)"
        << "\033[K\n";
    }

    std::time_t estimated_end_time = std::chrono::system_clock::to_time_t(
      start_time_point + std::chrono::duration_cast
      < std::chrono::system_clock::duration >( estimated_total_time ) );

    temp_oss << "MARLEY is estimated to terminate on "
      << put_time( std::localtime(&estimated_end_time), "%c %Z" ) << '\n';

    for ( size_t i = 0; i < output_files.size(); ++i ) temp_oss << "\033[F";
    temp_oss << "\033[F\033[F\033[F\033[F";

    return temp_oss.str();
  }

  class StatusInserter : public std::streambuf {
    public:
      StatusInserter( std::streambuf* dest, const long& ev_count,
        const long& num_events, const long& num_old_events,
        const std::chrono::system_clock::time_point& start_time_point,
        const std::vector< std::shared_ptr<marley::OutputFile> >& output_files )
        : std::streambuf(), myDest_( dest ), myIsAtStartOfLine_( true ),
        do_status_( true ), ev_count_( ev_count ), num_events_( num_events ),
        num_old_events_( num_old_events ),
        start_time_point_( start_time_point ), output_files_( output_files )
        {}

      inline void set_do_status( bool do_it ) { do_status_ = do_it; }

    protected:
      std::streambuf* myDest_;
      bool myIsAtStartOfLine_;
      bool do_status_;
      const long& ev_count_;
      const long& num_events_;
      const long& num_old_events_;
      const std::chrono::system_clock::time_point& start_time_point_;
      const std::vector< std::shared_ptr<marley::OutputFile> >& output_files_;

      int overflow( int ch ) override {
        int retval = 0;
        if ( ch != traits_type::eof() ) {
          if ( do_status_ && myIsAtStartOfLine_ ) {
            std::string status = makeStatusLines(ev_count_,
              num_events_, num_old_events_, start_time_point_, output_files_);
            myDest_->sputn( status.data(), status.size() );
          }
          myIsAtStartOfLine_ = ch == '\n';
          if ( myIsAtStartOfLine_ ) {
            std::string erase( "\033[K" );
            myDest_->sputn( erase.data(), erase.size() );
          }
          retval = myDest_->sputc( ch );
        }
        return retval;
      }
  };

}

bool marley::CommandHandler::cmd_generate( std::deque< std::string >& args ) {

  std::streambuf* cout_default_buf = std::cout.rdbuf();
  std::streambuf* cerr_default_buf = std::cerr.rdbuf();

  marley::Error::set_logging_status( false );

  try {

    marley::Logger::Instance().add_stream( std::cout,
      marley::Logger::LogLevel::INFO );

    std::string config_file_name;
    if ( !args.empty() ) config_file_name = args.front();

    bool cfn_empty = config_file_name.empty();

    if ( cfn_empty || config_file_name.front() == '-' )
    {
      if ( !cfn_empty && config_file_name != "-h"
        && config_file_name == "--help" )
      {
        std::cerr << "marley generate: unrecognized option '"
          << config_file_name << "'\n";
      }
      args.clear();
      args.push_front( "generate" );
      marley::CommandHandler::cmd_help( args );
      return false;
    }

    marley::JSON json = marley::JSON::load_file( config_file_name );
    marley::JSONConfig jc( json );

    std::chrono::system_clock::time_point start_time_point
      = std::chrono::system_clock::now();

    std::time_t start_time = std::chrono::system_clock::to_time_t(
      start_time_point );

    std::cout << "\nMARLEY started on "
      << put_time( std::localtime(&start_time), "%c %Z" ) << '\n';

    long num_old_events = 0;

    marley::JSON ex_set = json.get_object( "executable_settings" );
    long num_events = ex_set.get_long( "events", 1e3 );

    int status_update_interval = DEFAULT_STATUS_UPDATE_INTERVAL;
    if ( ex_set.has_key("status_update_interval") ) {
      const auto& sui = ex_set.at( "status_update_interval" );

      bool ok;
      int sui_value = sui.to_long( ok );

      if ( !ok || sui_value < 1 ) {
        throw marley::Error( "Invalid value " + sui.dump_string()
          + " given for the \"status_update_interval\" key in the"
          " job configuration file" );
      }
      else status_update_interval = sui_value;
    }

    std::vector< std::shared_ptr<marley::OutputFile> > output_files;

    if ( ex_set.has_key("output") ) {
      marley::JSON output_set = ex_set.at( "output" );
      if ( !output_set.is_array() ) throw marley::Error( "The"
        " \"output\" key in the executable settings must have a value that"
        " is a JSON array." );
      else for ( const auto& el : output_set.array_range() ) {
        output_files.push_back( marley::OutputFile::make_OutputFile(el) );
      }
    }
    else {
      std::string out_config_str = "{ format: \"ascii\","
        " file: \"events.hepmc3\", mode: \"overwrite\", force: false }";
      auto out_config = marley::JSON::load( out_config_str );

      output_files.push_back( marley::OutputFile::make_OutputFile(
        out_config ) );
    }

    std::unique_ptr< marley::Generator > gen;

    bool need_to_resume = false;
    for ( auto& file : output_files ) {
      if ( file->mode_is_resume() ) {
        if ( need_to_resume ) throw marley::Error( "Only one file may be used"
          " to resume a previous run." );
        else {
          need_to_resume = true;
          bool resume_ok = file->resume( gen, num_old_events );
          if ( !resume_ok ) throw marley::Error( "Failed to resume previous"
            " run from the file \"" + file->name() + '\"' );
        }
      }
    }

    if ( !need_to_resume ) gen = std::make_unique<marley::Generator>(
      jc.create_generator() );

    interrupted = false;
    std::signal( SIGINT, signal_handler );

    long ev_count = 1 + num_old_events;

    StatusInserter my_status_inserter( cout_default_buf, ev_count,
      num_events, num_old_events, start_time_point, output_files );
    std::cout.rdbuf( &my_status_inserter );
    std::cerr.rdbuf( &my_status_inserter );

    start_time_point = std::chrono::system_clock::now();
    start_time = std::chrono::system_clock::to_time_t( start_time_point );

    for (; ev_count <= num_events && !interrupted; ++ev_count) {

      auto event = gen->create_event();
      event->set_event_number( ev_count );

      for ( const auto& file : output_files ) {
        file->write_event( event.get() );
      }

      if ( (ev_count - num_old_events) % status_update_interval == 1
        || ev_count == num_events || status_update_interval == 1 )
      {
        my_status_inserter.set_do_status( false );

        std::cout << makeStatusLines( ev_count, num_events, num_old_events,
          start_time_point, output_files );

        my_status_inserter.set_do_status( true );
      }
    }

    std::cout.rdbuf( cout_default_buf );
    std::cerr.rdbuf( cerr_default_buf );

    std::cout << "\033[E";

    for (const auto& file : output_files) {
      std::cout << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string( file->bytes_written() )
        << "\033[K\n";
    }

    std::chrono::system_clock::time_point end_time_point
      = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(
      end_time_point);

    if (!interrupted) std::cout << "MARLEY terminated normally on ";
    else std::cout << "MARLEY was interrupted by the user on ";
    std::cout << put_time(std::localtime(&end_time), "%c %Z")
      << "\033[K\033[E\033[K\n\033[K";

    return true;
  }

  catch ( const std::exception& error ) {
    std::cout.rdbuf( cout_default_buf );
    std::cerr.rdbuf( cerr_default_buf );

    auto& log = marley::Logger::Instance();
    log.flush();
    MARLEY_LOG_ERROR() << error.what();
  }

  return false;
}
