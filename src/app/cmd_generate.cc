/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

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

// POSIX includes
#include <sys/ioctl.h>
#include <unistd.h>

// HepMC3 includes
#include "HepMC3/GenEvent.h"

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/Logger.hh"
#include "marley/OutputFile.hh"
#include "marley/marley_utils.hh"
#include "marley/hepmc3_utils.hh"

namespace {

  volatile static std::sig_atomic_t interrupted = false;
  volatile static std::sig_atomic_t terminal_resized = false;

  void signal_handler( int )
  {
    interrupted = true;
  }

  void sigwinch_handler( int )
  {
    terminal_resized = true;
  }

  constexpr int DEFAULT_STATUS_UPDATE_INTERVAL = 100;

  // Files above this count are reported as a single combined total line.
  // Adjust here; no other code needs to change.
  constexpr size_t MAX_FILE_STATUS_LINES = 2u;

  // Debounce period for terminal resize signals.
  constexpr auto RESIZE_DEBOUNCE = std::chrono::milliseconds( 150 );

  // True when the terminal is too small or stdout is not a TTY.
  // Reset to false at the start of each cmd_generate() call.
  static bool g_fallback_mode = false;

  struct TermSize {
    int rows;
    int cols;
  };

  // Query terminal dimensions via ioctl. Returns {0, 0} when stdout is not a
  // TTY (e.g., piped to a file), which callers treat as a signal to enter
  // fallback mode and suppress escape sequences.
  TermSize get_terminal_size() {
    struct winsize w;
    if ( ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0
         && w.ws_row > 0 && w.ws_col > 0 )
    {
      return { static_cast< int >( w.ws_row ),
        static_cast< int >( w.ws_col ) };
    }
    return { 0, 0 };
  }

  // Truncate s to at most max_len visible characters, appending the ellipsis
  // character (U+2026) if the string was cut. Counts bytes; acceptable because
  // all status content is ASCII or near-ASCII. Substitute a UTF-8 column
  // counter if non-ASCII filenames must be handled precisely.
  std::string truncate_for_terminal( const std::string& s, int max_len ) {
    if ( max_len <= 0 ) return {};
    if ( static_cast< int >( s.size() ) <= max_len ) return s;
    return s.substr( 0, (size_t)(max_len - 1) ) + "\xe2\x80\xa6"; // UTF-8 U+2026 …
  }

  std::string format_number( double number ) {
    static std::stringstream temp_stream;
    static bool configured = false;
    if ( !configured ) {
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

  // Declare (or re-declare) the scroll region and prepare the status zone.
  // Only the status zone rows are cleared; pre-existing log content above is
  // preserved. Safe to call both at startup and after a terminal resize.
  void setup_scroll_region( int num_status_lines ) {
    TermSize ts = get_terminal_size();
    int log_zone_end = ts.rows - num_status_lines;

    std::cout << "\033[1;" << log_zone_end << "r"; // declare scroll region

    for ( int r = log_zone_end + 1; r <= ts.rows; ++r )
      std::cout << "\033[" << r << ";1H\033[K";    // clear status zone only

    std::cout << "\033[" << log_zone_end << ";1H";  // park cursor at base of log zone
    std::flush( std::cout );
  }

  // Reset scroll region on every exit path.
  //
  // When num_status_lines > 0 and clear_status == true (normal / SIGINT exit):
  //   clears the status zone rows, resets the scroll region, and positions the
  //   cursor immediately below the last log content so the final summary prints
  //   without a gap.
  //
  // When num_status_lines > 0 and clear_status == false (exception exit):
  //   leaves the status zone intact so the user can see progress at the time of
  //   the throw, resets the scroll region, and positions the cursor just below
  //   the last status row so the error message appears directly beneath it.
  //
  // When num_status_lines == 0 (early exception or fallback mode):
  //   just resets the scroll region and leaves the cursor in place.
  void reset_terminal( int num_status_lines = 0, bool clear_status = true ) {
    TermSize ts = get_terminal_size();
    if ( ts.rows == 0 ) return; // non-TTY: no escape codes needed

    if ( num_status_lines > 0 ) {
      int log_zone_end = ts.rows - num_status_lines;
      if ( clear_status ) {
        // Clear status zone rows so they don't linger as blank lines
        for ( int r = log_zone_end + 1; r <= ts.rows; ++r )
          std::cout << "\033[" << r << ";1H\033[K";
        // Reset scroll region, position cursor right below last log content
        std::cout << "\033[r"
                  << "\033[" << log_zone_end << ";1H";
      } else {
        // Preserve status zone; position cursor just below it
        std::cout << "\033[r"
                  << "\033[" << ts.rows << ";1H\n";
      }
    } else {
      std::cout << "\033[r"; // reset scroll region; leave cursor in place
    }
    std::flush( std::cout );
  }

  // Switch to plain streaming output with no status display. Safe to call when
  // stdout is a pipe (ts.rows == 0); the notice is suppressed in that case.
  void enter_fallback_mode() {
    if ( g_fallback_mode ) return;
    g_fallback_mode = true;
    TermSize ts = get_terminal_size();
    if ( ts.rows > 0 ) {
      std::cout << "\033[r"  // reset scroll region (no-op if never set)
                << "[MARLEY] Terminal too small for status display."
                   " Running in log-only mode.\n";
    }
    std::flush( std::cout );
  }

  // Leave fallback mode and reinitialize the scroll region. Called when the
  // terminal is resized to a usable size while fallback mode is active.
  void exit_fallback_mode( int num_status_lines ) {
    if ( !g_fallback_mode ) return;
    g_fallback_mode = false;
    setup_scroll_region( num_status_lines );
    std::cout << "[MARLEY] Terminal size restored. Status display active.\n";
    std::flush( std::cout );
  }

  // Overwrite the status zone using absolute cursor positioning, then restore
  // the cursor to its saved position in the log zone. Never touches the log
  // zone during rendering; no line-counting required.
  void update_status_bars(
    long ev_count, long num_events, long num_old_events,
    std::chrono::system_clock::time_point start_time_point,
    const std::vector< std::shared_ptr<marley::OutputFile> >& output_files,
    int num_status_lines )
  {
    TermSize ts = get_terminal_size();
    int status_start_row = ts.rows - num_status_lines + 1;

    // Compute display values (same arithmetic as the original makeStatusLines())
    auto current_tp = std::chrono::system_clock::now();

    auto elapsed = std::chrono::duration_cast< marley_utils::seconds<double> >(
      current_tp - start_time_point );
    double avg_rate = (ev_count - num_old_events) / elapsed.count();
    double pct = static_cast<double>( ev_count ) / num_events * 100.;

    marley_utils::seconds<double> est_total =
      ( current_tp - start_time_point )
      * ( static_cast<double>(num_events - num_old_events)
          / (ev_count - num_old_events) );

    std::time_t est_end = std::chrono::system_clock::to_time_t(
      start_time_point + std::chrono::duration_cast<
      std::chrono::system_clock::duration>( est_total ) );

    std::ostringstream oss;
    oss << "\033[s"; // save cursor (current log zone position)

    // Jump to absolute row, clear it, write truncated content.
    auto write_line = [&]( int row, const std::string& content ) {
      oss << "\033[" << row << ";1H\033[K";
      oss << truncate_for_terminal( content, ts.cols );
    };

    int row = status_start_row;

    // Line 1: event count and rate
    { std::ostringstream l;
      l << "Event Count = " << ev_count << "/" << num_events
        << " (" << format_number(pct) << "% complete, "
        << format_number(avg_rate) << " events / s)";
      write_line( row++, l.str() ); }

    // Line 2: elapsed and estimated total run time
    { std::ostringstream l;
      l << "Elapsed time: "
        << marley_utils::elapsed_time_string(start_time_point, current_tp)
        << "  (Estimated total run time: "
        << marley_utils::duration_to_string< marley_utils::seconds<double> >(
           est_total ) << ")";
      write_line( row++, l.str() ); }

    // Line(s) 3+: output file status (individual or combined)
    if ( output_files.size() <= MAX_FILE_STATUS_LINES ) {
      for ( const auto& file : output_files ) {
        std::ostringstream l;
        l << "Data written to " << file->name() << "  "
          << marley_utils::num_bytes_to_string( file->bytes_written(), 2 )
          << " (estimate)";
        write_line( row++, l.str() );
      }
    } else {
      size_t total_bytes = 0;
      for ( const auto& file : output_files )
        total_bytes += file->bytes_written();
      std::ostringstream l;
      l << "Data written to " << output_files.size() << " output files  "
        << marley_utils::num_bytes_to_string( total_bytes, 2 ) << " (estimate)";
      write_line( row++, l.str() );
    }

    // Last line: estimated termination timestamp
    { std::ostringstream l;
      l << "MARLEY is estimated to terminate on "
        << put_time( std::localtime(&est_end), "%c %Z" );
      write_line( row, l.str() ); }

    oss << "\033[u"; // restore cursor to log zone
    std::cout << oss.str();
    std::flush( std::cout );
  }

} // anonymous namespace

bool marley::CommandHandler::cmd_generate( std::deque< std::string >& args ) {

  marley::Error::set_logging_status( false );

  // Declared here so the catch block can pass it to reset_terminal().
  // Remains 0 if an exception is thrown before output_files is populated.
  int num_status_lines = 0;

  // Flag that indicates whether caught exceptions need to call
  // reset_terminal(). After we set up the status line display, this becomes
  // important. Before that, it clears logging messages from the screen
  // unnecessarily.
  bool exception_needs_reset_terminal = false;

  // Create a pointer to store a buffered event and its corresponding
  // random number generator state string. We write out the buffered
  // event only when generation of the following event has successfully
  // completed. This saves space by allowing the generator state to be
  // saved only when needed at the end of a run (normal termination or
  // an early exit in response to an interruption/exception). The state can
  // be used to restart the MARLEY simulation job where it left off
  // without a drift in the random number state.
  std::shared_ptr< HepMC3::GenEvent > buffered_event;
  std::string cached_generator_state;

  std::vector< std::shared_ptr<marley::OutputFile> > output_files;

  try {

    std::string config_file_name;
    if ( !args.empty() ) config_file_name = args.front();

    bool cfn_empty = config_file_name.empty();

    if ( cfn_empty || config_file_name.front() == '-' )
    {
      if ( !cfn_empty && config_file_name != "-h"
        && config_file_name != "--help" )
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

    MARLEY_LOG( INFO, "app" ) << "Requested events: " << num_events
      << ", configuration file: \"" << config_file_name << "\"";

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

    // Fixed status line count for this run (computed once; used for scroll
    // region sizing and fallback threshold throughout).
    const int file_lines = ( output_files.size() <= MAX_FILE_STATUS_LINES )
      ? static_cast< int >( output_files.size() ) : 1;
    num_status_lines = 3 + file_lines;

    // Reset module-level state in case cmd_generate is called more than once.
    interrupted = false;
    terminal_resized = false;
    g_fallback_mode = false;

    std::signal( SIGINT,   signal_handler );
    std::signal( SIGWINCH, sigwinch_handler );

    // Initialize display before the event loop. Enter fallback mode if stdout
    // is not a TTY or the terminal is too small for the status zone.
    {
      TermSize ts = get_terminal_size();
      if ( ts.rows == 0 || ts.rows < num_status_lines + 2 ) {
        enter_fallback_mode();
      } else {
        setup_scroll_region( num_status_lines );
      }
    }

    exception_needs_reset_terminal = true;

    long ev_count = 1 + num_old_events;

    // Reset the start timestamp now that display setup is complete,
    // so rate and ETA calculations exclude initialization overhead.
    start_time_point = std::chrono::system_clock::now();
    start_time = std::chrono::system_clock::to_time_t( start_time_point );

    auto last_resize_time = std::chrono::steady_clock::time_point{};

    for (; ev_count <= num_events && !interrupted; ++ev_count) {

      // Handle terminal resize with clock-based debounce. The debounce
      // prevents thrashing during active window dragging; it is included
      // unconditionally because some physics configurations (e.g. CEvNS)
      // run at rates too fast to rely on event-loop latency alone.
      if ( terminal_resized ) {
        auto now = std::chrono::steady_clock::now();
        if ( now - last_resize_time >= RESIZE_DEBOUNCE ) {
          terminal_resized = false;
          last_resize_time = now;

          TermSize ts = get_terminal_size();
          if ( ts.rows == 0 || ts.rows < num_status_lines + 2 ) {
            enter_fallback_mode();
          } else if ( g_fallback_mode ) {
            exit_fallback_mode( num_status_lines );
          } else {
            setup_scroll_region( num_status_lines );
          }
        }
      }

      auto event = gen->create_event();
      event->set_event_number( ev_count );

      // If we have a buffered event from the previous iteration, we will
      // write it to the output file(s) now. The generator state string
      // will not be included here to save space. We can use implicit
      // conversion of the std::shared_ptr here since it default-constructs
      // to a nullptr and will therefore evaluate to false if we haven't
      // used it yet.
      if ( buffered_event ) {
        for ( const auto& file : output_files ) {
          file->write_event( buffered_event.get() );
        }
      }

      // Replace the buffered event with the current event. Also cache the
      // generator state string value corresponding to when the current
      // event was finished.
      buffered_event = event;
      cached_generator_state = gen->get_state_string();

      if ( !g_fallback_mode
        && ( (ev_count - num_old_events) % status_update_interval == 1
             || ev_count == num_events
             || status_update_interval == 1 ) )
      {
        update_status_bars( ev_count, num_events, num_old_events,
          start_time_point, output_files, num_status_lines );
      }

    } // event loop

    // We've exited the event loop, so write out the last completed event
    // (if any), which will be stored in the buffered_event pointer.
    // Attach the generator state string this time so that the MARLEY
    // job can be resumed from where it left off.
    // NOTE: This call to OutputFile::write_event() handles normal
    // termination and interruption via the SIGINT signal. The exception
    // exit path is handled separately below in the catch block.
    if ( buffered_event ) {
      marley::Generator::add_state_to_event( *buffered_event,
        cached_generator_state );
      for ( const auto& file : output_files ) {
        file->write_event( buffered_event.get() );
      }
    }

    reset_terminal( g_fallback_mode ? 0 : num_status_lines );

    for ( const auto& file : output_files ) {
      std::cout << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string( file->bytes_written() ) << '\n';
    }

    std::chrono::system_clock::time_point end_time_point
      = std::chrono::system_clock::now();
    std::time_t end_time
      = std::chrono::system_clock::to_time_t( end_time_point );

    if ( !interrupted ) {
      MARLEY_LOG( NOTICE, "app" ) << "Generated " << ( ev_count - 1
        - num_old_events ) << " event(s) successfully";
      std::cout << "MARLEY terminated normally on ";
      }
    else {
      MARLEY_LOG( NOTICE, "app" ) << "Generation interrupted after "
        << ( ev_count - 1 - num_old_events ) << " event(s)";
      std::cout << "MARLEY was interrupted by the user on ";
    }
    std::cout << put_time( std::localtime(&end_time), "%c %Z" ) << '\n';

    return true;
  }

  catch ( const std::exception& error ) {
    if ( exception_needs_reset_terminal ) {
      reset_terminal( g_fallback_mode ? 0 : num_status_lines,
        /*clear_status=*/false );
    }

    // Write out the buffered event (if any) that was successfully completed
    // before the exception occurred. Attach the generator state string
    // so that the MARLEY job can be restarted from the last successful
    // event for easier debugging.
    // NOTE: The buffered event was fully created before the exception
    // occurred, so no exceptions are expected to be thrown in this block.
    // Just in case, we wrap it with an additional try/catch to inform
    // the user if writing out the buffered event fails.
    try {
      if ( buffered_event ) {
        marley::Generator::add_state_to_event( *buffered_event,
          cached_generator_state );
        for ( const auto& file : output_files ) {
          file->write_event( buffered_event.get() );
        }
      }
    } catch ( const std::exception& except ) {
      MARLEY_LOG( WARN, "app" ) << std::flush << "Output of buffered"
        " MARLEY event and its generator state failed";
      MARLEY_LOG( ERROR, "app" ) << except.what();
    }

    MARLEY_LOG( ERROR, "app" ) << std::flush << error.what();
  }

  return false;
}
