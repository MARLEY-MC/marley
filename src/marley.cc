#include "marley/CommandHandler.hh"

int main( int argc, char* argv[] ) {

  // Initalize the CommandHandler which will manage calling
  // the functions requested on the command line
  marley::CommandHandler ch( argc, argv );

  bool ok = ch.execute();
  if ( !ok ) return 1;
  return 0;
}
