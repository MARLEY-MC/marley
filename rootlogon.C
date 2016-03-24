void rootlogon() {
  int return_code = gSystem->Load("libMARLEY.so");
  if (return_code == 0) std::cout << std::endl << "Successfully loaded MARLEY class dictionaries." << std::endl;
  else std::cout << "Error loading MARLEY class dictionaries." << std::endl;
}
