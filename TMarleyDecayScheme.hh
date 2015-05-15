#pragma once
#include <map>
#include <regex>
#include "TMarleyLevel.hh"

class TMarleyDecayScheme {
  public:
    TMarleyDecayScheme(std::string nucid, std::string filename);
    std::string get_nuc_id() const;
    void set_nuc_id(std::string id);
    void add_level(const TMarleyLevel level);
    std::map<std::string, TMarleyLevel>* get_levels();
    TMarleyLevel* get_level(std::string energy);
    std::vector<TMarleyLevel*>* get_sorted_level_pointers();
    void print_report(std::ostream& ostr = std::cout);
    TMarleyLevel* get_pointer_to_closest_level(double E_level);
    void do_cascade(std::string initial_energy);
    void do_cascade(double initial_energy);
    void do_cascade(TMarleyLevel* initial_level);
    void do_weisskopf(int i);
    void print_latex_table(std::ostream& ostr = std::cout);

  private:
    std::string nuc_id;
    std::string file_name;
    std::map<std::string, TMarleyLevel> levels;
    std::vector<TMarleyLevel*> pv_sorted_levels;
    std::vector<double> sorted_level_energies;
    static bool compare_level_energies(TMarleyLevel* first,
      TMarleyLevel* second);
    std::string process_continuation_records(std::ifstream &file_in,
      std::string &record, std::regex &rx_cont_record) const;

    // Functions for calculating Weisskopf estimates
    double doubleFact(int i); // Calculating double factorials x!!
    double calcBE(int i);
    double calcBM(int i);
    double calcf(int i);
    double calcTE(int i, double dE, double BE_i, double f_i);
    double calcTM(int i, double dE, double BM_i, double f_i);

};
