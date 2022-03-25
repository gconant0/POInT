//Gavin Conant  6/26/1999
#ifndef ___ANNEAL_TEMPLATE_H___
#define ___ANNEAL_TEMPLATE_H___

#include <iostream>
#include "libranlib.h"
#include "gen_dna_funcs.h"
#include "anneal_exchange.h"
#include <string>

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


#define LN25PER log(.25)

#ifndef FLOAT_TOL
#define FLOAT_TOL 1e-15
#endif

using namespace::std;

double gauss_prob_between (double low, double high, double std_dev);

struct list_avg_e {
  long double avg_val, var_val;
  list_avg_e *next, *last;
};


template <class SPACE_UNIT>
class  Space_point {
 public:
  Space_point () {cerr<<"Call to default constructor of class Space_point\n";};
  Space_point (Exchange *cexchange) {score=new_score=log_score=0;
                              point = new SPACE_UNIT(cexchange);};
  Space_point& operator=(Space_point & assign_from);
  double score, new_score, log_score;
  SPACE_UNIT *point;
};


template <class SPACE_UNIT>
class Space_point_array
{
 public:
  Space_point<SPACE_UNIT> **elements;

  Space_point_array()  {cerr<<"Call to default constructor of class Space_point_array\n";};
  Space_point_array(int num_points, Exchange *curr_exchange);
  Space_point<SPACE_UNIT> * operator[](int element);
  ~Space_point_array();
 protected:
  int size;
};

//Template class implementation of the Anneal class
//SPACE_UNIT class must provide operator= and Anneal derived
//class must provide constructor, destructor and protected functions move, initialize_params, set_params
template <class SPACE_UNIT>
class Anneal {
 public:
  Space_point<SPACE_UNIT> *best;

  Anneal ();
  Anneal (Exchange *cexchange, int nwalkers);
  void set_relax(double new_relax);
  void set_boltz(double new_boltz);
	double get_boltz() {return(k_boltz);};
  void set_temp_fact(double new_fact);
  void set_bail(BOOL blf) {bail_long_fails=blf;};
  void set_rate_stdev(double stdev);
	void set_temp_end(double stp)  {temp_end=stp;};
    void set_save_last(int val) {save_last=val;};
  void use_pop_replace() {population_replace=TRUE;};
	int get_num_walkers()  {return(num_walkers);};
	Space_point<SPACE_UNIT>* get_walker_point(int n)  {return((*current_walkers)[n]);};
  double begin_sim();
    void set_outfile(string out)                {outfile=out;};
  ~Anneal();  

 protected:
  //Variables
    int save_last, num_walkers, tree_swap_frac, time_avg, times_so_far, steps_since_replace;
    double temp, initial_temp,  orig, k_boltz,
        relax_time, last_relax, std_dev, rate_stdev, energy_delta, temp_fact, temp_end;
    long double average_e, varience_e, last_var, net_avg_e, net_var_delta, heat_capacity, last_heat_cap;
    string outfile;
    BOOL population_replace, bail_long_fails;
    Space_point_array<SPACE_UNIT>  *current_walkers;
    Space_point<SPACE_UNIT> *pre_move, *curr_move;
    Exchange *curr_exchange;
    list_avg_e *start_avge, *end_avge;
 

  //Anneal functions
  void update_eng(int start, int end, double eng[], double &tot_eng);
  double acceptance_prob(double delta_val);
  double ln_boltz_prob(double energy, double temperature);
  //void average_energies(double &heat_cap, double &relax, double &avg_e, double &var_e, double temperature);
  virtual void average_energies();
  void temp_decline();


  //Procedures used in annealing
  virtual void move(int walker)=0;
  virtual BOOL evaluate_move(int walker_num);
  virtual void after_move(int walker_num)   {};
  virtual void set_params(Space_point<SPACE_UNIT> *values)=0;
  virtual void initialize_params(Space_point<SPACE_UNIT> *curr_condition)=0;
  virtual void output() {};
};



#endif





