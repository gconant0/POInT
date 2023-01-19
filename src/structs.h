//This file contains all of the new enum and struct variables for
//the prosrch program
#ifndef ___STRUCTS_H___
#define ___STRUCTS_H___

#include "gen_dna_funcs.h"


struct pro_region {
  int end;
  double *rate_info;
};



struct pro_condition {
  double lnL_score, new_lnL_score, log_lnL_score;
  pro_region *regions;
};



struct list_pro_condition {
  int element_num;
  pro_condition current;
  list_pro_condition *last, *next;
};


//struct List_branch {
//  int element_num;
//  branch *current;
//  List_branch *last, *next;
//  BOOL end, empty;
//};


struct List_years {
  int current, element_num;
  List_years *last, *next;
  BOOL empty;
};


#endif

