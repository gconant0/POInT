#include "gen_dna_funcs.h"
#include <iostream>

#ifndef ___LINKED_LIST_H___
#define ___LINKED_LIST_H___
//This is a template class for a linked list.  The "atoms" of that 
//list need to have an = and and == operator defined


template <class LIST_UNIT>
class List_Element
{
  public:
  List_Element();
  LIST_UNIT* item() {return(data);};
  List_Element<LIST_UNIT>& operator=(LIST_UNIT & assign_from) {(*data)=assign_from; return(*this);}; 
  int operator==(LIST_UNIT test);
  void set_previous (List_Element<LIST_UNIT> *new_previous) {previous_element=new_previous;};
  void set_next (List_Element<LIST_UNIT> *new_next) {next_element=new_next;};
  List_Element<LIST_UNIT>* next() { return(next_element);};
  List_Element<LIST_UNIT>* previous() { return(previous_element);};  
  ~List_Element() {delete data;};
 private:
  LIST_UNIT *data;
  List_Element<LIST_UNIT> *next_element, *previous_element;
};


template <class LIST_UNIT>
class List 
{
 public:
  List ();
  List_Element<LIST_UNIT>* add_to_list(LIST_UNIT &new_unit);
  void insert_after(LIST_UNIT new_unit, int pos);
  void remove_from_list(List_Element<LIST_UNIT> *remove_unit);
  int get_list_length() { return(list_length);};
  List_Element<LIST_UNIT>* find_data(LIST_UNIT *find);
  List_Element<LIST_UNIT>* get_nth_element(int element_num);
  void return_to_start();
  LIST_UNIT* get_current_item();
  LIST_UNIT* get_first_item();
  List_Element<LIST_UNIT>* get_current() {return (list_position);};
  List_Element<LIST_UNIT>* get_first() {return (list_start);};
  List_Element<LIST_UNIT>* get_next();
  ~List ();
 private:
  int list_length;
  List_Element<LIST_UNIT> *list_start, *list_end, *list_position;
};

#endif
