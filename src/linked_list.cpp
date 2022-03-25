#include "linked_list.h"
using namespace std;


template <class LIST_UNIT>
List_Element<LIST_UNIT>::List_Element()
{
  next_element=0;
  previous_element=0;
  data=new LIST_UNIT;
}



template <class LIST_UNIT>
int List_Element<LIST_UNIT>::operator==(LIST_UNIT test)
{
  if ((*data)==test)
    return(0);
  else
    return(-1);
}


template <class LIST_UNIT>
List<LIST_UNIT>::List ()
{
 // list_start=new List_Element<LIST_UNIT>;
  list_end=list_start=0;
  list_position=list_start;
  list_length=0;
}





template <class LIST_UNIT>
List_Element<LIST_UNIT>* List<LIST_UNIT>::add_to_list(LIST_UNIT &new_unit)
{
	if (list_length==0) {  
		list_start=new List_Element<LIST_UNIT>;
		list_end=list_position=list_start;
		(*list_start)=new_unit;
	}
	else
    {
      list_end->set_next(new List_Element<LIST_UNIT>());
      list_end->next()->set_previous(list_end);
      list_end=list_end->next();
		(*list_end)=new_unit;
    }
  list_length++;
  return(list_end);
}


template <class LIST_UNIT>
void List<LIST_UNIT>::insert_after(LIST_UNIT new_unit, int pos)
{
	List_Element<LIST_UNIT> *temp, *prev_pos;

	if (pos != -1)
		prev_pos=get_nth_element(pos);
	else
		prev_pos=0;


	if (prev_pos == 0) {
		//cout<<"Inserting at start\n";
		temp=list_start;
		list_start=new List_Element<LIST_UNIT>;
		if (temp != 0)  {
			temp->set_previous(list_start);
			list_start->set_next(temp);
		}
		else 
			list_end=list_position=list_start;
		(*list_start)=new_unit;
	}
	else {
	//cout<<"Inserting in middle or end\n";
	  temp=prev_pos->next();
	  prev_pos->set_next(new List_Element<LIST_UNIT>());
	  prev_pos->next()->set_previous(prev_pos);
	  if (temp != 0) {
		temp->set_previous(list_end->next());
		prev_pos->next()->set_next(temp);
	  }
	  else 
		list_end=list_end->next();
      
	  (*prev_pos->next())=new_unit;	
	}
	list_length++;
}
 


template <class LIST_UNIT>
void List<LIST_UNIT>::remove_from_list(List_Element<LIST_UNIT> *remove_unit)
{
	//Make sure we're not deleting the place the pointer currently is
	if(list_position==remove_unit) {
		if (list_position->next() !=0)
			list_position=list_position->next();
		else
			list_position=list_position->previous();
	}
  
	if (remove_unit==list_start)
	{
		if (remove_unit->next() != 0)
			remove_unit->next()->set_previous(0);
		list_start=remove_unit->next();
    }

	else
    {
      remove_unit->previous()->set_next(remove_unit->next());
      if (list_position==remove_unit)
	  list_position=remove_unit->previous();

      if (list_end==remove_unit)
	  list_end=remove_unit->previous();
      else
	  remove_unit->next()->set_previous(remove_unit->previous());
    }

  delete remove_unit;
  list_length--;
}




template <class LIST_UNIT>
List_Element<LIST_UNIT>* List<LIST_UNIT>::find_data(LIST_UNIT *find)
{
  List_Element<LIST_UNIT> *orig, *value=0;

  orig=list_position;

  list_position=list_start;

  while(!((*list_position)==(*find)) && list_position->next()!=0)
    list_position=list_position->next();
  

  if ((*list_position)==(*find))
      value=list_position;

  list_position=orig;
  return(value);
}




template <class LIST_UNIT>
List_Element<LIST_UNIT>* List<LIST_UNIT>::get_nth_element(int element_num)
{
  int i;
  List_Element<LIST_UNIT> *orig, *value;

  orig=list_position;

  if (element_num>list_length)
    cerr<<"Requested element larger than size of list\n";
  else
    {
      list_position=list_start;
    
      for (i=0; i<element_num; i++)
		list_position=list_position->next();

      value=list_position;
    }
  list_position=orig;
  return(value);
}




template <class LIST_UNIT>
void List<LIST_UNIT>::return_to_start()
{
  list_position=list_start;
}




template <class LIST_UNIT>
LIST_UNIT* List<LIST_UNIT>::get_current_item()
{
	if (list_position != 0)
		return(list_position->item());
	else
		return(0);
}
 
template <class LIST_UNIT>
LIST_UNIT* List<LIST_UNIT>::get_first_item()
{
	return(list_start->item());
}

template <class LIST_UNIT>
List_Element<LIST_UNIT>* List<LIST_UNIT>::get_next()
{
  list_position=list_position->next();
  return(list_position);
}



template <class LIST_UNIT>
List<LIST_UNIT>::~List ()
{
  int i;
  List_Element<LIST_UNIT> *the_next;

  list_position=list_start;

  for (i=0; i<list_length; i++)
    {
		the_next=list_position->next();
		delete list_position;
		list_position=the_next;
    }

} //End List<LIST_UNIT>::~List

