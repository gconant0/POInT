//Gavin Conant 6/26/1999
#include "anneal_template.h"
#include <math.h>
#include <iostream>
#include <fstream>



#define DEBUG
#define RECIP_SQRT_TWO 0.707106781
#define SMALLER(a,b) (b>a ? a :b)
#define ONE_MINUS_RECIP_E (1.0-1/exp(1))


using namespace::std;

template <class SPACE_UNIT>
Space_point<SPACE_UNIT>& Space_point<SPACE_UNIT>::operator=(Space_point<SPACE_UNIT> & assign_from)
{
  score=assign_from.score;

  log_score=log(assign_from.score);
  (*point)=(*assign_from.point);
  return(*this);
}  //End Space_point<SPACE_UNIT>::operator=


template <class SPACE_UNIT>
Anneal<SPACE_UNIT>::Anneal()
{
  cerr<<"Default constructor for class Anneal called\n";
}


template <class SPACE_UNIT>
Space_point_array<SPACE_UNIT>::Space_point_array(int num_points, Exchange *curr_exchange)
{
  int i;
  size=num_points;
  
  elements=new Space_point<SPACE_UNIT> * [size];

  for(i=0; i<size; i++)
    elements[i]= new Space_point<SPACE_UNIT> (curr_exchange);

}


template <class SPACE_UNIT>
Space_point<SPACE_UNIT> * Space_point_array<SPACE_UNIT>::operator[](int element)
{
  if ((element >=0) && (element < size))
    return(elements[element]);
  else {
    cerr<<"Invalid index for Space_point array\n";
    return(0);
  }
}



template <class SPACE_UNIT>
Space_point_array<SPACE_UNIT>::~Space_point_array()
{
  int i;
  
  for(i=0; i<size; i++)
    delete elements[i];
  delete[] elements;
}


template <class SPACE_UNIT>
Anneal<SPACE_UNIT>::Anneal(Exchange *cexchange,int nwalkers)
{
  int i, j;
 
  population_replace=FALSE;
  curr_exchange=cexchange;
  num_walkers=nwalkers;
  current_walkers=new Space_point_array<SPACE_UNIT>(num_walkers, cexchange);

  best=new Space_point<SPACE_UNIT> (cexchange);
  pre_move= new Space_point<SPACE_UNIT> (cexchange);


  time_avg=50;
  times_so_far=0;

  end_avge=new  list_avg_e;
  start_avge=end_avge;
  for (i=1; i<time_avg; i++)
    {
      end_avge->avg_val=0;
      end_avge->next=new list_avg_e;
      end_avge->next->last=end_avge;
      end_avge=end_avge->next;
    }

  end_avge->avg_val=0;
  end_avge->next=start_avge;
  start_avge->last=end_avge;


  //save_last defines how many iterations will be made before cooling begins
  save_last=20;

  std_dev=(3*num_walkers);
  rate_stdev=0.05;
  last_heat_cap=0.1;
  relax_time=last_relax=1;
  k_boltz=0.003;
  energy_delta=0.00002;
  temp_fact=2.0;
	last_var=0;
	temp_end=0.15;
 
  net_var_delta=0;
  //Determines, on average, how many iterations occur between TBR actions (2 means every 2 iterations)
  tree_swap_frac=1;

  bail_long_fails=TRUE;
	

} //End Anneal::Anneal


template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::set_relax(double new_relax)
{
	relax_time=last_relax=new_relax;
}


template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::set_boltz(double new_boltz)
{
	k_boltz=new_boltz;
}


template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::set_temp_fact(double new_fact) 
{
	temp_fact=new_fact;
}

template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::set_rate_stdev(double stdev)
{
  rate_stdev=stdev;
}


template <class SPACE_UNIT>
double Anneal<SPACE_UNIT>::begin_sim()
  //This procedure actually does the simlulated annealing.
  //The system is first 'simmered' at the initial temperature for Anneal::save_last steps
  //Then a while loop kicks in, declining the temperature with each step.  When the temperature
  //is so low that the kinetic energ of the system only allows moves of less than 3 sites, or
  //when 1000 moves have be attemped without a successful one, the program considers the minima
  //found and reports the results

{
  int i, j, k, steps_since_accept=0, worst_loc;
	long long int net_its=0;
  double worst;
  BOOL step_cont=TRUE;
  


  temp=initial_temp=1;
  initialize_params((*current_walkers)[0]);
  initialize_params(best);
  (*current_walkers)[0]->score=(*current_walkers)[0]->new_score;
 
    if (save_last !=0) {
        move(0);
        (*current_walkers)[0]->score=(*current_walkers)[0]->new_score;
    }

  (*best)=(*(*current_walkers)[0]);
      
  for(i=1; i<num_walkers; i++)
    {  
      initialize_params((*current_walkers)[i]);
      (*current_walkers)[i]->score=(*current_walkers)[i]->new_score;
      move(i);
      (*current_walkers)[i]->score=(*current_walkers)[i]->new_score;


    
      if ((*current_walkers)[i]->score<best->score)
	*best=*(*current_walkers)[i];
    }
  
  
  orig=best->score;
  temp=initial_temp=2*fabs(best->score);
  cout<<"Starting score: "<<orig<<endl;
  //Filling the array that calculates the heat capacity and relaxation time
  for (i=1; i<save_last; i++)
    for (j=0; j<num_walkers; j++)
      {
		move(j);
		evaluate_move(j);
      }
  
	temp=initial_temp=sqrt(fabs(best->score))/temp_fact;
	average_energies();    
	net_avg_e=average_e;

	cout<<"Anneal::begin_sim--Begining Annealing.  Temperature: "<<temp
      <<" Current best score: "<<best->score
      <<" Relax time: "<<relax_time<<endl;   
	steps_since_replace=0;	

  //Begin simulated Annealing
  while (step_cont && temp>(temp_end*(initial_temp)))
    {
		for (i=0; i<num_walkers; i++)
		{
			move(i);
			if(evaluate_move(i)==TRUE)
				steps_since_accept=0; 
		}
      

		steps_since_accept++;
		average_energies();     
		temp_decline();
		if (net_its%40000==0) {
			cout<<"Temp: "<<temp<<" RT: "<<relax_time<<" HC: "<<heat_capacity<<"\t";
		for(i=0; i<num_walkers; i++) cout<<(*current_walkers)[i]->score<<"\t";
			cout<<endl;
		}
		++net_its;

		if (population_replace==TRUE) {
		  if(steps_since_replace <4000)
		    steps_since_replace++;
		  else {
		    steps_since_replace=0;
		    worst=(*current_walkers)[0]->score;
		    worst_loc=0;
		    for (i=1; i<num_walkers; i++)
		      {
			if((*current_walkers)[i]->score>worst) {
			  worst=(*current_walkers)[i]->score;
			  worst_loc=i;
			}
		      }
		    cout<<"Replacing "<<(*current_walkers)[worst_loc]->score<<" with "<<best->score<<endl;
		    (*(*current_walkers)[worst_loc])=(*best);
		    
		  }

		  if (bail_long_fails == TRUE) {
			if (steps_since_accept>=(1500000/num_walkers))
				step_cont=FALSE;
		  }
		}
    } //End Anneal While loop 
  
 
  cout<<"Total iterations: "<<net_its<<endl;
  
  // cout<<"Annealing complete: Final score: "<<best->score<<endl;

  
  cout<<"Annealing complete: Final temp: "<<temp<<endl;
  set_params(best);

  output();
  return(best->score);


}  //End Anneal::begin_sim




template <class SPACE_UNIT>
Anneal<SPACE_UNIT>::~Anneal()
{
  delete current_walkers;
}





//Anneal Private functions
template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::update_eng(int start, int end, double eng[], double &tot_eng)
{
  for(int i=start; i<end; i++)
    tot_eng-=eng[i]; 
}  //End Anneal::update_eng



template <class SPACE_UNIT>
double Anneal<SPACE_UNIT>::acceptance_prob(double delta_val)
{
  if(delta_val<=0.0)
    return(1.0);
  else
    return(exp(-(delta_val/(temp*k_boltz))));
} //End Anneal::acceptance_prob
 


template <class SPACE_UNIT>
double Anneal<SPACE_UNIT>::ln_boltz_prob(double energy, double temperature)
{
  return(-energy/(k_boltz*temperature));
} //End Anneal::acceptance_prob



template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::average_energies()
{
  int i;
  long double ln_boltz_p, z=1.0;
  
  average_e=varience_e=1.0;

  for (i=0; i<num_walkers; i++)
    {
		(*current_walkers)[i]->log_score=log(fabs((*current_walkers)[i]->score));
      ln_boltz_p=ln_boltz_prob(fabs((*current_walkers)[i]->score), temp);
      z=logadd(ln_boltz_p, z);
      average_e=logadd((*current_walkers)[i]->log_score+ln_boltz_p, average_e);
      varience_e=logadd((2*(*current_walkers)[i]->log_score)+ln_boltz_p, varience_e);
    }

	//cout<<"AVGE: "<<average_e<<" Z: "<<z<<" VAR: "<<varience_e<<" LBP: "<<ln_boltz_p<<endl;
  average_e=average_e-z;
  varience_e=varience_e-z;
  //heat_capacity= exp(logsubtract(varience_e, 2*(average_e)))/(k_boltz*SQR(temp));

  if (end_avge->avg_val!=0) {
    net_avg_e-=end_avge->last->avg_val-end_avge->avg_val;
    net_var_delta-=(end_avge->var_val-(end_avge->avg_val*end_avge->avg_val));
  }
  else
    times_so_far++;
  
  end_avge->avg_val=exp(average_e);
  end_avge->var_val=exp(varience_e); 
  start_avge=end_avge;
  end_avge=end_avge->last;
  


  if (fabs(start_avge->next->avg_val-0)>0.001)
    { 
      net_avg_e+=start_avge->avg_val-start_avge->next->avg_val;
    }

  net_var_delta+=start_avge->var_val-(start_avge->avg_val*start_avge->avg_val);
 heat_capacity= ((double)net_var_delta/times_so_far)/((double)k_boltz*(temp*temp));
 //cout<<"Avg e: "<<average_e<<" Var e: "<<varience_e<<" Z: "<<z<<" net d: "<<net_var_delta<<" Heat cap: "<<heat_capacity<<" ";

 //relax_time=-1/(net_avg_e/times_so_far);
  
  //cout<<" relax time: "<<relax_time<< " last: "<<last_relax<<endl;

  // if (relax_time>0 && relax_time<1e20)  
    // last_relax=relax_time;
  //  else
      relax_time=last_relax;

  if ((heat_capacity>1e-8) && (heat_capacity<=1e5))  last_heat_cap=heat_capacity;
  else heat_capacity=last_heat_cap;
    
}  //End Anneal::average_energies




template <class SPACE_UNIT>
void Anneal<SPACE_UNIT>::temp_decline()
{ 
  temp-=(energy_delta*temp)/(relax_time*sqrt(heat_capacity)); 
  //cout<<" New temp: "<<temp<<endl;
}  //End Anneal::temp_decline



template <class SPACE_UNIT>
BOOL Anneal<SPACE_UNIT>::evaluate_move(int walker_num)
{
  int i,j;
  double delta;
  BOOL accept=FALSE, good=FALSE;
    
    //cout<<"Evaluating: "<<(*current_walkers)[walker_num]->new_score<<" verses old :"<<(*current_walkers)[walker_num]->score<<endl;
 
  delta=(*current_walkers)[walker_num]->new_score-(*current_walkers)[walker_num]->score;

//Evaluate the quality of that move
  if (delta<=0.0)
    {
      accept=TRUE;
    }
  else
    if( ranf()<=acceptance_prob(delta))
      accept=TRUE;

  if (accept==TRUE)
    {
       // cout<<"Accepted: "<<delta<<"\n";
	  after_move(walker_num);

	  if ((*current_walkers)[walker_num]->new_score<best->score)
			good=TRUE;

      (*current_walkers)[walker_num]->score=(*current_walkers)[walker_num]->new_score;
      if (good==TRUE)
		{
			cout<<"Found a new best score: "<<(*current_walkers)[walker_num]->score<<"\n"<<flush;
			(*best)=(*(*current_walkers)[walker_num]);  
		}    
     // else
			//cout<<"Accepted: "<<(*current_walkers)[walker_num]->score<<" (prob was "<<acceptance_prob(delta)<<"\n";
      
	  (*current_walkers)[walker_num]->log_score=log((*current_walkers)[walker_num]->score);

     }   
  
  else
    {
      set_params(pre_move);
      *(*current_walkers)[walker_num]=(*pre_move);
    }
 return(good);
}  //End Anneal::evaluate_move


















double gauss_prob_between (double low, double high, double std_dev)
{
  low=(low/std_dev)*RECIP_SQRT_TWO;
  high=(high/std_dev)*RECIP_SQRT_TWO;
  if (low<0 &&high>0)
    return(0.5*high+0.5*low);
  else
    return(0.5*high-0.5*low);
}
























