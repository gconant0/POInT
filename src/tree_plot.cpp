#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "stdio.h"
#include "math.h"
#include "plot.h"
#include "gen_dna_funcs.h"
#include "read_tree.h"
#include "tree.h"

using namespace::std;

int draw_rooted_tree(Exchange *curr_exchange, Tree *the_tree, string plotfile, double lnL, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss);
int draw_unrooted_tree(Exchange *curr_exchange, Tree *the_tree, string plotfile, double lnL);
LKMODEL get_model(string model_string);
void draw_branch(Branch *the_brn, double scalefact, double startx, double old_ycenter, double lowy, double highy);
int count_tips(Branch *curr_branch);




LKMODEL get_model(string model_string)
{
    cout<<"Matching "<<model_string<<endl;
	if (model_string.find("Jukes-Cantor 1969 model") != std::string::npos)
		return(JC_NUCLEOTIDE);
	if (model_string.find("Kimura 2-Parameter (1980) model") != std::string::npos)
		return(K2P_NUCLEOTIDE);
	if (model_string.find("Hasegawa-Kishino-Yano (1985) model") != std::string::npos)
		return(HKY_NUCLEOTIDE);
	if (model_string.find("Muse-Gaut/Goldman-Yang 1994 model") != std::string::npos)
		return(MG_94_HKY);
	if (model_string.find("Amino acid groups model") != std::string::npos)
		return(C_00_HKY);
	if (model_string.find("Linear Comb. of Amino Acid Props. model") != std::string::npos)
		return(LCAP_HKY);
	if (model_string.find("Three-state duplicate-loss model") != std::string::npos)
		return(DUPL);
	if (model_string.find("Four-state duplicate loss model with fixed duplicates") != std::string::npos)
			return(DUPL_FIX);
	if (model_string.find("Five-state duplicate loss model with parallel losses") != std::string::npos)
			return(DUPL_PARALLEL);
	if (model_string.find("Six-state duplicate loss model with parallel losses and fixed duplicates") != std::string::npos)
			return(DUPL_PARALLEL_FIX);
	if (model_string.find("Six-state duplicate loss model with subfunctionalizing losses and fixed duplicates") != std::string::npos)
			return(DUPL_PARALLEL_FIX_SUBF);
	if (model_string.find("Six-state duplicate loss model with only subfunctionaliztion") != std::string::npos);
	if(model_string.find("Six-state duplicate loss model with subfunction. losses, fixed dupls, and all rates") != std::string::npos)
	    return(DUPL_SUBF_3_RATE);
    if (model_string.find("Six-state duplicate loss model with subfunctionalizing losses and fixed duplicates and inferred tracking") != std::string::npos)
	  return(DUPL_PARALLEL_FIX_SUBF_NOSTATE);
	if (model_string.find("Six-state duplicate loss model without subfunction. losses and with fixed dupls, all rates, and inferred tracking") != std::string::npos)
		return(DUPL_2_RATE_NOSUBF_NOSTATE);
 }





int draw_rooted_tree(Exchange *curr_exchange, Tree *the_tree, string plotfile, double lnL, BOOL bitmap, BOOL IPC_call, std::stringstream *plot_ss)
{
  int i, j, k, thandle, size, num_size_index, xcoord1, xcoord2, ycoord1, ycoord2, realx, realy, mypos, name_len=0;
  double spacex, spacey, border=30, xcenter, ycenter, label_size, maxlen=0, sumlen,
		branch_scale, width, up_ratio, down_ratio, sum_ratio, centery;
  char dummy_string[50], write_string[100],size_string[25];
    char* buf_pt;
  Branch *next;
  FILE *outfile;
    std::size_t msize;
  
    if (IPC_call==FALSE)
        outfile=fopen(plotfile.c_str(), "w");
    else {
        msize=0;
        //printf ("At open buf = `%s', size = %zu\n", buf_pt, size);
        outfile= open_memstream (&buf_pt, &msize);
        
        //fprintf(outfile, "TESTING");
        fflush(outfile);
    }
   
    spacex=1000.0;
    spacey=500.0;
	
       
  xcenter=ycenter=spacex/2.0;
  
  num_size_index=6;
	
 
  /* set a Plotter parameter */
  /*pl_parampl ("PAGESIZE", "letter"); */ 

  /* create a Postscript Plotter that writes to standard output */
    if (bitmap == FALSE) {
        if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0) {
            fprintf (stderr, "Couldn't create Plotter\n");
            return 1;
        }
        pl_selectpl (thandle);       /* select the Plotter for use */
    }
    else {
        strcpy(size_string, "2000x1000");
        realx=2000;
        realy=1000;
                
        /* set a Plotter parameter */
        pl_parampl ("BITMAPSIZE", size_string);
       
        if ((thandle = pl_newpl ("png", stdin, outfile, stderr)) < 0) {
            fprintf (stderr, "Couldn't create Plotter\n");
            return 1;
        }
       
        pl_selectpl (thandle);
    }

  if (pl_openpl () < 0)       /* open Plotter */
    {
      fprintf (stderr, "Couldn't open Plotter\n");
      return 1;
    }
  pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
  pl_flinewidth (1.0);       /* line thickness in user coordinates */
  pl_pencolorname ("black");    /* path will be drawn in red */
  pl_erase ();                /* erase Plotter's graphics display */
  /*pl_fmove (600.0, 300.0);*/    /* position the graphics cursor */
    
    if (bitmap == TRUE) {
        pl_fontname ("HersheySans-Oblique"); /* choose a Postscript font */
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (18.0);
    }
    else {
        pl_fontname("Helvetica-Oblique");
        pl_ftextangle (0);         /* text inclination angle (degrees) */
        pl_ffontsize (18);
    }
  

    for(i=0; i<curr_exchange->get_num_taxa(); i++) {
        sumlen=0;
        next=(*the_tree)[i];
        sumlen+=next->expect_subs_site();
        while(next->get_parent() !=0) {
            next=next->get_parent();
            sumlen+=next->expect_subs_site();
            cout<<(*the_tree)[i]->get_name()<<": Updated sum to "<<sumlen<<endl;
        }

        if (sumlen > maxlen) {
            maxlen=sumlen;
        }
    
        pl_ffontsize (21.0);
        //cout<<"PRinting "<<the_brn->get_name()<<endl;
        width = pl_flabelwidth ((*the_tree)[i]->get_name());
        if (width>name_len) name_len=width;
  }


  branch_scale=(spacex-2*border-name_len)/maxlen;
	
  cout<<"Scaling branches by "<<branch_scale<<endl;

		
  draw_branch(the_tree->find_root(), branch_scale, 0+border, 0, 0, (spacey-2.0*border));
  

  double_to_string(dummy_string, 49, 4, lnL);
  strcpy(write_string, "Log-likelihood of this tree:  ");
  strcat(write_string, dummy_string);
 
  width = pl_flabelwidth (write_string);
  pl_fmove(border+width, spacey-border);
  pl_alabel ('c', 'c', write_string);       

#if 0
    
  if((curr_exchange->get_model() == DUPL_PARALLEL_FIX) || (curr_exchange->get_model() == DUPL_FIX) ||
 (curr_exchange->get_model() == DUPL_PARALLEL_FIX_SUBF) || (curr_exchange->get_model() == DUPL_SUBF_3_RATE) || 
(curr_exchange->get_model() == DUPL_SUBF_ONLY)|| (DUPL_PARALLEL_FIX_SUBF_NOSTATE)) {
      double_to_string(dummy_string, 49, 4, curr_exchange->get_dupl_fix_rate());
      strcpy(write_string, "Relative rate of duplicate fixation:  ");
      strcat(write_string, dummy_string);
 
      width = pl_flabelwidth (write_string);
      pl_fmove(border+width, border+width+10);
      pl_alabel ('c', 'c', write_string);       

  }
  if((curr_exchange->get_model() == DUPL_PARALLEL_FIX) || (curr_exchange->get_model() == DUPL_PARALLEL) || 
(curr_exchange->get_model() == DUPL_PARALLEL_FIX_SUBF)|| (curr_exchange->get_model() == DUPL_SUBF_3_RATE) ||
 (curr_exchange->get_model() == DUPL_SUBF_ONLY) || (DUPL_PARALLEL_FIX_SUBF_NOSTATE)) {
      double_to_string(dummy_string, 49, 4, curr_exchange->get_dupl_parallel_rate());
      strcpy(write_string, "Relative rate of duplicate convergence:  ");
      strcat(write_string, dummy_string);
 
      width = pl_flabelwidth (write_string);
      pl_fmove(border+width, border+width+20);
      pl_alabel ('c', 'c', write_string);       

  }

  if( (curr_exchange->get_model() == DUPL_SUBF_3_RATE) || (curr_exchange->get_model() == DUPL_SUBF_ONLY)) {
      double_to_string(dummy_string, 49, 4, curr_exchange->get_loss_rate_scale());
      strcpy(write_string, "Relative rate of duplicate loss after  convergence:  ");
      strcat(write_string, dummy_string);
 
      width = pl_flabelwidth (write_string);
      pl_fmove(border+width, border+width+20);
      pl_alabel ('c', 'c', write_string);       

  }
if( (curr_exchange->get_model() == DUPL_SUBF_3_RATE)) {
      double_to_string(dummy_string, 49, 4, curr_exchange->get_fix_rate_scale());
      strcpy(write_string, "Relative rate of duplicate fixation after convergence:  ");
      strcat(write_string, dummy_string);
 
      width = pl_flabelwidth (write_string);
      pl_fmove(border+width, border+width+20);
      pl_alabel ('c', 'c', write_string);       

  }
if( (curr_exchange->get_model() == DUPL_PARALLEL_FIX_SUBF_NOSTATE)) {
      double_to_string(dummy_string, 49, 4, curr_exchange->get_strand_switch_prob());
      strcpy(write_string, "Probability of a track switch between two adjacent genes:  ");
      strcat(write_string, dummy_string);
 
      width = pl_flabelwidth (write_string);
      pl_fmove(border+width, border+width+20);
      pl_alabel ('c', 'c', write_string);       

  }

#endif
	

  if (pl_closepl () < 0)      /* close Plotter */
    {
      fprintf (stderr, "Couldn't close Plotter\n");
      return 1;
    }

  pl_selectpl (0);            /* select default Plotter */
  if (pl_deletepl (thandle) < 0) /* delete Plotter we used */
    {
      fprintf (stderr, "Couldn't delete Plotter\n");
      return 1;
    }

    if (IPC_call==TRUE) {
        std::cout<<"Extracting from buffer\n";
        fflush(outfile);
        
        //printf ("At close buf = `%s', size = %zu\n", buf_pt, size);
        mypos=(size_t)0;
        
       
        while(mypos < msize) {
            plot_ss->put((char)buf_pt[mypos]);
            //cout<<"P: "<<mypos<<" = "<<(int)buf_pt[mypos]<<endl;
            ++mypos;
        }
        
    }

    fclose(outfile);
    
  return 0;
}


void draw_branch(Branch *the_brn, double scalefact, double startx, double old_ycenter, double lowy, double highy)
{

    char brnlen_string[20];
    string name, new_name;
    double width, endx, up_ratio, topy, bottomy, down_ratio, sum_ratio, centery;
    std::size_t found;
	
	if(the_brn->is_tip() ==FALSE) {
		up_ratio=count_tips(the_brn->get_child(0));
		down_ratio=count_tips(the_brn->get_child(1));
		sum_ratio=up_ratio+down_ratio;
		up_ratio=up_ratio/sum_ratio;
		down_ratio=down_ratio/sum_ratio;
		
		centery=(down_ratio*(highy-lowy))+lowy;
	}
	else centery=(0.5*(highy-lowy))+lowy;
	
    if (old_ycenter!=0) pl_fline(startx, old_ycenter, startx, centery);
	
    endx=startx+(the_brn->expect_subs_site()*scalefact);
    
	cout<<"Frame "<<lowy<<": "<<highy<<" X loc: "<<startx<<endl;
	
    cout<<"Drawing "<<startx<<" "<<centery<<" "<<endx<<" "<<centery<<" Up down ratio "<<up_ratio<<": "<<down_ratio<<endl;
    pl_fline(startx, centery, endx, centery);

    pl_ffontsize (17.0);
    double_to_string(brnlen_string, 19, 2, the_brn->expect_subs_site());
   
    width = pl_flabelwidth (brnlen_string);	
    pl_fmove(startx+(the_brn->expect_subs_site()*scalefact*0.5), centery+10);
    pl_alabel ('c', 'c', brnlen_string); 
    
    if(the_brn->is_tip() ==FALSE) {
		topy=centery+(0.5*(highy-lowy)*up_ratio);
		bottomy=centery-(0.5*(highy-lowy)*down_ratio);
		
		//pl_fline(endx, bottomy, endx, topy);
		
		
		draw_branch(the_brn->get_child(0), scalefact, endx, centery, centery, highy);
		draw_branch(the_brn->get_child(1), scalefact, endx, centery, lowy, centery);
    }
    else {
        pl_ffontsize (21.0);
        name = the_brn->get_name();
        
        found=name.find("_V3");
        if (found!=std::string::npos)
            name=name.substr(0, found);
        
        found=name.find("_V2");
        if (found!=std::string::npos)
            name=name.substr(0, found);
        
        found=name.find("_Col-0_v10");
        if (found!=std::string::npos)
            name=name.substr(0, found);
        
        found=name.find("_rerun");
        if (found!=std::string::npos)
            name=name.substr(0, found);
        
                
        std::size_t found = name.find("_");
        while (found!=std::string::npos) {
            new_name = name.substr(0, found);
            new_name+=" ";
            
            cout<<"HEad is "<<new_name<<"| found is "<<found<<endl;
            
            new_name += name.substr(found+1, name.length()-found);
            cout<<"Created name "<<new_name<<endl;
            name=new_name;
            cout<<"Current name is |"<<name<<endl;
            found=std::string::npos;
            std::size_t found = name.find("_");
        }
        
        name[0] = toupper(name[0]);
        
		//cout<<"PRinting "<<the_brn->get_name()<<endl;
		//width = pl_flabelwidth (the_brn->get_name());
        width = pl_flabelwidth (name.c_str());
		pl_fmove(endx+0.6*width, centery);
		//pl_alabel ('c', 'c', the_brn->get_name());
        pl_alabel ('c', 'c', name.c_str());
    }
}



int draw_unrooted_tree(Exchange *curr_exchange, Tree *the_tree, string plotfile, double lnL)
{
	cerr<<"ERROR: Cannot draw unrooted trees yet\n";
}



int count_tips(Branch *curr_branch)
{
	if (curr_branch->is_tip() == FALSE) return(count_tips(curr_branch->get_child(0))+count_tips(curr_branch->get_child(1)));
	else return(1);
}

