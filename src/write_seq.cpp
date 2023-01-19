//Copyright 1999-2002 Gavin Conant

#include <iostream>
#include <fstream>
#include "string.h"
#include "write_seq.h"

using namespace::std;
Write_Sequence::Write_Sequence()
{
  cerr<<"Called default constructor for class Write_Sequence\n";
}


Write_Sequence::Write_Sequence(const char *output_file, DATATYPE cdata)
{
  failed=FALSE;
  set_output_type(cdata);
  outfile.open(output_file);
  if (outfile.fail())
    {
      cerr<<"Cannot create file "<<output_file<<endl;
      failed=TRUE;
    }
}


Write_Sequence::Write_Sequence(const char *output_file, DATATYPE cdata, BOOL append)
{
	failed=FALSE;
	set_output_type(cdata);
	if (append == FALSE)
		outfile.open(output_file);
	else
		outfile.open(output_file, ios::app);
    
    //cout<<"opened output file "<<output_file<<" Append is "<<append<<endl;
	
	if (outfile.fail())
    {
		cerr<<"Cannot create file "<<output_file<<endl;
		failed=TRUE;
    }
}


void Write_Sequence::write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset)
{ 
  int i;

 
  for (i=0; i<ntaxa; i++)
      write_sequence(nchars, (*dataset)[i]);
    
    //cout<<"Wrote sequence with lengths. CLosing\n";
 
  outfile.close();	
}


void Write_Sequence::write_dataset(int ntaxa, Sequence_dataset *dataset)
{
    int i;

 
    for (i=0; i<ntaxa; i++)
      write_sequence((*dataset)[i].Sequence_size(), (*dataset)[i]);
    
    //cout<<"Wrote sequence taxa only. CLosing\n";
    outfile.close();
}

void Write_Sequence::set_output_type(DATATYPE cdata)
{
  curr_data=cdata;
  switch (cdata)
    {
    case PROTEIN:
      convert_int=&num_to_aa;
      break;
    case NUCLEIC:
      convert_int=&num_to_base;
      break;
	case DUPL_STATUS:
		convert_int=&num_to_dupl_data;
    default:
        convert_int=&num_to_base;
        break;
        
    }
} //End Write_Sequence::set_output_type

BOOL Write_Sequence::fail()
{
  return(failed);
}




//Write_PIR functions
Write_PIR::Write_PIR() : Write_Sequence()
{
}

Write_PIR::Write_PIR(const char *output_file, DATATYPE cdata) : Write_Sequence(output_file, cdata)
{
}


void Write_PIR::write_sequence(int nchars, Molecule_Sequence &seq_data)
{
  int i, so_far, linelen;
  so_far=0;
  linelen=60;


  switch (curr_data) {
    case PROTEIN:    
       outfile<<">P1;"<<seq_data.Sequence_name()<<endl<<endl;
       break;
    case NUCLEIC:
      outfile<<">DL;"<<seq_data.Sequence_name()<<endl<<endl;
      break;
	case DUPL_STATUS:
		outfile<<">DS;"<<seq_data.Sequence_name()<<endl<<endl;
    default:
        outfile<<">UK;"<<seq_data.Sequence_name()<<endl<<endl;
    }
 
  
  while (nchars-so_far>linelen)
    {
      for(i=0; i<linelen; i++)
	outfile<<convert_int(seq_data[i+so_far]);
      outfile<<endl;
      so_far+=linelen;
    }
  
  for(i=0; i<nchars-so_far; i++)
      outfile<<convert_int(seq_data[i+so_far]);
  
  outfile<<endl<<"*"<<endl;
}  //End Write_PIR::write_sequence





//Write_Nexus functions
Write_Nexus::Write_Nexus() : Write_Sequence()
{
}

Write_Nexus::Write_Nexus(const char *output_file, DATATYPE cdata) : Write_Sequence(output_file, cdata)
{
  has_comments=FALSE;
}


Write_Nexus::Write_Nexus(const char *output_file, DATATYPE cdata, char** data_comments, int c_lines) : Write_Sequence(output_file, cdata)
{
  int i;

  has_comments=TRUE;
  comment_lines=c_lines;
  comments=new char *[comment_lines];
  for(i=0; i<comment_lines; i++)
    {
     comments[i]=new char [200];
     strcpy(comments[i], data_comments[i]);
    }
}

void Write_Nexus::write_dataset(int ntaxa, Sequence_dataset *dataset)
{
  write_dataset(ntaxa, (*dataset)[0].Sequence_size(), dataset);
}




void Write_Nexus::write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset)
{
  int i,j, namelen, stop, local_stop;
  char scan_name[700];
  BOOL has_under, first_line;

  if (has_comments==FALSE)
    outfile<<"#NEXUS"<<endl<<endl<<"Begin data;"<<endl<<"\tDimensions ntax="
	   <<ntaxa<<" nchar="<<nchars<<";"<<endl;
  
  else
    {
      outfile<<"#NEXUS"<<endl<<endl;
      for (i=0; i<comment_lines; i++)
	outfile<<"[!"<<comments[i]<<"]\n";
      outfile<<endl<<"Begin data;"<<endl<<"\tDimensions ntax="
	     <<ntaxa<<" nchar="<<nchars<<";"<<endl;
  
    }

  switch (curr_data)
    {
    case PROTEIN:    
      outfile<<"\tFormat datatype=protein";
      break;
    case NUCLEIC:
      outfile<<"\tFormat datatype=nucleotide";    
      break;
	case DUPL_STATUS:
		//Note that PAUP won't understand this
		outfile<<"\tFormat datatype=dupl_status";
    default:
        //Note that PAUP won't understand this
        outfile<<"\tFormat datatype=UNKNONWN";
    }
  outfile<<" gap=? missing=- matchchar=.;"<<endl<<"\tMatrix"<<endl;

  stop=2;

  for (i=0; i<ntaxa; i++)
    {
      namelen=strlen((*dataset)[i].Sequence_name());
      strcpy(scan_name, (*dataset)[i].Sequence_name());

      has_under=FALSE;
      
      for (j=0; j<namelen; j++)
	if ((scan_name[j]=='_')  || (scan_name[j]=='-'))
	  has_under=TRUE;
	
      if (namelen>=stop && has_under==FALSE)
	 stop=namelen;
      else if (namelen>=stop-2 && has_under==TRUE)
	    stop=namelen+2;
	  
	
    }

 
    for (i=0; i<ntaxa; i++)
      {
	namelen=strlen((*dataset)[i].Sequence_name());
	strcpy(scan_name, (*dataset)[i].Sequence_name());

	has_under=FALSE;
       first_line=TRUE;

       for (j=0; j<namelen; j++)
           if ((scan_name[j]=='_')  || (scan_name[j]=='-'))
               has_under=TRUE;
       
       if (has_under==TRUE)
	 {
	   Write_Nexus::outfile<<'\'';
	   local_stop=stop-2;
	 }
       else
	 local_stop=stop;

	 for(j=0; j<local_stop; j++)
	   {
	     if (j<namelen)
	       outfile<<scan_name[j];
	     else
	       outfile<<' ';
	     if (has_under==TRUE && j==namelen-1)
	       outfile<<'\''; 
	   }

	 outfile<<' ';

	 for(j=0; j<nchars; j++)
	   {
	     outfile<<convert_int((*dataset)[i][j]);
	     if((j+1-(60-stop-1))%60==0 && first_line==FALSE)
	       outfile<<endl;
	     else if ((j+1)%(60-stop-1)==0 && first_line==TRUE)
	       {
		 outfile<<endl;
		 first_line=FALSE;
	       }
	   }
	 outfile<<endl;
     }

   outfile<<"\t;"<<endl;
   outfile<<"End;"<<endl;

   outfile.close();
   
} //End Write_Nexus::write_dataset



void Write_Nexus::write_sequence(int nchars, Molecule_Sequence &seq_data)
{
  cerr<<"Cannot write single sequences in Nexus format\n";  
}


Write_Nexus::~Write_Nexus()
{
  int i;
  if (has_comments==TRUE)
    {
      for(i=0; i<comment_lines; i++)
	delete[] comments[i];
      delete[] comments;
    }
}

//Write_Phylip_interleave functions
Write_Phylip_interleave::Write_Phylip_interleave() : Write_Sequence()
{
}


Write_Phylip_interleave::Write_Phylip_interleave(const char *output_file, DATATYPE cdata) : Write_Sequence(output_file, cdata)
{
}


void Write_Phylip_interleave::write_dataset(int ntaxa, Sequence_dataset *dataset)
{
  write_dataset(ntaxa, (*dataset)[0].Sequence_size(), dataset);
}



void Write_Phylip_interleave::write_dataset(int ntaxa, int nchars, Sequence_dataset *dataset)
{
 int i,j, k, linelen, so_far=0, namelen, remain_group;
 char scan_name[700];

 outfile<<"     "<<ntaxa<<"    "<<nchars<<endl;

 if (nchars-so_far<50)
   linelen=nchars-so_far;
 else
   linelen=50;
 
 for (i=0; i<ntaxa; i++)
   {
     namelen=strlen((*dataset)[i].Sequence_name());
     strcpy(scan_name, (*dataset)[i].Sequence_name());
     
     for(j=0; j<10; j++)
       {
	 if (j<namelen)
	   outfile<<scan_name[j];
	 else
	   outfile<<' ';
       }
     outfile<<' ';
	    
	    
     for(j=so_far; j<linelen+so_far; j++)
       {
	 outfile<<convert_int((*dataset)[i][j]);
	 if((j+1)%10==0)
	   outfile<<' ';

       }
     outfile<<endl;
	    
   }	

 outfile<<endl;
 
 so_far+=linelen;
 
 if(nchars-so_far%50==0)
   remain_group=(int)((nchars-so_far)/50);
 else  
   remain_group=(int)((nchars-so_far)/50)+1;
 
       
 for (k=0; k<remain_group; k++)
   {
     if (nchars-so_far<50)
       linelen=nchars-so_far;
     else
       linelen=50;
     
     for (i=0; i<ntaxa; i++)
       {
	 outfile<<"           ";
	 
	 for(j=so_far; j<linelen+so_far; j++)
	   {
	     outfile<<convert_int((*dataset)[i][j]);
	     if((j+1)%10==0)
	       outfile<<' '; 
	   }
	 outfile<<endl;
       }	
     
     so_far+=linelen;
     outfile<<endl;
   }
 outfile.close();	

}


void Write_Phylip_interleave::write_sequence(int nchars, Molecule_Sequence &seq_data)
{
  cerr<<"Cannot write single sequences in interleave format\n";
}



Write_FASTA::Write_FASTA() : Write_Sequence()
{
}

Write_FASTA::Write_FASTA(const char *output_file, DATATYPE cdata) : Write_Sequence(output_file, cdata)
{
}


Write_FASTA::Write_FASTA(const char *output_file, DATATYPE cdata, BOOL append) : Write_Sequence(output_file, cdata, append)
{
}

void Write_FASTA::write_sequence(int nchars, Molecule_Sequence &seq_data)
{
    int i, linelen=60;


    outfile<<'>'<<seq_data.Sequence_name()<<endl;

    for(i=0; i<nchars; i++) {
        if ( ((i%linelen) == 0) &&( i !=0) && (i!= nchars)) outfile<<endl;
        outfile<<convert_int(seq_data[i]);
        
    }
  
  
    outfile<<endl<<endl;
    
}  //End Write_FASTA::write_sequence
