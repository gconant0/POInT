//Copyright 1999-2002 Gavin Conant

#include "read_seq.h"
#include "gen_dna_funcs.h"
#include <string>
#include <iostream>
#include <fstream>



BOOL test_interleave(char *filename)
{
  int num_taxa, num_chars;
  char dum_name[800], line[3000], test;
  ifstream infile;

  infile.open(filename);
  if(infile.fail()) {
    cerr<<"Error in Read_Sequence test_interleave: file "<<filename<<" not found\n";
    return(FALSE);
  }

  infile>>num_taxa>>num_chars;
  infile.get(dum_name, 799);
  infile.get(test);
  infile.get(line, 2999);
  infile.get(test);
  infile.get(dum_name, 799);
  infile.get(test);
  
  infile.close();

if (test == '\t')
    return(TRUE);
  else
    return(FALSE);
  

}



Read_Sequence::Read_Sequence()
{
	std_in=FALSE;
	finished = FALSE;
	var_stream=0;
	in_seq=0;
	sizes=0;
	build_taxa= new List_Molecule_Sequence;
	taxa_start=build_taxa;
	build_taxa->next=0;
	ntaxa=0;
	size=0;
}



int Read_Sequence::get_num_taxa()
{
  return(ntaxa);
}



int Read_Sequence::get_num_chars()
{
  return(size);
}




BOOL Read_Sequence::protein_data()
{
  return(aa_data);
}

#if 0
Read_Sequence::~Read_Sequence()
{
    if (build_taxa !=0) {
       build_taxa=taxa_start;
      
      while (build_taxa->next!=0)
        {
          temp_taxa=build_taxa->next;
            delete build_taxa->current;
          delete build_taxa;
          build_taxa=temp_taxa;
        }
      delete build_taxa;
    }

} //End Read_Sequence::~Read_Sequence
#endif
//End Read_Sequence Public functions





//Read_Sequence Private Functions
void Read_Sequence::get_sequence_array_size_known(Molecule_Sequence *new_sequence)
{ 
	int i;
	char inbase;
      
	inbase=next_char();
        
  while(valid_input(inbase)!=TRUE)
	  inbase=next_char();

    //infile.get(inbase);
  new_sequence->Assign_site(0, convert_char(inbase));
 
  i=1;

  //while((i<size) && (!infile.eof()))
	while((i<size) && (!read_finished()))
    {	
		inbase=next_char();
      
      if (valid_input(inbase)==TRUE)
	{
	  new_sequence->Assign_site(i, convert_char(inbase));
	  i++;
	} 
    }
  if (i<size-1)
    while (i<size) 
      {
	  new_sequence->Assign_site(i, convert_char('-'));
	  i++;
		  }
} //End Read_Sequence::get_sequence_array_size_known



void Read_Sequence::get_sequence_array(Molecule_Sequence *&new_sequence)
{
   int size_so_far=0, sitecount;
	char inbase;
	sequence *temp;

	inbase=next_char();
	
	in_seq_start=new sequence;	
	in_seq=in_seq_start;
	in_seq->next_site=0;

	while((valid_input(inbase)!=TRUE) && (inbase != '>') && !read_finished())
		inbase=next_char();

 
	if ((inbase == '>') || (read_finished())) {
		new_sequence=new Molecule_Sequence(0);
		return;
	}

	in_seq->site=convert_char(inbase);
 
	inbase=next_char();
	size_so_far++;
   
	while((inbase!= '>') && (inbase !='*') && (!read_finished()))
	{	
			if (valid_input(inbase)==TRUE)
			{
				get_next_sequence_element(in_seq);
				in_seq->site=convert_char(inbase);
				size_so_far++;
			}
		inbase=next_char();
	}

	size=size_so_far;

	in_seq=in_seq_start;
  
	new_sequence=new Molecule_Sequence(size);

	for(sitecount=0; sitecount<size; sitecount++)
		{
		  new_sequence->Assign_site(sitecount, in_seq->site);
		  in_seq=in_seq->next_site;   
		}    
	 
	in_seq=in_seq_start;
	
	while(in_seq != 0) {
		temp=in_seq->next_site;
		delete in_seq;
		in_seq=temp;
	}
	
	last_char=inbase;
} //End Read_Sequence::get_sequence_array
  


void Read_Sequence::set_datatype(DATATYPE cdata)
{
  switch (cdata)
    {
    case NUCLEIC:
      aa_data=FALSE;
	  ctype=cdata;
      valid_input=&is_base;
      convert_char=&readchar_to_base;
      break;
    case PROTEIN:
      aa_data=TRUE;
	  ctype=cdata;
      valid_input=&is_aa;
      convert_char=&readchar_to_aa;
      break;
	case DUPL_STATUS:
		aa_data=FALSE;
		ctype=cdata;
		valid_input=&is_dupl_data;
		convert_char=&readchar_to_dupl;
		break;
	case SNP_DATA:
		aa_data=FALSE;
		ctype=cdata;
		valid_input=&is_snpstate;
		convert_char=&readsnp_to_snpstate;
        break;
    default:
        cerr<<"ERROR: Trying to read datafile of non-readable type: e.g., codons or arbitrary files\n";
        aa_data=FALSE;
        ctype=cdata;
        break;
    }

} //End Read_Sequence::set_datatype



void Read_Sequence::get_next_list_element(List_Molecule_Sequence *&thelist)
{
  if (thelist->next==0)
    {
      thelist->next= new List_Molecule_Sequence;
      thelist->last=thelist;
      thelist->next->next=0;
    }
 
  thelist=thelist->next;
}



void Read_Sequence::get_next_sequence_element(sequence *&thelist)
{
  if (thelist->next_site==0)
    {
      in_seq->next_site=new sequence;
      in_seq->next_site->next_site=0;
    }

  in_seq=in_seq->next_site;
} //End Read_Sequence::get_next_sequence_element
//End of Read_Sequence protected functions


BOOL Read_Sequence::read_finished()
{
	if (std_in == FALSE) { 
		if (var_stream->eof()) finished=TRUE;	
	}
	else {
		if (inbase == '#') finished=TRUE;
	}
	return(finished);
}


char Read_Sequence::next_char()
{
	var_stream->get(inbase);
	read_finished();
	return(inbase);
}

//Read_PIR Public functions

Sequence_dataset * Read_PIR::get_dataset(int &num_taxa, int &num_sites, const char *filename, BOOL is_aa_seq)
{
	if(is_aa_seq == TRUE)
		return(get_dataset(num_taxa, num_sites, filename, PROTEIN));
	else
		return(get_dataset(num_taxa, num_sites, filename, NUCLEIC));
}



Sequence_dataset * Read_PIR::get_dataset(int &num_taxa, int &num_sites, const char *filename, DATATYPE type)
{
	var_stream=&file_stream;
	file_stream.open(filename);
	
	if (file_stream.fail())
	{
        cerr<<"Read_PIR::get_dataset: Cannot find file "<<filename<<endl;
		return(0);
	}
	else
		return(get_dataset(num_taxa, num_sites, type));
}


Sequence_dataset * Read_PIR::get_dataset(int &num_taxa, int &num_sites, DATATYPE type)
{
	int i;
	char names[700];
	Sequence_dataset *the_data;
		
	 
	set_datatype(type);

	if (var_stream == 0) {
		var_stream = &cin;
		std_in=TRUE;
	}

	build_taxa= taxa_start;

	
	num_taxa=0;
	  
	inbase=next_char();
	  
	while (inbase!='>')
		inbase=next_char();
	while (inbase!=';')
		inbase=next_char();
		  
	var_stream->getline(names, 698);
	//infile.getline(names, 698);
		  
	get_sequence_array( build_taxa->current);
		
	build_taxa->current->Assign_name(names);
		  
	num_sites=size;
	get_next_list_element(build_taxa);
		  
	num_taxa++;
		  
	//while (!infile.eof())
	while (!read_finished())
	{
		inbase=next_char();
		
	  while (inbase!='>' && !read_finished())
		  inbase=next_char();
		  
	  while (inbase!=';' && !read_finished())
		 inbase=next_char();
		  
		  
	  if (!read_finished())
		{
			var_stream->getline(names, 698);
		  //infile.getline(names, 698);
			  
		  build_taxa->current= new Molecule_Sequence(size);
	  
		  get_sequence_array_size_known(build_taxa->current);
		  build_taxa->current->Assign_name(names);
			  
		  get_next_list_element(build_taxa);  
		  num_taxa++;
		}
	}
		  
	ntaxa=num_taxa;
		  
	 the_data= new Sequence_dataset(num_taxa, size, ctype);
		  
	  build_taxa=taxa_start;
		  
	  for(i=0; i<num_taxa; i++)
		{
		  (*the_data)[i]=*build_taxa->current;
		  
		  if (i!=num_taxa-1)
			build_taxa=build_taxa->next;
		}
		  
		if (std_in == FALSE)
			file_stream.close();
		  //infile.close();
		  return(the_data);
} //End Read_PIR::get_dataset

Read_PIR::~Read_PIR ()
{
    if (build_taxa !=0) {
       build_taxa=taxa_start;
      
      while (build_taxa->next!=0)
        {
          temp_taxa=build_taxa->next;
            delete build_taxa->current;
          delete build_taxa;
          build_taxa=temp_taxa;
        }
      delete build_taxa;
    }
}

//End of Read_PIR functions




//Read_Nexus functions
Sequence_dataset * Read_Nexus::get_dataset(int &num_taxa, int &num_sites, const char *filename, BOOL is_aa_seq)
{
	if (is_aa_seq == TRUE)
		return(get_dataset(num_taxa, num_sites, filename, PROTEIN));
	else
		return(get_dataset(num_taxa, num_sites, filename, NUCLEIC));
}


Sequence_dataset * Read_Nexus::get_dataset(int &num_taxa, int &num_sites, const char *filename, DATATYPE type)
{	
	
	var_stream = &file_stream;
	file_stream.open(filename);
	
	if (file_stream.fail())
    {
		cerr<<"Read_Nexus::get_dataset:Cannot find file "<<filename<<endl;
		return(0);
    }
	else
		return(get_dataset(num_taxa, num_sites, type));
		
}


Sequence_dataset * Read_Nexus::get_dataset(int &num_taxa, int &num_sites, DATATYPE type)
{
	int i,j, namepos;
    string match, line, taxa, chars, new_name;
	Sequence_dataset *the_data;

	set_datatype(type);
	if (var_stream == 0) {
		var_stream=&cin;
		std_in=TRUE;
	}

	
    std::getline((*var_stream), line);
    
    match="Begin Data";
	while (word_match(line, match)!=1)
		std::getline((*var_stream), line);
      
    match="dimensions";
	while (word_match(line, match)!=1)
        std::getline((*var_stream), line);
      
      
    i=line.find("ntax");
    taxa= line.substr(i+5);
    //cout<<"Found ntax in "<<line<<" at "<<i<<" giving "<<taxa<<endl;
    //i=loc_word_match(line.c_str(), 201, "ntax", 4)+5;
    num_taxa=string_to_int(taxa.c_str());
      
    
    ntaxa=num_taxa;
      
    i=line.find("nchar");
    chars=line.substr(i+6);
    //cout<<"Found nchars in "<<line<<" at "<<i<<" giving "<<chars<<endl;
    //i=loc_word_match(line.c_str(), 201, "nchar", 5)+6;
    num_sites=string_to_int(chars.c_str());
      
    
    size=num_sites;
    
    the_data = new Sequence_dataset(num_taxa, size, ctype);
      
    std::getline((*var_stream), line);
      
    match="Matrix";
    while (word_match(line, match)!=1)
        std::getline((*var_stream), line);
      
      
	for ( i=0; i<num_taxa; i++) {
        (*var_stream)>>new_name;

	new_name.erase(std::remove(new_name.begin(), new_name.end(), '\''), new_name.end());
	  
        get_sequence_array_size_known(&(*the_data)[i]);
        (*the_data)[i].Assign_name(new_name.c_str());
	}
		
	  if (std_in == FALSE)
		  file_stream.close();
	
      return(the_data);

} //End Read_Nexus::get_dataset

Read_Nexus::~Read_Nexus ()
{
    if (build_taxa !=0) {
       build_taxa=taxa_start;
      
      while (build_taxa->next!=0)
        {
          temp_taxa=build_taxa->next;
            delete build_taxa->current;
          delete build_taxa;
          build_taxa=temp_taxa;
        }
      delete build_taxa;
    }
}


//End of Read_Nexus functions
  



//Read_Phylip_interleave functions
Sequence_dataset * Read_Phylip_interleave::get_dataset(int &num_taxa, int &num_sites, const char *filename, BOOL is_aa_seq)
{
	if(is_aa_seq == TRUE)
		return(get_dataset(num_taxa, num_sites, filename, PROTEIN));
	else
		return(get_dataset(num_taxa, num_sites, filename, NUCLEIC));
}



Sequence_dataset * Read_Phylip_interleave::get_dataset(int &num_taxa, int &num_sites, const char *filename, DATATYPE type)
{
	var_stream=&file_stream;
	
	file_stream.open(filename);
	if (file_stream.fail())
    {
		cerr<<"Read_Phylip_interleave::get_dataset:Cannot find file "<<filename<<endl;
		return(0);
    }
	else
		return(get_dataset(num_taxa, num_sites, type));
		
}


Sequence_dataset * Read_Phylip_interleave::get_dataset(int &num_taxa, int &num_sites, DATATYPE type)
{
	int i, j, k, l, total_count=0, end, blocks, subblocks, name_offset, name_len;
	char dump, inbase, new_name[50], real_name[50];
	BOOL first_line=TRUE;
	Sequence_dataset *the_data;

	set_datatype(type);
	
	if (var_stream == 0) {
		var_stream=&cin;
		std_in=TRUE;
	}
       
	(*var_stream)>>num_taxa;
	(*var_stream)>>num_sites;
      

	ntaxa=num_taxa;
    size=num_sites;
      
    the_data=new Sequence_dataset(num_taxa, size, ctype);
      
    if (size%50==0)
		blocks=(int)size/50;
	else
		blocks=(int)size/50+1;
      
      
	dump=next_char();      
	while(dump!='\n')
		dump=next_char();
      
	for(l=0; l<blocks; l++) {
	  if (l==blocks-1)
	    if ((size-l*50)%10==0)
	      subblocks=(int)((size-l*50)/10);
	    else
	      subblocks=(int)((size-l*50)/10)+1;
	  
	  
	  else
	    subblocks=5;
	  
		
	  for (i=0; i<num_taxa; i++)   {
			if ((i==0) && (first_line == TRUE)) {
				name_len=0;
				for(j=0; j<11; j++) {
					new_name[j]=next_char();
					if (new_name[j] != ' ') name_len++;
				}
				for(j=0; j<name_len; j++) real_name[j]=new_name[j];
				real_name[name_len]='\0';
				
				(*the_data)[i].Assign_name(real_name);
				
				
				name_offset=11;
				dump=next_char();
				while (dump == ' ') {
					name_offset++;
					dump=next_char();
				}
				
				
				if (l==(blocks-1) && j == (subblocks-1)) end=(int)(size-((l*50)+10*(subblocks-1)));
				else	end=10;
				
				total_count=50*(l);
				inbase=dump;
				
				(*the_data)[i].Assign_site(total_count, convert_char(inbase));
				total_count++;
				
				for(k=1; k<end; k++) {
					inbase=next_char();
					(*the_data)[i].Assign_site(total_count, convert_char(inbase));
					total_count++;
				}
				dump=next_char();
				
				
				for(j=1; j<subblocks; j++) {
					if (l==(blocks-1) && j == (subblocks-1))
						end=(int)(size-((l*50)+10*(subblocks-1)));
					else
						end=10;
					
					for(k=0; k<end; k++) {
						inbase=next_char();
						(*the_data)[i].Assign_site(total_count, convert_char(inbase));
						total_count++;
					}
					dump=next_char();
				}
				while (dump!='\n')
					dump=next_char();
			
			}
			else {
				name_len=0;
				for(j=0; j<11; j++) {
					if(first_line==TRUE) {
						new_name[j]=next_char();
						if (new_name[j] != ' ') name_len++;
					}
					else
						dump=next_char();
				}
				
				for(j=0; j<name_len; j++) real_name[j]=new_name[j];
				real_name[name_len]='\0';
				
				for(j=11; j<name_offset; j++)	dump=next_char();

				if (first_line == TRUE)
					(*the_data)[i].Assign_name(real_name);
				
				total_count=50*(l);
			  
				for(j=0; j<subblocks; j++) {
					if (l==(blocks-1) && j == (subblocks-1))
						end=(int)(size-((l*50)+10*(subblocks-1)));
					else
						end=10;
			  
					for(k=0; k<end; k++) {
						inbase=next_char();
						(*the_data)[i].Assign_site(total_count, convert_char(inbase));
						total_count++;
						
						}
					dump=next_char();
				}
				while (dump!='\n')
				  dump=next_char();
				

			}
		}
		dump=next_char();
		while (dump != '\n') dump=next_char();
		
		first_line=FALSE;
	}
    
	if (std_in == FALSE)
		file_stream.close();
    return(the_data);
    
}

Read_Phylip_interleave::~Read_Phylip_interleave ()
{
    if (build_taxa !=0) {
       build_taxa=taxa_start;
      
      while (build_taxa->next!=0)
        {
          temp_taxa=build_taxa->next;
            delete build_taxa->current;
          delete build_taxa;
          build_taxa=temp_taxa;
        }
      delete build_taxa;
    }
}




Sequence_dataset * Read_Phylip_noninterleave::get_dataset(int &num_taxa, int &num_sites, const char *filename, BOOL is_aa_seq)
{
	if(is_aa_seq == TRUE)
		return(get_dataset(num_taxa, num_sites, filename, PROTEIN));
	else
		return(get_dataset(num_taxa, num_sites, filename, NUCLEIC));
}

Sequence_dataset * Read_Phylip_noninterleave::get_dataset(int &num_taxa, int &num_sites, const char *filename, DATATYPE type)
{
	var_stream=&file_stream;
	
	file_stream.open(filename);
	if (file_stream.fail())
    {
		cerr<<"Read_Phylip_noninterleave::get_dataset:Cannot find file "<<filename<<endl;
		return(0);
    }
	else
		return(get_dataset(num_taxa, num_sites, type));
	
}


Sequence_dataset * Read_Phylip_noninterleave::get_dataset(int &num_taxa, int &num_sites, DATATYPE type)
{
	int i;
	char dump, new_name[700];
	Sequence_dataset *the_data;
  
	set_datatype(type);
  
	if (var_stream == 0) {
		var_stream=&cin;
		std_in=TRUE;
	}

    (*var_stream)>>num_taxa>>num_sites;
      
	ntaxa=num_taxa;
      
	size=num_sites;
    
    the_data = new Sequence_dataset(num_taxa, size, ctype);
      
    dump=next_char();
      
	for (i=0; i<num_taxa; i++)
	{
	 
	  (*var_stream)>>new_name; 
	  
	  get_sequence_array_size_known(&(*the_data)[i]);
	  (*the_data)[i].Assign_name(new_name);
	}
	
	if (std_in == FALSE)
		file_stream.close();
	
	return(the_data);
}


//Read_FASTA functions
Sequence_dataset * Read_FASTA::get_dataset(int &num_taxa, int &num_sites, const char *filename, BOOL is_aa_seq)
{
	if (is_aa_seq==TRUE)
		return(get_dataset(num_taxa, num_sites, filename, PROTEIN));
	else
		return(get_dataset(num_taxa, num_sites, filename, NUCLEIC));
}


Sequence_dataset * Read_FASTA::get_dataset(int &num_taxa, int &num_sites, const char *filename, DATATYPE type)
{
	var_stream=&file_stream;
	file_stream.clear();
	file_stream.open(filename);
	
	
	if (file_stream.fail())
    {
		cerr<<"Read_FASTA::get_dataset:Cannot find file "<<filename<<endl;
		return(0);
    }
	else
		return(get_dataset(num_taxa, num_sites, type));
		
}


Read_FASTA::~Read_FASTA()
{
	if (sizes !=0) delete[] sizes;
}


Sequence_dataset * Read_FASTA::get_dataset(int &num_taxa, int &num_sites, DATATYPE type)
{
	int i;
    char inbase;
    string names;
	//char inbase, names[700];
	Sequence_dataset *the_data;
  
 
	set_datatype(type);
	
	if (var_stream == 0) {
		var_stream=&cin;
		std_in = TRUE;
	}
		

	build_taxa= taxa_start;
 
    num_taxa=0;

    inbase=next_char();
	last_char=inbase;
	//cout<<"Initial read: "<<inbase<<endl;
    while (!read_finished())
	{
		inbase=last_char;
		//cout<<"Current inbase: "<<inbase<<endl; 
		while ((inbase!='>') && !read_finished())
			inbase=next_char();
		//cout<<"SEq end base: "<<inbase<<endl;
		if (!read_finished())
		{
            getline((*var_stream), names);
			    //var_stream->getline(names, 699);
            names=names.substr(0, 699);
			//cout<<"Name line: "<<names<<endl;
			get_sequence_array(build_taxa->current);
			build_taxa->current->Assign_name(names.c_str());
			get_next_list_element(build_taxa);
			num_taxa++;
		}
	 
	}
      
    ntaxa=num_taxa;
    sizes=new int [ntaxa];
      
    build_taxa=taxa_start;
      
    for (i=0; i<ntaxa; i++)
	{
	  //  cout<<"Size of sequence "<<i<<": "<<build_taxa->current->Sequence_size()<<endl;
	  sizes[i]=build_taxa->current->Sequence_size();
	  if (i!=ntaxa-1)
	    build_taxa=build_taxa->next;
	}

      the_data= new Sequence_dataset(num_taxa, sizes, ctype);
      
      build_taxa=taxa_start;
      
      for(i=0; i<num_taxa; i++)
	{
	  (*the_data)[i]=(*build_taxa->current);
	  
	  if (i!=num_taxa-1)
	    build_taxa=build_taxa->next;
	}
      
      num_sites=(*the_data)[0].Sequence_size();
      
    if (std_in == FALSE)
		file_stream.close();
	
    return(the_data);
    
}



Sequence_dataset * Read_FASTA::get_data_subset(int size, string filename, int &num_taxa, DATATYPE type)
{
    int i;
    char inbase;
    string names;
    Sequence_dataset *the_data;
    
    
    
    if (started == FALSE) {
        //cout<<"First read: openning "<<filename<<flush<<endl;
        set_datatype(type);
    
        file_stream.clear();
        file_stream.open(filename.c_str());
    
        var_stream=&file_stream;
        
        if (file_stream.fail()) {
            cerr<<"Cannot find file "<<filename<<endl;
            return(0);
        }
        else
            started=TRUE;
    }
    
    build_taxa= taxa_start;
    
    num_taxa=0;
    
    inbase=next_char();
    last_char=inbase;
        //cout<<"Initial read: "<<inbase<<flush<<endl;
    while ((!read_finished()) && (num_taxa < size)) {
        inbase=last_char;
        //cout<<"Current inbase: "<<inbase<<flush<<endl;
        while ((inbase!='>') && (!read_finished()))
            inbase=next_char();
        //cout<<"SEq end base: "<<inbase<<flush<<endl;
        if (!read_finished())
        {
            getline(file_stream, names);
           //cout<<"Name line: "<<names<<flush<<endl;
            get_sequence_array(build_taxa->current);

            build_taxa->current->Assign_name(names.c_str());
            get_next_list_element(build_taxa);
            
            num_taxa++;
        }
        
    }
    
    if (num_taxa > 0) {
        sizes=new int [num_taxa];
        
        build_taxa=taxa_start;
        
        for (i=0; i<num_taxa; i++) {
             //cout<<"Size of sequence "<<i<<": "<<build_taxa->current->Sequence_size()<<flush<<endl;
            sizes[i]=build_taxa->current->Sequence_size();
            if (i!=num_taxa-1)
                build_taxa=build_taxa->next;
        }

        the_data= new Sequence_dataset(num_taxa, sizes, ctype);
        
        build_taxa=taxa_start;
        
        for(i=0; i<num_taxa; i++) {
            (*the_data)[i]=(*build_taxa->current);
            
            if (i!=num_taxa-1)
                build_taxa=build_taxa->next;
        }
        
        
        build_taxa=taxa_start;
        while(build_taxa !=0) {
            delete build_taxa->current;
            temp_taxa=build_taxa->next;
            delete build_taxa;
            build_taxa=temp_taxa;
        }
        taxa_start=new List_Molecule_Sequence;
        build_taxa=taxa_start;
        
        if (read_finished()) {
            file_stream.close();
            finished=TRUE;
        }
        return(the_data);
    }
    else
        return(0);
    
}


Read_Phylip_noninterleave::~Read_Phylip_noninterleave ()
{
    if (build_taxa !=0) {
       build_taxa=taxa_start;
      
      while (build_taxa->next!=0)
        {
          temp_taxa=build_taxa->next;
            delete build_taxa->current;
          delete build_taxa;
          build_taxa=temp_taxa;
        }
      delete build_taxa;
    }
}

