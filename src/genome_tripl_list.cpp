#include "genome_tripl_list.h"



WGX_Locus::WGX_Locus()
{
    dupl_level=2;
    alloc_arrays();
	species_num=0;
}



WGX_Locus::WGX_Locus(int sp_num, int d_level,  string *nms, Genome *genome)
{
    int i;
    
	species_num=sp_num;
    the_genome=genome;
    dupl_count=0;

    dupl_level=d_level;
    names=new string[dupl_level];
    has_ith=new BOOL [dupl_level];
    contig_index=new int[dupl_level];
    gene_index=new int [dupl_level];
    
    for (i=0; i<dupl_level; i++) {
        //cout<<"For genome "<<genome->get_name_string()<<": level "<<dupl_level<<" name: "<<nms[i]<<endl;
        if (nms[i] != "NONE") {
            names[i]=nms[i];
            has_ith[i]=TRUE;
            find_id(names[i], contig_index[i], gene_index[i]);
            dupl_count++;
        }
        else {
            has_ith[i]=FALSE;
            names[i]="";
            contig_index[i]=-1;
            gene_index[i]=-1;
        }
    }
	
}


Gene * WGX_Locus::get_gene_obj(int index)
{
	if (has_ith[index] == TRUE)
		return(&(*the_genome)[contig_index[index]][gene_index[index]]);
	else
		return(0);

}

BOOL WGX_Locus::has_all_duplicates()
{
    int i;
    BOOL retval=TRUE;
    
    for (i=0; i<dupl_level; i++) {
        if (has_ith[i] == FALSE) retval=FALSE;
    }
    return(retval);
}


WGX_Locus::~WGX_Locus()
{
    delete[] names;
    delete[] has_ith;
    delete[] contig_index;
    delete[] gene_index;
}

void WGX_Locus::alloc_arrays()
{
    int i;
    
    contig_index=new int [dupl_level];
    gene_index=new int [dupl_level];
    names=new string [dupl_level];
    has_ith=new BOOL [dupl_level];
    
    for(i=0; i<dupl_level; i++) {
        contig_index[i]=0;
        gene_index[i]=0;
        names[i]="NONE";
        has_ith[i]=FALSE;
    }
}


void WGX_Locus::find_id(string name, int &contig_id, int &gene_id)
{
	int i, j;
	BOOL found=FALSE;


	i=0;

	while((i<the_genome->get_num_contigs()) && (found == FALSE)) {
		j=0;
		while((j<(*the_genome)[i].get_num_genes()) && (found == FALSE)) {
			if (strcmp(name.c_str(), (*the_genome)[i][j].get_name()) == 0) {
				found =TRUE;
				contig_id=i;
				gene_id=j;
			}
			j++;
		}
		i++;

	}

	if (found == FALSE) 
		cerr<<"Error: Cannot find gene "<<name<<" in contigs\n";
}




Homologs_DX::Homologs_DX()
{
	the_loci=0;
    dupl_level=2;
	null_locus=new WGX_Locus();
	the_genomes=0;
}



Homologs_DX::Homologs_DX(int d_level, string **orthos_n, Clade *genomes)
{
	int i;

	the_genomes=genomes;
    dupl_level=d_level;
	null_locus=new WGX_Locus();

	the_loci =new WGX_Locus * [the_genomes->get_num_genomes()];
	
	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_loci[i]=new WGX_Locus(i, dupl_level, orthos_n[i], &(*the_genomes)[i]);

}


WGX_Locus& Homologs_DX::operator[] (int element)
{
	if (the_genomes ==0)
		return(*null_locus);

	if ((0<=element) && (element < the_genomes->get_num_genomes()))
		return(*the_loci[element]);
	else
		return(*null_locus);

}



Homologs_DX::~Homologs_DX()
{
	int i;

	delete null_locus;

	if (the_loci != 0)
	{
		for(i=0; i<the_genomes->get_num_genomes(); i++)
			delete the_loci[i];
		delete[] the_loci;
	}

}



WGX_Data::WGX_Data()
{
	num_homologs=0;
	null_homolog=new Homologs_DX();
    dupl_level=2;

	the_homologs=0;
}
	

WGX_Data::WGX_Data(int nhomologs, int d_level, string ***orthos, Clade *genomes)
{
	int i;

	num_homologs=nhomologs;
    dupl_level=d_level;
	null_homolog=new Homologs_DX();

	the_genomes=genomes;

	the_homologs=new Homologs_DX* [num_homologs];

	for(i=0; i<num_homologs; i++)
		the_homologs[i] = new Homologs_DX(dupl_level, orthos[i], the_genomes);

}



Homologs_DX & WGX_Data::operator[] (int element)
{
	if ((0<=element) && (num_homologs > element))
		return(*the_homologs[element]);
	else
		return(*null_homolog);
}


WGX_Data::~WGX_Data()
{
	int i;

	delete null_homolog;

	if (num_homologs != 0) {
		for(i=0; i<num_homologs; i++)
			delete the_homologs[i];
		delete[] the_homologs;
	}
}



Read_WGX_Data::Read_WGX_Data()
{
	species_indexes=0;
	genome_names=0;
    orthos=0;
}
	


WGX_Data * Read_WGX_Data::get_data(string filename, int d_level, int &num_sites, Clade *genomes)
{
	int i, j, k;
	string line;
	WGX_Data *return_data;
    
    dupl_level=d_level;

	the_genomes=genomes;
	num_genomes=the_genomes->get_num_genomes();

	num_sites=0;

	infile.open(filename.c_str());

	if(!infile.fail()) {
		genome_names=new string [the_genomes->get_num_genomes()];
		species_indexes = new int [the_genomes->get_num_genomes()];
		

		for(i=0; i<the_genomes->get_num_genomes(); i++)
			infile>>genome_names[i];

		set_indexes();
		while(!infile.eof()) {
			getline(infile, line);
			num_sites++;
		}

		//Should have tried one failed read
		num_sites -=2;
		sites=num_sites;

		//Reset the file to actually read the data--lazy
		infile.close();
		infile.clear();
		infile.open(filename.c_str());

        cout<<"Reading "<<num_sites<<" orthologous sites\n";
		orthos=new string ** [num_sites];
		
        for(i=0; i<num_sites; i++) {
			orthos[i]=new string * [the_genomes->get_num_genomes()];
            for(j=0; j<the_genomes->get_num_genomes(); j++) {
                orthos[i][j]=new string [dupl_level];
            }
            
        }

		//Clear the Genome Name line
		getline(infile, line);

		for(i=0; i<num_sites; i++) {
			//Read the pillars
            for(k=0; k<dupl_level; k++) {
                for(j=0; j<the_genomes->get_num_genomes(); j++)
                    infile>>orthos[i][species_indexes[j]][k];
            }

		}
		infile.close();
		return_data = new WGX_Data(num_sites, dupl_level, orthos, the_genomes);
		return(return_data);
	}
    else return(0);


}



Read_WGX_Data::~Read_WGX_Data()
{
	int i, j;


	if (genome_names !=0)
		delete[] genome_names;


	if (species_indexes != 0)
		delete[] species_indexes;


	if (orthos != 0) {
		for(i=0; i<sites; i++) {
			for(j=0; j<num_genomes; j++) {
				delete[] orthos[i][j];
			}
			delete[] orthos[i];
		}
		delete[] orthos;
	}
}



void Read_WGX_Data::set_indexes()
{
	int i, j;
	

	for(i=0; i<the_genomes->get_num_genomes(); i++) {
		j=0;
		while((j<the_genomes->get_num_genomes()) && 
			((*the_genomes)[j].get_name_string ()!=genome_names[i]))	j++;

		

		if (j == the_genomes->get_num_genomes())
			cerr<<"Error: could not find genome "<<genome_names[i]<<" in data\n";
		else
			species_indexes[i]=j;
	
	}

}


Gene_Track_List_DX::Gene_Track_List_DX()
{
	my_locus=0;
	index_num=0;
	next=last=0;
	partition_num=0;
    in_track_line=FALSE;
}


Gene_Track_List_DX& Gene_Track_List_DX::operator= (Gene_Track_List_DX &assign_from)
{
	my_locus=assign_from.my_locus;
	next=assign_from.next;
	last=assign_from.last;
	index_num=assign_from.index_num;
	return(*this);
}




Tracking_List_DX::Tracking_List_DX()
{
	cerr<<"Error: Call to default constructor of class Tracking_List\n";
	last=next=0;
    dupl_level=2;
	element = new Gene_Track_List_DX*[2];
	element[0]=element[1]=0;
}
	

Tracking_List_DX::Tracking_List_DX(int d_level, class Gene_Track_List_DX **ele, Tracking_List_DX *lst, int n)
{
    int i;
    
	Tracking_List_DX *old_next=0;
    dupl_level=d_level;

	element = new Gene_Track_List_DX*[dupl_level];

    for (i=0; i<dupl_level; i++)
        element[i]=ele[i];

	if (lst !=0)
		old_next=lst->next;

	last=lst;

	if (lst != 0)
		lst->next=this;

	next=old_next;

	if(old_next != 0)
		old_next->last=this;

	num=n;

}


Track_List_List_DX::Track_List_List_DX()
{
	last=next=0;
	element=0;
}

Track_List_List_DX::Track_List_List_DX(Gene_Track_List_DX *new_ele, Track_List_List_DX *old_top)
{
	next=old_top;
	last=0;
	element=new_ele;
}

Tracking_List_DX::~Tracking_List_DX()
{
	delete[] element;
}


Track_Stack_DX::Track_Stack_DX()
{
	top=bottom=0;
	size=0;
}



Gene_Track_List_DX * Track_Stack_DX::get_bottom()
{
	return(bottom->element);
}



Gene_Track_List_DX * Track_Stack_DX::get_next()
{
	Gene_Track_List_DX *retval=0;

	if ((size > 0) && (pos_num < size)) {
		retval=curr_pos->element;
		curr_pos=curr_pos->last;
		pos_num++;
	}
	return(retval);
}



void Track_Stack_DX::delete_element(Gene_Track_List_DX *element)
{
	Track_List_List_DX *temp, *temp2;
	temp=top;

	while((temp !=0) && (temp->element != element)) temp=temp->next;

	if (temp!=0) {
		if (temp == curr_pos) curr_pos=curr_pos->last;
		if (temp->last != 0) {
			temp2=temp->last;
			if (temp->next !=0) {
				//Element is neither top nor bottom
				temp->next->last=temp2;
				temp2->next=temp->next;
			}
			else
			{
				//Element is bottom
				bottom=temp2;
				bottom->next=0;
			}
		}
		else
		{
			//Element is top
			top=top->next;
			if (top != 0)
				top->last=0;
			else
				bottom=0;
			
		}
		
        delete temp;
		size--;
		
	}

}




void Track_Stack_DX::pop()
{
	Track_List_List_DX *temp;

	temp=bottom->last;


    delete bottom;
	bottom=temp;
	if (size > 1) 
		bottom->next=0;
	else {
		top=0;
		bottom=0;
	}
	size--;
}



void Track_Stack_DX::push(Gene_Track_List_DX *entry)
{
	Track_List_List_DX *new_item=0;

    new_item=new Track_List_List_DX(entry, top);
	
	if (top != 0) {
		top->last=new_item;
		new_item->next=top;
		top=new_item;
	}
	else
	{
		//This is the first element
		bottom=new_item;
		top=bottom;
		pos_num=0;
	}
	
	size++;
}





Track_Stack_DX::~Track_Stack_DX()
{
	int i;
	Track_List_List_DX *temp;

	for(i=0; i<size; i++)
	{
		temp=bottom->last;
		delete bottom;
		bottom=temp;
	}

}


Genome_Track_DX::Genome_Track_DX()
{
	cerr<<"Error: Call to default constructor of class Genome_Track_DX\n";
	hold_reversal_list=0;
    new_inferred_tracking=0;
    inferred_tracking=0;
    to_ends=0;
    broken_lefts=0;
    broken_rights=0;

    tracks=0;
}


Genome_Track_DX::Genome_Track_DX(int id, Genome *genome, WGX_Data *homologs)
{
	int i, j;
	Gene_Track_List_DX **ele;

	taxa_id=id;
    
    
	the_genome=genome;
	the_homologs=homologs;
    
    ele=new Gene_Track_List_DX*[the_homologs->get_dupl_level()];
    tracks=new int [the_homologs->get_dupl_level()];

    to_ends=new Gene_Track_List_DX * [the_homologs->get_dupl_level()];

    new_inferred_tracking  =0;
    
	inferred_tracking = new Gene_Track_List_DX* [homologs->get_num_homologs()];
	
	for(i=0; i<homologs->get_num_homologs(); i++)
		inferred_tracking[i]=new Gene_Track_List_DX [the_homologs->get_dupl_level()];
	
    for(i=0; i<homologs->get_num_homologs(); i++) {
        for(j=0; j<homologs->get_dupl_level(); j++) {
            inferred_tracking[i][j].my_locus=0;
            inferred_tracking[i][j].index_num=-1;
            inferred_tracking[i][j].to_track_num=-1;
            inferred_tracking[i][j].track_pos=i;
            inferred_tracking[i][j].next=0;
            inferred_tracking[i][j].last=0;
        }
    }
    
    broken_lefts=new Gene_Track_List_DX* [the_homologs->get_dupl_level()*the_homologs->get_dupl_level()];
    broken_rights=new Gene_Track_List_DX* [the_homologs->get_dupl_level()*the_homologs->get_dupl_level()];
    
    

	pow_N[0]=1;
	for(i=1; i<32; i++)
		pow_N[i] =pow_N[i-1]*((long long int)the_homologs->get_dupl_level());

	partial_track_start=partial_track_pos=partial_track_end=0;
	list_len=last_list_pos=0;

	hold_reversal_list =new Tracking_List_DX* [the_homologs->get_num_homologs()];
	reversal_pos=0;

	for(i=0; i<the_homologs->get_num_homologs(); i++) {
        for(j=0; j<the_homologs->get_dupl_level(); j++)
            ele[j]=new Gene_Track_List_DX;
		hold_reversal_list[i]=new Tracking_List_DX(the_homologs->get_dupl_level(), ele, 0, i);
	}
}




void Genome_Track_DX::number_contigs(int *order)
{
	number_contigs(order, the_homologs->get_num_homologs());
}



void Genome_Track_DX::number_contigs(int *order, int size)
{
	int i, j;
    Gene_Track_List_DX *dummy;
	Track_Stack_DX *lefts=0, *rights=0;
	
	for(j=0; j<size; j++) {
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            
            if ((*the_homologs)[order[j]][taxa_id].has_duplicate(i) == TRUE) {
                inferred_tracking[j][i].my_locus=&(*the_homologs)[order[j]][taxa_id];
                inferred_tracking[j][i].index_num=i;
            }
            else 	inferred_tracking[j][i].my_locus=0;
            
            
        }
        
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            inferred_tracking[j][i].last=inferred_tracking[j][i].next=0;
            inferred_tracking[j][i].dist_to_last=0;
        }
			
    }
 
#if 0
    for(j=0; j<size; j++) {
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            if (inferred_tracking[j][i].my_locus!=0) cout<<inferred_tracking[j][i].my_locus->get_gene_obj(i)->get_name()<<" (index "<<inferred_tracking[j][i].index_num<<")\t";
            else cout<<"NONE\t";
        }
        cout<<endl;
    }
#endif
    for(j=0; j<size; j++) {
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            inferred_tracking[j][i].to_track_num=-1;
            inferred_tracking[j][i].in_track_line=FALSE;
        }
    }
	//recurse_assemble(lefts, rights, 0, size);
    recurse_assemble_V2(0, size);
	delete lefts;
	delete rights;


	//ID each locus
	for(j=0; j<size; j++) {
        for(i=0; i<the_homologs->get_dupl_level(); i++)
            inferred_tracking[j][i].partition_num=j;
	}

	//Recheck all breaks to see if we can gain any breaks
	//for(j=1; j<size; j++)
	//	break_and_rejoin(j, size);
	

	//Puts the above list into an easy-to-understand 2Xn table
	store_tracking();
}



void Genome_Track_DX::print_tracking()
{
	print_tracking(TRUE);
}


void Genome_Track_DX::print_tracking(BOOL use_file)
{
	int i, j;
    char backlink;
	string filename;
	ofstream fout;
	BOOL link;

	if (use_file == TRUE) {
		filename =the_genome->get_name();
		filename = filename + "_tracking.txt";
		fout.open(filename.c_str());
	}

	for(j=0; j<the_homologs->get_dupl_level(); j++) {
		for(i=0; i<the_homologs->get_num_homologs(); i++) {
            if (inferred_tracking[i][j].last ==0) backlink='X';
            else backlink='<';
			if (inferred_tracking[i][j].my_locus != 0) {
				if (use_file == TRUE) 
					fout<<backlink<<" "<<(*the_genome)[inferred_tracking[i][j].my_locus->get_contig(
							inferred_tracking[i][j].index_num)][inferred_tracking[i][j].my_locus->get_gene(
							inferred_tracking[i][j].index_num)].get_name();
				else
					cout<<backlink<<" "<<(*the_genome)[inferred_tracking[i][j].my_locus->get_contig(
							inferred_tracking[i][j].index_num)][inferred_tracking[i][j].my_locus->get_gene(
							inferred_tracking[i][j].index_num)].get_name();

				if (inferred_tracking[i][j].next == 0)
					link=FALSE;
				else
					link=TRUE;
				
				if (link == TRUE) {
					if (use_file == TRUE)
						fout<<"\t->\t";
					else
						cout<<"\t->\t";
				}
				else {
					if (use_file == TRUE)
						fout<<"\t|\t";
					else
						cout<<"\t|\t";
				}
			}
			else {
				if (use_file == TRUE)
					fout<<backlink<<" "<<"NONE\t*\t";
				else
					cout<<backlink<<" "<<"NONE\t*\t";
			}

			
		
		}

		if (use_file == TRUE)
			fout<<endl<<endl;
		else
			cout<<endl<<endl;
	}

	if (use_file == TRUE)
		fout.close();	
}




Gene_Track_List_DX* Genome_Track_DX::get_gene_track(int locus_num, int track_num)
{
	if ((locus_num >=0) && (locus_num<the_homologs->get_num_homologs())) {
        return(&inferred_tracking[locus_num][track_num]);
	}
	else {
		cerr<<"Request for invalid index "<<locus_num<<" in Genome_Track_DX\n";
 		return(0);
	}
}
	


BOOL Genome_Track_DX::has_back_link(int locus_num, int track_num)
{
    if ((locus_num >=1) && (locus_num<the_homologs->get_num_homologs())) {
        if (inferred_tracking[locus_num][track_num].last != 0)
            return(TRUE);
        else
            return(FALSE);
    }
	else return(FALSE);

}

BOOL Genome_Track_DX::has_forward_link(int locus_num, int track_num)
{
    if ((locus_num >=0) && (locus_num<the_homologs->get_num_homologs())) {
        if (inferred_tracking[locus_num][track_num].next != 0)
            return(TRUE);
        else
            return(FALSE);
    }
    else return(FALSE);
    
}



//This could still be wrong

BOOL Genome_Track_DX::has_tracked_full_break(int n)
{
	int i, j;
	BOOL stop, retval=TRUE;
	Gene_Track_List_DX *from_previous;

	i=n-1;
	stop=FALSE;

		//Find all the possible left endpoints

	while((stop==FALSE) && (retval == TRUE)) {
		to_ends[0]=inferred_tracking[i][0].next;
			
		if (i != n-1)
			from_previous=inferred_tracking[i+1][0].last;
		else
			from_previous=0;

        for(j=1; j<the_homologs->get_dupl_level(); j++) {
            if (inferred_tracking[i][j].my_locus != 0)
                to_ends[j]=inferred_tracking[i][j].next;
            else
                to_ends[j] = 0;
        }
	

		if (to_ends[0] == 0) {
			if (i != n-1) {
				//If the previous entry has a join but not to this entry, then
				//this one must be on the other track.  Hence, stop
				if (from_previous !=0)
					stop=TRUE;
			}
		}
		else if (to_ends[0]->partition_num > inferred_tracking[n-1][0].partition_num)
			retval=FALSE;
        for(j=1; j<the_homologs->get_dupl_level(); j++) {
            if ((inferred_tracking[i][j].my_locus != 0) && (to_ends[j] ==0)) {
                //This is a left end and a stopping point
                stop=TRUE;
            }
            else if ((to_ends[j] != 0) && (to_ends[j]->partition_num > inferred_tracking[n-1][0].partition_num))
                    retval=FALSE;
        }
		

		if (i == 0)
			stop=TRUE;

		i--;
	}
	
	return(retval);
}



int Genome_Track_DX::get_dist_to_last(int locus_num, int track_num)
{
	if (has_back_link(locus_num, track_num) == FALSE)
		//Error condition
		return(-1);
	else
		return(inferred_tracking[locus_num][track_num].dist_to_last);
}



int Genome_Track_DX::count_num_breaks()
{
	return(count_num_breaks(the_homologs->get_num_homologs()));
}



int Genome_Track_DX::count_num_breaks(int size)
{
	int i, j, retval=0;

	for(i=0; i<size; i++) {
        for(j=0; j<the_homologs->get_dupl_level(); j++) {
			if ((inferred_tracking[i][j].my_locus !=0) && (inferred_tracking[i][j].last == 0))
				retval++;
        }
	}
	return(retval);
}


int Genome_Track_DX::count_num_full_breaks()
{
    int i, j, breaks, retval=0;
    
    for(i=0; i<the_homologs->get_num_homologs(); i++) {
        breaks=0;
        for(j=0; j<the_homologs->get_dupl_level(); j++) {
            //if ((inferred_tracking[i][j].my_locus !=0) && (inferred_tracking[i][j].last == 0))
            if (inferred_tracking[i][j].last == 0)                breaks++;
        }
        if (breaks == the_homologs->get_dupl_level()) retval++;
    }
    return(retval);
}

int Genome_Track_DX::count_num_list_track_breaks()
{
	int i, retval=0;
	Tracking_List_DX *curr_pos;

	curr_pos=partial_track_start;

	while(curr_pos != 0) {
        for (i=0; i<the_homologs->get_dupl_level(); i++) {
            if ((curr_pos->element[i]->my_locus != 0) &&(curr_pos->element[i]->next == 0))
                retval++;
        }
		curr_pos=curr_pos->next;
	}

	return(retval);
	
}

#if 0
void Genome_Track_DX::start_track_list(int n)
{
	Gene_Track_List_DX **locus;
	Tracking_List_DX *new_member;
	
	locus=create_new_track_locus(&(*the_homologs)[n][taxa_id]);

	new_member=new Tracking_List_DX(the_homologs->get_dupl_level(), locus, 0, 0);

	partial_track_start=partial_track_pos=partial_track_end=new_member;
	last_list_pos=0;
	list_len=1;
}


Tracking_List_DX * Genome_Track_DX::create_list_element(int locus_num)
{
	Gene_Track_List_DX **locus;
	Tracking_List_DX *new_member;

	locus=create_new_track_locus(&(*the_homologs)[locus_num][taxa_id]);
	
	new_member=new Tracking_List_DX(the_homologs->get_dupl_level(), locus, 0, locus_num);

	return(new_member);
}
#endif

int Genome_Track_DX::joins_after(int n)
{
	int i, num_joins=0;

    for(i=0; i<the_homologs->get_dupl_level(); i++)
        tracks[i]=n;

	for(i=0; i<the_homologs->get_dupl_level(); i++) {
		while((tracks[i] >= 0)  && (inferred_tracking[tracks[i]][i].my_locus == 0) ) {
			tracks[i]--;
		}
		if ((tracks[i]>=0) && (inferred_tracking[tracks[i]][i].next != 0))
			num_joins++;
	}

	return(num_joins);
	
}



BOOL Genome_Track_DX::all_poss_joins_after(int n)
{
    int i;
    BOOL retval=TRUE;
    for(i=0; i<the_homologs->get_dupl_level(); i++) {
        if ((inferred_tracking[n][i].my_locus != 0) && (inferred_tracking[n][i].next == 0))
            retval=FALSE;
    }
    return(retval);
}


BOOL Genome_Track_DX::diagnose_list()
{
	int i, j, k;
	BOOL retval=TRUE;
	Tracking_List_DX *curr_pos;

	curr_pos = partial_track_start;

	while(curr_pos != 0) {
		for(i=0; i<the_homologs->get_dupl_level(); i++) {
			if ((curr_pos->element[i]->my_locus != 0) && (curr_pos->element[i]->last != 0)) {
                for(j=0; j<the_homologs->get_num_homologs(); j++) {
                    for (k=0; k<the_homologs->get_dupl_level(); k++)
                        if (curr_pos->element[i]->last == hold_reversal_list[j]->element[k])
						retval = FALSE;
                }
				if ((curr_pos->element[i]->last->my_locus->get_gene_obj(curr_pos->element[i]->last->index_num)->get_neighbor(0) !=
					curr_pos->element[i]->my_locus->get_gene_obj(curr_pos->element[i]->index_num) ) &&
					(curr_pos->element[i]->last->my_locus->get_gene_obj(curr_pos->element[i]->last->index_num)->get_neighbor(1) != 
					curr_pos->element[i]->my_locus->get_gene_obj(curr_pos->element[i]->index_num)))
					retval=FALSE;
			}

			if ((curr_pos->element[i]->my_locus != 0) && (curr_pos->element[i]->next != 0)) {
                for(j=0; j<the_homologs->get_num_homologs(); j++) {
                     for (k=0; k<the_homologs->get_dupl_level(); k++)
                         if (curr_pos->element[i]->next == hold_reversal_list[j]->element[k] )
                             retval = FALSE;
                }
				
				if ((curr_pos->element[i]->next->my_locus->get_gene_obj(curr_pos->element[i]->next->index_num)->get_neighbor(0) != 
					curr_pos->element[i]->my_locus->get_gene_obj(curr_pos->element[i]->index_num) ) &&
					(curr_pos->element[i]->next->my_locus->get_gene_obj(curr_pos->element[i]->next->index_num)->get_neighbor(1) != 
					curr_pos->element[i]->my_locus->get_gene_obj(curr_pos->element[i]->index_num)))
					retval=FALSE;
			}

				
			
		}

		for(i=0; i<the_homologs->get_num_homologs(); i++)
		{
			if (hold_reversal_list[i] == curr_pos)
				retval=TRUE;
		}
		curr_pos=curr_pos->next;
	}

	return(retval);

}


int Genome_Track_DX::get_nth_tracking_list_pos(int n)
{
	int i;

	if (n<last_list_pos)  {
		partial_track_pos=partial_track_start;
		last_list_pos=0;
	}



	while (last_list_pos<n) {
		partial_track_pos = partial_track_pos->next;
		last_list_pos++;
	}
	
	return(partial_track_pos->num);

}




int Genome_Track_DX::get_list_position_number(int n)
{
	return(get_list_pos(n)->num);
}






Genome_Track_DX::~Genome_Track_DX()
{
	int i, j;

    
    cout<<"in track DX descructor: "<<new_inferred_tracking<<"\n";
    if (inferred_tracking!= 0) {
        //for(i=0; i<the_homologs->get_num_homologs(); i++)
        //    delete[] inferred_tracking[i];
        delete[] inferred_tracking;
        cout<<"Deleted inferred tracking\n"<<flush;
    }
    
    if (tracks !=0) delete[] tracks;
    cout<<"Deleted tracks\n"<<flush;
	
    

    delete[] to_ends;
    cout<<"Deleted to ends\n"<<flush;

	if (hold_reversal_list != 0)
	{
		for(i=0; i<the_homologs->get_num_homologs(); i++)
		{
            for(j=0; j<the_homologs->get_dupl_level(); j++)
                delete hold_reversal_list[i]->element[j];
			delete hold_reversal_list[i];
		}

		delete[] hold_reversal_list;
	}

    if (broken_lefts !=0) delete[] broken_lefts;
    if (broken_rights !=0) delete[] broken_rights;
    
	if (partial_track_start != 0) {
		partial_track_pos=partial_track_start;

		while(partial_track_pos != 0) {
			partial_track_end=partial_track_pos->next;

			for(j=0; j<the_homologs->get_dupl_level(); j++)
                delete partial_track_pos->element[j];
			delete[] partial_track_pos->element;

			delete partial_track_pos;
			partial_track_pos=partial_track_end;

		}

	}
}


void Genome_Track_DX::mark_used(int left_pos, int right_pos)
{
    int mark_pos, empty_index;
    
    for(mark_pos=left_pos+1; mark_pos<right_pos; mark_pos++) {
        empty_index=0;
        while(((inferred_tracking[mark_pos][empty_index].my_locus !=0) ||
               (inferred_tracking[mark_pos][empty_index].in_track_line == TRUE)) &&
              (empty_index<the_homologs->get_dupl_level())) {
            empty_index++;
        }
        if (empty_index==the_homologs->get_dupl_level())
            cerr<<"ERROR: Made a join between "<<left_pos<<" and "<<right_pos<<" without an avaliable path\n";
        inferred_tracking[mark_pos][empty_index].in_track_line=TRUE;
    }
    
}




int Genome_Track_DX::num_pass_throughs (int pos)
{
    int index, ret_val;
    
    ret_val = the_homologs->get_dupl_level();
    
    for(index=0; index<the_homologs->get_dupl_level(); index++) {
        if ((inferred_tracking[pos][index].my_locus !=0) ||
            (inferred_tracking[pos][index].in_track_line ==TRUE))
            ret_val--;
    }
    return(ret_val);
}



void Genome_Track_DX::recurse_assemble(Track_Stack_DX *&lefts, Track_Stack_DX *&rights, int left_end, int size)
{
    int i, j,  size_left=0, size_right=0, depth, num_joins=0, join, left_pos, right_pos, curr_used, first_non_empty, allowed_left_joins, allowed_right_joins, num_deep_left_joins, num_deep_right_joins, num_paths;
    Gene_Track_List_DX *curr_left, *curr_right, *valid_left, *valid_right, *stop_track;
    Track_Stack_DX *left_lefts, *left_rights, *right_rights, *right_lefts;
    BOOL stop;
    
    if (size > 2) {
        size_left=size_right=depth=0;
        while (pow_N[depth] < size) depth++;
        size_left=pow_N[depth-1];
        size_right=size-size_left;
        
        
        //cout<<"Joining "<<inferred_tracking[left_end+size_left-1][0].my_locus->get_gene_obj(0)->get_name()<<" to "<<inferred_tracking[left_end+size_left][0].my_locus->get_gene_obj(0)->get_name()<<"Recurse...."<<endl;
        recurse_assemble(left_lefts, left_rights, left_end, size_left);
        recurse_assemble(right_lefts, right_rights, left_end+size_left, size_right);
        
        
        left_pos=left_end+size_left-1;
        right_pos=left_end+size_left;
//cout<<"Starting assembly of size "<<size<<" from "<<left_pos<<" ("<<left_end<<") to "<<right_pos<<" ("<<right_pos+size_right<<")"<<endl;
        
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            if (inferred_tracking[left_pos][i].my_locus !=0) {
                for(j=0; j<the_homologs->get_dupl_level(); j++) {
                    if (inferred_tracking[right_pos][j].my_locus!=0) {
                        join=assign_adjacency(&inferred_tracking[left_pos][i], &inferred_tracking[right_pos][j]);
                        
                      //  if (join == 1) cout<<"Local Joined: "<<inferred_tracking[left_pos][i].my_locus->get_gene_obj(inferred_tracking[left_pos][i].index_num)->get_name()<<" at "<<left_pos<<" to R: "<<inferred_tracking[right_pos][j].my_locus->get_gene_obj(inferred_tracking[right_pos][j].index_num)->get_name()<<" at "<<right_pos<<endl;
                    }
                }
            }
        }
        
        allowed_left_joins=num_pass_throughs(left_pos);
        
        allowed_right_joins=num_pass_throughs(right_pos);
        
        //cout<<"Can make "<<allowed_left_joins<<" going left and "<<allowed_right_joins<<" going right\n";
        
        num_deep_left_joins=0;
        num_deep_right_joins=0;
        
        //Join
        while(left_rights->get_size()>0) {
            curr_left=left_rights->get_bottom();
            join=0;
            if (curr_left->next == 0) {   //It's possible this end was joined in the loop over the two ends--if so, ignore
                right_lefts->reset();
                
                curr_right=right_lefts->get_next();
                while((curr_right != 0) && (join == 0)) { //It's possible this end was joined in the loop over the two ends--if so, ignore
                     //cout<<"End: "<<left_end<<": Test joining "<<curr_left->track_pos<<" to "
                       // <<curr_right->track_pos<<endl;
                    if (curr_right->last == 0) {
                        join=assign_adjacency(curr_left, curr_right);
                        if (join == 1) {
                            //We made a join
                            if (curr_right->track_pos != left_end+size_left) num_deep_right_joins++;
                            if (curr_left->track_pos != left_end+size_left-1) num_deep_left_joins++;
                            //right_lefts->delete_element(curr_right);
                            //left_rights->pop();
                            mark_used(curr_left->track_pos, curr_right->track_pos);
                            num_joins++;
                            
                           // cout<<"Joined L "<<size<<" "<<curr_left->my_locus->get_gene_obj(curr_left->index_num)->get_name()<<" at "<<curr_left->track_pos<<" to R: "<<curr_right->my_locus->get_gene_obj(curr_right->index_num)->get_name()<<" at "<<curr_right->track_pos<<endl;
                            
                           
                            //if (num_deep_left_joins == allowed_left_joins) {
                                //We can't allow further pass-through joins. At best, use any open genes at the right end of the left
                                // cout<<"Resetting stacks\n";
                                delete left_rights;
                                left_rights=new Track_Stack_DX();
                                
                                
                                left_pos=left_end+size_left-1;
                                stop=FALSE;
                            while ((stop ==FALSE) && (left_pos>=left_end)) {
                                for (j=0; j<the_homologs->get_dupl_level(); j++) {
                                    if ((inferred_tracking[left_pos][j].my_locus != 0) && (inferred_tracking[left_pos][j].next == 0)) {
                                            left_rights->push(&inferred_tracking[left_pos][j]);
                                            //cout<<"Reset: Pushed L"<<inferred_tracking[left_pos][j].my_locus->get_gene_obj(j)->get_name()<<endl;
                                    }
                                }
                                num_paths=num_pass_throughs (left_pos);
                                
                                
                                if (num_paths < 1) stop=TRUE;
                                left_pos--;
                            }
                                
                          // }
                           // if (num_deep_right_joins == allowed_right_joins) {
                                //We can't allow further pass-through joins. At best, use any open genes at the left end of the right
                                delete right_lefts;
                                right_lefts=new Track_Stack_DX();
                                    
                                right_pos=left_end+size_left;
                                stop=FALSE;
                            
                            while((stop ==FALSE) && (right_pos<(left_end+size))) {
                                for (j=0; j<the_homologs->get_dupl_level(); j++) {
                                    if ((inferred_tracking[right_pos][j].my_locus != 0) && (inferred_tracking[right_pos][j].last == 0)) {
                                        right_lefts->push(&inferred_tracking[right_pos][j]);
                                       // cout<<"Reset: Pushed R"<<inferred_tracking[right_pos][j].my_locus->get_gene_obj(j)->get_name()<<endl;
                                    }
                                }
                                right_pos++;
                                
                                num_paths=num_pass_throughs (right_pos);
                                
                                
                                if (num_paths < 1) stop=TRUE;
                            }
                            //}
                        }
                        else
                            curr_right=right_lefts->get_next();
                        
                    }
                    else curr_right=right_lefts->get_next(); //We ignored this end because it had been joined--keep going.
                }
            }
            
            if (join == 0) {
                //We've been through all the rights for this left and can't join it--get rid of it
                left_rights->pop();
            }
        }
        
        delete left_rights;
        delete right_lefts;
        delete left_lefts;
        delete right_rights;
        
        //Now find the ends to pass off to the next level
        lefts=new Track_Stack_DX();
        rights=new Track_Stack_DX();
        stop=FALSE;
        
        //Find all the possible left endpoints
        left_pos=left_end;
        
        
        while((stop==FALSE) && (left_pos<left_end+size)) {
            for (j=0; j<the_homologs->get_dupl_level(); j++) {
                if ((inferred_tracking[left_pos][j].my_locus != 0) && (inferred_tracking[left_pos][j].last == 0)) {
                    lefts->push(&inferred_tracking[left_pos][j]);
                    //cout<<"From Recurse: Pushed L"<<inferred_tracking[left_pos][j].my_locus->get_gene_obj(j)->get_name()<<endl;
                }
            }
            num_paths=num_pass_throughs (left_pos);
            
            
            if (num_paths < 1) stop=TRUE;
            left_pos++;
        }
        
        stop=FALSE;
        right_pos=left_end+size-1;
        
        while((stop==FALSE) && (right_pos>=left_end)) {
            for (j=0; j<the_homologs->get_dupl_level(); j++) {
                if ((inferred_tracking[right_pos][j].my_locus != 0) && (inferred_tracking[right_pos][j].next == 0)) {
                    rights->push(&inferred_tracking[right_pos][j]);
                    //cout<<"From Recurse: Pushed R"<<inferred_tracking[right_pos][j].my_locus->get_gene_obj(j)->get_name()<<endl;
                }
            }
            num_paths=num_pass_throughs (right_pos);
            
            if (num_paths < 1) stop=TRUE;
            right_pos--;
        }
        
    }
    else {
        //Base recursion
        lefts=new Track_Stack_DX();
        rights=new Track_Stack_DX();
        if (size == 1) {
            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                if (inferred_tracking[left_end][i].my_locus !=0)
                    lefts->push(&inferred_tracking[left_end][i]);
            }
            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                if (inferred_tracking[left_end][i].my_locus !=0)
                    rights->push(&inferred_tracking[left_end][i]);
            }
        }
        
        else {
            //cout<<"Base Joining "<<inferred_tracking[left_end][0].my_locus->get_gene_obj(0)->get_name()<<" to "<<inferred_tracking[left_end+1][0].my_locus->get_gene_obj(0)->get_name()<<endl;
            //Check for all possible 2 gene joins
            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                for(j=0; j<the_homologs->get_dupl_level(); j++)
                    assign_adjacency(&inferred_tracking[left_end][i], &inferred_tracking[left_end+1][j]);
            }
            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                if (inferred_tracking[left_end][i].my_locus != 0)
                    lefts->push(&inferred_tracking[left_end][i]);
            }
            
            first_non_empty=0;
            while (inferred_tracking[left_end][first_non_empty].my_locus ==0) first_non_empty++;
            
            if (inferred_tracking[left_end][first_non_empty].my_locus->has_all_duplicates() == FALSE) {
                for (i=0; i<the_homologs->get_dupl_level(); i++) {
                    if ((inferred_tracking[left_end+1][i].my_locus != 0) && (inferred_tracking[left_end+1][i].last ==0))
                        lefts->push(&inferred_tracking[left_end+1][i]);
                }
            }  //Done with Lefts
            
            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                if (inferred_tracking[left_end+1][i].my_locus != 0)
                    rights->push(&inferred_tracking[left_end+1][i]);
            }
            first_non_empty=0;
            while (inferred_tracking[left_end+1][first_non_empty].my_locus ==0) first_non_empty++;
            
            if (inferred_tracking[left_end+1][first_non_empty].my_locus->has_all_duplicates() == FALSE) {
                for (i=0; i<the_homologs->get_dupl_level(); i++) {
                    if ((inferred_tracking[left_end][i].my_locus != 0) && (inferred_tracking[left_end][i].next ==0))
                        rights->push(&inferred_tracking[left_end][i]);
                }
            }   //Done with Rights
            
        }	//Done with 2-gene join
    }					//Done with <3 gene joins (1 or 2)
    
    
}

void Genome_Track_DX::recurse_assemble_V2(int left_end, int size)
{
    int i, j,  size_left=0, size_right=0, depth, join, left_pos, right_pos, curr_used, first_non_empty, num_paths;
    Gene_Track_List_DX *curr_left, *curr_right, *valid_left, *valid_right, *stop_track;
    Track_Stack_DX *right_lefts, *left_rights;
    BOOL stop;
    
    if (size > 2) {
        size_left=size_right=depth=0;
        while (pow_N[depth] < size) depth++;
        size_left=pow_N[depth-1];
        size_right=size-size_left;
        
        
        //cout<<"Joining "<<inferred_tracking[left_end+size_left-1][0].my_locus->get_gene_obj(0)->get_name()<<" to "<<inferred_tracking[left_end+size_left][0].my_locus->get_gene_obj(0)->get_name()<<"Recurse...."<<endl;
        recurse_assemble_V2(left_end, size_left);
        recurse_assemble_V2(left_end+size_left, size_right);
        
        
        left_pos=left_end+size_left-1;
        right_pos=left_end+size_left;
        //cout<<"Starting assembly of size "<<size<<" from "<<left_pos<<" ("<<left_end<<") to "<<right_pos<<" ("<<right_pos+size_right<<")"<<endl;
        
        for(i=0; i<the_homologs->get_dupl_level(); i++) {
            if (inferred_tracking[left_pos][i].my_locus !=0) {
                for(j=0; j<the_homologs->get_dupl_level(); j++) {
                    if (inferred_tracking[right_pos][j].my_locus!=0) {
                        join=assign_adjacency(&inferred_tracking[left_pos][i], &inferred_tracking[right_pos][j]);
                        
                        //  if (join == 1) cout<<"Local Joined: "<<inferred_tracking[left_pos][i].my_locus->get_gene_obj(inferred_tracking[left_pos][i].index_num)->get_name()<<" at "<<left_pos<<" to R: "<<inferred_tracking[right_pos][j].my_locus->get_gene_obj(inferred_tracking[right_pos][j].index_num)->get_name()<<" at "<<right_pos<<endl;
                    }
                }
            }
        }
        
        set_stacks(left_rights, left_pos, left_end-1, TRUE);
        set_stacks(right_lefts, left_end+size_left, left_end+size, FALSE);
        
        //Join
        while(left_rights->get_size()>0) {
            curr_left=left_rights->get_bottom();
            join=0;
            
            right_lefts->reset();
            curr_right=right_lefts->get_next();
            
            while((curr_right != 0) && (join == 0)) {
                    //cout<<"End: "<<left_end<<": Test joining "<<curr_left->track_pos<<" to "
                    // <<curr_right->track_pos<<endl;
                
                        join=assign_adjacency(curr_left, curr_right);
                        if (join == 1) {
                            //We made a join
                            mark_used(curr_left->track_pos, curr_right->track_pos);
                            
                            // cout<<"Joined L "<<size<<" "<<curr_left->my_locus->get_gene_obj(curr_left->index_num)->get_name()<<" at "<<curr_left->track_pos<<" to R: "<<curr_right->my_locus->get_gene_obj(curr_right->index_num)->get_name()<<" at "<<curr_right->track_pos<<endl;
                            
                            
                            //if (num_deep_left_joins == allowed_left_joins) {
                            //We can't allow further pass-through joins. At best, use any open genes at the right end of the left
                            // cout<<"Resetting stacks\n";
                            delete left_rights;
                         
                            delete right_lefts;
                            set_stacks(left_rights, left_pos, left_end-1, TRUE);
                            set_stacks(right_lefts, left_end+size_left, left_end+size, FALSE);
                        }
                        else
                            curr_right=right_lefts->get_next();
            }
            
            if (join == 0)  left_rights->pop();              //We've been through all the rights for this left and can't join it--get rid of it
 
        }
        
        delete left_rights;
        delete right_lefts;
        
    }
    else {
        //Base recursion
        if (size > 1) {
            //cout<<"Base Joining "<<inferred_tracking[left_end][0].my_locus->get_gene_obj(0)->get_name()<<" to "<<inferred_tracking[left_end+1][0].my_locus->get_gene_obj(0)->get_name()<<endl;
            //Check for all possible 2 gene joins
            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                for(j=0; j<the_homologs->get_dupl_level(); j++)
                    assign_adjacency(&inferred_tracking[left_end][i], &inferred_tracking[left_end+1][j]);
            }
            
        }   //Done with 2-gene join
    }                    //Done with <3 gene joins (1 or 2)
    
}




#if 0
int Genome_Track_DX::break_and_rejoin(int break_before, int size)
{
	int num_orig, num_new, save=-1;
	
	num_orig=find_partition_breaks(break_before);
	set_stacks(&new_right, break_before, size, FALSE);
	set_stacks(&new_left, break_before-1, -1, TRUE);
	num_new=make_joins(&new_left, &new_right, inferred_tracking[break_before-1],
		inferred_tracking[break_before], TRUE, save);
	
	new_left.clear();
	new_right.clear();

	return(num_new-num_orig);	
}
#endif



Tracking_List_DX * Genome_Track_DX::get_list_pos(int pos)
{
	int i;
	Tracking_List_DX *curr_pos;

	if ((last_list_pos > pos) || (last_list_pos == list_len)) {
		partial_track_pos=partial_track_start;
		last_list_pos=0;
	}

	curr_pos=partial_track_pos;

	if (pos > list_len) {
		cerr<<"Error: Trying to read past end of list\n";
		pos=list_len;
	}
	
	i=last_list_pos;

	while(i<pos-1) {
		curr_pos=curr_pos->next;
		i++;
	}
	return(curr_pos);
}


void Genome_Track_DX::set_stacks(Track_Stack_DX *&new_stack, int end, int stop_point, BOOL left)
//Note that stop point should be 1 bigger or smaller than the block end
{
	int position, j, num_paths;
	BOOL stop=FALSE;


    new_stack=new Track_Stack_DX();
    
	//Find all the possible left/right endpoints
    
	position=end;
    
	while((stop==FALSE) && (position != stop_point)) {
		if (left == FALSE) {
            for (j=0; j<the_homologs->get_dupl_level(); j++) {
                if ((inferred_tracking[position][j].my_locus != 0) && (inferred_tracking[position][j].last == 0))
                    new_stack->push(&inferred_tracking[position][j]);
            }
            
            num_paths=num_pass_throughs (position);
            
            if (num_paths < 1) stop=TRUE;
            position++;
           
		}
		else {
            for (j=0; j<the_homologs->get_dupl_level(); j++) {
                if ((inferred_tracking[position][j].my_locus != 0) && (inferred_tracking[position][j].next == 0))
                    new_stack->push(&inferred_tracking[position][j]);
            }
            
            num_paths=num_pass_throughs (position);
            
            if (num_paths < 1) stop=TRUE;
            position--;
		}
    }
}





void Genome_Track_DX::set_link_counts(Track_Stack_DX *lefts, Track_Stack_DX *rights)
{
	int join;
	Gene_Track_List_DX *curr_left, *curr_right;


	lefts->reset();
	curr_left=lefts->get_next();
	while(curr_left != 0) {
		curr_left->num_joins=0;
		curr_left=lefts->get_next();
	}

	rights->reset();
	curr_right=rights->get_next();
	while(curr_right != 0) {
		curr_right->num_joins=0;
		curr_right=rights->get_next();
	}
	
	//Do a preliminary check to find any ends that have more than one possible join
	lefts->reset();
	curr_left=lefts->get_next();
	while(curr_left != 0) {
			rights->reset();
			curr_right=rights->get_next();
			while(curr_right != 0) {			
				join=assign_adjacency(curr_left, curr_right, FALSE);
				curr_left->num_joins+=join;
				curr_right->num_joins+=join;
				curr_right=rights->get_next();
			}
			curr_left=lefts->get_next();
		}

}



#if 0
int Genome_Track_DX::make_joins(Track_Stack_DX *lefts, Track_Stack_DX *rights, Gene_Track_List_DX *last_left, Gene_Track_List_DX *new_right, BOOL do_link, int &save_index)
{
	int i, num_joins=0, join, runlen;
	Gene_Track_List_DX *curr_left, *curr_right, *valid_left, *valid_right;
    BOOL unjoined_at_locus;

	set_link_counts(lefts, rights);
		
	lefts->reset();
	//Join
	while((lefts->get_size()>0) && (num_joins<the_homologs->get_dupl_level())) {
			curr_left=lefts->get_bottom();
			rights->reset();
			join=0;
			curr_right=rights->get_next();
			while((curr_right != 0) && (join == 0)) {			
				//Catches the possiblity of more that one join for two ends--we don't want to join these two
				//because it precludes a second join elsewhere
				if ((curr_left->num_joins == the_homologs->get_dupl_level()) && (curr_right->num_joins == the_homologs->get_dupl_level()))
					join=0;
				else 
					join=assign_adjacency(curr_left, curr_right, do_link, save_index);

				if (join == 1) {
				//We made a join
					rights->delete_element(curr_right);
					lefts->pop();
					num_joins++;
					
                    if (num_joins == (the_homologs->get_dupl_level()-1))
                    {
                        unjoined_at_locus=FALSE;
                        for (i=0; i<the_homologs->get_dupl_level(); i++)
                            if (last_left[i].next ==0) unjoined_at_locus=TRUE;
                        
                        
                        if (unjoined_at_locus == TRUE) {
                            valid_left=last_left;
                            lefts->reset();
                            runlen=lefts->get_size();
                            for(i=0; i<runlen; i++) lefts->pop();
                            
                            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                                valid_left=&last_left[i];
                                
                                if (valid_left->next ==0)
                                    lefts->push(valid_left);
                            }
                        }
                        
                        unjoined_at_locus=FALSE;
                        for (i=0; i<the_homologs->get_dupl_level(); i++)
                            if (new_right[i].next ==0) unjoined_at_locus=TRUE;
                        
                        if (unjoined_at_locus == TRUE) {
                            valid_right=new_right;
                            rights->reset();
                            runlen=rights->get_size();
                            for(i=0; i<runlen; i++) rights->pop();
                            
                            for (i=0; i<the_homologs->get_dupl_level(); i++) {
                                valid_right=&new_right[i];
                                
                                if (valid_right->last ==0)
                                    rights->push(valid_right);
                            }

                        }
                        
                    }
                    
					
					set_link_counts(lefts, rights);

				}
				else
					curr_right=rights->get_next();

			}

			if (join == 0) {
			//We've been through all the rights for this left and can't join it--get rid of it
				lefts->pop(); }
		}

	lefts->restore_orig();
	rights->restore_orig();

	return(num_joins);
}
#endif




int Genome_Track_DX::find_partition_breaks(int break_before)
{
	int i, j, num_breaks=0;
	i=break_before-1;
	while((i>=0) &&(num_breaks < the_homologs->get_dupl_level())) {
		for(j=0; j<the_homologs->get_dupl_level(); j++) {
			if ((inferred_tracking[i][j].next != 0) && (inferred_tracking[i][j].next->partition_num >= break_before)) {
				broken_lefts[num_breaks]=&inferred_tracking[i][j];
				broken_rights[num_breaks]=inferred_tracking[i][j].next;
				inferred_tracking[i][j].next->last=0;
				inferred_tracking[i][j].next=0;
				num_breaks++;
			}
		}
		i--;
	}

	return(num_breaks);
}


void Genome_Track_DX::reset_partition_numbers()
{
	int i, j;

	for(i=0; i<the_homologs->get_num_homologs(); i++)
	{
        for (j=0; j<the_homologs->get_dupl_level(); j++)
            inferred_tracking[i][j].partition_num=i;
	}
}


void Genome_Track_DX::reset_list_partition_numbers()
{
    int i;
	Tracking_List_DX *use_pos;
	use_pos=partial_track_start;

    for(i=0; i<the_homologs->get_dupl_level(); i++)
        use_pos->element[i]->partition_num=0;
	
	use_pos=use_pos->next;
	while (use_pos != 0) {
        for(i=0; i<the_homologs->get_dupl_level(); i++)
            use_pos->element[i]->partition_num=use_pos->last->element[i]->partition_num+1;
		use_pos=use_pos->next;
	}

}




int Genome_Track_DX::assign_adjacency(Gene_Track_List_DX *locus1, Gene_Track_List_DX *locus2)
{
	int temp=0;
	return(assign_adjacency(locus1, locus2, TRUE, temp));
}


int Genome_Track_DX::assign_adjacency(Gene_Track_List_DX *locus1, Gene_Track_List_DX *locus2, BOOL do_link)
{
	int temp=0;
	return(assign_adjacency(locus1, locus2, do_link, temp));
}



int Genome_Track_DX::assign_adjacency(Gene_Track_List_DX *locus1, Gene_Track_List_DX *locus2, BOOL do_link, int &save_index)
{
	int retval=0;
	if ((locus1->my_locus != 0) && (locus2->my_locus != 0)) {
        if ((locus1->next == 0) && (locus2->last == 0)) {
            if (locus1->my_locus->get_contig(locus1->index_num) == locus2->my_locus->get_contig(locus2->index_num)) {
                if ((locus1->my_locus->get_gene(locus1->index_num) == locus2->my_locus->get_gene(locus2->index_num)+1) ||
                    (locus1->my_locus->get_gene(locus1->index_num) == locus2->my_locus->get_gene(locus2->index_num)-1)) {
                        {
                            //if (strcmp(locus1->my_locus->get_gene_obj(locus1->index_num)->get_name(),"Sbay_19.117") == 0){
                            //    cout<<"Looking at "<<locus1->my_locus->get_gene_obj(locus1->index_num)->get_name()<<" and "<<locus2->my_locus->get_gene_obj(locus2->index_num)->get_name()<<endl;
                            //}
                            
                            if (do_link == TRUE) {
                                locus1->next=locus2;
                                if (save_index != -1) {
                                    new_lasts[save_index]=locus1;
                                    new_nexts[save_index++]=locus2;
                                }
                                locus2->last=locus1;
                            }
                            retval=1;
                        }
                }
			}
		}
    }

	return(retval);

}




void Genome_Track_DX::store_tracking()
{
    int i, j, track_num, num_break, num_all_break;
    Gene_Track_List_DX *dummy, **lasts;
    
    lasts =new Gene_Track_List_DX * [the_homologs->get_dupl_level()];
    
    for(j=0; j<the_homologs->get_dupl_level(); j++) lasts[j]=0;
    
    add_track_starts();
    set_track_nums();
    
    //for(j=0; j<the_homologs->get_dupl_level(); j++) {
   //     cout<<"Site 9073: level "<<j<<" to track: "<<inferred_tracking[9073][j].to_track_num<<"\t";
    //    if (inferred_tracking[9073][j].my_locus!=0) cout<<inferred_tracking[9073][j].my_locus->get_gene_obj(inferred_tracking[9073][j].index_num)->get_name()<<endl;
   //     else cout<<"NONE\n";
   // }
    
    if (new_inferred_tracking==0){
        new_inferred_tracking = new Gene_Track_List_DX* [the_homologs->get_num_homologs()];
        
        for(i=0; i<the_homologs->get_num_homologs(); i++)
            new_inferred_tracking[i]=new Gene_Track_List_DX [the_homologs->get_dupl_level()];
        
        for(i=0; i<the_homologs->get_num_homologs(); i++) {
            for(j=0; j<the_homologs->get_dupl_level(); j++) {
                new_inferred_tracking[i][j].my_locus=0;
                new_inferred_tracking[i][j].to_track_num=j;
                new_inferred_tracking[i][j].index_num=-1;
                new_inferred_tracking[i][j].track_pos=i;
                new_inferred_tracking[i][j].last=0;
                new_inferred_tracking[i][j].next=0;
            }
        }
        
        for(i=0; i<the_homologs->get_num_homologs(); i++) {
            for(j=0; j<the_homologs->get_dupl_level(); j++) {
                if (inferred_tracking[i][j].my_locus != 0) {
                   // cout<<"Old tracking for "<<i<<" track "<<j<<" is "<<inferred_tracking[i][j].my_locus->get_gene_obj(inferred_tracking[i][j].index_num)->get_name()<<" headed for "<<inferred_tracking[i][j].to_track_num<<endl;
                    
                     dummy=&inferred_tracking[i][j];
                     track_num=inferred_tracking[i][j].to_track_num;
                    if ((track_num <0) || (track_num >= the_homologs->get_dupl_level())) {
                        dummy=&inferred_tracking[i][j];
                        cout<<"ERROR: track number from assembly for point "<<i<<" is invalid: "<<track_num<<" in "<<the_genome->get_name()<<endl;
                    }
                   // cout<<"Setting track "<<track_num<<" for point "<<i<<endl;
                    new_inferred_tracking[i][track_num].my_locus = inferred_tracking[i][j].my_locus;
                    new_inferred_tracking[i][track_num].index_num=inferred_tracking[i][j].index_num;
                    new_inferred_tracking[i][track_num].to_track_num=track_num;
                    
                    if (inferred_tracking[i][j].last != 0) {
                        new_inferred_tracking[i][track_num].last = &new_inferred_tracking[inferred_tracking[i][j].last->track_pos][track_num];
                        new_inferred_tracking[inferred_tracking[i][j].last->track_pos][track_num].next=&new_inferred_tracking[i][track_num];
                    }
                    if (inferred_tracking[i][j].next !=0) lasts[track_num]=&new_inferred_tracking[i][track_num];
                    else lasts[track_num]=0;
                    
                    //cout<<"Position "<<i<<": "<<track_num<<": "<<new_inferred_tracking[i][track_num].my_locus->get_gene_obj(new_inferred_tracking[i][track_num].index_num)->get_name()<<" has next: "<<inferred_tracking[i][j].next;
                    //if (inferred_tracking[i][j].next!=0)
                    //    cout<<"= "<<inferred_tracking[i][j].next->my_locus->get_gene_obj(inferred_tracking[i][j].next->index_num)->get_name()<<endl;
                    //else cout<<endl;
                    
                }
            }
            
            for(j=0; j<the_homologs->get_dupl_level(); j++) {
                if (new_inferred_tracking[i][j].my_locus == 0) {
                    new_inferred_tracking[i][j].last=lasts[j];
                   // cout<<"Pointing locus "<<i<<","<<j<<" to "<<lasts[j];
                   // if (lasts[j]!=0) cout<<lasts[j]->my_locus->get_gene_obj(lasts[j]->index_num)->get_name()<<endl;
                   // else cout<<endl;
                }
            }
            
        }
        
        //num_all_break=0;
        
        for(i=0; i<the_homologs->get_num_homologs(); i++)
        {
            //num_break=0;
            for (j=0; j<the_homologs->get_dupl_level(); j++) {
                if (new_inferred_tracking[i][j].last != 0) {
                    new_inferred_tracking[i][j].dist_to_last=new_inferred_tracking[i][j].track_pos-new_inferred_tracking[i][j].last->track_pos;
                }
               // else num_break++;
            }
            
            //if (num_break == the_homologs->get_dupl_level()) num_all_break++;
        }
        for(i=0; i<the_homologs->get_num_homologs(); i++)
            delete[] inferred_tracking[i];
        delete[] inferred_tracking;
        
        inferred_tracking=new_inferred_tracking;
        //cout<<"Assigned inferred+Tracking to "<<inferred_tracking<<endl;
        new_inferred_tracking=0;
        
        //cout<<"Completed tracking store for "<<the_genome->get_name()<<". Found "<<num_all_break<<" full breaks. Break count retrieves: "<<count_num_full_breaks()<<"\n";
        
       // for(j=0; j<the_homologs->get_dupl_level(); j++) {
       //     cout<<"AFTER Site 9073: level "<<j;
        //    if (inferred_tracking[9073][j].my_locus!=0) cout<<inferred_tracking[9073][j].my_locus->get_gene_obj(inferred_tracking[9073][j].index_num)->get_name()<<endl;
        //    else cout<<"NONE\n";
       // }
        delete[] lasts;
    }
    else {
        cerr<<"ERROR: there is an existing tracking stored at new_inferred_tracking\n";
    }
}






void Genome_Track_DX::set_null_lasts()
{
	int i, j, end;
	Gene_Track_List_DX *curr=0;

	for(j=0; j<the_homologs->get_dupl_level(); j++) {
		end=the_homologs->get_num_homologs()-1;
		while(inferred_tracking[end][j].my_locus==0) end--;
		curr=&inferred_tracking[end][j];

		for(i=end-1; i>=0; i--) {
			if (inferred_tracking[i][j].my_locus == 0) {
				inferred_tracking[i][j].last=curr->last;
			}
			else
				curr=&inferred_tracking[i][j];
			
		}
	}

}


#if 0
Gene_Track_List_DX ** Genome_Track_DX::create_new_track_locus(WGX_Locus *the_locus)
{
    int i;
	Gene_Track_List_DX **new_locus_data;

	new_locus_data=new Gene_Track_List_DX *[the_homologs->get_dupl_level()];
    
    for(i=0; i<the_homologs->get_dupl_level(); i++) {
        new_locus_data[i]=new Gene_Track_List_DX();
        
        if (the_locus->has_duplicate(i) ==TRUE) {
            new_locus_data[i]->my_locus=the_locus;
            new_locus_data[i]->index_num=i;
        }
    }
	return(new_locus_data);
}
#endif


BOOL Genome_Track_DX::check_list()
{
	int i;
	BOOL valid=TRUE;
	Tracking_List_DX *curr;
	curr=partial_track_start;
	

	if (partial_track_start->last !=0)
		valid=FALSE;
	if(partial_track_end->next != 0)
		valid=FALSE;

	for(i=1; i<list_len; i++) {
		curr=curr->next;

		if ((curr->next !=0) && (curr->next->last != curr))
			valid=FALSE;

		//if (curr->element[0]->my_locus ==0)
		//	valid=FALSE;
	}

	if (curr != partial_track_end)
		valid=FALSE;

	return(valid);
}


void Genome_Track_DX::add_track_starts()
{
    int i, j;
    
    
    for (i=0; i<the_homologs->get_num_homologs(); i++) {
        for (j=0; j<the_homologs->get_dupl_level(); j++) {
            if ((inferred_tracking[i][j].my_locus!=0) && (inferred_tracking[i][j].last ==0))
                track_starts.push(&inferred_tracking[i][j]);
        }
    }
}


void Genome_Track_DX::set_track_nums()
{
    int i, my_track;
    BOOL used;
    Gene_Track_List_DX *my_track_pos;
    
    while (track_starts.get_size()>0) {
        my_track_pos=track_starts.get_bottom();
        my_track=0;
        //if (my_track_pos->track_pos==9073) cout<<"Placing 9073 as a start "<<my_track_pos->track_pos<<endl;
        do {
            used = FALSE;
            for(i=0; i<the_homologs->get_dupl_level(); i++) {
                if (inferred_tracking[my_track_pos->track_pos][i].to_track_num == my_track) used=TRUE;
            }
            if (used == TRUE) my_track++;
        } while((my_track < the_homologs->get_dupl_level()) && (used==TRUE));
       
        if (my_track == the_homologs->get_dupl_level()) cerr<<"ERROR: at location "<<my_track_pos->track_pos<<" all possible tracking positions are in use\n";
        
        my_track_pos->to_track_num=my_track;
        //cout<<"Placing "<<my_track_pos->my_locus->get_gene_obj(my_track_pos->index_num)->get_name()<<" to "<<my_track<<endl;
        
        if (my_track_pos->next !=0) extend_track(my_track_pos, my_track_pos->next, my_track);
        
        track_starts.pop();
    }
}


void Genome_Track_DX::extend_track(Gene_Track_List_DX *my_pos, Gene_Track_List_DX *next_pos, int track_num)
{
    int i, next_track, curr_pos;
    BOOL stop;
    
    curr_pos = my_pos->track_pos+1;
 
    while(curr_pos != next_pos->track_pos) {
        next_track=0;
        
        for (i=0; i<the_homologs->get_dupl_level(); i++) {
            if (inferred_tracking[curr_pos][i].to_track_num == track_num) cerr<<"ERROR at position "<<curr_pos<<": track number "<<track_num<<" is already assigned to "<<i<<endl;
        }
        stop=FALSE;
        do  {
            if ((inferred_tracking[curr_pos][next_track].my_locus == 0) && (inferred_tracking[curr_pos][next_track].to_track_num == -1)) stop = TRUE;
            else next_track++;
            
            if (next_track == the_homologs->get_dupl_level()) stop=TRUE;
            
        } while ( stop == FALSE);
        
        if (next_track < the_homologs->get_dupl_level())
            inferred_tracking[curr_pos][next_track].to_track_num=track_num;
        else cerr<<"ERROR: No empty tracks avalaible at position "<<curr_pos<<endl;
        curr_pos++;
    }
    
    for (i=0; i<the_homologs->get_dupl_level(); i++) {
        if (inferred_tracking[next_pos->track_pos][i].to_track_num == track_num) cerr<<"ERROR at Next site position "<<next_pos->track_pos<<": track number "<<track_num<<" is already assigned to "<<i<<endl;
    }
    
    next_pos->to_track_num=track_num;
    
    if (next_pos->next !=0) extend_track(next_pos, next_pos->next, track_num);
}

WGX_Tracks::WGX_Tracks(WGX_Data *homologs, Clade *genomes)
{
	int i;

	the_homologs=homologs;
	the_genomes=genomes;
	tracking_correct=FALSE;

	order=new int[the_homologs->get_num_homologs()];
	fully_connected=new BOOL[the_homologs->get_num_homologs()+1];
	move_locations = new int [the_homologs->get_num_homologs()];
	move_section_starts = new int [the_homologs->get_num_homologs()];
	move_section_ends = new int [the_homologs->get_num_homologs()];
	real_to_move_locs = new int [the_homologs->get_num_homologs()+1];

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=i;

	the_trackings=new Genome_Track_DX * [the_genomes->get_num_genomes()];
    for(i=0; i<the_genomes->get_num_genomes(); i++)
        the_trackings[i]=new Genome_Track_DX (i, &(*the_genomes)[i], the_homologs);

    //for(i=0; i<the_genomes->get_num_genomes(); i++) cout<<"Tracking "<<i<<" is genome "<<the_trackings[i]->get_genome()->get_name()<<endl;
    
    single_genome_num_trackings=recurse_factorial(the_homologs->get_dupl_level());
    num_orders=single_genome_num_trackings;
    
    for(i=1;i<the_genomes->get_num_genomes(); i++) num_orders *=single_genome_num_trackings;
    
    //std::cout<<"There are "<<num_orders<<" total track orderings to check\n";
}


WGX_Tracks::WGX_Tracks(WGX_Data *homologs, Clade *genomes, int *new_order)
{
    int i;
    
    the_homologs=homologs;
    the_genomes=genomes;
    tracking_correct=FALSE;
    
    order=new int[the_homologs->get_num_homologs()];
    fully_connected=new BOOL[the_homologs->get_num_homologs()+1];
    move_locations = new int [the_homologs->get_num_homologs()];
    move_section_starts = new int [the_homologs->get_num_homologs()];
    move_section_ends = new int [the_homologs->get_num_homologs()];
    real_to_move_locs = new int [the_homologs->get_num_homologs()+1];
    
    for(i=0; i<the_homologs->get_num_homologs(); i++)
        order[i]=new_order[i];
    
    the_trackings=new Genome_Track_DX * [the_genomes->get_num_genomes()];
    for(i=0; i<the_genomes->get_num_genomes(); i++)
        the_trackings[i]=new Genome_Track_DX (i, &(*the_genomes)[i], the_homologs);
    
    //for(i=0; i<the_genomes->get_num_genomes(); i++) cout<<"Tracking "<<i<<" is genome "<<the_trackings[i]->get_genome()->get_name()<<endl;
    
    single_genome_num_trackings=recurse_factorial(the_homologs->get_dupl_level());
    num_orders=single_genome_num_trackings;
    
    for(i=1;i<the_genomes->get_num_genomes(); i++) num_orders *=single_genome_num_trackings;
    
    //std::cout<<"There are "<<num_orders<<" total track orderings to check\n";
}


 int WGX_Tracks::get_num_positions_w_breaks()
{
    int pillar, taxa, dupl_level, num_w_breaks=0;
    BOOL has_break;
    for (pillar=0; pillar<the_homologs->get_num_homologs(); pillar++) {

    //for (pillar=1; pillar<the_homologs->get_num_homologs(); pillar++) {
        has_break=FALSE;
        for(taxa=0; taxa<the_genomes->get_num_genomes(); taxa++) {
            for (dupl_level =0; dupl_level<the_homologs->get_dupl_level(); dupl_level++) {
                if ((the_trackings[taxa]->get_gene_track(pillar, dupl_level)->my_locus !=0) && (the_trackings[taxa]->has_back_link(pillar, dupl_level) == FALSE))
                //if ((the_trackings[taxa]->has_back_link(pillar, dupl_level) == FALSE))
                    has_break=TRUE;
            }
        }
        if (has_break ==TRUE) num_w_breaks++;
    }
    return(num_w_breaks);
}

int WGX_Tracks::get_num_breaks()
{
	int i, retval=0;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		retval+=the_trackings[i]->count_num_breaks();

	return(retval);

	
}


int WGX_Tracks::get_num_full_breaks()
{
    int i, retval=0;
    
    for(i=0; i<the_genomes->get_num_genomes(); i++)
        retval+=the_trackings[i]->count_num_full_breaks();
    
    return(retval);
    
    
}



int WGX_Tracks::get_num_list_breaks()
{
	int i, retval=0;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		retval+=the_trackings[i]->count_num_list_track_breaks();

	return(retval);
}



void WGX_Tracks::change_order(int *new_order)
{
	int i;

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=new_order[i];

	tracking_correct=FALSE;
}



void WGX_Tracks::change_order()
{
	int i;

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=the_trackings[0]->get_list_position_number(i+1);

	tracking_correct=FALSE;
}



void WGX_Tracks::update_tracking()
{
	int i;
	
	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_trackings[i]->number_contigs(this->order);

	tracking_correct=TRUE;
}



void WGX_Tracks::print_all_tracks()
{
	print_all_tracks(TRUE);
}



void WGX_Tracks::print_all_tracks(BOOL use_file)
{
	int i;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_trackings[i]->print_tracking(use_file);
}
	

void WGX_Tracks::print_homolog_file(string outfile)
{
    int i, j, k;
    ofstream fout;
    
    fout.open(outfile.c_str());
    if( fout.fail()) {
        cerr<<"ERROR: Cannot open output file "<<outfile<<endl;
        return;
    }
    else {
        fout<<(*the_genomes)[0].get_name_string();
        for(i=1; i<the_genomes->get_num_genomes(); i++)
            fout<<"\t"<<(*the_genomes)[i].get_name_string();
        fout<<endl;
        
        for(i=0; i<the_homologs->get_num_homologs(); i++) {
            for(k=0; k<the_homologs->get_dupl_level(); k++) {
                for (j=0; j<the_genomes->get_num_genomes(); j++) {
                    if ((k==0) && (j==0)) {
                        if ((*the_homologs)[order[i]][j].has_duplicate(k)==TRUE)
                            fout<<(*the_homologs)[order[i]][j].get_gene_obj(k)->get_name_string();
                        else
                            fout<<"NONE";
                    }
                    else {
                        if ((*the_homologs)[order[i]][j].has_duplicate(k)==TRUE)
                            fout<<"\t"<<(*the_homologs)[order[i]][j].get_gene_obj(k)->get_name_string();
                        else
                            fout<<"\tNONE";
                    }
                }
            }
            fout<<endl;
        }
    }
    fout.close();
}



void WGX_Tracks::set_all_null_lasts()
{
	int i;

	for(i=0; i<the_genomes->get_num_genomes(); i++) 
		the_trackings[i]->set_null_lasts();
}



Genome_Track_DX & WGX_Tracks::operator[] (int index)
{
	return((*the_trackings[index]));
}
	


void WGX_Tracks::swap_elements (int index1, int index2)
{
	int temp;

	temp=order[index1];
	order[index1]=order[index2];
	order[index2]=temp;

	tracking_correct=FALSE;

}

void WGX_Tracks::get_move_section_start_end(int move_pos, int &start, int &end)
{
	int i;

	i=0;
	while(move_section_ends[i] < move_locations[move_pos]) i++;

	start=real_to_move_locs[move_section_starts[i]];
	end=real_to_move_locs[move_section_ends[i]];
}








WGX_Tracks::~WGX_Tracks()
{
    int i;
    
    delete[] order;
    if (the_trackings != 0) {
        for(i=0; i<the_genomes->get_num_genomes(); i++)
            delete[] the_trackings[i];
        delete[] the_trackings;
    }
	delete[] move_locations;
	delete[] move_section_starts;
	delete[] move_section_ends;
	delete[] real_to_move_locs;
    delete[] fully_connected;
}



void WGX_Tracks::set_array_for_pos(int pos, int val, int size, int *array)
{
	int i;

	for(i=size-1; i>pos; i--) 
		array[i]=array[i-1];

	array[pos]=val;
}



void WGX_Tracks::check_for_fully_connected_positions()
{
	int i, j, num_joins, num_fully_connected=0, pos;
	BOOL all_break;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_trackings[i]->reset_partition_numbers();

	fully_connected[0]=FALSE;
	for(i=1; i<the_homologs->get_num_homologs(); i++)
		cout<<i<<"\t";
	cout<<endl;

	for(i=1; i<the_homologs->get_num_homologs(); i++) {
	
		fully_connected[i] = TRUE;
		for(j=0; j<the_genomes->get_num_genomes(); j++) 
			if (the_trackings[j]->all_poss_joins_after(i-1) == FALSE)
				fully_connected[i] = FALSE;
		cout<<fully_connected[i]<<"\t";
		if (fully_connected[i] == TRUE)
			 num_fully_connected++;
	
	}
	cout<<endl;

	fully_connected[the_homologs->get_num_homologs()]=FALSE;

	num_pseudo_chroms=0;
	move_section_starts[0]=0;
	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		all_break=TRUE;
		for(j=0; j<the_genomes->get_num_genomes(); j++) {
			if (the_trackings[j]->has_tracked_full_break(i) == FALSE)
				all_break=FALSE;
		}
		if (all_break == TRUE) {
			move_section_ends[num_pseudo_chroms++]=i-1;
			move_section_starts[num_pseudo_chroms]=i;
		}

	}
	move_section_ends[num_pseudo_chroms++]=the_homologs->get_num_homologs();

	move_locations[0]=0;
	real_to_move_locs[0]=0;
	pos = 0;

	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		if (fully_connected[i] == FALSE) 
			move_locations[++pos]=i;
		real_to_move_locs[i]=pos;

	}
	//Account for a possible move at the end of the tracking
	move_locations[++pos]=the_homologs->get_num_homologs();
	real_to_move_locs[the_homologs->get_num_homologs()]=pos;
	num_move_locations=pos+1;

	cout<<"Tracking has "<<num_fully_connected<<" pillar joins\n";
}
