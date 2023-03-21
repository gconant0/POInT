#include "genome_list.h"



Gene::Gene()
{
	name="";
	cerr<<"Error: call to default constructor of class Gene\n";
	neighbors=new Gene*[2];
    neighbors[0]=0;
    neighbors[1]=0;
	num_neighbors=0;
	keep=TRUE;
    start_pos=-1;
    end_pos=-1;
    chrom_string="NONE";
}



Gene::Gene(int num, char *new_name)
{
	gene_num=num;
	name=new_name;
	neighbors=new Gene* [2];
    neighbors[0]=0;
    neighbors[1]=0;
	keep = TRUE;
    num_neighbors=0;
    start_pos=-1;
    end_pos=-1;
    chrom_string="NONE";
}

Gene::Gene(int num,  string new_name)
{
    gene_num=num;
    name=new_name;
    neighbors=new Gene* [2];
    neighbors[0]=0;
    neighbors[1]=0;
    keep = TRUE;
    num_neighbors=0;
    start_pos=-1;
    end_pos=-1;
    chrom_string="NONE";
}


Gene& Gene::operator= (Gene &assign_from)
{
	gene_num=assign_from.get_gene_num();
	name=assign_from.get_name();
    //num_neighbors=assign_from.get_num_neighbors();
	return(*this);
}


Gene * Gene::get_neighbor(int index)
{
	if (index == 0)
		return(neighbors[0]);
	else
		return(neighbors[1]);
}
	

void Gene::set_neighbor(Gene * new_neighbor, int index)
{	
    if (index == 0) {
        if ((neighbors[0] == 0) && (new_neighbor !=0)) num_neighbors++;
		neighbors[0]=new_neighbor;
    }
    else {
        if ((neighbors[1] == 0) && (new_neighbor !=0)) num_neighbors++;
		neighbors[1]=new_neighbor;
	
    }
}
	
void Gene::set_location(string ch_string, int s, int e)
{
    chrom_string=ch_string;
    start_pos=abs(s);
    end_pos=abs(e);
    if (start_pos>end_pos) {
        end_pos=abs(s);
        start_pos=abs(e);
    }
    
}

Gene::~Gene()
{
	delete[] neighbors;
}


Contig::Contig()
{
    string dummy_name;
	num_genes=0;
	the_genes=0;
	
    dummy_name="Dummy Gene";
    null_gene=new Gene(0, dummy_name.c_str());
	cerr<<"Error: call to default constructor of class Contig\n";

}



Contig::Contig(int ngenes, char **gene_names)
{
	int i;
    string dummy_name;
    
	num_genes=ngenes;
	the_genes=0;
	

	if (num_genes > 0) {
		the_genes = new Gene* [num_genes];

		for(i=0; i<num_genes; i++)
			the_genes[i]=new Gene(i, gene_names[i]);
	}

    dummy_name="DummyGene";
	null_gene=new Gene(0, dummy_name.c_str());

	assign_neighbors();
}


Contig::Contig(int ngenes, string *gene_names)
{
    int i;
    string dummy_name;
    
    num_genes=ngenes;
    the_genes=0;
    
    
    if (num_genes > 0) {
        the_genes = new Gene* [num_genes];
        
        for(i=0; i<num_genes; i++)
        the_genes[i]=new Gene(i, gene_names[i]);
    }
    
    dummy_name="DummyGene";
    null_gene=new Gene(0, dummy_name.c_str());
    
    assign_neighbors();
}




Gene& Contig::operator[] (int element)
{
	if ((0<=element) &&(element<num_genes))
		return(*the_genes[element]);
	else
		return(*null_gene);
}


Contig& Contig::operator=(Contig & assign_from)
{
	int i;
    string temp_name;
    
    temp_name="";

	if (num_genes != assign_from.get_num_genes()) {
		for(i=0; i<num_genes; i++)
			delete the_genes[i];
		if (the_genes != 0)
			delete[] the_genes;		
		num_genes=assign_from.get_num_genes();
		the_genes=new Gene * [num_genes];
		for(i=0; i<num_genes; i++)
			the_genes[i]=new Gene(i, temp_name.c_str());
	}


	for(i=0; i<num_genes; i++)
		(*the_genes[i])=assign_from[i];

	assign_neighbors();

	return(*this);

}


Contig::~Contig()
{
	int i;

	delete null_gene;

	if (the_genes != 0)
	{
		for(i=0; i<num_genes; i++)
			delete the_genes[i];
		delete[] the_genes;
	}

}


void Contig::assign_neighbors()
{
	int i;

	for(i=0; i<num_genes; i++) {
		if (i!=0)
			the_genes[i]->set_neighbor(the_genes[i-1], 0);
		else
			the_genes[i]->set_neighbor(0, 0);

		if (i!= num_genes-1)
			the_genes[i]->set_neighbor(the_genes[i+1], 1);
		else
			the_genes[i]->set_neighbor(0,1);
	}
}


List_Contig::List_Contig()
{
	last=next=0;
	the_contig=0;
	cerr<<"Error: Call to default constructor of class List_Contig\n";
}



List_Contig::List_Contig(List_Contig *lst, Contig *contig)
{
	next=0;
	last=lst;
	the_contig=contig;
}



Genome::Genome()
{
    string *empty=0;
    
    web_link="NONE";
	num_contigs=0;
	the_contigs=0;
	
	null_contig = new Contig(0, empty);
	cerr<<"Error: Call to default constructor of class Genome\n";
}



Genome::Genome(int ncontigs, char *name, List_Contig *start_contigs)
{
	int i;
    string *empty=0;
    
    web_link="NONE";
	num_contigs=ncontigs;
	genome_name=name;

	if (num_contigs != 0) {
		the_contigs = new Contig* [num_contigs];
		for(i=0; i<num_contigs; i++)
		{
			the_contigs[i] = new Contig(0,empty);
			(*the_contigs[i])=(*start_contigs->the_contig);
			start_contigs=start_contigs->next;
		}
	}	
	else
		the_contigs=0;

	null_contig = new Contig(0, empty);
}

Genome::Genome(int ncontigs, string name, List_Contig *start_contigs)
{
    int i;
    string *empty=0;
    
    web_link="NONE";
    num_contigs=ncontigs;
    genome_name=name;
    
    if (num_contigs != 0) {
        the_contigs = new Contig* [num_contigs];
        for(i=0; i<num_contigs; i++)
        {
            the_contigs[i] = new Contig(0,empty);
            (*the_contigs[i])=(*start_contigs->the_contig);
            start_contigs=start_contigs->next;
        }
    }	
    else
    the_contigs=0;
    
    null_contig = new Contig(0, empty);
}


Contig& Genome::operator [] (int element)
{
	if ((0<=element) && (element < num_contigs))
		return(*(the_contigs[element]));
	else
		return(*null_contig);
}



Genome& Genome::operator= (Genome &assign_from)
{
	int i;
    string *empty=0;
   

    
	if (num_contigs != assign_from.get_num_contigs() ) {
		for(i=0; i<num_contigs; i++)
			delete the_contigs[i];
		if (the_contigs != 0)
		delete[] the_contigs;

		num_contigs=assign_from.get_num_contigs();

		the_contigs=new Contig*[num_contigs];
		for(i=0; i<num_contigs; i++)
			the_contigs[i]=new Contig(0,empty);
	}

	genome_name=assign_from.get_name_string();
    web_link=assign_from.get_web_link();
    
	for(i=0; i<num_contigs; i++)
		(*the_contigs[i])=assign_from[i];

	return(*this);

}


void Genome::print_genome(string filename)
{
    int i,j;
    ofstream fout;
    
    fout.open(filename.c_str());
    
    if (fout.fail()) {cerr<<"ERROR: Could not open output file "<<filename<<". Exiting\n"<<endl;}
    else {
        fout<<genome_name<<endl;
        
        for(i=0; i<num_contigs; i++) {
            cout<<"Genome "<<genome_name<<", Contig "<<i<<" has "<<(*the_contigs)[i].get_num_genes()<<" genes\n";
            for(j=0; j<(*the_contigs)[i].get_num_genes(); j++) {
                fout<<i+1<<"\t"<<(*the_contigs)[i][j].get_name_string()<<endl;
            }
        }
        fout.close();
    }
}

Genome::~Genome()
{
	int i;

	delete null_contig;

	if (num_contigs != 0)
	{
		for(i=0; i<num_contigs; i++)
			delete the_contigs[i];
		delete[] the_contigs;

	}
}



Clade::Clade()
{
    string genome_name;
    
    genome_name="NONE";
	num_genomes=0;
	the_genomes=0;
	null_genome=new Genome(0, genome_name.c_str(), 0);
	cerr<<"Error: Call to default constructor of class Clade\n";

}



Clade::Clade(int ngenomes, Genome **genomes)
{
	int i;
    string genome_name;
    
    genome_name="NONE";
    
	num_genomes=ngenomes;
   
	if (num_genomes != 0) {
		the_genomes=new Genome * [num_genomes];

		for(i=0; i<num_genomes; i++) {
			the_genomes[i]=new Genome(0, genome_name.c_str(), 0);
			(*the_genomes[i])=(*genomes[i]);
		}
	}
	else 
		the_genomes=0;

	null_genome=new Genome(0, genome_name.c_str(), 0);
}



Genome& Clade::operator [] (int element)
{
	if ((0<=element) && (element < num_genomes))
		return(*the_genomes[element]);
	else
		return(*null_genome);
}



Clade::~Clade()
{
	int i;

	delete null_genome;

	if(num_genomes > 0) {
		for(i=0; i<num_genomes; i++)
			delete the_genomes[i];
		delete[] the_genomes;
	}



}




List_Names::List_Names()
{
	next=last=0;
	name="";
	cerr<<"Error: Call to default constructor of class List_Names\n";
}



List_Names::List_Names(List_Names *lst, char *new_name)
{
	next=0;
	last=lst;
	name=new_name;
}

List_Names::List_Names(List_Names *lst, string new_name)
{
    next=0;
    last=lst;
    name=new_name;
}


Read_Genome::Read_Genome()
{
	start_names=list_names=0;
	start_contigs=list_contigs=0;

	//Possible bug--code assumes contigs are numbered starting at a
	//number greater than -1
	last_contig_num=-1;
}



Genome * Read_Genome::get_genome(const char *filename)
{
	string genome_name;
	Contig *next_contig;
	Genome *the_genome;

	infile.open(filename);

	if (!infile.fail()) {
		num_contigs=0;
		infile>>genome_name;

		while(!infile.eof()) {
			next_contig=get_contig();
			
			if (start_contigs == 0) {
				start_contigs=new List_Contig(0, next_contig);
				list_contigs=start_contigs;
			}
			else {
				list_contigs->next=new List_Contig(list_contigs, next_contig);
				list_contigs=list_contigs->next;
			}
			num_contigs++;
		}

		the_genome=new Genome(num_contigs, genome_name, start_contigs);
		infile.close();
		list_contigs=start_contigs;

		while(list_contigs != 0)
		{
			delete list_contigs->the_contig;
			start_contigs=list_contigs->next;
			delete list_contigs;
			list_contigs=start_contigs;
		}
		return(the_genome);

	}
	else 
		return(0);

}


Genome * Read_Genome::get_genome(string filename)
{
    string genome_name;
    Contig *next_contig;
    Genome *the_genome;
    
    infile.open(filename.c_str());
    
    if (!infile.fail()) {
        num_contigs=0;
        infile>>genome_name;
        //cout<<"Reading genome "<<genome_name<<endl;
        while(!infile.eof()) {
            next_contig=get_contig();
            
            if (start_contigs == 0) {
                start_contigs=new List_Contig(0, next_contig);
                list_contigs=start_contigs;
            }
            else {
                list_contigs->next=new List_Contig(list_contigs, next_contig);
                list_contigs=list_contigs->next;
            }
            num_contigs++;
            //cout<<"Read contig "<<num_contigs<<endl;
        }
        
        the_genome=new Genome(num_contigs, genome_name, start_contigs);
        cout<<"Creating genome for "<<genome_name<<" with "<<num_contigs<<" contigs"<<endl;
        infile.close();
        list_contigs=start_contigs;
        
        while(list_contigs != 0)
        {
            delete list_contigs->the_contig;
            start_contigs=list_contigs->next;
            delete list_contigs;
            list_contigs=start_contigs;
        }
        return(the_genome);
        
    }
    else 
    return(0);
    
}


Contig* Read_Genome::get_contig()
{
	int i, num_genes;
    char dump;
	string new_name;
	Contig *new_contig;

	

	contig_num=-1;
	
	if(!infile.eof()) {
		infile>>contig_num>>new_name;
		infile.get(dump);
		if (last_contig_num == -1) {
			last_contig_num=contig_num;
			start_names=new List_Names(0, new_name);
			list_names=start_names;
			num_genes=1;
			contig_num=-1;
			infile>>contig_num>>new_name;
			infile.get(dump);
		}
		else
		{
			start_names = new List_Names(0, next_contig_name1);
			list_names=start_names;
			num_genes=1;
		}

	
		while( (contig_num == last_contig_num) && ( (!infile.eof()) || (contig_num != -1)) )
		{
			list_names->next=new List_Names(list_names, new_name);
			list_names=list_names->next;
		
			num_genes++;
			contig_num=-1;
			infile>>contig_num>>new_name;
            //cout<<"Read "<<contig_num<<" and name "<<new_name<<"("<<num_genes<<")"<<endl;
			infile.get(dump);
		}
		
		if ((!infile.eof()) || (contig_num != -1))
			next_contig_name1=new_name;
		else {
			//Save the last gene in the last contig
			if ((new_name != next_contig_name1) && (contig_num != -1)) {
				list_names->next=new List_Names(list_names, new_name);
				list_names=list_names->next;
		
				num_genes++;
			}

		}
		last_contig_num=contig_num;

		name_list=new string [num_genes];



		list_names=start_names;

		for(i=0; i<num_genes; i++)  {
			name_list[i]=list_names->name;
			list_names=list_names->next;
		}

		list_names=start_names;
		while(list_names != 0) {
			start_names=list_names->next;
			delete list_names;
			list_names=start_names;
		}

		new_contig = new Contig (num_genes, name_list);
		return(new_contig);
	}
    else return(0);


}



WGD_Locus::WGD_Locus()
{
	species_num=contig_index[0]=contig_index[1]=gene_index[0]=gene_index[1]=0;
	strcpy(name1, "NONE");
	strcpy(name2, "NONE");
	has_second=FALSE;

}



WGD_Locus::WGD_Locus(int sp_num, char *n1, char *n2, Genome *genome)
{
	species_num=sp_num;

	strcpy(name1, n1);
	has_second=FALSE;
	the_genome=genome;

	find_id(name1, contig_index[0], gene_index[0]);
	

	if (strcmp(n2, "NONE") != 0) {
		strcpy(name2, n2);
		has_second=TRUE;
		find_id(name2, contig_index[1], gene_index[1]);
	}

}

Gene * WGD_Locus::get_gene_obj(int index)					
{
	if ((index ==0) || ((index ==1) && (has_second == TRUE)))
		return(&(*the_genome)[contig_index[index]][gene_index[index]]);
	else
		return(0);

}




void WGD_Locus::find_id(char *name, int &contig_id, int &gene_id)
{
	int i, j;
	BOOL found=FALSE;


	i=0;

	while((i<the_genome->get_num_contigs()) && (found == FALSE)) {
		j=0;
		while((j<(*the_genome)[i].get_num_genes()) && (found == FALSE)) {
			if (strcmp(name, (*the_genome)[i][j].get_name()) == 0) {
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




Homologs::Homologs()
{
	the_loci=0;
	null_locus=new WGD_Locus();
	the_genomes=0;
}



Homologs::Homologs(char **first_orthos, char **second_orthos, Clade *genomes)
{
	int i;

	the_genomes=genomes;
	null_locus=new WGD_Locus();

	the_loci =new WGD_Locus * [the_genomes->get_num_genomes()];
	
	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_loci[i]=new WGD_Locus(i, first_orthos[i], second_orthos[i], &(*the_genomes)[i]);

}


WGD_Locus& Homologs::operator[] (int element)
{
	if (the_genomes ==0)
		return(*null_locus);

	if ((0<=element) && (element < the_genomes->get_num_genomes()))
		return(*the_loci[element]);
	else
		return(*null_locus);

}



Homologs::~Homologs()
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



WGD_Data::WGD_Data()
{
	num_homologs=0;
	null_homolog=new Homologs();

	the_homologs=0;
}
	

WGD_Data::WGD_Data(int nhomologs, char ***first_orthos, char ***second_orthos, Clade *genomes)
{
	int i;

	num_homologs=nhomologs;
	null_homolog=new Homologs();

	the_genomes=genomes;

	the_homologs=new Homologs* [num_homologs];

	for(i=0; i<num_homologs; i++)
		the_homologs[i] = new Homologs(first_orthos[i], second_orthos[i], the_genomes);

}
	

Homologs & WGD_Data::operator[] (int element)
{
	if ((0<=element) && (num_homologs > element))
		return(*the_homologs[element]);
	else
		return(*null_homolog);
}


WGD_Data::~WGD_Data()
{
	int i;

	delete null_homolog;

	if (num_homologs != 0) {
		for(i=0; i<num_homologs; i++)
			delete the_homologs[i];
		delete[] the_homologs;
	}
}



Read_WGD_Data::Read_WGD_Data()
{
	species_indexes=0;
	genome_names=0;
	first_orthos=second_orthos=0;	
}
	


WGD_Data * Read_WGD_Data::get_data(char *filename, int &num_sites, Clade *genomes)
{
	int i, j;
	char line[2000];
	WGD_Data *return_data;

	the_genomes=genomes;
	num_genomes=the_genomes->get_num_genomes();

	num_sites=0;

	infile.open(filename);

	if(!infile.fail()) {
		genome_names=new char * [the_genomes->get_num_genomes()];
		species_indexes = new int [the_genomes->get_num_genomes()];
		

		for(i=0; i<the_genomes->get_num_genomes(); i++)
		{
			genome_names[i]=new char [50];
			infile>>genome_names[i];
		}

		set_indexes();
		while(!infile.eof()) {
			infile.getline(line, 1999);
			num_sites++;
		}

		//Should have tried one failed read
		num_sites -=2;
		sites=num_sites;

		//Reset the file to actually read the data--lazy
		infile.close();
		infile.clear();
		infile.open(filename);


		first_orthos=new char ** [num_sites];
		second_orthos=new char ** [num_sites];

		for(i=0; i<num_sites; i++) {
			first_orthos[i]=new char * [the_genomes->get_num_genomes()];
			second_orthos[i]=new char * [the_genomes->get_num_genomes()];

			for(j=0; j<the_genomes->get_num_genomes(); j++) {
				first_orthos[i][j]=new char[50];
				second_orthos[i][j]=new char[50];
			}
		}

		//Clear the Genome Name line
		infile.getline(line, 1999);

		for(i=0; i<num_sites; i++) {
			//Read the pillars
			for(j=0; j<the_genomes->get_num_genomes(); j++)
				infile>>first_orthos[i][species_indexes[j]];
			for(j=0; j<the_genomes->get_num_genomes(); j++)
				infile>>second_orthos[i][species_indexes[j]];

		}
		infile.close();
		return_data = new WGD_Data(num_sites, first_orthos, second_orthos, the_genomes);
		return(return_data);
	}
    else return(0);


}



Read_WGD_Data::~Read_WGD_Data()
{
	int i, j;


	if (genome_names !=0) {
		for(i=0; i<num_genomes; i++)
			delete[] genome_names[i];
		delete[] genome_names;
	}


	if (species_indexes != 0)
		delete[] species_indexes;


	if (first_orthos != 0) {
		for(i=0; i<sites; i++) {
			for(j=0; j<num_genomes; j++) {
				delete[] first_orthos[i][j];
				delete[] second_orthos[i][j];
			}
			delete[] first_orthos[i];
			delete[] second_orthos[i];
		}
		delete[] first_orthos;
		delete[] second_orthos;
	}
}



void Read_WGD_Data::set_indexes()
{
	int i, j;
	

	for(i=0; i<the_genomes->get_num_genomes(); i++) {
		j=0;
		while((j<the_genomes->get_num_genomes()) && 
			(strcmp((*the_genomes)[j].get_name(), genome_names[i])!=0))	j++;

		

		if (j == the_genomes->get_num_genomes())
			cerr<<"Error: could not find genome "<<genome_names[i]<<" in data\n";
		else
			species_indexes[i]=j;
	
	}

}


Gene_Track_List::Gene_Track_List()
{
	my_locus=0;
	index_num=0;
	next=last=0;
	partition_num=0;
}


Gene_Track_List& Gene_Track_List::operator= (Gene_Track_List &assign_from)
{
	my_locus=assign_from.my_locus;
	next=assign_from.next;
	last=assign_from.last;
	index_num=assign_from.index_num;
	return(*this);
}


Reorg_Pattern::Reorg_Pattern()
{
	last_join_same[0]=last_join_same[1]=FALSE;
	next_join_same[0]=next_join_same[1]=FALSE;
	valid_without_add=TRUE;
}

int	Reorg_Pattern::get_piece_num(int n)
{
	//Don't error check for speed
	return(piece_num[n]);
}
	

int Reorg_Pattern::get_reversal(int n)
{	
	//Don't error check for speed
	return(reversal[n]);
}
	

void Reorg_Pattern::set_piece_num(int n, int val)
{
	if ((n>=0) && (n<3))
		piece_num[n]=val;
	else
		piece_num[0]=val;
}


void Reorg_Pattern::set_reversal(int n, int val)
{
	if ((n>=0) && (n<3))
		reversal[n]=val;
	else
		reversal[0]=val;
}



All_Reorg_Patterns::All_Reorg_Patterns()
{
	int i, j;

	num_patterns=26;

	the_patterns=new Reorg_Pattern [num_patterns];
	//Null all reversals so we only have to set real reversals
	for (i=0; i<26; i++) {
		for(j=0; j<3; j++)
			the_patterns[i].set_reversal(j, 0);
	}

	//How many sections?
	for (i=0; i<2; i++)
		the_patterns[i].set_num_pieces(1);

	for(i=2; i<6; i++)
		the_patterns[i].set_num_pieces(2);

	for(i=6; i<26; i++)
		the_patterns[i].set_num_pieces(3);


	//Where is the new pillar?
	the_patterns[0].set_insert_loc(0);
	the_patterns[1].set_insert_loc(1);

	for(i=2; i<12; i++)
		the_patterns[i].set_insert_loc(1);

	for(i=12; i<20; i++)
		the_patterns[i].set_insert_loc(2);

	for(i=20; i<24; i++)
		the_patterns[i].set_insert_loc(1);

	for(i=24; i<26; i++)
		the_patterns[i].set_insert_loc(2);


	//Where is the first piece?
	for(i=0; i<20; i++)
		the_patterns[i].set_piece_num(0, 0);

	for(i=20; i<26; i++)
		the_patterns[i].set_piece_num(0,1);

	//Where are the second and third pieces?
	for(i=2; i<8; i++)
		the_patterns[i].set_piece_num(1,1);

	the_patterns[6].set_piece_num(2,2);
	the_patterns[7].set_piece_num(2,2);

	the_patterns[8].set_piece_num(1,2);
	the_patterns[9].set_piece_num(1,2);
	the_patterns[8].set_piece_num(2,1);
	the_patterns[9].set_piece_num(2,1);


	for(i=10; i<16; i++) {
		the_patterns[i].set_piece_num(1,1);
		the_patterns[i].set_piece_num(2,2);
	}
	
	for(i=16; i<20; i++) {
		the_patterns[i].set_piece_num(1,2);
		the_patterns[i].set_piece_num(2,1);
	}
	
	for(i=20; i<26; i++) {
		the_patterns[i].set_piece_num(1,0);
		the_patterns[i].set_piece_num(2,2);
	}


	//Set reversals
	the_patterns[3].set_reversal(1, 1);
	the_patterns[4].set_reversal(0, 1);
	the_patterns[5].set_reversal(0, 1);
	the_patterns[5].set_reversal(1, 1);
	
	the_patterns[6].set_reversal(1, 1);
	the_patterns[7].set_reversal(1, 1);
	the_patterns[7].set_reversal(2, 1);

	the_patterns[9].set_reversal(2, 1);
	
	the_patterns[10].set_reversal(0, 1);
	the_patterns[10].set_reversal(1, 1);
	the_patterns[11].set_reversal(0, 1);
	the_patterns[11].set_reversal(1, 1);
	the_patterns[11].set_reversal(2, 1);
	
	the_patterns[12].set_reversal(1, 1);
	the_patterns[13].set_reversal(1, 1);
	the_patterns[13].set_reversal(2, 1);
	
	
	the_patterns[14].set_reversal(0, 1);
	the_patterns[14].set_reversal(1, 1);
	the_patterns[15].set_reversal(0, 1);
	the_patterns[15].set_reversal(1, 1);
	the_patterns[15].set_reversal(2, 1);
	
	the_patterns[17].set_reversal(0, 1);
	the_patterns[18].set_reversal(1, 1);
	the_patterns[19].set_reversal(0, 1);
	the_patterns[19].set_reversal(1, 1);

	the_patterns[21].set_reversal(2, 1);
	
	the_patterns[22].set_reversal(1, 1);
	the_patterns[23].set_reversal(1, 1);
	the_patterns[23].set_reversal(2, 1);
	
	the_patterns[25].set_reversal(0, 1);	
		
	//Set BOOLs which tell us what to recalculate
	the_patterns[2].set_next_same(0);
	the_patterns[3].set_last_same(0);

	the_patterns[4].set_next_same(0);
	the_patterns[5].set_last_same(0);

	the_patterns[6].set_next_same(0);
	the_patterns[6].set_next_same(1);

	the_patterns[7].set_next_same(0);
	the_patterns[7].set_last_same(0);
	the_patterns[7].set_last_same(1);

	the_patterns[8].set_last_same(0);
	the_patterns[8].set_next_same(0);
	the_patterns[8].set_next_same(1);

	the_patterns[9].set_last_same(0);
	the_patterns[9].set_last_same(1);

	the_patterns[10].set_next_same(0);
	the_patterns[10].set_next_same(1);

	the_patterns[11].set_last_same(0);
	the_patterns[11].set_last_same(1);

	the_patterns[12].set_next_same(0);
	the_patterns[12].set_next_same(1);

	the_patterns[13].set_last_same(0);
	the_patterns[13].set_last_same(1);

	the_patterns[14].set_next_same(0);
	the_patterns[14].set_next_same(0);

	the_patterns[15].set_last_same(0);
	the_patterns[15].set_last_same(1);

	the_patterns[20].set_next_same(0);
	the_patterns[20].set_next_same(1);

	the_patterns[21].set_last_same(0);
	the_patterns[21].set_last_same(1);
	the_patterns[21].set_next_same(0);

	the_patterns[22].set_last_same(0);
	the_patterns[22].set_next_same(0);
	the_patterns[22].set_next_same(1);

	the_patterns[23].set_last_same(0);
	the_patterns[23].set_last_same(1);

}
	

All_Reorg_Patterns::All_Reorg_Patterns(BOOL noadds)
{
	int i, j;

	num_patterns=20;
	the_patterns=new Reorg_Pattern [num_patterns];

	//Null all reversals so we only have to set real reversals
	for (i=0; i<num_patterns; i++) {
		for(j=0; j<3; j++)
			the_patterns[i].set_reversal(j, 0);
	}

	//How many sections?


	for(i=0; i<3; i++)
		the_patterns[i].set_num_pieces(2);

	for(i=3; i<num_patterns; i++)
		the_patterns[i].set_num_pieces(3);


	//Where is the first piece?
	for(i=0; i<14; i++)
		the_patterns[i].set_piece_num(0, 0);

	for(i=14; i<num_patterns; i++)
		the_patterns[i].set_piece_num(0,1);

	//Where are the second and third pieces?
	for(i=0; i<5; i++) {
		the_patterns[i].set_piece_num(1,1);
		the_patterns[i].set_piece_num(2,2);
	}

	for(i=5; i<8; i++) {
		the_patterns[i].set_piece_num(1,2);
		the_patterns[i].set_piece_num(2,1);
	}

	for(i=8; i<11; i++) {
		the_patterns[i].set_piece_num(1,1);
		the_patterns[i].set_piece_num(2,2);
	}
	
	for(i=11; i<14; i++) {
		the_patterns[i].set_piece_num(1,2);
		the_patterns[i].set_piece_num(2,1);
	}

	for(i=14; i<20; i++) {
		the_patterns[i].set_piece_num(1,0);
		the_patterns[i].set_piece_num(2,2);
	}
	

	//Set reversals
	the_patterns[0].set_reversal(1, 1);
	the_patterns[1].set_reversal(0, 1);
	the_patterns[2].set_reversal(0, 1);
	the_patterns[2].set_reversal(1, 1);
	
	the_patterns[3].set_reversal(1, 1);
	the_patterns[4].set_reversal(1, 1);
	the_patterns[4].set_reversal(2, 1);

	the_patterns[6].set_reversal(2, 1);
	the_patterns[7].set_reversal(1, 1);
	
	the_patterns[8].set_reversal(0, 1);
	the_patterns[8].set_reversal(1, 1);
	the_patterns[9].set_reversal(0, 1);
	the_patterns[9].set_reversal(1, 1);
	the_patterns[9].set_reversal(2, 1);
	
	the_patterns[10].set_reversal(0, 1);
	the_patterns[10].set_reversal(2, 1);
	
	
	the_patterns[11].set_reversal(0, 1);
	the_patterns[11].set_reversal(1, 1);
	the_patterns[12].set_reversal(0, 1);
	the_patterns[13].set_reversal(0, 1);
	the_patterns[13].set_reversal(2, 1);
	
	the_patterns[15].set_reversal(2, 1);
	the_patterns[16].set_reversal(1, 1);
	the_patterns[17].set_reversal(1, 1);
	the_patterns[17].set_reversal(2, 1);

	the_patterns[18].set_reversal(0, 1);
	
	the_patterns[19].set_reversal(0, 1);
	the_patterns[19].set_reversal(2, 1);

		
	//Set BOOLs which tell us what to recalculate
	the_patterns[3].set_next_same(0);
	the_patterns[4].set_last_same(0);

	the_patterns[5].set_next_same(0);
	the_patterns[6].set_last_same(0);

	the_patterns[8].set_next_same(0);
	the_patterns[9].set_last_same(0);

	the_patterns[12].set_next_same(0);
	the_patterns[13].set_last_same(0);

	the_patterns[14].set_next_same(0);
	the_patterns[15].set_last_same(0);

	the_patterns[16].set_next_same(0);
	the_patterns[17].set_last_same(0);

	the_patterns[18].set_next_same(0);
	the_patterns[19].set_last_same(0);

}

Reorg_Pattern & All_Reorg_Patterns::operator[] (int element)
{
	return(the_patterns[element]);
}

All_Reorg_Patterns::~All_Reorg_Patterns()
{
	delete[] the_patterns;
}

Tracking_List::Tracking_List()
{
	cerr<<"Error: Call to default constructor of class Tracking_List\n";
	last=next=0;
	element = new Gene_Track_List*[2];
	element[0]=element[1]=0;
}
	

Tracking_List::Tracking_List(Gene_Track_List **ele, Tracking_List *lst, int n)
{	
	Tracking_List *old_next=0;

	element = new Gene_Track_List*[2];

	element[0]=ele[0];
	element[1]=ele[1];

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


Track_List_List::Track_List_List()
{
	last=next=0;
	element=0;
}

Track_List_List::Track_List_List(Gene_Track_List *new_ele, Track_List_List *old_top)
{
	next=old_top;
	last=0;
	element=new_ele;
}

Tracking_List::~Tracking_List()
{
	delete[] element;
}


Track_Stack::Track_Stack()
{
	top=bottom=full_top=0;
	orig_top=orig_bottom=0;
	size=full_size=orig_size=0;
}



Gene_Track_List * Track_Stack::get_bottom()
{
	return(bottom->element);
}



Gene_Track_List * Track_Stack::get_next()
{
	Gene_Track_List *retval=0;

	if ((size > 0) && (pos_num < size)) {
		retval=curr_pos->element;
		curr_pos=curr_pos->last;
		pos_num++;
	}
	return(retval);
}



void Track_Stack::delete_element(Gene_Track_List *element)
{
	Track_List_List *temp, *temp2;
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
		
		temp->last=orig_top;

		if(orig_top != 0) {
			orig_top->last=temp;
			temp->next=orig_top;
			orig_top=temp;
		}
		else {
			orig_top=temp;
			orig_top->last=orig_top->next=0;
			orig_bottom=orig_top;
		}
		orig_top->last=0;
		size--;
		
	}

}


void Track_Stack::restore_orig()
{
	if (orig_size > size) {
		size=orig_size;
	
		if (top != 0) {
			orig_bottom->next=top;
			top->last=orig_bottom;
			top=orig_top;
		}
		else {
			bottom=orig_bottom;
			top=orig_top;
		}
	}
	orig_top=orig_bottom=0;
	reset();
}
	




void Track_Stack::pop()
{
	Track_List_List *temp;

	temp=bottom->last;

	if (orig_top != 0) {
		orig_top->last=bottom;
		bottom->next=orig_top;
		orig_top=bottom;
	}
	else {
		orig_top=bottom;
		orig_bottom=orig_top;
		
	}	
	
	orig_top->last=0;

	bottom=temp;
	if (size > 1) 
		bottom->next=0;
	else {
		top=0;
		bottom=0;
	}
	size--;
}



void Track_Stack::push(Gene_Track_List *entry)
{
	Track_List_List *new_item;

	if (full_size == orig_size) {
		new_item=new Track_List_List(entry, top);
		full_size++;
	}
	else {
		new_item=full_top;
		full_top=full_top->next;

		new_item->next=new_item->last=0;
		
		if(full_top != 0)
			full_top->last=0;

		new_item->element=entry;
	}


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
		if (full_size == 1) 
			full_top=0;
		pos_num=0;
	}
	
	size++;
	orig_size++;
}


void Track_Stack::add_back(Gene_Track_List *entry)
{
	Track_List_List *temp;

	temp=orig_top;

	while((temp != 0 ) && (temp->element != entry)) temp=temp->next;

	if (temp != 0) {
		if (temp == orig_top) {
			orig_top=orig_top->next;
			orig_top->last=0;
			if (temp == orig_bottom)
				orig_bottom=0;
		}
		else if (temp == orig_bottom) {
			orig_bottom=orig_bottom->last;
			orig_bottom->next=0;
		}
		else {
			temp->last->next=temp->next;
			temp->next->last=temp->last;
		}
		
		if (top != 0) {
			top->last=temp;
			temp->next=top;
			top=top->last;
		}
		else {
			top=bottom=temp;
			top->next=top->last=0;
		}
		size++;

	}

}


void Track_Stack::clear()
{
	restore_orig();
	if (size > 0){
		if (full_top != 0) {
			full_top->last=bottom;
			bottom->next=full_top;
		}
		full_top=top;
	}

	top=bottom=orig_top=orig_bottom=0;
	size=orig_size=0;

}




Track_Stack::~Track_Stack()
{
	int i;
	Track_List_List *temp;

	for(i=0; i<size; i++)
	{
		temp=bottom->last;
		delete bottom;
		bottom=temp;
	}



	while(full_top != 0) {
		temp=full_top->last;
		delete full_top;
		full_top=temp;
	}

}


Genome_Track::Genome_Track()
{
	cerr<<"Error: Call to default constructor of class Genome_Track\n"; 
	hold_reversal_list=0;
	reversal_pointers=0;
}

Genome_Track::Genome_Track(int id, Genome *genome, WGD_Data *homologs)
{
	int i, j;
	Gene_Track_List *ele[2];

	taxa_id=id;

	the_genome=genome;
	the_homologs=homologs;

	reorg_patterns=new All_Reorg_Patterns(TRUE);

	inferred_tracking = new Gene_Track_List* [homologs->get_num_homologs()];
	
	for(i=0; i<homologs->get_num_homologs(); i++) {
		inferred_tracking[i]=new Gene_Track_List [2];
	
	}

	pow2[0]=1;
	for(i=1; i<32; i++)
		pow2[i] =pow2[i-1]*2;

	partial_track_start=partial_track_pos=partial_track_end=0;
	list_len=last_list_pos=0;

	hold_reversal_list =new Tracking_List* [the_homologs->get_num_homologs()];
	reversal_pos=0;

	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		ele[0]=new Gene_Track_List;
		ele[1]=new Gene_Track_List;
		hold_reversal_list[i]=new Tracking_List(ele, 0, i);
	}
	reversal_pointers = new Gene_Track_List * [2*the_homologs->get_num_homologs()];
}


void Genome_Track::number_contigs(int *order)
{
	number_contigs(order, the_homologs->get_num_homologs());
}



void Genome_Track::number_contigs(int *order, int size)
{
	int i, j, contig1, contig2, temp;
	Track_Stack *lefts, *rights;
	
	for(j=0; j<size; j++) {
		inferred_tracking[j][0].my_locus=&(*the_homologs)[order[j]][taxa_id];
		inferred_tracking[j][0].index_num=0;

		if ((*the_homologs)[order[j]][taxa_id].has_duplicate() == TRUE) {
			inferred_tracking[j][1].my_locus=&(*the_homologs)[order[j]][taxa_id];
			inferred_tracking[j][1].index_num=1;
		}
		else 	inferred_tracking[j][1].my_locus=0;

		inferred_tracking[j][0].last=inferred_tracking[j][0].next=0;
		inferred_tracking[j][1].last=inferred_tracking[j][1].next=0;
		inferred_tracking[j][0].dist_to_last=0;
		inferred_tracking[j][1].dist_to_last=0;
			
	}

	recurse_assemble(lefts, rights, 0, size);
	delete lefts;
	delete rights;


	//ID each locus
	for(j=0; j<size; j++) {
		inferred_tracking[j][0].partition_num=j;
		inferred_tracking[j][1].partition_num=j;

	}

	//Recheck all breaks to see if we can gain any breaks
	for(j=1; j<size; j++) 
		break_and_rejoin(j, size);
	

	//Puts the above list into an easy-to-understand 2Xn table
	store_tracking(size);


}

void Genome_Track::print_tracking()
{
	print_tracking(TRUE);
}


void Genome_Track::print_tracking(BOOL use_file)
{
	int i, j;
	char filename[100];
	ofstream fout;
	BOOL link;

	if (use_file == TRUE) {
		strcpy(filename, the_genome->get_name());
		strcat(filename, "_tracking.txt");
		fout.open(filename);
	}

	for(j=0; j<2; j++) {
		for(i=0; i<the_homologs->get_num_homologs(); i++) {
			if (inferred_tracking[i][j].my_locus != 0) {
				if (use_file == TRUE) 
					fout<<(*the_genome)[inferred_tracking[i][j].my_locus->get_contig(
							inferred_tracking[i][j].index_num)][inferred_tracking[i][j].my_locus->get_gene(
							inferred_tracking[i][j].index_num)].get_name();
				else
					cout<<(*the_genome)[inferred_tracking[i][j].my_locus->get_contig(
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
					fout<<"NONE\t*\t";
				else
					cout<<"NONE\t*\t";
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




Gene_Track_List* Genome_Track::get_gene_track(int locus_num, int track_num)
{
	if ((locus_num >=0) && (locus_num<the_homologs->get_num_homologs())) {
		if (track_num == 0)
			return(&inferred_tracking[locus_num][0]);
		else
			return(&inferred_tracking[locus_num][1]);
	}
	else {
		cerr<<"Request for invalid index "<<locus_num<<" in Genome_Track\n";
 		return(0);
	}
}
	


BOOL Genome_Track::has_back_link(int locus_num, int track_num)
{
	if ((locus_num >=1) && (locus_num<the_homologs->get_num_homologs())) {
		if (track_num == 0) {
			if (inferred_tracking[locus_num][0].last != 0)
				return(TRUE);
			else
				return(FALSE);
		}
		else {
			if (inferred_tracking[locus_num][1].last != 0)
				return(TRUE);
			else
				return(FALSE);

		}
	}
	else return(FALSE);

}


BOOL Genome_Track::has_double_break(int n)
{
	BOOL stop, retval=TRUE;
	Gene_Track_List *to_end, *from_previous, *to_end_2nd;
	Tracking_List *curr_pos, *end;


	stop=FALSE;

		//Find all the possible left endpoints
	curr_pos=get_list_pos(n);
	end=curr_pos;
	while((stop==FALSE) && (retval == TRUE)) {
		to_end=curr_pos->element[0]->next;
			
		if (curr_pos != end)
			from_previous=curr_pos->next->element[0]->last;
		else
			from_previous=0;

		if (curr_pos->element[1]->my_locus != 0)
			to_end_2nd=curr_pos->element[1]->next;
		else
			to_end_2nd = 0;
	

		if (to_end == 0) {
			if (curr_pos != end) {
				//If the previous entry has a join but not to this entry, then
				//this one must be on the other track.  Hence, stop
				if (from_previous !=0)
					stop=TRUE;
			}
		}
		else if (to_end->partition_num > end->element[0]->partition_num)
			retval=FALSE;

		if ((curr_pos->element[1]->my_locus != 0) && (to_end_2nd ==0)) {
			//This is a left end and a stopping point
			stop=TRUE;
		}
		else if ((to_end_2nd != 0) && (to_end_2nd->partition_num > end->element[0]->partition_num))
				retval=FALSE;
		

		if (curr_pos->last == 0)
			stop=TRUE;

		curr_pos=curr_pos->last;
	}
	
	return(retval);
}

BOOL Genome_Track::has_tracked_double_break(int n)
{
	int i;
	BOOL stop, retval=TRUE;
	Gene_Track_List *to_end, *from_previous, *to_end_2nd;

	i=n-1;
	stop=FALSE;

		//Find all the possible left endpoints

	while((stop==FALSE) && (retval == TRUE)) {
		to_end=inferred_tracking[i][0].next;
			
		if (i != n-1)
			from_previous=inferred_tracking[i+1][0].last;
		else
			from_previous=0;

		if (inferred_tracking[i][1].my_locus != 0)
			to_end_2nd=inferred_tracking[i][1].next;
		else
			to_end_2nd = 0;
	

		if (to_end == 0) {
			if (i != n-1) {
				//If the previous entry has a join but not to this entry, then
				//this one must be on the other track.  Hence, stop
				if (from_previous !=0)
					stop=TRUE;
			}
		}
		else if (to_end->partition_num > inferred_tracking[n-1][0].partition_num)
			retval=FALSE;

		if ((inferred_tracking[i][1].my_locus != 0) && (to_end_2nd ==0)) {
			//This is a left end and a stopping point
			stop=TRUE;
		}
		else if ((to_end_2nd != 0) && (to_end_2nd->partition_num > inferred_tracking[n-1][0].partition_num))
				retval=FALSE;
		

		if (i == 0)
			stop=TRUE;

		i--;
	}
	
	return(retval);
}



int Genome_Track::get_dist_to_last(int locus_num, int track_num)
{
	if (has_back_link(locus_num, track_num) == FALSE)
		//Error condition
		return(-1);
	else
		return(inferred_tracking[locus_num][track_num].dist_to_last);
}



int Genome_Track::count_num_breaks()
{
	return(count_num_breaks(the_homologs->get_num_homologs()));
}

int Genome_Track::count_num_breaks(int size)
{
	int i, retval=0;

	for(i=0; i<size; i++) {
			if ((inferred_tracking[i][0].my_locus !=0) && (inferred_tracking[i][0].last == 0))
				retval++;
			if ((inferred_tracking[i][1].my_locus !=0) && (inferred_tracking[i][1].last == 0))
				retval++;


	}
	return(retval);

}


int Genome_Track::count_num_list_track_breaks()
{
	int retval=0;
	Tracking_List *curr_pos;

	curr_pos=partial_track_start;

	while(curr_pos != 0) {
		if (curr_pos->element[0]->next == 0)
			retval++;

		if ((curr_pos->element[1]->my_locus != 0) &&(curr_pos->element[1]->next == 0))
			retval++;

		curr_pos=curr_pos->next;

	}

	return(retval);
	
}

void Genome_Track::start_track_list(int n)
{
	Gene_Track_List **locus;
	Tracking_List *new_member;
	
	locus=create_new_track_locus(&(*the_homologs)[n][taxa_id]);

	new_member=new Tracking_List(locus, 0, 0);

	partial_track_start=partial_track_pos=partial_track_end=new_member;
	last_list_pos=0;
	list_len=1;
}


Tracking_List * Genome_Track::create_list_element(int locus_num)
{
	Gene_Track_List **locus;
	Tracking_List *new_member;

	locus=create_new_track_locus(&(*the_homologs)[locus_num][taxa_id]);
	
	new_member=new Tracking_List(locus, 0, locus_num);

	return(new_member);
}


int Genome_Track::joins_after(int n)
{
	int i, track[2], num_joins=0;

	track[0]=track[1]=n;

	for(i=0; i<2; i++) {
		while((track[i] >= 0)  && (inferred_tracking[track[i]][i].my_locus == 0) ) {
			track[i]--;
		}
		if ((track[i]>=0) && (inferred_tracking[track[i]][i].next != 0))
			num_joins++;
	}

	return(num_joins);
	
}


int Genome_Track::joins_after_list(int n)
{
	//COMPLETELY NON_FUNCTIONAL FUNCTION!
	int num_joins;
	Tracking_List *curr_pos;

	curr_pos=get_list_pos(n);
	if (curr_pos->element[0]->next != 0)
		num_joins++;

	curr_pos=curr_pos->last;
	while((num_joins) && (curr_pos)) {
		num_joins++;

	}
	return(num_joins);

}


BOOL Genome_Track::all_poss_joins_after(int n)
{
	BOOL retval=TRUE;
	if ((inferred_tracking[n][0].my_locus != 0) && (inferred_tracking[n][0].next == 0))
		retval=FALSE;
	if ((inferred_tracking[n][1].my_locus != 0) && (inferred_tracking[n][1].next == 0))
		retval=FALSE;

	return(retval);
}


BOOL Genome_Track::diagnose_list()
{
	int i, j;
	BOOL retval=TRUE;
	Tracking_List *curr_pos;

	curr_pos = partial_track_start;

	while(curr_pos != 0) {
		for(i=0; i<2; i++) {
			if ((curr_pos->element[i]->my_locus != 0) && (curr_pos->element[i]->last != 0)) {
				for(j=0; j<the_homologs->get_num_homologs(); j++)
					if ((curr_pos->element[i]->last == hold_reversal_list[j]->element[0] ) || 
						(curr_pos->element[i]->last == hold_reversal_list[j]->element[1]))
						retval = FALSE;
				if ((curr_pos->element[i]->last->my_locus->get_gene_obj(curr_pos->element[i]->last->index_num)->get_neighbor(0) != 
					curr_pos->element[i]->my_locus->get_gene_obj(curr_pos->element[i]->index_num) ) &&
					(curr_pos->element[i]->last->my_locus->get_gene_obj(curr_pos->element[i]->last->index_num)->get_neighbor(1) != 
					curr_pos->element[i]->my_locus->get_gene_obj(curr_pos->element[i]->index_num)))
					retval=FALSE;
			}

			if ((curr_pos->element[i]->my_locus != 0) && (curr_pos->element[i]->next != 0)) {
				for(j=0; j<the_homologs->get_num_homologs(); j++)
					if ((curr_pos->element[i]->next == hold_reversal_list[j]->element[0] ) ||
						(curr_pos->element[i]->next == hold_reversal_list[j]->element[1]))
						retval = FALSE;
				
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


int Genome_Track::get_nth_tracking_list_pos(int n)
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




int Genome_Track::get_list_position_number(int n)
{
	return(get_list_pos(n)->num);
}




int Genome_Track::add_to_position_n(int n, Tracking_List *new_locus_start, Tracking_List *new_locus_end, int size, BOOL do_link)
{
	int i, num_joins=0, lost_breaks=0, save_index = 0;
	Tracking_List *curr_pos, *temp, *right_pos;
	
	curr_pos=get_list_pos(n);

	//cout<<"n: "<<n<<"\t"<<curr_pos<<endl;
	right_pos=curr_pos->next;

	if (n > 0) {
		temp=curr_pos->next;
		curr_pos->next=new_locus_start;
		new_locus_start->last=curr_pos;

		if (temp != 0) {
			temp->last=new_locus_end;
			new_locus_end->next=temp;
		}
	}
	else {
		partial_track_start->last=new_locus_end;
		new_locus_end->next=partial_track_start;
	}


	if ((n != 0) && (n<list_len)) {
		section_mins[0]=partial_track_start->element[0]->partition_num;
		section_maxs[0]=curr_pos->element[0]->partition_num;
		lost_breaks = find_partition_breaks(curr_pos->next, partial_track_start, 1, section_mins, section_maxs);
	}

	if(n!=0) {
		set_stacks(&new_left, curr_pos, partial_track_start, TRUE);
		set_stacks(&new_middle_left, new_locus_start, new_locus_end, FALSE);
		num_joins = make_joins(&new_left, &new_middle_left, curr_pos->element[0], new_locus_start->element[0], TRUE, save_index);
		new_left.clear();
		new_middle_left.clear();
		
		if (n != list_len) {
			set_stacks(&new_middle_right, new_locus_end, partial_track_start, TRUE);
			set_stacks(&new_right, right_pos, partial_track_end, FALSE);
			num_joins += make_joins(&new_middle_right, &new_right, new_locus_end->element[0], right_pos->element[0], do_link);

			new_middle_right.clear();
			new_right.clear();
		}

	}
	else {
		set_stacks(&new_middle_left, new_locus_end, new_locus_start, TRUE);
		set_stacks(&new_left, partial_track_start, partial_track_end, FALSE);
		num_joins += make_joins(&new_middle_left, &new_left, new_locus_end->element[0], partial_track_start->element[0], do_link);

		new_middle_left.clear();
		new_left.clear();
	}
	

	

	if ( do_link == FALSE) {
		//Fix the mess we just made
		if (n > 0) {
			curr_pos->next=new_locus_end->next;
			if (new_locus_end->next != 0 )
				new_locus_end->next->last=curr_pos;
		}
		else
			partial_track_start->last=0;

		for(i=0; i<2; i++) 
			new_locus_start->element[i]->last=new_locus_end->element[i]->next=0;
		new_locus_end->next=new_locus_start->last=0;

		for(i=0; i<save_index; i++) {
			new_lasts[i]->next=0;
			new_nexts[i]->last=0;
			
		}

		for(i=0; i<lost_breaks; i++) {
				broken_lefts[i]->next=broken_rights[i];
				broken_rights[i]->last=broken_lefts[i];
		}

		if (n != 0) {
			last_list_pos++;
			partial_track_pos=partial_track_pos->next;
		}
	}
	else {
		if (n ==0)
			partial_track_start=new_locus_start;
		else if (n == list_len)
			partial_track_end=new_locus_end;
		list_len+=size;

		partial_track_pos=partial_track_start;
		last_list_pos=0;
		reset_list_partition_numbers();
	}

	return(num_joins-lost_breaks);
}



void Genome_Track::try_position_reorgs(int n, int **scores, BOOL *fully_connected)
{
		//Note that we never try n=list len, since this possibility is implicit in n=0
	int i, j, k, lost_breaks, new_lost_breaks, pos[3], rev[3], section_pos[3], 
		save_index=0, last_score, last_score2, old_rev_pos, min, 
		num_check_breaks[2][3], break_sec_num[2][3], sec_cnt, num_sections;
	BOOL do_2nd_rejoin;
	Tracking_List *curr_pos, *curr_pos2, *find_pos, *dum;
	
	curr_pos=get_list_pos(n);
	dum=get_list_pos(6);

	section_mins[0]=partial_track_start->element[0]->partition_num;
	section_maxs[0]=curr_pos->element[0]->partition_num;
	lost_breaks=find_partition_breaks(curr_pos->next, partial_track_start, 1 , section_mins, section_maxs);
	
	num_old=lost_breaks;

	for(i=0; i<lost_breaks; i++) {
		old_lasts[i]=broken_lefts[i];
		old_nexts[i]=broken_rights[i];
	}
	
	
	find_possible_internal_rejoins(lost_breaks, 0, 0, break_sec_num, num_check_breaks, 
		curr_pos);

	section_starts[0][0]=partial_track_start;
	section_ends[0][0]=curr_pos;

	section_starts[1][0]=curr_pos->next;
	section_ends[1][0]=partial_track_end;

	reverse_section(partial_track_start, curr_pos, section_starts[0][1], section_ends[0][1]);
	
	old_rev_pos=reversal_pos;

	reverse_section(curr_pos->next, partial_track_end, section_starts[1][1], section_ends[1][1]);

	find_reversal_internal_rejoins(0, num_check_breaks, section_ends[0][1], section_starts[1][1]);

	//Only need to initialize the stacks once
	for(j=0; j<2; j++) {
		for(k=0; k<2; k++) {
			set_stacks(&section_stacks[j][k][0], section_ends[j][k], section_starts[j][k], TRUE);
			set_stacks(&section_stacks[j][k][1], section_starts[j][k], section_ends[j][k], FALSE);
		}
	}
		

	for(j=0; j<3; j++) {
		for(k=0; k<2; k++) {
				section_pos[(*reorg_patterns)[j].get_piece_num(k)]=k;
				pos[k]=(*reorg_patterns)[j].get_piece_num(k);
				rev[k]=(*reorg_patterns)[j].get_reversal(k);
		}

				
		save_index=0;

		section_starts[pos[0]][rev[0]]->last=0;
		section_ends[pos[1]][rev[1]]->next=0;


		section_ends[pos[0]][rev[0]]->next=section_starts[pos[1]][rev[1]];
		section_starts[pos[1]][rev[1]]->last=section_ends[pos[0]][rev[0]];
			
		
		scores[0][j]+=make_joins(&section_stacks[pos[0]][rev[0]][0], &section_stacks[pos[1]][rev[1]][1], 
			section_ends[pos[0]][rev[0]]->element[0], section_starts[pos[1]][rev[1]]->element[0], TRUE, save_index)-lost_breaks;

		if (num_check_breaks[0][0] != 0)  {
			num_sections = set_section_min_maxs(j, rev, section_pos, pos, 0, 0, break_sec_num);
			scores[0][j]+=break_and_rejoin(check_break_pos[0][rev[section_pos[break_sec_num[0][0]]]][0], 
				section_starts[pos[0]][rev[0]],	section_ends[pos[1]][rev[1]], save_index, 0, num_sections, section_mins, section_maxs);	
		}

		if (num_check_breaks[1][0] !=0 ) {
			num_sections = set_section_min_maxs(j, rev, section_pos, pos, 0, 1, break_sec_num);
			scores[0][j]+=break_and_rejoin(check_break_pos[1][rev[section_pos[break_sec_num[1][0]]]][0], section_starts[pos[0]][rev[0]],
				section_ends[pos[1]][rev[1]], save_index, 0, num_sections, section_mins, section_maxs);	
		}

		for(k=num_old-1; k>=lost_breaks; k--) {
				old_lasts[k]->next=old_nexts[k];
				old_nexts[k]->last=old_lasts[k];
		}
	
		if (save_index != 0) {
			for(k=0; k<save_index; k++) {
				if (new_lasts[k]->next == new_nexts[k]) 
					new_lasts[k]->next=0;
				if (new_nexts[k]->last == new_lasts[k]) 
					new_nexts[k]->last=0;
			}
		}		
			
		num_old=lost_breaks;
			
	}

	section_ends[0][0]->next=section_starts[1][0];
	section_starts[1][0]->last=section_ends[0][0];


	for(i=0; i<2; i++) {
		section_stacks[1][i][0].clear();
		section_stacks[1][i][1].clear();
	}

	section_ends[2][0]=partial_track_end;

	for(i=n+1; i<list_len; i++) {
		curr_pos2=get_list_pos(i);

		section_ends[1][0]=curr_pos2;
		section_starts[2][0]=curr_pos2->next;

		if (fully_connected[i] == FALSE) {
		

			if (num_check_breaks[1][0] != 0) {
				find_pos=section_starts[1][0];
				break_sec_num[1][0]=2;
				while(find_pos != section_ends[1][0]->next) {
					if (find_pos == check_break_pos[1][0][0]->next)
						break_sec_num[1][0]=1;
					find_pos=find_pos->next;
				}

				if ((find_pos == section_ends[1][0]) || (find_pos == section_starts[2][0]))
					do_2nd_rejoin=FALSE;
				else
					do_2nd_rejoin=TRUE;
			
			}
			
			reversal_pos=old_rev_pos;
		
			for(j=0; j<2; j++) {
				section_mins[j]=section_starts[j][0]->element[0]->partition_num;
				section_maxs[j]=section_ends[j][0]->element[0]->partition_num;	
			}
			new_lost_breaks=find_partition_breaks(section_starts[2][0], partial_track_start, 2, section_mins, section_maxs);
		
			num_old=lost_breaks+new_lost_breaks;

			for(j=0; j<new_lost_breaks; j++) {
				old_lasts[j+lost_breaks]=broken_lefts[j];
				old_nexts[j+lost_breaks]=broken_rights[j];
			}

			find_possible_internal_rejoins(new_lost_breaks, lost_breaks, 1, break_sec_num, num_check_breaks, 
				curr_pos2);

		

			reverse_section(section_starts[1][0], section_ends[1][0], section_starts[1][1], section_ends[1][1]);
			reverse_section(section_starts[2][0], section_ends[2][0], section_starts[2][1], section_ends[2][1]);

			section_ends[2][1]->next=section_starts[1][1];
			section_starts[1][1]->last=section_ends[2][1];
			section_ends[1][1]->next=section_starts[0][1];
			section_starts[0][1]->last=section_ends[1][1];
			section_starts[2][1]->last=0;
			section_ends[0][1]->next=0;

			find_reversal_internal_rejoins(1, num_check_breaks, section_ends[1][1], section_starts[2][1]);

			/*	section_starts[0][1]->next=section_starts[1][1]->next=section_starts[2][1]->next=0;
			section_ends[0][1]->next=section_ends[1][1]->next=section_ends[2][1]->next=0;
			section_starts[0][1]->last=section_starts[1][1]->last=section_starts[2][1]->last=0;
			section_ends[0][1]->last=section_ends[1][1]->last=section_ends[2][1]->last=0;*/


			for(j=1; j<3; j++) {
				for(k=0; k<2; k++) {
					set_stacks(&section_stacks[j][k][0], section_ends[j][k], section_starts[j][k], TRUE);
					set_stacks(&section_stacks[j][k][1], section_starts[j][k], section_ends[j][k], FALSE);
				}
			}


			for (j=3; j<reorg_patterns->get_num_patterns(); j++) {
				for(k=0; k<3; k++) {
					pos[k]=(*reorg_patterns)[j].get_piece_num(k);
					rev[k]=(*reorg_patterns)[j].get_reversal(k);
					section_pos[(*reorg_patterns)[j].get_piece_num(k)]=k;
				}

				save_index=0;
						
				section_starts[pos[0]][rev[0]]->last=0;
				section_ends[pos[2]][rev[2]]->next=0;
						
				//Make 2 joins 
				section_ends[pos[0]][rev[0]]->next = section_starts[pos[1]][rev[1]];
				section_starts[pos[1]][rev[1]]->last = section_ends[pos[0]][rev[0]];
						
				section_ends[pos[1]][rev[1]]->next = section_starts[pos[2]][rev[2]];
				section_starts[pos[2]][rev[2]]->last = section_ends[pos[1]][rev[1]];

				scores[i][j]+=make_joins(&section_stacks[pos[0]][rev[0]][0], &section_stacks[pos[1]][rev[1]][1], 
						section_ends[pos[0]][rev[0]]->element[0], section_starts[pos[1]][rev[1]]->element[0], TRUE, save_index)-lost_breaks-new_lost_breaks;
						
						
				set_stacks(&new_left, section_ends[pos[1]][rev[1]], section_starts[pos[0]][rev[0]], TRUE);			
				scores[i][j]+=make_joins(&new_left, &section_stacks[pos[2]][rev[2]][1], section_ends[pos[1]][rev[1]]->element[0], 
					section_starts[pos[2]][rev[2]]->element[0], TRUE, save_index);
				new_left.clear();


				if (rev[0] == 0) {
					section_mins[0]=section_starts[pos[0]][rev[0]]->element[0]->partition_num;
					section_maxs[0]=section_ends[pos[0]][rev[0]]->element[0]->partition_num;
				}
				else {
					section_maxs[0]=section_starts[pos[0]][rev[0]]->element[0]->partition_num;
					section_mins[0]=section_ends[pos[0]][rev[0]]->element[0]->partition_num;
				}

				scores[i][j]+=break_and_rejoin(section_ends[pos[0]][rev[0]], section_starts[pos[0]][rev[0]],
					section_ends[pos[2]][rev[2]], save_index, 0, 1, section_mins, section_maxs);	
				
			
				if (num_check_breaks[0][0] != 0) {
					num_sections = set_section_min_maxs(j, rev, section_pos, pos, 0, 0, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[0][rev[section_pos[break_sec_num[0][0]]]][0], section_starts[pos[0]][rev[0]],
						section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);
				}
				if ((num_check_breaks[1][0] != 0) && (do_2nd_rejoin == TRUE)) {
					num_sections = set_section_min_maxs(j, rev, section_pos, pos, 0, 1, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[1][rev[section_pos[break_sec_num[1][0]]]][0], section_starts[pos[0]][rev[0]],
						section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);
				}

				

				if (num_check_breaks[0][1] != 0) {
					num_sections = set_section_min_maxs(j, rev, section_pos, pos, 1, 0, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[0][rev[section_pos[break_sec_num[0][1]]]][1], section_starts[pos[0]][rev[0]],
						section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);
				}

				if (num_check_breaks[1][1] != 0) {
					num_sections = set_section_min_maxs(j, rev, section_pos, pos, 1, 1, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[1][rev[section_pos[break_sec_num[1][1]]]][1], section_starts[pos[0]][rev[0]],
						section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);
				}


				for (k=lost_breaks+new_lost_breaks; k<num_old; k++) {
					old_lasts[k]->next=old_nexts[k];
					old_nexts[k]->last=old_lasts[k];
				}
				num_old=lost_breaks+new_lost_breaks;

				if (save_index != 0) {
					for(k=0; k<save_index; k++) {
						if (new_lasts[k]->next == new_nexts[k]) 
							new_lasts[k]->next=0;
						if (new_nexts[k]->last == new_lasts[k])
							new_nexts[k]->last=0;
						
					}
				}	
			
				
				
			}		
			
				
			//Reset joins
			section_ends[0][0]->next=section_starts[1][0];
			section_starts[1][0]->last=section_ends[0][0];
			section_ends[1][0]->next=section_starts[2][0];
			section_starts[2][0]->last=section_ends[1][0];

			section_ends[2][1]->next=section_starts[1][1];
			section_starts[1][1]->last=section_ends[2][1];
			section_ends[1][1]->next=section_starts[0][1];
			section_starts[0][1]->last=section_ends[1][1];

			for(j=0; j<new_lost_breaks; j++) {
				old_lasts[j+lost_breaks]->next=old_nexts[j+lost_breaks];
				old_nexts[j+lost_breaks]->last=old_lasts[j+lost_breaks];
			}
			num_old -= new_lost_breaks;

			for(j=1; j<3; j++) {
				for(k=0; k<2; k++) {
					section_stacks[j][k][0].clear();
					section_stacks[j][k][1].clear();
				}
			}
		}  //if (fully_connected[i] == FALSE
	} //for (i=n+1...

	//Reform the list as before
	section_starts[0][0]->last=0;
	section_ends[2][0]->next=0;
	section_ends[0][0]->next=section_starts[1][0];
	section_starts[1][0]->last=section_ends[0][0];

	for(j=0; j<3; j++) {
		for(k=0; k<2; k++) {
			section_stacks[j][k][0].clear();
			section_stacks[j][k][1].clear();
		}
	}

	for(i=0; i<lost_breaks; i++) {
		old_lasts[i]->next=old_nexts[i];
		old_nexts[i]->last=old_lasts[i];
	}

	if (n+1 < list_len) {

		//There are three sections
		section_ends[1][0]->next=section_starts[2][0];
		section_starts[2][0]->last=section_ends[1][0];
	
		for(j=0; j<new_lost_breaks; j++) {
			old_lasts[j+lost_breaks]->next=old_nexts[j+lost_breaks];
			old_nexts[j+lost_breaks]->last=old_lasts[j+lost_breaks];
		}
	}
	


	reversal_pos=0;
	partial_track_pos=partial_track_start;
	last_list_pos=0;
	

}



void Genome_Track::try_position_inserts(int start_pos, int end_pos, int **scores, BOOL *fully_connected)
{
	int i, j, k, rejoins, lost_breaks[3], save_index=0, save_index2=2, 
		last_score, last_score2, reused_losts=0, break_sec_num[2][3], num_check_breaks[2][3], num_sections;
	BOOL do_2nd_rejoin;
	Tracking_List *curr_pos, *save_pos, *find_pos, *dum;
	

	dum=get_list_pos(4);
	section_starts[1][0]=get_list_pos(start_pos)->next;
	section_ends[1][0]=get_list_pos(end_pos);
	save_pos=section_starts[1][0]->last;

	precut_max_partition=save_pos->element[0]->partition_num;
	postcut_min_partition=section_ends[1][0]->next->element[0]->partition_num;

	section_mins[0]=partial_track_start->element[0]->partition_num;
	section_maxs[0]=save_pos->element[0]->partition_num;
	lost_breaks[0]=find_partition_breaks(section_starts[1][0], partial_track_start, 1, section_mins, section_maxs);

	for(i=0; i<lost_breaks[0]; i++) {
		old_lasts[i]=broken_lefts[i];
		old_nexts[i]=broken_rights[i];
	}

	num_old=lost_breaks[0];

	find_possible_internal_rejoins(lost_breaks[0], 0, 0, break_sec_num, num_check_breaks, 
		save_pos);


	if (num_check_breaks[1][0] != 0) {
		if (check_break_pos[1][0][0] == section_ends[1][0]) {
			check_break_pos[1][0][0]=save_pos;
			break_sec_num[1][0]=2;
		}
		else {
			find_pos=section_starts[1][0];
			break_sec_num[1][0]=2;
			while(find_pos != section_ends[1][0]->next) {
				if (find_pos == check_break_pos[1][0][0]->next)
					break_sec_num[1][0]=1;
				find_pos=find_pos->next;
			}

			if ((check_break_pos[1][0][0] == section_ends[1][0]) || (check_break_pos[1][0][0] == section_ends[1][0]->next))
				do_2nd_rejoin=FALSE;
			else
				do_2nd_rejoin=TRUE;
		}
		
	}
		
	section_mins[0]=partial_track_start->element[0]->partition_num;
	section_maxs[0]=save_pos->element[0]->partition_num;
	section_mins[1]=section_starts[1][0]->element[0]->partition_num;
	section_maxs[1]=section_ends[1][0]->element[0]->partition_num;
	lost_breaks[1]=find_partition_breaks(section_ends[1][0]->next, section_starts[1][0], 2, section_mins, section_maxs);

	for(i=0; i<lost_breaks[1]; i++) {
		old_lasts[i+num_old]=broken_lefts[i];
		old_nexts[i+num_old]=broken_rights[i];
	}
	
	num_old += lost_breaks[1];

	find_possible_internal_rejoins(lost_breaks[1], lost_breaks[0], 1, break_sec_num, num_check_breaks, 
		section_ends[1][0]);
		
	list_len-=(end_pos-start_pos);

	section_starts[1][0]->last->next=section_ends[1][0]->next;
	
	section_ends[1][0]->next->last=section_starts[1][0]->last;


	section_starts[1][0]->last=0;
	section_ends[1][0]->next=0;
	reverse_section(section_starts[1][0], section_ends[1][0], section_starts[1][1], section_ends[1][1]);

	if ((num_check_breaks[1][0] == 1) && (break_sec_num[1][0] == 1)) 
		find_reversal_internal_rejoins(0, num_check_breaks, 0, section_starts[1][1]);
		
	if (num_check_breaks[0][1] == 1) 
		find_reversal_internal_rejoins(1, num_check_breaks, section_ends[1][1], 0);			


	set_stacks(&section_stacks[1][0][0], section_starts[1][0], section_ends[1][0], FALSE);
	set_stacks(&section_stacks[1][0][1], section_starts[1][1], section_ends[1][1], FALSE);
	

	set_stacks(&new_left, save_pos, partial_track_start, TRUE);
	set_stacks(&new_right, save_pos->next, partial_track_end, FALSE);
	rejoins = make_joins(&new_left, &new_right, save_pos->element[0], save_pos->next->element[0], TRUE, save_index);

	new_right.clear();
	new_left.clear();
	new_middle_left.clear();



	section_starts[0][0]=partial_track_start;
	section_ends[0][0]=partial_track_start;
	section_starts[2][0]=partial_track_start;
	section_ends[2][0]=partial_track_end;
	
	if (num_check_breaks[0][0] != 0) {
			if (check_break_pos[0][0][0] == section_starts[0][0])
				break_sec_num[0][0]=0;
			else
				break_sec_num[0][0]=2;
	}



	for(i=0; i<list_len-1; i++) {
		section_ends[0][0]=section_starts[2][0];
		section_starts[2][0]=section_starts[2][0]->next;

		if ((num_check_breaks[0][0] != 0) && (check_break_pos[0][0][0] == section_ends[0][0]))
			break_sec_num[0][0]=0;

		if ((num_check_breaks[1][0] != 0) &&(break_sec_num[1][0] != 1) && (check_break_pos[1][0][0] == section_ends[0][0])) 
			break_sec_num[1][0]=0;

		if ((num_check_breaks[1][1] != 0) && (check_break_pos[1][0][1] == section_ends[0][0]))
			break_sec_num[1][1]=0;
	
		if (((i<start_pos) && (fully_connected[i] == FALSE)) || (fully_connected[i+(end_pos-start_pos)]== FALSE)) {

			section_mins[0]=section_starts[0][0]->element[0]->partition_num;
			section_maxs[0]=section_ends[0][0]->element[0]->partition_num;
			lost_breaks[2]=find_partition_breaks(section_starts[2][0], section_starts[0][0], 1, section_mins, section_maxs);

			for(j=0; j<lost_breaks[2]; j++){
				old_lasts[j+num_old]=broken_lefts[j];
				old_nexts[j+num_old]=broken_rights[j];
			}

			num_old=lost_breaks[2]+lost_breaks[0]+lost_breaks[1];
			find_possible_internal_rejoins(lost_breaks[2], lost_breaks[0]+lost_breaks[1], 2, break_sec_num, num_check_breaks, 
				section_ends[0][0]);
			
			//Reset the section numbers to account for the clipped section
			if (num_check_breaks[0][2] == 1)
				break_sec_num[0][2] = 0;
			if (num_check_breaks[1][2] == 1)
				break_sec_num[1][2] = 2;

			set_stacks(&section_stacks[0][0][0], section_ends[0][0], section_starts[0][0], TRUE);
			set_stacks(&section_stacks[2][1][0], section_starts[2][0], section_ends[2][0], FALSE);

			
			for(j=0; j<2; j++) {
				save_index2=save_index;
				section_ends[0][0]->next=section_starts[1][j];
				section_starts[1][j]->last=section_ends[0][0];
				section_starts[2][0]->last=section_ends[1][j];
				section_ends[1][j]->next=section_starts[2][0];

			
		
				scores[i][j]+=make_joins(&section_stacks[0][0][0], &section_stacks[1][0][j], 
					section_ends[0][0]->element[0], section_starts[1][j]->element[0], 
					TRUE, save_index2)-lost_breaks[0]-lost_breaks[1]-lost_breaks[2]+rejoins;
				
				set_stacks(&new_middle_left, section_ends[1][j], section_starts[0][0], TRUE);
				scores[i][j]+=make_joins(&new_middle_left, &section_stacks[2][1][0], section_ends[1][j]->element[0],
					section_starts[2][0]->element[0], TRUE, save_index2);
				
				new_middle_left.clear();	

			
				if ((num_check_breaks[0][0] != 0) && 
					(check_break_pos[0][0][0] != section_ends[0][0]) && (check_break_pos[0][0][0] != section_starts[2][0])) 
				{
					num_sections=set_section_min_maxs(0, 0, 0, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[0][0][0], section_starts[0][0],
							section_ends[2][0], save_index2, save_index, num_sections, section_mins, section_maxs);	
				}
				
				if (num_check_breaks[0][1] != 0) {
					num_sections=set_section_min_maxs(1, 0, j, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[0][j][1], section_starts[0][0],
							section_ends[2][0], save_index2, save_index,  num_sections ,section_mins, section_maxs);
				}

				if (num_check_breaks[0][2] != 0) {
					num_sections=set_section_min_maxs(2, 0, 0, break_sec_num);
					scores[i][j]+=break_and_rejoin(check_break_pos[0][0][2], section_starts[0][0],
							section_ends[2][0], save_index2, save_index, num_sections, section_mins, section_maxs);
				}

				if (num_check_breaks[1][0] !=0) {
					if (break_sec_num[1][0] == 1)  {
						num_sections=set_section_min_maxs(0, 1, j, break_sec_num);
						scores[i][j]+=break_and_rejoin(check_break_pos[1][j][0], section_starts[0][0],
							section_ends[2][0], save_index2, save_index, num_sections, section_mins, section_maxs);
					}
					else {
						if ((check_break_pos[1][0][0] != section_ends[0][0]) && (check_break_pos[1][0][0] != 
							section_starts[2][0])) {
								num_sections=set_section_min_maxs(0, 1, 0, break_sec_num);
								scores[i][j]+=break_and_rejoin(check_break_pos[1][0][0], section_starts[0][0],
									section_ends[2][0], save_index2, save_index,  num_sections, section_mins, section_maxs);
						}
					}
						
				}

				if ((num_check_breaks[1][1] != 0) && (check_break_pos[1][0][1] != section_ends[0][0]) 
					&& (check_break_pos[1][0][1] != section_starts[2][0])) {
						num_sections=set_section_min_maxs(1, 1, 0, break_sec_num);
						scores[i][j]+=break_and_rejoin(check_break_pos[1][0][1], section_starts[0][0],
								section_ends[2][0], save_index2, save_index, num_sections, section_mins, section_maxs);
				}

				if (num_check_breaks[1][2] != 0) {
					num_sections=set_section_min_maxs(2, 1, 0, break_sec_num);	
					scores[i][j]+=break_and_rejoin(check_break_pos[1][0][2], section_starts[0][0],
								section_ends[2][0], save_index2, save_index, num_sections, section_mins, section_maxs);
				}


				for(k=num_old-1; k>=lost_breaks[2]+lost_breaks[1]+lost_breaks[0]; k--) {
					old_lasts[k]->next=old_nexts[k];
					old_nexts[k]->last=old_lasts[k];
				}
				
				num_old = lost_breaks[1]+lost_breaks[0]+lost_breaks[2];

				for(k=save_index; k<save_index2; k++) {
					if (new_lasts[k]->next == new_nexts[k]) 
						new_lasts[k]->next=0;
					if (new_nexts[k]->last == new_lasts[k])
						new_nexts[k]->last=0;
				}
				
			}

			section_starts[1][0]->last=section_starts[1][1]->last=section_ends[1][0]->next=section_ends[1][1]->next=0;
			section_ends[0][0]->next=section_starts[2][0];
			section_starts[2][0]->last=section_ends[0][0];

			section_stacks[0][0][0].clear();
			section_stacks[2][1][0].clear();

			for(j=lost_breaks[1]+lost_breaks[0]; j<num_old; j++) {
				old_lasts[j]->next=old_nexts[j];
				old_nexts[j]->last=old_lasts[j];
			}
			num_old -= lost_breaks[2];
		}

	}	



	curr_pos=save_pos->next;
	save_pos->next=section_starts[1][0];
	section_starts[1][0]->last=save_pos;
	curr_pos->last=section_ends[1][0];
	section_ends[1][0]->next=curr_pos;

	section_stacks[1][0][0].clear();
	section_stacks[1][0][1].clear();

	//Clear any joins we made when we removed the section
	for(k=0; k<save_index; k++) {
		new_lasts[k]->next=0;
		new_nexts[k]->last=0;
	}

	for(i=0; i<num_old; i++) {
		old_lasts[i]->next=old_nexts[i];
		old_nexts[i]->last=old_lasts[i];
	}


	list_len += (end_pos-start_pos);
	reversal_pos=0;
	partial_track_pos=partial_track_start;
	last_list_pos=0;
}







int Genome_Track::do_reorg(int n, int m, int pattern_num, BOOL keep)
{
		//Note that we never try n=list len, since this possibility is implicit in n=0
	int i, j, k, lost_breaks, new_lost_breaks=0, pos[3], rev[3], section_pos[3], save_index=0, break_sec_num[2][3], 
		num_check_breaks[2][3], new_joins=0, num_sections;
	BOOL do_2nd_rejoin;
	Tracking_List *curr_pos, *curr_pos2, *temp,  *find_pos, *dum;

	
	dum=get_list_pos(11);
	curr_pos=get_list_pos(n);
	section_mins[0]=partial_track_start->element[0]->partition_num;
	section_maxs[0]=curr_pos->element[0]->partition_num;
	lost_breaks=find_partition_breaks(curr_pos->next, partial_track_start, 1, section_mins, section_maxs);
			
	num_old=lost_breaks;

	for(i=0; i<lost_breaks; i++) {
		old_lasts[i]=broken_lefts[i];
		old_nexts[i]=broken_rights[i];
	}
	
	
	find_possible_internal_rejoins(lost_breaks, 0, 0, break_sec_num, num_check_breaks, 
		curr_pos);

	section_starts[0][0]=partial_track_start;
	section_ends[0][0]=curr_pos;

	section_starts[1][0]=curr_pos->next;
	section_ends[1][0]=partial_track_end;

	for(i=0; i<(*reorg_patterns)[pattern_num].get_num_pieces(); i++) {
			section_pos[(*reorg_patterns)[pattern_num].get_piece_num(i)]=i;
			pos[i]=(*reorg_patterns)[pattern_num].get_piece_num(i);
			rev[i]=(*reorg_patterns)[pattern_num].get_reversal(i);
	}


	if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[0]) == 1) {
		reverse_section(partial_track_start, curr_pos, section_starts[0][1], section_ends[0][1]);
	
		if (keep == TRUE) {
			//Move the old list to the reversal section
			temp=partial_track_start;
			i=0;
			while(temp != curr_pos->next) 
			{
				hold_reversal_list[i++]=temp;
				temp=temp->next;
			}
		}
	}
		

	if (pattern_num < 3) {
		//Only need to initialize the stacks once
		if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 1) {
			reverse_section(curr_pos->next, partial_track_end, section_starts[1][1], section_ends[1][1]);
		
			if (keep ==TRUE) {
				//Move the old list to the reversal section
				temp=curr_pos->next;
				i=0;
				while(hold_reversal_list[i] != section_starts[1][1]) i++;


				while(temp != partial_track_end->next) 
				{
					hold_reversal_list[i++]=temp;
					temp=temp->next;
				}
			}
		}

		if (((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 1) ||
			((*reorg_patterns)[pattern_num].get_reversal(section_pos[0]) == 1)) {
			if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[0]) == 0) 
				section_ends[0][1]=0;
			if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 0)
				section_starts[1][1] = 0;

			find_reversal_internal_rejoins(0, num_check_breaks, section_ends[0][1], section_starts[1][1]);

		}

		save_index =0;
		section_starts[pos[0]][rev[0]]->last=0;
		section_ends[pos[1]][rev[1]]->next=0;

		partial_track_start=section_starts[pos[0]][rev[0]];
		partial_track_end=section_ends[pos[1]][rev[1]];

		section_ends[pos[0]][rev[0]]->next=section_starts[pos[1]][rev[1]];
		section_starts[pos[1]][rev[1]]->last=section_ends[pos[0]][rev[0]];
			
		set_stacks(&new_left, section_ends[pos[0]][rev[0]],	section_starts[pos[0]][rev[0]], TRUE);
		set_stacks(&new_right, section_starts[pos[1]][rev[1]], section_ends[pos[1]][rev[1]], FALSE);			
		new_joins+=	make_joins(&new_left, &new_right, section_ends[pos[0]][rev[0]]->element[0], 
			section_starts[pos[1]][rev[1]]->element[0], 
			TRUE, save_index);
		
		new_left.clear();
		new_right.clear();

		if (num_check_breaks[0][0] != 0) {
			num_sections = set_section_min_maxs(pattern_num, rev, section_pos, pos, 0, 0, break_sec_num);
			new_joins+=	break_and_rejoin(check_break_pos[0][rev[section_pos[break_sec_num[0][0]]]][0], section_starts[pos[0]][rev[0]],
				section_ends[pos[1]][rev[1]], save_index, 0, num_sections, section_mins, section_maxs);	
		}

		if (num_check_breaks[1][0] !=0 ) {
			num_sections = set_section_min_maxs(pattern_num, rev, section_pos, pos, 0, 1, break_sec_num);
			new_joins+=	break_and_rejoin(check_break_pos[1][rev[section_pos[break_sec_num[1][0]]]][0], section_starts[pos[0]][rev[0]],
				section_ends[pos[1]][rev[1]], save_index, 0, num_sections, section_mins, section_maxs);	
		}

		
	}
	else {
		curr_pos2=get_list_pos(m);
		section_ends[1][0]=curr_pos2;
		section_starts[2][0]=curr_pos2->next;
		section_ends[2][0]=partial_track_end;

		if (num_check_breaks[1][0] != 0) {
			find_pos=section_starts[1][0];
			break_sec_num[1][0]=2;
			while(find_pos != section_ends[1][0]->next) {
				if (find_pos == check_break_pos[1][0][0]->next)
					break_sec_num[1][0]=1;
				find_pos=find_pos->next;
			}

			if ((find_pos == section_ends[1][0]) || (find_pos == section_starts[2][0]))
				do_2nd_rejoin=FALSE;
			else
				do_2nd_rejoin=TRUE;
		
		}

		for(j=0; j<2; j++) {
			section_mins[j]=section_starts[j][0]->element[0]->partition_num;
			section_maxs[j]=section_ends[j][0]->element[0]->partition_num;
		}

		new_lost_breaks=find_partition_breaks(section_starts[2][0], section_starts[0][0], 2, section_mins, section_maxs);

		num_old=lost_breaks+new_lost_breaks;

		for(j=0; j<new_lost_breaks; j++) {
			old_lasts[j+lost_breaks]=broken_lefts[j];
			old_nexts[j+lost_breaks]=broken_rights[j];
		}

		find_possible_internal_rejoins(new_lost_breaks, lost_breaks, 1, break_sec_num, num_check_breaks, 
			curr_pos2);


		if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 1) {
			reverse_section(curr_pos->next, curr_pos2, section_starts[1][1], section_ends[1][1]);
			if (keep == TRUE) {
				//Move the old list to the reversal section
				temp=curr_pos->next;
				
				i=0;
				while(hold_reversal_list[i] != section_starts[1][1]) i++;

				while(temp != curr_pos2->next) 
				{
					hold_reversal_list[i++]=temp;
					temp=temp->next;
				}
			}
		}
		
		if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[2]) == 1) {
			reverse_section(curr_pos2->next, partial_track_end, section_starts[2][1], section_ends[2][1]);
		
			if (keep == TRUE) {
				//Move the old list to the reversal section
				temp=curr_pos2->next;
			
				i=0;
				while(hold_reversal_list[i] != section_starts[2][1]) i++;

				while(temp != partial_track_end->next) 
				{
					hold_reversal_list[i++]=temp;
					temp=temp->next;
				}
			}
		}

		if (((*reorg_patterns)[pattern_num].get_reversal(section_pos[0]) == 1) ||
			((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 1) ||
			((*reorg_patterns)[pattern_num].get_reversal(section_pos[2]) == 1)) {
			
			if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[0]) == 0)
				section_ends[0][1]=0;
			if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 0) 
				section_starts[1][1]=section_ends[1][1]=0;
			if ((*reorg_patterns)[pattern_num].get_reversal(section_pos[2]) == 0)
				section_starts[2][1]=0;
	
			if (((*reorg_patterns)[pattern_num].get_reversal(section_pos[0]) == 1) ||
				((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 1)) {
				if ((do_2nd_rejoin == TRUE) && (break_sec_num[1][0] == 1))
					find_reversal_internal_rejoins(0, num_check_breaks, 
							section_ends[0][1], section_starts[1][1]);
				else
					find_reversal_internal_rejoins(0, num_check_breaks, section_ends[0][1], 0);
			}


			if (((*reorg_patterns)[pattern_num].get_reversal(section_pos[1]) == 1) ||
				((*reorg_patterns)[pattern_num].get_reversal(section_pos[2]) == 1))
					find_reversal_internal_rejoins(1, num_check_breaks, 
						section_ends[1][1], section_starts[2][1]);

			
		}

		section_starts[pos[0]][rev[0]]->last=0;
		section_ends[pos[2]][rev[2]]->next=0;
		
		partial_track_start=section_starts[pos[0]][rev[0]];
		partial_track_end=section_ends[pos[2]][rev[2]];

		section_ends[pos[0]][rev[0]]->next = section_starts[pos[1]][rev[1]];
		section_starts[pos[1]][rev[1]]->last = section_ends[pos[0]][rev[0]];
					
		section_ends[pos[1]][rev[1]]->next = section_starts[pos[2]][rev[2]];
		section_starts[pos[2]][rev[2]]->last = section_ends[pos[1]][rev[1]];

		set_stacks(&new_left, section_ends[pos[0]][rev[0]],	section_starts[pos[0]][rev[0]], TRUE);
		set_stacks(&new_middle_right, section_starts[pos[1]][rev[1]], section_ends[pos[1]][rev[1]], FALSE);	
		new_joins+=		make_joins(&new_left, &new_middle_right, section_ends[pos[0]][rev[0]]->element[0], 
			section_starts[pos[1]][rev[1]]->element[0], TRUE, save_index);
	
			
		set_stacks(&new_middle_left, section_ends[pos[1]][rev[1]], section_starts[pos[0]][rev[0]], TRUE);
		set_stacks(&new_right, section_starts[pos[2]][rev[2]], section_ends[pos[2]][rev[2]], FALSE);
		new_joins+=		make_joins(&new_middle_left, &new_right, section_ends[pos[1]][rev[1]]->element[0], 
				section_starts[pos[2]][rev[2]]->element[0], TRUE, save_index);
		
		new_left.clear();
		new_middle_right.clear();
		new_middle_left.clear();
		new_right.clear();

	
		if (rev[0] == 0) {
			section_mins[0]=section_starts[pos[0]][rev[0]]->element[0]->partition_num;
			section_maxs[0]=section_ends[pos[0]][rev[0]]->element[0]->partition_num;
		}
		else {
			section_maxs[0]=section_starts[pos[0]][rev[0]]->element[0]->partition_num;
			section_mins[0]=section_ends[pos[0]][rev[0]]->element[0]->partition_num;
		}
		new_joins+=	break_and_rejoin(section_ends[pos[0]][rev[0]], section_starts[pos[0]][rev[0]],
			section_ends[pos[2]][rev[2]], save_index, 0, 1, section_mins, section_maxs);	
			
		
		if (num_check_breaks[0][0] != 0) {
			num_sections = set_section_min_maxs(pattern_num, rev, section_pos, pos, 0, 0, break_sec_num);
			new_joins+=	break_and_rejoin(check_break_pos[0][rev[section_pos[break_sec_num[0][0]]]][0], section_starts[pos[0]][rev[0]],
				section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);
			
		}
		if ((num_check_breaks[1][0] != 0) && (do_2nd_rejoin == TRUE)) {
			num_sections = set_section_min_maxs(pattern_num, rev, section_pos, pos, 0, 1, break_sec_num);
			new_joins+=	break_and_rejoin(check_break_pos[1][rev[section_pos[break_sec_num[1][0]]]][0], section_starts[pos[0]][rev[0]],
				section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);
		}

			

		if (num_check_breaks[0][1] != 0) {
			num_sections = set_section_min_maxs(pattern_num, rev, section_pos, pos, 1, 0, break_sec_num);
			new_joins+=	break_and_rejoin(check_break_pos[0][rev[section_pos[break_sec_num[0][1]]]][1], section_starts[pos[0]][rev[0]],
				section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);	
		}

		if (num_check_breaks[1][1] != 0) {
			num_sections = set_section_min_maxs(pattern_num, rev, section_pos, pos, 1, 1, break_sec_num);
			new_joins+=	break_and_rejoin(check_break_pos[1][rev[section_pos[break_sec_num[1][1]]]][1], section_starts[pos[0]][rev[0]],
				section_ends[pos[2]][rev[2]], save_index, 0, num_sections, section_mins, section_maxs);	
		}

			
	}		
		
	if (keep == FALSE) {
		partial_track_start=section_starts[0][0];
		section_ends[0][0]->next=section_starts[1][0];
		section_starts[1][0]->last=section_ends[0][0];
		if (pattern_num >= 3) {
			section_ends[1][0]->next=section_starts[2][0];
			section_starts[2][0]->last=section_ends[1][0];
			partial_track_end=section_ends[2][0];
		}
		else
			partial_track_end=section_ends[1][0];
			partial_track_start->last=partial_track_end->next=0;
	
			for(k=num_old-1; k>=lost_breaks+new_lost_breaks; k--) {
				old_lasts[k]->next=old_nexts[k];
				old_nexts[k]->last=old_lasts[k];
			}

			num_old=lost_breaks+new_lost_breaks;

			if (save_index != 0) {
				for(k=0; k<save_index; k++) {
					if (new_lasts[k]->next == new_nexts[k]) 
						new_lasts[k]->next=0;
					if (new_nexts[k]->last == new_lasts[k])
						new_nexts[k]->last=0;
					
				}
			}	

			for (k=0; k<num_old; k++) {
				old_lasts[k]->next=old_nexts[k];
				old_nexts[k]->last=old_lasts[k];
			}

			num_old=0;

	}


	reversal_pos=0;
	partial_track_pos=partial_track_start;
	reset_list_partition_numbers();
	last_list_pos=0;
	return(new_joins-lost_breaks-new_lost_breaks);

}



int Genome_Track::do_section_insert(int start_pos, int end_pos, int loc, int rev, BOOL do_link)
{
	int i, j, k, rejoins, num_joins=0, lost_breaks[3], save_index=0, save_index2=2, 
		last_score, last_score2, reused_losts=0, section1_max, num_sections;
	BOOL do_3rd_rejoin;
	Tracking_List *curr_pos, *save_pos, *save_pos2, *find_pos, *dum, *temp,
		*old_start, *old_end, *use_end, *use_start;
	
	dumele=get_list_pos(11)->element[1];

	/*cout<<"Before "<<start_pos<<"\t"<<end_pos<<"\tloc: "<<loc<<endl;
	old_start=partial_track_start;
	while(old_start != 0) {
		cout<<"("<<old_start->element[0]<<", "<<old_start->element[1]<<")\t";
		old_start=old_start->next;
	}
	cout<<endl;
	for(i=0; i<the_homologs->get_num_homologs(); i++)
		cout<<"("<<hold_reversal_list[i]->element[0]<<", "<<hold_reversal_list[i]->element[1]<<")\t";
	cout<<endl;*/

	old_start=partial_track_start;
	old_end=partial_track_end;
	partial_track_pos=partial_track_start;

	if (start_pos != 0)
		section_starts[1][0]=get_list_pos(start_pos)->next;
	else
		section_starts[1][0]=partial_track_start;

	if (end_pos != the_homologs->get_num_homologs())
		section_ends[1][0]=get_list_pos(end_pos);
	else
		section_ends[1][0]=partial_track_end;

	save_pos=section_starts[1][0]->last;
	save_pos2=section_ends[1][0]->next;

	if (start_pos != 0)
		precut_max_partition=save_pos->element[0]->partition_num;
	else
		precut_max_partition=save_pos2->element[0]->partition_num;

	if (end_pos != the_homologs->get_num_homologs())
		postcut_min_partition=save_pos2->element[0]->partition_num;
	else
		postcut_min_partition=precut_max_partition;

	if (start_pos != 0) {
		section_mins[0]=partial_track_start->element[0]->partition_num;
		section_maxs[0]=section_starts[1][0]->last->element[0]->partition_num;
		lost_breaks[0]=find_partition_breaks(section_starts[1][0], partial_track_start, 1, section_mins, section_maxs);

		for(i=0; i<lost_breaks[0]; i++) {
			old_lasts[i]=broken_lefts[i];
			old_nexts[i]=broken_rights[i];
		}

	}
	else {
		partial_track_start=save_pos2;
		lost_breaks[0]=0;
	}

	num_old=lost_breaks[0];


	if (end_pos != the_homologs->get_num_homologs()) {
		if (start_pos != 0) {
			section_mins[0]=partial_track_start->element[0]->partition_num;
			section_maxs[0]=save_pos->element[0]->partition_num;
			section_mins[1]=section_starts[1][0]->element[0]->partition_num;
			section_maxs[1]=section_ends[1][0]->element[0]->partition_num;
			lost_breaks[1]=find_partition_breaks(section_ends[1][0]->next, section_starts[1][0], 2, section_mins, section_maxs);
		}
		else {
			section_mins[0]=section_starts[1][0]->element[0]->partition_num;
			section_maxs[0]=section_ends[1][0]->element[0]->partition_num;
			lost_breaks[1]=find_partition_breaks(section_ends[1][0]->next, section_starts[1][0], 1, section_mins, section_maxs);
		}
	}
	else {
		lost_breaks[1] = 0;
		partial_track_end=section_starts[1][0]->last;
	}

	for(i=0; i<lost_breaks[1]; i++) {
		old_lasts[i+num_old]=broken_lefts[i];
		old_nexts[i+num_old]=broken_rights[i];
	}
	
	num_old += lost_breaks[1];


	list_len-=(end_pos-start_pos);

	if ((start_pos != 0) && (end_pos != the_homologs->get_num_homologs())) {
		section_starts[1][0]->last->next=section_ends[1][0]->next;
		section_ends[1][0]->next->last=section_starts[1][0]->last;
	}
	else {
		if (start_pos == 0)
			section_ends[1][0]->next->last=0;
		else
			section_starts[1][0]->last->next=0;
	}

	section_starts[1][0]->last=0;
	section_ends[1][0]->next=0;
	if (rev == 1) {
		reverse_section(section_starts[1][0], section_ends[1][0], section_starts[1][1], section_ends[1][1]);
		if (do_link == TRUE) {
				//Move the old list to the reversal section
				temp=section_starts[1][0];
			
				i=0;
				while(hold_reversal_list[i] != section_starts[1][1]) i++;

				while(temp != section_ends[1][1]->next) 
				{
					hold_reversal_list[i++]=temp;
					temp=temp->next;
				}
		}

	}
				


	if (rev == 0)
		set_stacks(&section_stacks[1][0][0], section_starts[1][0], section_ends[1][0], FALSE);
	else
		set_stacks(&section_stacks[1][0][1], section_starts[1][1], section_ends[1][1], FALSE);
	
	if ((start_pos != 0) && (end_pos != the_homologs->get_num_homologs())) {
		set_stacks(&new_left, save_pos, partial_track_start, TRUE);
		set_stacks(&new_right, save_pos->next, partial_track_end, FALSE);
		num_joins += make_joins(&new_left, &new_right, save_pos->element[0], save_pos->next->element[0], TRUE, save_index);
	}


	new_right.clear();
	new_left.clear();
	new_middle_left.clear();
	
	last_list_pos=0;
	partial_track_pos=partial_track_start;

	section_starts[0][0]=partial_track_start;

	if (loc != 0) {
		curr_pos=get_list_pos(loc);
		section_ends[0][0]=curr_pos;
		if (curr_pos->next != 0) {
			section_starts[2][0]=curr_pos->next;
			section_ends[2][0]=partial_track_end;
		}
		else
			section_starts[2][0]=section_ends[2][0]=0;
	}
	else {
		section_ends[0][0]=partial_track_end;
		section_starts[2][0]=section_ends[2][0]=0;
	}


	if (section_starts[2][0] != 0) {
		section_mins[0]=section_starts[0][0]->element[0]->partition_num;
		section_maxs[0]=section_ends[0][0]->element[0]->partition_num;
		lost_breaks[2]=find_partition_breaks(section_starts[2][0], section_starts[0][0], 1, section_mins, section_maxs);
	}
	else
		lost_breaks[2]=0;

	for(j=0; j<lost_breaks[2]; j++){
		old_lasts[j+num_old]=broken_lefts[j];
		old_nexts[j+num_old]=broken_rights[j];
	}

	num_old=lost_breaks[2]+lost_breaks[0]+lost_breaks[1];

	set_stacks(&section_stacks[0][0][0], section_ends[0][0], section_starts[0][0], TRUE);
	
	if (section_starts[2][0] != 0) 
		set_stacks(&section_stacks[2][1][0], section_starts[2][0], section_ends[2][0], FALSE);


	save_index2=save_index;
	if (loc != 0) {
		section_ends[0][0]->next=section_starts[1][rev];
		section_starts[1][rev]->last=section_ends[0][0];
	}
	else {
		section_ends[1][rev]->next=section_starts[0][0];
		section_starts[0][0]->last=section_ends[1][rev];
	}

	if (section_starts[2][0] != 0) {
		section_starts[2][0]->last=section_ends[1][rev];
		section_ends[1][rev]->next=section_starts[2][0];
	}
	
	if (loc != 0) {
		num_joins+=make_joins(&section_stacks[0][0][0], &section_stacks[1][0][rev], 
				section_ends[0][0]->element[0], section_starts[1][rev]->element[0], 
				TRUE, save_index2)-lost_breaks[0]-lost_breaks[1]-lost_breaks[2];
	
		if (section_starts[2][0] != 0) {
			set_stacks(&new_middle_left, section_ends[1][rev], section_starts[0][0], TRUE);
			num_joins+=make_joins(&new_middle_left, &section_stacks[2][1][0], section_ends[1][rev]->element[0],
						section_starts[2][0]->element[0], TRUE, save_index2);

			new_middle_left.clear();
		}
		
	}
	else {
			set_stacks(&new_left, section_ends[1][rev], section_starts[1][rev], TRUE);
			set_stacks(&new_right, section_starts[0][0], section_ends[0][0], FALSE);
			num_joins+=make_joins(&new_left, &new_right, section_ends[1][rev]->element[0],
						section_starts[1][0]->element[0], TRUE, save_index2)-
						lost_breaks[0]-lost_breaks[1]-lost_breaks[2];
					
			new_left.clear();
			new_right.clear();	
	}

	if (section_ends[2][0] == 0) {
		if (loc == 0) {
			use_start=section_starts[1][rev];
			use_end=section_ends[0][0];
		}
		else {
			use_end=section_ends[1][rev];
			use_start=section_starts[0][0];
		}
	}
	else {
		use_start=section_starts[0][0];
		use_end=section_ends[2][0];
	}
		


	section1_max=section_ends[1][0]->element[0]->partition_num;
	do_3rd_rejoin = TRUE;
		
	if (loc != 0) {
		//Rejoin from the current join out to the next dupl
		find_pos=section_ends[0][0]->last;
	
		if (find_pos != 0) {
			do {
				if ((section1_max < find_pos->element[0]->partition_num) && (start_pos != 0)) {
					num_sections=2;
					section_mins[0]=section_starts[0][0]->element[0]->partition_num;
					section_maxs[0]=section_starts[1][0]->element[0]->partition_num-1;
					section_mins[1]=section_ends[1][0]->element[0]->partition_num+1;
					section_maxs[1]=find_pos->element[0]->partition_num;
				}
				else {
					num_sections=1;
					section_mins[0]=section_starts[0][0]->element[0]->partition_num;
					section_maxs[0]=find_pos->element[0]->partition_num;
				}
				
				num_joins+=break_and_rejoin(find_pos, use_start, use_end, save_index2, save_index,
					num_sections, section_mins, section_maxs, do_link);	
				find_pos=find_pos->last;
			} while ((find_pos != 0) && (find_pos->element[1]->my_locus == 0));
		}

		find_pos=section_starts[1][rev];

		if ((find_pos != 0) && (find_pos != section_ends[1][rev])) {
			do {
				if ((section1_max < section_ends[0][0]->element[0]->partition_num) && (start_pos !=0)) {
						num_sections=3;
						section_mins[0]=section_starts[0][0]->element[0]->partition_num;
						section_maxs[0]=section_starts[1][0]->element[0]->partition_num-1;
						section_mins[1]=section_ends[1][0]->element[0]->partition_num+1;
						section_maxs[1]=section_ends[0][0]->element[0]->partition_num;
							
						if (rev == 0) {
							section_mins[2]=section_starts[1][0]->element[0]->partition_num;
							section_maxs[2]=find_pos->element[0]->partition_num;
						}
						else {
							section_mins[2]=find_pos->element[0]->partition_num;
							section_maxs[2]=section_starts[1][1]->element[0]->partition_num;
						}
					}
					else {
						num_sections=2;
						section_mins[0]=section_starts[0][0]->element[0]->partition_num;
						section_maxs[0]=section_ends[0][0]->element[0]->partition_num;

						if (rev == 0) {
							section_mins[1]=section_starts[1][0]->element[0]->partition_num;
							section_maxs[1]=find_pos->element[0]->partition_num;
						}
						else {
							section_mins[1]=find_pos->element[0]->partition_num;
							section_maxs[1]=section_starts[1][1]->element[0]->partition_num;
						}
					}
				
				num_joins+=break_and_rejoin(find_pos, use_start, use_end, save_index2, save_index,
					num_sections, section_mins, section_maxs, do_link);	
				find_pos=find_pos->next;
			} while((find_pos != section_ends[1][rev]) &&(find_pos->element[1]->my_locus == 0));

			if (find_pos == section_ends[1][rev])
				do_3rd_rejoin = FALSE;

		}
	}

	find_pos=section_ends[1][rev]->last;

	if ((find_pos != 0) && (do_3rd_rejoin == TRUE)  && (find_pos != section_starts[1][rev]->last)) {
			do {
				if (loc != 0) {
					if ((section1_max < section_ends[0][0]->element[0]->partition_num) && (start_pos != 0)) {
						num_sections=3;
						section_mins[0]=section_starts[0][0]->element[0]->partition_num;
						section_maxs[0]=section_starts[1][0]->element[0]->partition_num-1;
						section_mins[1]=section_ends[1][0]->element[0]->partition_num+1;
						section_maxs[1]=section_ends[0][0]->element[0]->partition_num;

						if (rev == 0) {
							section_mins[2]=section_starts[1][0]->element[0]->partition_num;
							section_maxs[2]=find_pos->element[0]->partition_num;
						}
						else {
							section_mins[2]=find_pos->element[0]->partition_num;
							section_maxs[2]=section_starts[1][1]->element[0]->partition_num;
						}
					}
					else {
						num_sections=2;
						section_mins[0]=section_starts[0][0]->element[0]->partition_num;
						section_maxs[0]=section_ends[0][0]->element[0]->partition_num;
						if (rev == 0) {
							section_mins[1]=section_starts[1][rev]->element[0]->partition_num;
							section_maxs[1]=find_pos->element[0]->partition_num;
						}
						else {
							section_mins[1]=find_pos->element[0]->partition_num;
							section_maxs[1]=section_starts[1][1]->element[0]->partition_num;
						}
					}
				}
				else {
					num_sections =1;
					if (rev == 0) {
						section_mins[0]=section_starts[1][0]->element[0]->partition_num;
						section_maxs[0]=find_pos->element[0]->partition_num;
					}
					else {
						section_mins[0]=find_pos->element[0]->partition_num;
						section_maxs[0]=section_starts[1][1]->element[0]->partition_num;
					}
				}
				
				num_joins+=break_and_rejoin(find_pos, use_start, use_end, save_index2, save_index,
					num_sections, section_mins, section_maxs, do_link);	
				find_pos=find_pos->last;
			} while ((find_pos != section_starts[1][rev]->last) && (find_pos->element[1]->my_locus == 0) );

	}

	if (!((loc != 0) &&(section_starts[2][0] == 0))) {
		if (section_starts[2][0] == 0) 
			find_pos=section_starts[0][0];
		else
			find_pos=section_starts[2][0];

		if (find_pos->next != 0) {
			do {
				if (loc != 0) {
					if (section1_max > section_starts[2][0]->element[0]->partition_num) {
						num_sections=3;
						section_mins[0]=section_starts[0][0]->element[0]->partition_num;
						section_maxs[0]=section_ends[0][0]->element[0]->partition_num;
						section_mins[1]=section_starts[1][0]->element[0]->partition_num;
						section_maxs[1]=section_ends[1][0]->element[0]->partition_num;
						section_mins[2]=section_starts[2][0]->element[0]->partition_num;
						section_maxs[2]=find_pos->element[0]->partition_num;
					}
					else  {
						num_sections=1;
						section_mins[0]=0;
						section_maxs[0]=find_pos->element[0]->partition_num;
					}
				}
				else {
					if (find_pos->element[0]->partition_num < section_starts[1][0]->element[0]->partition_num) {
						num_sections = 2;
						section_mins[0]=section_starts[1][0]->element[0]->partition_num;
						section_maxs[0]=section_ends[1][0]->element[0]->partition_num;
						section_mins[1]=section_starts[0][0]->element[0]->partition_num;
						section_maxs[1]=find_pos->element[0]->partition_num;
					}
					else {
						num_sections = 1;
						section_mins[0]=0;
						section_maxs[0]=find_pos->element[0]->partition_num;
					}
				}
				num_joins+=break_and_rejoin(find_pos, use_start, use_end, save_index2, save_index,
					num_sections, section_mins, section_maxs, do_link);
				find_pos=find_pos->next;
			} while ((find_pos->next != 0) && (find_pos->element[1]->my_locus == 0));
		}

		if (save_pos != 0)
			find_pos=save_pos->last;
		else
			find_pos =0;

		if ((find_pos != 0) && (find_pos != section_ends[1][rev])) {
			do {
				if (((section_starts[2][0] != 0) && 
					(find_pos->element[0]->partition_num < section_starts[2][0]->element[0]->partition_num)) || 
					((loc != 0) && (section_starts[2][0] == 0))) {
					num_sections = 1;
					section_mins[0]=section_starts[0][0]->element[0]->partition_num;
					section_maxs[0]=find_pos->element[0]->partition_num;
				}
				else {
					num_sections=2;
					section_mins[0]=section_starts[0][0]->element[0]->partition_num;
					section_maxs[0]=find_pos->element[0]->partition_num;
					section_mins[1]=section_starts[1][0]->element[0]->partition_num;
					section_maxs[1]=section_ends[1][0]->element[0]->partition_num;
				}

				num_joins+=break_and_rejoin(find_pos, partial_track_start, partial_track_end, save_index2, save_index,
					num_sections, section_mins, section_maxs, do_link);	
				find_pos=find_pos->last;
			}	while ((find_pos != 0) && (find_pos->element[1]->my_locus == 0) && (find_pos != section_ends[1][rev]));

		}
		
		find_pos=save_pos2;

		if ((save_pos != 0) && (find_pos !=0) && (find_pos->next != 0)  && (find_pos->next != section_starts[1][rev])) {
			do {
				if (((section_starts[2][0] != 0) && 
					(find_pos->element[0]->partition_num < section_starts[2][0]->element[0]->partition_num)) || 
					((loc != 0) && (section_starts[2][0] == 0))) {
					num_sections = 2;
					section_mins[0]=section_starts[0][0]->element[0]->partition_num;
					section_maxs[0]=save_pos->element[0]->partition_num;
					section_mins[1]=save_pos2->element[0]->partition_num;
					section_maxs[1]=find_pos->element[0]->partition_num;
				}
				else {
					num_sections=3;
					section_mins[0]=section_starts[0][0]->element[0]->partition_num;
					section_maxs[0]=save_pos->element[0]->partition_num;
					section_mins[1]=section_starts[1][0]->element[0]->partition_num;
					section_maxs[1]=section_ends[1][0]->element[0]->partition_num;
					section_mins[2]=save_pos2->element[0]->partition_num;
					section_maxs[2]=find_pos->element[0]->partition_num;
				}
				
				num_joins+=break_and_rejoin(find_pos, partial_track_start, partial_track_end, save_index2, save_index,
					num_sections, section_mins, section_maxs, do_link);	
				find_pos=find_pos->next;

			} while((find_pos->next != 0) && (find_pos->element[1]->my_locus == 0) && (find_pos->next != section_starts[1][rev]));
		}
	

	}

	if (do_link == FALSE) {
		if (section_starts[2][0] != 0) {
			section_ends[0][0]->next=section_starts[2][0];
			section_starts[2][0]->last=section_ends[0][0];
		}
		else
			partial_track_end->next=0;

		if (start_pos != 0) {
			save_pos->next=section_starts[1][0];
			section_starts[1][0]->last=save_pos;
		}
		else {
			partial_track_start=section_starts[1][0];
			
		}
		
		if (end_pos != the_homologs->get_num_homologs()) {
			save_pos2->last=section_ends[1][0];
			section_ends[1][0]->next=save_pos2;
			save_pos2->last=section_ends[1][0];
		}
		else {
			partial_track_end=section_ends[1][0];
		}

		partial_track_start->last = 0;
		partial_track_end->next = 0;

		for(k=num_old-1; k>=lost_breaks[0]+lost_breaks[1]+lost_breaks[2]; k--) {
				old_lasts[k]->next=old_nexts[k];
				old_nexts[k]->last=old_lasts[k];
		}

		num_old=lost_breaks[0]+lost_breaks[1]+lost_breaks[2];

		if (save_index2 != 0) {
			for(k=save_index; k<save_index2; k++) {
				if (new_lasts[k]->next == new_nexts[k]) 
					new_lasts[k]->next=0;
				if (new_nexts[k]->last == new_lasts[k])
					new_nexts[k]->last=0;
				
			}
		}	

		for (k=lost_breaks[0]+lost_breaks[1]; k<num_old; k++) {
			old_lasts[k]->next=old_nexts[k];
			old_nexts[k]->last=old_lasts[k];
		}
	
		num_old=lost_breaks[0]+lost_breaks[1];

		for(k=0; k<save_index; k++) {
				if (new_lasts[k]->next == new_nexts[k]) 
					new_lasts[k]->next=0;
				if (new_nexts[k]->last == new_lasts[k])
					new_nexts[k]->last=0;
				
			}

		for (k=0; k<num_old; k++) {
			old_lasts[k]->next=old_nexts[k];
			old_nexts[k]->last=old_lasts[k];
		}
	
		num_old=0;	

		if (old_start != partial_track_start)
			cerr<<"Error: wrong start\n";
		if (old_end != partial_track_end)
			cerr<<"Error: wrong end\n";
		if ((partial_track_start->last != 0) || (partial_track_start->next == 0))
			cerr<<"Error: bad list start\n";
		if((partial_track_end->next != 0) ||(partial_track_end->last ==0))
			cerr<<"Error: bad list end\n";
	}	
	else {
		if (loc == 0)
			partial_track_start=section_starts[1][rev];
		else
			partial_track_start=section_starts[0][0];

		if (section_ends[2][0] == 0) {
			if (loc == 0)
				partial_track_end=section_ends[0][0];
			else
				partial_track_end=section_ends[1][rev];
		}
		else
			partial_track_end=section_ends[2][0];
	}


	
	section_stacks[0][0][0].clear();
	section_stacks[2][1][0].clear();

	section_stacks[1][0][0].clear();
	section_stacks[1][0][1].clear();


	list_len += (end_pos-start_pos);
	reversal_pos=0;
	partial_track_pos=partial_track_start;
	reset_list_partition_numbers();
	last_list_pos=0;
	return(num_joins);
}



Genome_Track::~Genome_Track()
{
	int i;
	Tracking_List *temp;

	for(i=0; i<the_homologs->get_num_homologs(); i++) 
		delete[] inferred_tracking[i];
	
	if (reversal_pointers != 0)
		delete[] reversal_pointers;
	
	delete[] inferred_tracking;


	if (hold_reversal_list != 0)
	{
		for(i=0; i<the_homologs->get_num_homologs(); i++)
		{
			delete hold_reversal_list[i]->element[0];
			delete hold_reversal_list[i]->element[1];
			delete hold_reversal_list[i];
		}

		delete[] hold_reversal_list;
	}

	if (partial_track_start != 0) {
		partial_track_pos=partial_track_start;

		while(partial_track_pos != 0) {
			partial_track_end=partial_track_pos->next;

			delete partial_track_pos->element[0];
			delete partial_track_pos->element[1];
			delete[] partial_track_pos->element;

			delete partial_track_pos;
			partial_track_pos=partial_track_end;

		}

	}
}




void Genome_Track::recurse_assemble(Track_Stack *&lefts, Track_Stack *&rights,
												int left_end, int size)

{
	int i, j,  size_left, size_right, depth, num_joins=0, join;
	Gene_Track_List *curr_left, *curr_right, *valid_left, *valid_right;
	Track_Stack *left_lefts, *left_rights, *right_rights, *right_lefts;
	BOOL stop;

	if (size > 2) {
		size_left=size_right=depth=0;
		while (pow2[depth] < size) depth++;
		size_left=pow2[depth-1];
		size_right=size-size_left;

		recurse_assemble(left_lefts, left_rights, left_end, size_left);
		recurse_assemble(right_lefts, right_rights, left_end+size_left, size_right);

		//Join
		while((left_rights->get_size()>0) && (num_joins<2)) {
			curr_left=left_rights->get_bottom();
			right_lefts->reset();
			join=0;
			curr_right=right_lefts->get_next();
			while((curr_right != 0) && (join == 0)) {			
				join=assign_adjacency(curr_left, curr_right);
				if (join == 1) {
				//We made a join
					right_lefts->delete_element(curr_right);
					left_rights->pop();
					num_joins++;
					
					//Check left stack--If our join "crossed" a locus that will be the only valid remaining end
					if (inferred_tracking[left_end+size_left-1][0].my_locus->has_duplicate() == FALSE)
					{
						if (curr_left != &inferred_tracking[left_end+size_left-1][0]) {
							valid_left=&inferred_tracking[left_end+size_left-1][0];
							delete left_rights;
							left_rights=new Track_Stack;
							left_rights->push(valid_left);
						}

					}
					//Check right stack
					if (inferred_tracking[left_end+size_left][0].my_locus->has_duplicate() == FALSE)
					{
						if (curr_right != &inferred_tracking[left_end+size_left][0]) {
							valid_right=&inferred_tracking[left_end+size_left][0];
							delete right_lefts;
							right_lefts=new Track_Stack;
							right_lefts->push(valid_right);
						}

					}

				}
				else
					curr_right=right_lefts->get_next();

			}

			if (join == 0) {
			//We've been through all the rights for this left and can't join it--get rid of it
				left_rights->pop(); }
		}

		delete left_lefts;
		delete right_rights;
		delete left_rights;
		delete right_lefts;

		//Now find the ends to pass off to the next level
		lefts=new Track_Stack();
		rights=new Track_Stack();
		stop=FALSE;

		//Find all the possible left endpoints
		i=left_end;
		while((stop==FALSE) && (i<left_end+size)) {
			if (inferred_tracking[i][0].last == 0) {
				//This is a left end
				lefts->push(&inferred_tracking[i][0]);

				if (i>left_end) {
					//If the previous entry has a join but not to this entry, then
					//this one must be on the other track.  Hence, stop
					if (inferred_tracking[i-1][0].next !=0)
						stop=TRUE;
				}
			}
			if ((inferred_tracking[i][1].my_locus != 0) && (inferred_tracking[i][1].last ==0)) {
				//This is a left end and a stopping point
				lefts->push(&inferred_tracking[i][1]);
				stop=TRUE;
			}

			i++;
		}

		stop=FALSE;
		//Find all the possible right endpoints
		i=left_end+size-1;
		while((stop==FALSE) && (i>=left_end)) {
			if (inferred_tracking[i][0].next == 0) {
				//This is a left end
				rights->push(&inferred_tracking[i][0]);

				if (i<left_end+size-1) {
					//If the previous entry has a join but not to this entry, then
					//this one must be on the other track.  Hence, stop
					if (inferred_tracking[i+1][0].last !=0)
						stop=TRUE;
				}
			}
			if ((inferred_tracking[i][1].my_locus != 0) && (inferred_tracking[i][1].next == 0)) {
				//This is a left end and a stopping point
				rights->push(&inferred_tracking[i][1]);
				stop=TRUE;
			}

			i--;
		}		
		
	}

	else {                
		//Base recursion level
		lefts=new Track_Stack();
		rights=new Track_Stack();
		if (size == 1) {

			lefts->push(&inferred_tracking[left_end][0]);
			if (inferred_tracking[left_end][1].my_locus !=0)
				lefts->push(&inferred_tracking[left_end][1]);

			rights->push(&inferred_tracking[left_end][0]);
			if (inferred_tracking[left_end][1].my_locus !=0)
				rights->push(&inferred_tracking[left_end][1]);
		}

		else {
			//Check for all possible 2 gene joins
			assign_adjacency(&inferred_tracking[left_end][0], &inferred_tracking[left_end+1][0]);
			assign_adjacency(&inferred_tracking[left_end][0], &inferred_tracking[left_end+1][1]);
			assign_adjacency(&inferred_tracking[left_end][1], &inferred_tracking[left_end+1][0]);
			assign_adjacency(&inferred_tracking[left_end][1], &inferred_tracking[left_end+1][1]);
			

			//Do left ends
			lefts->push(&inferred_tracking[left_end][0]);

			if (inferred_tracking[left_end][1].my_locus != 0)
				lefts->push(&inferred_tracking[left_end][1]);
			else {
				//We have a hole in the second left position:  we could have as many as 3 left ends
				if (inferred_tracking[left_end][0].next == 0) {
					//No connections between first and second--still 3 possible ends
					lefts->push(&inferred_tracking[left_end+1][0]);
					if (inferred_tracking[left_end+1][1].my_locus != 0)
						lefts->push(&inferred_tracking[left_end+1][1]);
				}
				else {
					if (inferred_tracking[left_end][0].next == &inferred_tracking[left_end+1][0]) {
						if (inferred_tracking[left_end+1][1].my_locus != 0)
							lefts->push(&inferred_tracking[left_end+1][1]);
					}
					else 
						lefts->push(&inferred_tracking[left_end+1][0]);
				}

			}

			//Do rights
			rights->push(&inferred_tracking[left_end+1][0]);

			if (inferred_tracking[left_end+1][1].my_locus != 0)
				rights->push(&inferred_tracking[left_end+1][1]);
			else {
				//We have a hole in the second right position:  we could have as many as 3 left ends
				if (inferred_tracking[left_end+1][0].last == 0) {
					//No connections between first and second--still 3 possible ends
					rights->push(&inferred_tracking[left_end][0]);
					if (inferred_tracking[left_end][1].my_locus != 0)
						rights->push(&inferred_tracking[left_end][1]);
				}
				else {
					if (inferred_tracking[left_end+1][0].last == &inferred_tracking[left_end][0]) {
						if (inferred_tracking[left_end][1].my_locus != 0)
							rights->push(&inferred_tracking[left_end][1]);
					}
					else 
						rights->push(&inferred_tracking[left_end][0]);
				}

			} 			//Done with rights

		}				//Done with 2-gene join

	}					//Done with <3 gene joins (1 or 2)


}



int Genome_Track::break_and_rejoin(int break_before, int size)
{
	int num_orig, num_new;
	
	num_orig=find_partition_breaks(break_before);
	set_stacks(&new_right, break_before, size, FALSE);
	set_stacks(&new_left, break_before-1, -1, TRUE);
	num_new=make_joins(&new_left, &new_right, &inferred_tracking[break_before-1][0], 
		&inferred_tracking[break_before][0], TRUE);
	
	new_left.clear();
	new_right.clear();

	return(num_new-num_orig);	
}

int Genome_Track::break_and_rejoin(int break_before)
{
	int num_orig, num_new, temp=0;
	Tracking_List *curr_pos;

	curr_pos=get_list_pos(break_before);
	num_old=0;
	section_mins[0]=partial_track_start->element[0]->partition_num;
	section_maxs[0]=curr_pos->element[0]->partition_num;

	return(break_and_rejoin(curr_pos, partial_track_start, partial_track_end, temp, 0, 1, section_mins, section_maxs));
}



int Genome_Track::break_and_rejoin(Tracking_List *right_end, Tracking_List *start, Tracking_List *end, 
								   int &save_index, int start_index, int num_left_sections, int *left_mins, int *left_maxs)
{
	return(break_and_rejoin(right_end, start, end, save_index, start_index, num_left_sections, left_mins, left_maxs, TRUE));
}



int Genome_Track::break_and_rejoin(Tracking_List *right_end, Tracking_List *start, Tracking_List *end, 
								   int &save_index, int start_index, int num_left_sections, int *left_mins, int *left_maxs, BOOL do_link)
{
	//This is now broken for the code which tests all possible reorgs--it will only work for the section_insert code
	int i, j, k, num_orig, num_new, save_index_orig;
	BOOL broken_old[2];

	broken_old[0]=broken_old[1]=TRUE;
	save_index_orig=save_index;

	num_orig=find_partition_breaks(right_end->next, start, num_left_sections, left_mins, left_maxs);

	//Clear the broken ends from the saved list
	if (save_index != 0) {
		for(i=0; i<num_orig; i++) {
			j=start_index;
			while(!((broken_lefts[i] == new_lasts[j]) && (broken_rights[i] == new_nexts[j])) && (j<save_index)) j++;
			
			if(j == save_index )
				broken_old[i] = TRUE;
			else
				broken_old[i] = FALSE;
		}
	}

	if (do_link == FALSE) {
		for(i=0; i<num_orig; i++) {
			old_lasts[i+num_old]=broken_lefts[i];
			old_nexts[i+num_old]=broken_rights[i];
			}
		num_old+=num_orig;
	}

	set_stacks(&new_left, right_end, start, TRUE);
	set_stacks(&new_right, right_end->next, end, FALSE);
	num_new=make_joins(&new_left, &new_right, right_end->element[0], right_end->next->element[0], TRUE, save_index);

	if (do_link == FALSE) {
		i=save_index-num_new;
		while (i<save_index) {
			j=0;
			while((j<num_orig) && !((broken_lefts[j] == new_lasts[i]) && (broken_rights[j] == new_nexts[i]))) j++;
			if (j != num_orig) {
			//This was an existing join
				if (broken_old[j] == TRUE) {
				//We don't want to erase this join
					if (i != save_index-1) {
						new_nexts[i]=new_nexts[i+1];
						new_lasts[i]=new_lasts[i+1];
					}
					save_index--;
					i--;
				}

			}
			i++;
		}
	}
	
	new_right.clear();
	new_left.clear();

	return(num_new-num_orig);
}


Tracking_List * Genome_Track::get_list_pos(int pos)
{
	int i;
	Tracking_List *curr_pos;

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

void Genome_Track::set_stacks(Track_Stack *new_stack, int end, int stop_point, BOOL left)
{
	int i;
	BOOL stop;
	Gene_Track_List *to_end, *to_end_previous, *to_end_2nd;
	stop=FALSE;

	//Find all the possible left endpoints
	i=end;
	while((stop==FALSE) && (i != stop_point)) {
		if (left == FALSE) {
			to_end=inferred_tracking[i][0].last;
			to_end_previous=inferred_tracking[i-1][0].next;
			if (inferred_tracking[i][1].my_locus != 0)
				to_end_2nd=inferred_tracking[i][1].last;
		}
		else {
			to_end=inferred_tracking[i][0].next;
			to_end_previous=inferred_tracking[i+1][0].last;
			if (inferred_tracking[i][1].my_locus != 0)
				to_end_2nd=inferred_tracking[i][1].next;
		}

		if (to_end == 0) {
			//This is a left end
			new_stack->push(&inferred_tracking[i][0]);

			if (i != end) {
				//If the previous entry has a join but not to this entry, then
				//this one must be on the other track.  Hence, stop
				if (to_end_previous !=0)
					stop=TRUE;
			}
		}
		if ((inferred_tracking[i][1].my_locus != 0) && (to_end_2nd ==0)) {
			//This is a left end and a stopping point
			new_stack->push(&inferred_tracking[i][1]);
			stop=TRUE;
		}

		if (left == FALSE)
			i++;
		else
			i--;

			
	}


}

void Genome_Track::set_stacks(Track_Stack *new_stack, Tracking_List *end, Tracking_List *stop_point, BOOL left)
{
	BOOL stop;
	Gene_Track_List *to_end, *from_previous, *to_end_2nd;
	Tracking_List *curr_pos;


	stop=FALSE;

		//Find all the possible left endpoints
	curr_pos=end;
	while(stop==FALSE) {
		if (left == TRUE) {
			to_end=curr_pos->element[0]->next;
			
			if (curr_pos != end)
				from_previous=curr_pos->next->element[0]->last;
			else
				from_previous=0;

			if (curr_pos->element[1]->my_locus != 0)
				to_end_2nd=curr_pos->element[1]->next;
		}
		else {
			to_end=curr_pos->element[0]->last;
			
			if (curr_pos != end)
				from_previous=curr_pos->last->element[0]->next;
			else
				from_previous=0;

			if (curr_pos->element[1]->my_locus != 0)
				to_end_2nd=curr_pos->element[1]->last;
		}

		if (to_end == 0) {
			//This is an end
			new_stack->push(curr_pos->element[0]);

			if (curr_pos != end) {
				//If the previous entry has a join but not to this entry, then
				//this one must be on the other track.  Hence, stop
				if (from_previous !=0)
					stop=TRUE;
			}
		}
		if ((curr_pos->element[1]->my_locus != 0) && (to_end_2nd ==0)) {
			//This is a left end and a stopping point
			new_stack->push(curr_pos->element[1]);
			stop=TRUE;
		}

		if (curr_pos == stop_point)
			stop=TRUE;

		if (left == TRUE)
			curr_pos=curr_pos->last;
		else
			curr_pos=curr_pos->next;
	}
	

}

void Genome_Track::set_link_counts(Track_Stack *lefts, Track_Stack *rights)
{
	int join;
	Gene_Track_List *curr_left, *curr_right;


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


int Genome_Track::make_joins(Track_Stack *lefts, Track_Stack *rights, Gene_Track_List *last_left, Gene_Track_List *new_right, BOOL do_link)
{
	int temp=-1;

	return(make_joins(lefts, rights, last_left, new_right, do_link, temp));
}



int Genome_Track::make_joins(Track_Stack *lefts, Track_Stack *rights, Gene_Track_List *last_left, Gene_Track_List *new_right, BOOL do_link, int &save_index)
{
	int i, num_joins=0, join, runlen;
	Gene_Track_List *curr_left, *curr_right, *valid_left, *valid_right;

	set_link_counts(lefts, rights);
		
	lefts->reset();
	//Join
	while((lefts->get_size()>0) && (num_joins<2)) {
			curr_left=lefts->get_bottom();
			rights->reset();
			join=0;
			curr_right=rights->get_next();
			while((curr_right != 0) && (join == 0)) {			
				//Catches the possiblity of more that one join for two ends--we don't want to join these two
				//because it precludes a second join elsewhere
				if ((curr_left->num_joins == 2) && (curr_right->num_joins == 2))
					join=0;
				else 
					join=assign_adjacency(curr_left, curr_right, do_link, save_index);

				if (join == 1) {
				//We made a join
					rights->delete_element(curr_right);
					lefts->pop();
					num_joins++;
					
					//Check left stack--If our join "crossed" a locus that will be the only valid remaining end
					if (last_left->my_locus->has_duplicate() == FALSE)
					{
						if (curr_left != last_left) {
							valid_left=last_left;
							lefts->reset();
							runlen=lefts->get_size();
							for(i=0; i<runlen; i++) lefts->pop();
							lefts->add_back(valid_left);
							//lefts->clear();
							//lefts->push(valid_left);
						}

					}
					//Check right stack
					if (new_right->my_locus->has_duplicate() == FALSE)
					{
						if (curr_right != new_right) {
							valid_right=new_right;
							rights->reset();
							runlen=rights->get_size();
							for(i=0; i<runlen; i++)  rights->pop();
							rights->add_back(valid_right);
							//rights->clear();
							//rights->push(valid_right);
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



int Genome_Track::find_partition_breaks(Tracking_List *left_end, Tracking_List *start, int num_sections_left, int *left_mins, int *left_maxs)
{
	int i, num_breaks=0;;
	BOOL stop=FALSE, do_break;
	Tracking_List *curr_pos;

	curr_pos=left_end->last;

	while (!stop) {
		
		if (curr_pos != left_end->last) {
			if ((curr_pos->next->element[0]->last != curr_pos->element[0]) && (curr_pos->next->element[0]->last != 0))
				stop=TRUE;
		}

		if (curr_pos->element[0]->next != 0) {
			do_break = TRUE;
			for(i=0; i<num_sections_left; i++) { 
				if ((curr_pos->element[0]->next->partition_num >= left_mins[i]) && 
					(curr_pos->element[0]->next->partition_num <= left_maxs[i])) 
						do_break=FALSE;
			}
			if (do_break == TRUE) {
				broken_lefts[num_breaks]=curr_pos->element[0];
				broken_rights[num_breaks]=curr_pos->element[0]->next;
				curr_pos->element[0]->next->last=0;
				curr_pos->element[0]->next=0;
				num_breaks++;
			}
		}
		if ((curr_pos->element[1]->my_locus != 0) && (curr_pos->element[1]->next != 0)) {
			do_break=TRUE;
			for(i=0; i<num_sections_left; i++) {
				if ((curr_pos->element[1]->next->partition_num >= left_mins[i]) && 
					(curr_pos->element[1]->next->partition_num <= left_maxs[i]))
						do_break=FALSE;
			}
			if (do_break == TRUE) {
				broken_lefts[num_breaks]=curr_pos->element[1];
				broken_rights[num_breaks]=curr_pos->element[1]->next;
				curr_pos->element[1]->next->last=0;
				curr_pos->element[1]->next=0;
				num_breaks++;
			}
		}

		if (curr_pos->element[0]->my_locus->has_duplicate() == TRUE)
			stop=TRUE;

		curr_pos=curr_pos->last;
		
		if ((num_breaks==2) || (curr_pos == start->last))
			stop=TRUE;
	}
	return(num_breaks);
}


int Genome_Track::find_partition_breaks(int break_before)
{
	int i, j, num_breaks=0;
	i=break_before-1;
	while((i>=0) &&(num_breaks < 2)) {
		for(j=0; j<2; j++) {
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


void Genome_Track::reset_partition_numbers()
{
	int i;

	for(i=0; i<the_homologs->get_num_homologs(); i++)
	{
		inferred_tracking[i][0].partition_num=i;
		inferred_tracking[i][1].partition_num=i;
	}
}

void Genome_Track::reset_list_partition_numbers()
{
	Tracking_List *use_pos;
	use_pos=partial_track_start;

	use_pos->element[0]->partition_num=0;
	use_pos->element[1]->partition_num=0;

	use_pos=use_pos->next;
	while (use_pos != 0) {
		use_pos->element[0]->partition_num=use_pos->last->element[0]->partition_num+1;
		use_pos->element[1]->partition_num=use_pos->last->element[1]->partition_num+1;
		use_pos=use_pos->next;
	}

}

void Genome_Track::find_possible_internal_rejoins(int num_lost_breaks, int old_lost_breaks, int break_num, int break_sec_num[2][3], 
	  int num_check_breaks[2][3], Tracking_List *start_pos)
{
	int i,j;
	Tracking_List *find_pos;


	num_check_breaks[0][break_num]=0;
	num_check_breaks[1][break_num]=0;
	for(i=0; i<num_lost_breaks; i++) {
		find_pos=start_pos;
		while((find_pos->element[0] != old_lasts[i+old_lost_breaks]) && 
			(find_pos->element[1] != old_lasts[i+old_lost_breaks])) find_pos=find_pos->last;

		if (find_pos != start_pos) {
			num_check_breaks[0][break_num]=1;
			check_break_pos[0][0][break_num]=find_pos;
			break_sec_num[0][break_num]=break_num;
			
		}

		find_pos=start_pos->next;
		while((find_pos->element[0] != old_nexts[i+old_lost_breaks]) && 
			(find_pos->element[1] != old_nexts[i+old_lost_breaks])) find_pos=find_pos->next;

		if (find_pos != start_pos->next) {
			check_break_pos[1][0][break_num]=find_pos->last;
			break_sec_num[1][break_num]=break_num+1;
			num_check_breaks[1][break_num]=1;
		}
			
	}

}


void Genome_Track::find_reversal_internal_rejoins(int break_num, int num_check_breaks[2][3], 
		Tracking_List *sectionA_end, Tracking_List *sectionB_start)
{
	Tracking_List *find_pos;

	if ((num_check_breaks[0][break_num] != 0) && (sectionA_end != 0)) {
		find_pos=sectionA_end;
		while(find_pos->element[0]->my_locus != check_break_pos[0][0][break_num]->element[0]->my_locus) find_pos=find_pos->last;

		check_break_pos[0][1][break_num]=find_pos->last;
	}

	if ((num_check_breaks[1][break_num] != 0) && (sectionB_start != 0)) {
		find_pos=sectionB_start;
		while(find_pos->element[0]->my_locus != check_break_pos[1][0][break_num]->next->element[0]->my_locus) 
			find_pos=find_pos->next;

		check_break_pos[1][1][break_num]=find_pos;
	}
	
}


int Genome_Track::set_section_min_maxs(int pattern, int *rev, int *section_pos, int *pos, int break_num, int break_side,
						  int break_sec_num[2][3])
{
	int i, num_sections;
	
	num_sections=0;
	for(i=0; i<=section_pos[break_sec_num[break_side][break_num]]; i++) {
		if (rev[i] == 0) {
			section_mins[i]=section_starts[pos[i]][0]->element[0]->partition_num;
			section_maxs[i]=section_ends[pos[i]][0]->element[0]->partition_num;
		}
		else {
			section_maxs[i]=section_starts[pos[i]][1]->element[0]->partition_num;
			section_mins[i]=section_ends[pos[i]][1]->element[0]->partition_num;
		}
		num_sections++;
	}

	if (rev[section_pos[break_sec_num[break_side][break_num]]]== 0) 
		section_maxs[section_pos[break_sec_num[break_side][break_num]]]
			=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
	else 
		section_mins[section_pos[break_sec_num[break_side][break_num]]]
			=check_break_pos[break_side][1][break_num]->element[0]->partition_num;
	

	return(num_sections);

}


int Genome_Track::set_section_min_maxs(int break_num, int break_side, int rev, int break_sec_num[2][3])
{
	int num_sections, section1_end;

	if (break_sec_num[break_side][break_num] == 0) 
		section1_end=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
	else
		section1_end=section_ends[0][0]->element[0]->partition_num;

	if (section1_end > precut_max_partition) {
		num_sections=break_sec_num[break_side][break_num]+2;
			section_mins[0]=section_starts[0][0]->element[0]->partition_num;
			section_maxs[0]=precut_max_partition;
			section_mins[1]=postcut_min_partition;
		
		if (break_sec_num[break_side][break_num] == 0) 
			section_maxs[1]=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
		else {
			section_maxs[1]=section_ends[0][0]->element[0]->partition_num;
			if (rev ==0) {
				section_mins[2]=section_starts[1][rev]->element[0]->partition_num;
				if (break_sec_num[break_side][break_num] == 1) 
					section_maxs[2]=check_break_pos[break_side][rev][break_num]->element[0]->partition_num;
				else {
					//There are actually 3 sections ahead of the break
					section_maxs[2]=section_ends[1][rev]->element[0]->partition_num;
					section_mins[3]=section_starts[2][0]->element[0]->partition_num;
					section_maxs[3]=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
				}
			}
			else {
				section_maxs[2]=section_starts[1][rev]->element[0]->partition_num;
				if (break_sec_num[break_side][break_num] == 1) 
					section_mins[2]=check_break_pos[break_side][rev][break_num]->element[0]->partition_num;
				else {
					//There are actually 3 sections ahead of the break
					section_mins[2]=section_ends[1][rev]->element[0]->partition_num;
					section_mins[3]=section_starts[2][0]->element[0]->partition_num;
					section_maxs[3]=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
				}

			}
		}
	}
	else {
		num_sections=break_sec_num[break_side][break_num]+1;

		section_mins[0]=section_starts[0][0]->element[0]->partition_num;
		
		if (break_sec_num[break_side][break_num] == 0) 
			section_maxs[0]=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
		else {
			section_maxs[0]=section_ends[0][0]->element[0]->partition_num;

			if (rev == 0) {
				section_mins[1]=section_starts[1][rev]->element[0]->partition_num;
				if (break_sec_num[break_side][break_num] == 1) 
					section_maxs[1]=check_break_pos[break_side][rev][break_num]->element[0]->partition_num;
				else {
					//There are actually 3 sections ahead of the break
					section_maxs[1]=section_ends[1][rev]->element[0]->partition_num;
					section_mins[2]=section_starts[2][0]->element[0]->partition_num;
					section_maxs[2]=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
				}
			}
			else {
				section_maxs[1]=section_starts[1][rev]->element[0]->partition_num;
				if (break_sec_num[break_side][break_num] == 1) 
					section_mins[1]=check_break_pos[break_side][rev][break_num]->element[0]->partition_num;
				else {
					//There are actually 3 sections ahead of the break
					section_mins[1]=section_ends[1][rev]->element[0]->partition_num;
					section_mins[2]=section_starts[2][0]->element[0]->partition_num;
					section_maxs[2]=check_break_pos[break_side][0][break_num]->element[0]->partition_num;
				}

			}
		}
	}
	return(num_sections);
}



void Genome_Track::reverse_section(Tracking_List *start, Tracking_List *end, 
								   Tracking_List *&new_start, Tracking_List *&new_end)
{
	int i, max_size, min_partition, max_partition;
	Tracking_List *curr, *new_list;

	curr=start;
	min_partition=start->element[0]->partition_num;
	max_partition=end->element[0]->partition_num;

	//We need to index the pointers in the original list so we can set the reversed list
	i=0;
		
	do {
		curr->element[0]->mark_num=i++;
		curr->element[1]->mark_num=i++;

		curr=curr->next;
	} while (curr != end->next);

	max_size=i;
	i=max_size-1;
	new_start=hold_reversal_list[reversal_pos++];
	new_list=new_start;
	new_list->last=0;
	reversal_pointers[i--]=new_list->element[1];
	reversal_pointers[i--]=new_list->element[0];

	new_list->element[0]->next=new_list->element[0]->last=0;
	new_list->element[1]->next=new_list->element[1]->last=0;
	
	while(i > 0) {
		new_list->next=hold_reversal_list[reversal_pos++];
		new_list->next->last=new_list;
		new_list=new_list->next;
		reversal_pointers[i--]=new_list->element[1];
		reversal_pointers[i--]=new_list->element[0];
		new_list->element[0]->next=new_list->element[0]->last=0;
		new_list->element[1]->next=new_list->element[1]->last=0;
	}

	new_list=new_start;

	curr=end;

	do {
		for(i=0; i<2; i++) {
			new_list->element[i]->my_locus=curr->element[i]->my_locus;
			new_list->element[i]->index_num=curr->element[i]->index_num;
	
	
			if ((curr->element[i]->next != 0) && 
				(curr->element[i]->next->partition_num >= min_partition) && 
				(curr->element[i]->next->partition_num <= max_partition)) 
				new_list->element[i]->last=reversal_pointers[curr->element[i]->next->mark_num];
			else
				new_list->element[i]->last=0;

			if ((curr->element[i]->last != 0) && 
				(curr->element[i]->last->partition_num >= min_partition) &&
				(curr->element[i]->last->partition_num <= max_partition)) 
				new_list->element[i]->next=reversal_pointers[curr->element[i]->last->mark_num];
			else
				new_list->element[i]->next=0;

			new_list->element[i]->partition_num=curr->element[i]->partition_num;
		}
		new_list->num=curr->num;
		
		if (curr != start) 
			new_list=new_list->next;
		curr=curr->last;
		//}
	} while (curr != start->last);
	new_end=new_list;
	new_end->next=0;

}


int Genome_Track::assign_adjacency(Gene_Track_List *locus1, Gene_Track_List *locus2) 
{
	int temp=0;
	return(assign_adjacency(locus1, locus2, TRUE, temp));
}

int Genome_Track::assign_adjacency(Gene_Track_List *locus1, Gene_Track_List *locus2, BOOL do_link)
{
	int temp=0;
	return(assign_adjacency(locus1, locus2, do_link, temp));
}


int Genome_Track::assign_adjacency(Gene_Track_List *locus1, Gene_Track_List *locus2, BOOL do_link, int &save_index)
{
	int retval=0;
	if ((locus1->my_locus != 0) && (locus2->my_locus != 0)) {
		if (locus1->my_locus->get_contig(locus1->index_num) == locus2->my_locus->get_contig(locus2->index_num)) {
			if ((locus1->my_locus->get_gene(locus1->index_num) == locus2->my_locus->get_gene(locus2->index_num)+1) ||
				(locus1->my_locus->get_gene(locus1->index_num) == locus2->my_locus->get_gene(locus2->index_num)-1)) {
					{
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

	return(retval);

}







void Genome_Track::store_tracking( int size)
{
	int i, j, track_cnt[2], track_num, other_track_num;
	Gene_Track_List *curr_track, *last;
	Track_Stack the_runs;

	new_inferred_tracking = new Gene_Track_List* [the_homologs->get_num_homologs()];
	
	for(i=0; i<the_homologs->get_num_homologs(); i++) 
		new_inferred_tracking[i]=new Gene_Track_List [2];
	
	for(i=0; i<the_homologs->get_num_homologs(); i++) {
		new_inferred_tracking[i][0].my_locus=0;
		new_inferred_tracking[i][1].my_locus=0;
	}

	
	track_cnt[0]=track_cnt[1]=0;
	inferred_tracking[0][0].to_track_num=0;
	the_runs.push(&inferred_tracking[0][0]);
	curr_track=0;


	while(the_runs.get_size() > 0) {
		if (curr_track == 0)
		{
			curr_track=the_runs.get_bottom();
			track_num=curr_track->to_track_num;

			//cout<<track_cnt[track_num]<<endl<<flush;
			while ((&inferred_tracking[track_cnt[track_num]][0] != curr_track) && 
				(&inferred_tracking[track_cnt[track_num]][1] != curr_track)) track_cnt[track_num]++;

			if (track_num ==0)
				other_track_num=1;
			else
				other_track_num=0;
			last=0;
		}
	
		new_inferred_tracking[track_cnt[track_num]][track_num].my_locus=curr_track->my_locus;
		new_inferred_tracking[track_cnt[track_num]][track_num].index_num=curr_track->index_num;
		if (last !=0)
		{
			new_inferred_tracking[track_cnt[track_num]][track_num].last=last;
			last->next=&new_inferred_tracking[track_cnt[track_num]][track_num];
		}
		
		if ((curr_track->my_locus->has_duplicate() == TRUE) &&
			((track_cnt[track_num] > track_cnt[other_track_num]) || ((track_num == 0)&&(track_cnt[other_track_num]==0)))) {
			//If this locus is duplicated we may have a second track
			if ((curr_track->index_num == 0) && (inferred_tracking[track_cnt[track_num]][1].last == 0)) {
				inferred_tracking[track_cnt[track_num]][1].to_track_num=other_track_num;

				the_runs.push(&inferred_tracking[track_cnt[track_num]][1]);
			}
			else
				if (inferred_tracking[track_cnt[track_num]][0].last == 0) {
					inferred_tracking[track_cnt[track_num]][0].to_track_num=other_track_num;

					the_runs.push(&inferred_tracking[track_cnt[track_num]][0]);
				}
		}

		if (curr_track->next != 0) {
			last=&new_inferred_tracking[track_cnt[track_num]][track_num];
			track_cnt[track_num]++;	

			while((&inferred_tracking[track_cnt[track_num]][0] != curr_track->next) && 
				(&inferred_tracking[track_cnt[track_num]][1] != curr_track->next)) {
				if ((inferred_tracking[track_cnt[track_num]][0].last == 0) &&
					(track_cnt[track_num] >= track_cnt[other_track_num])) {
					inferred_tracking[track_cnt[track_num]][0].to_track_num=other_track_num;
					//If we don't have one, make this gene the start of the other track
					the_runs.push(&inferred_tracking[track_cnt[track_num]][0]);
				}
				track_cnt[track_num]++;						
			}
			
		}
		else {
			//We've finished with this run--kill it
			the_runs.pop();

			if (the_runs.get_size()==0) {
				//We have a double strand break--get the next position
				if (track_cnt[0] < track_cnt[1]) 
					track_cnt[0]=track_cnt[1];
				track_cnt[0]++;

				if (track_cnt[0] < size) {
					inferred_tracking[track_cnt[0]][0].to_track_num=0;
					the_runs.push(&inferred_tracking[track_cnt[0]][0]);
				}
			}
			

		}

		curr_track=curr_track->next;
	
	}

	last = 0;
	for(i=size-1; i>=0; i--) {
		if (new_inferred_tracking[i][0].my_locus != 0) {
			if (new_inferred_tracking[i][0].last != 0)
				last=new_inferred_tracking[i][0].last;
			else
				last = 0;
		}
		else if (last !=0) {
			new_inferred_tracking[i][0].last=last;
		}

	}


	last = 0;
	for(i=size-1; i>=0; i--) {
		if (new_inferred_tracking[i][1].my_locus != 0) {
			if (new_inferred_tracking[i][1].last != 0)
				last=new_inferred_tracking[i][1].last;
			else
				last = 0;
		}
		else if (last !=0) {
			new_inferred_tracking[i][1].last=last;
		}

	}

	for(i=0; i<size; i++)
	{
		if (new_inferred_tracking[i][0].last != 0)
		{
			j=i-1;
			while(new_inferred_tracking[i][0].last != &new_inferred_tracking[j][0]) j--;
			new_inferred_tracking[i][0].dist_to_last=i-j;
		}

		if (new_inferred_tracking[i][1].last != 0)
		{
			j=i-1;
			while(new_inferred_tracking[i][1].last != &new_inferred_tracking[j][1]) j--;
			new_inferred_tracking[i][1].dist_to_last=i-j;
		}
	}

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		delete[] inferred_tracking[i];
	delete[] inferred_tracking;

	inferred_tracking=new_inferred_tracking;
	


}


void Genome_Track::set_null_lasts()
{
	int i, j, end;
	Gene_Track_List *curr=0;

	for(j=0; j<2; j++) {
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

Gene_Track_List ** Genome_Track::create_new_track_locus(WGD_Locus *the_locus)
{
	Gene_Track_List **new_locus_data;

	new_locus_data=new Gene_Track_List *[2];
	new_locus_data[0]=new Gene_Track_List();
	new_locus_data[1]=new Gene_Track_List();

	new_locus_data[0]->my_locus=the_locus;
	new_locus_data[0]->index_num=0;

	if (the_locus->has_duplicate() == TRUE) {
		new_locus_data[1]->my_locus=the_locus;
		new_locus_data[1]->index_num=1;
	}
	return(new_locus_data);
}


BOOL Genome_Track::check_list()
{
	int i;
	BOOL valid=TRUE;
	Tracking_List *curr;
	curr=partial_track_start;
	

	if (partial_track_start->last !=0)
		valid=FALSE;
	if(partial_track_end->next != 0)
		valid=FALSE;

	for(i=1; i<list_len; i++) {
		curr=curr->next;

		if ((curr->next !=0) && (curr->next->last != curr))
			valid=FALSE;

		if (curr->element[0]->my_locus ==0)
			valid=FALSE;
	}

	if (curr != partial_track_end)
		valid=FALSE;

	return(valid);
}


WGD_Tracks::WGD_Tracks(WGD_Data *homologs, Clade *genomes)
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

	the_trackings=new Genome_Track * [the_genomes->get_num_genomes()];

	for(i=0; i<the_genomes->get_num_genomes(); i++) 
		the_trackings[i]=new Genome_Track (i, &(*the_genomes)[i], the_homologs);
	
}
	

int WGD_Tracks::get_num_breaks()
{
	int i, retval=0;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		retval+=the_trackings[i]->count_num_breaks();

	return(retval);

	
}

int WGD_Tracks::get_num_list_breaks()
{
	int i, retval=0;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		retval+=the_trackings[i]->count_num_list_track_breaks();

	return(retval);
}

void WGD_Tracks::change_order(int *new_order)
{
	int i;

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=new_order[i];

	tracking_correct=FALSE;
}


void WGD_Tracks::change_order()
{
	int i;

	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=the_trackings[0]->get_list_position_number(i+1);

	tracking_correct=FALSE;
}


void WGD_Tracks::update_tracking()
{
	int i;
	
	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_trackings[i]->number_contigs(order);

	tracking_correct=TRUE;
}



void WGD_Tracks::print_all_tracks()
{
	print_all_tracks(TRUE);
}

void WGD_Tracks::print_all_tracks(BOOL use_file)
{
	int i;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_trackings[i]->print_tracking(use_file);
}
	


void WGD_Tracks::set_all_null_lasts()
{
	int i;

	for(i=0; i<the_genomes->get_num_genomes(); i++) 
		the_trackings[i]->set_null_lasts();
}



Genome_Track & WGD_Tracks::operator[] (int index)
{
	if (tracking_correct==FALSE)
		update_tracking();
	
	return((*the_trackings[index]));
}
	


void WGD_Tracks::swap_elements (int index1, int index2)
{
	int temp;

	temp=order[index1];
	order[index1]=order[index2];
	order[index2]=temp;

	tracking_correct=FALSE;

}

void WGD_Tracks::get_move_section_start_end(int move_pos, int &start, int &end)
{
	int i;

	i=0;
	while(move_section_ends[i] < move_locations[move_pos]) i++;

	start=real_to_move_locs[move_section_starts[i]];
	end=real_to_move_locs[move_section_ends[i]];
}



void WGD_Tracks::search_for_best()
{
	int i, j, k, l, best_score, best_pos, best_pos_m, old_best, best_rev, 
		best_pattern, best_loc, **new_scores, max_size, stop_count,
		start, end, zero_level=0, temp;
	BOOL made_improve=TRUE, *full_break;
	Gene_Track_List **new_locus_data;
	Tracking_List **new_loci;

	new_loci =new Tracking_List * [the_genomes->get_num_genomes()];
	new_scores=new int*[the_homologs->get_num_homologs()];
	full_break=new BOOL [the_genomes->get_num_genomes()];
	

	for(i=0; i<the_homologs->get_num_homologs(); i++) 
		new_scores[i]=new int [20];

	for(j=0; j<the_genomes->get_num_genomes(); j++)  
		the_trackings[j]->start_track_list(order[0]);
	

	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		best_score=0;
		best_pos=-1;
		cout<<"Pos "<<i<<"\t";
		
		for(j=0; j<the_genomes->get_num_genomes(); j++)
			new_loci[j]=the_trackings[j]->create_list_element(order[i]);

		for(j=0; j<=i; j++)
			new_scores[j][0] =0;

		for(k=0; k<the_genomes->get_num_genomes(); k++) 
			for(j=0; j<=i; j++) 
				new_scores[j][0]+=the_trackings[k]->add_to_position_n(j, new_loci[k], new_loci[k], 1, FALSE);

		best_pos=0;
		best_score=new_scores[0][0];

		for(j=1; j<=i; j++) {
			if (new_scores[j][0] > best_score) {
			best_score=new_scores[j][0];
			best_pos=j;
			}
		}
		
	//	if (best_score != 0)
			cout<<i<<" Best: "<<best_score<<" pos: "<<best_pos<<endl<<flush;
		
		for(k=0; k<the_genomes->get_num_genomes(); k++) 
			the_trackings[k]->add_to_position_n(best_pos, new_loci[k], new_loci[k], 1, TRUE);

		//for(j=0; j<=i; j++)
		//	cout<<the_trackings[0]->get_nth_tracking_list_pos(j)<<"\t";
		//cout<<endl;
		
	}

	for(k=0; k<the_genomes->get_num_genomes(); k++) 
		the_trackings[k]->reset_list_partition_numbers();
	for(i=1; i<the_homologs->get_num_homologs(); i++) {
		for(k=0; k<the_genomes->get_num_genomes(); k++) 
			the_trackings[k]->break_and_rejoin(i);
	}


	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=the_trackings[0]->get_nth_tracking_list_pos(i);

	update_tracking();
	check_for_fully_connected_positions();
	cout<<"Tracking has "<<get_num_breaks()<<" breaks"<<endl<<flush;
	cout<<"List tracking has "<<get_num_list_breaks()<<" breaks"<<endl;


//	print_all_tracks(FALSE);

#if 0
	while (made_improve == TRUE ) {
		made_improve=FALSE;

		do {	
			best_score = -1;
	
			for(i=1; i<the_homologs->get_num_homologs(); i++) {
				if (fully_connected[i] == FALSE) {
					for(j=0; j<the_homologs->get_num_homologs(); j++)
						for(k=0; k<20; k++)
							new_scores[j][k]=0;

						for(k=0; k<the_genomes->get_num_genomes(); k++) {
							the_trackings[k]->try_position_reorgs(i, new_scores, fully_connected);
						//	cout<<i<<"\t"<<k<<"\t"<<get_num_list_breaks()<<endl;
						//	if (i == 5)
						//		cout<<new_scores[10][5]<<endl;
						}
						
					


					for(k=0; k<3; k++) {
						if (new_scores[0][k] > best_score) {
							best_pos=i;
							best_pos_m=-1;
							best_pattern=k;
							best_score=new_scores[0][k];
						}
			
					}
				
					for(j=2; j<the_homologs->get_num_homologs(); j++) {
						for(k=3; k<20; k++){
							if (new_scores[j][k] > best_score) {
								best_pos=i;
								best_pos_m=j;
								best_pattern=k;
								best_score=new_scores[j][k];
							}
						}
					}
				}

				cout<<i<<": "<<best_score<<endl;
			//	cout<<i<<"\tList breaks: "<<get_num_list_breaks()<<endl;
			}

			if (best_score > 0) {
				made_improve=TRUE;
				cout<<"Prior list breaks: "<<get_num_list_breaks()<<endl;
				zero_level=0;
				for(k=0; k<the_genomes->get_num_genomes(); k++)
					the_trackings[k]->do_reorg(best_pos, best_pos_m, best_pattern, FALSE);
				cout<<"List breaks after fake: "<<get_num_list_breaks()<<endl;

				for(k=0; k<the_genomes->get_num_genomes(); k++)
					the_trackings[k]->do_reorg(best_pos, best_pos_m, best_pattern, TRUE);

				for(k=0; k<the_genomes->get_num_genomes(); k++) 
					the_trackings[k]->reset_list_partition_numbers();
				for(i=1; i<the_homologs->get_num_homologs(); i++) {
					for(k=0; k<the_genomes->get_num_genomes(); k++) 
						the_trackings[k]->break_and_rejoin(i);
				}

				cout<<"Found rearrangement "<<best_pattern<<"\t"<<best_pos<<"\t"<<best_pos_m<<" which improves score by "<<best_score<<endl;
				//for(j=0; j<the_homologs->get_num_homologs(); j++)
				//	cout<<the_trackings[0]->get_nth_tracking_list_pos(j)<<"\t";
				//cout<<endl;

				for(i=0; i<the_homologs->get_num_homologs(); i++)
					order[i]=the_trackings[0]->get_nth_tracking_list_pos(i);

				update_tracking();
				check_for_fully_connected_positions();

				cout<<"Tracking has "<<get_num_breaks()<<" breaks"<<endl<<flush;
				cout<<"List tracking has "<<get_num_list_breaks()<<" breaks"<<endl;

			//	print_all_tracks(FALSE);
			}

		} while (best_score > 0);


		cout<<"Starting inserts: List tracking has "<<get_num_list_breaks()<<" breaks"<<endl;


		
		do {
			best_score=-1;
		
		
			for(i=1; i<the_homologs->get_num_homologs()-2; i++) {
				j=i+1;
				if (fully_connected[i] == FALSE) {
					for(k=0; k<the_genomes->get_num_genomes(); k++)
						full_break[k]=FALSE;
					stop_count = 0;
					while ((j<the_homologs->get_num_homologs()-1) && (stop_count < the_genomes->get_num_genomes())) {
						if (fully_connected[j] == FALSE) {
							for (k=0; k<the_homologs->get_num_homologs(); k++)
								new_scores[k][0]=new_scores[k][1]=0;

						//	cout<<"Old list#: "<<get_num_list_breaks()<<endl;
							for(k=0; k<the_genomes->get_num_genomes(); k++) {
								the_trackings[k]->try_position_inserts(i,j,new_scores, fully_connected);
							//	if (i==6 && j==10)
								//	cout<<i<<"\t"<<j<<"\t"<<new_scores[0][0]<<"\t"<<new_scores[0][1]<<endl;
							//	cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<get_num_list_breaks()<<endl;
							}
							//cout<<"New list#: "<<get_num_list_breaks()<<endl;
							for (k=0; k<the_homologs->get_num_homologs(); k++) {
								if ((new_scores[k][0] > best_score) || (best_score == -1)) {
									best_score=new_scores[k][0];
									best_rev=0;
									best_pos=i;
									best_pos_m=j;
									best_loc=k;
								}
								if ((new_scores[k][1] > best_score) || (best_score == -1)) {
									best_score=new_scores[k][1];
									best_rev=1;
									best_pos=i;
									best_pos_m=j;
									best_loc=k;
								}
							}
						}
						for(k=0; k<the_genomes->get_num_genomes(); k++) {
							if (full_break[k] == FALSE) {
								if (the_trackings[k]->has_double_break(j)==TRUE) {
									stop_count++;
									full_break[k] = TRUE;
								}
							}
						}
						j++;
					}
					
				}
				cout<<"i: "<<i<<"jmax: "<<j<<" Best: "<<best_score<<endl;
			}
				
			if (best_score > 0) {
				made_improve=TRUE;
				
				cout<<"Prior list tracking has "<<get_num_list_breaks()<< " breaks\n";
				for(k=0; k<the_genomes->get_num_genomes(); k++) {
					the_trackings[k]->do_section_insert(best_pos, best_pos_m, best_loc+1, best_rev, FALSE);
					cout<<k<<": "<<get_num_list_breaks()<<endl;
				}
				cout<<"After test list breaks: "<<get_num_list_breaks()<<endl;

				for(k=0; k<the_genomes->get_num_genomes(); k++)
					the_trackings[k]->do_section_insert(best_pos, best_pos_m, best_loc+1, best_rev, TRUE);

				for(k=0; k<the_genomes->get_num_genomes(); k++) 
					the_trackings[k]->reset_list_partition_numbers();
				for(i=1; i<the_homologs->get_num_homologs(); i++) {
					for(k=0; k<the_genomes->get_num_genomes(); k++) 
						the_trackings[k]->break_and_rejoin(i);
				}


				cout<<"Found insert "<<best_loc<<"\t"<<best_pos<<"\t"<<best_pos_m<<"\t"<<best_rev<<" which improves score by "<<best_score<<endl;
				for(i=0; i<the_homologs->get_num_homologs(); i++)
					order[i]=the_trackings[0]->get_nth_tracking_list_pos(i);

				update_tracking();
				check_for_fully_connected_positions();

				cout<<"Tracking has "<<get_num_breaks()<<" breaks"<<endl<<flush;
				cout<<"List tracking has "<<get_num_list_breaks()<< " breaks\n";
			//	print_all_tracks(FALSE);	
			}
		}while(best_score > 0);

	

	}
#endif

	if (the_homologs->get_num_homologs()/2 ==0)
		max_size=the_homologs->get_num_homologs()/2;
	else
		max_size=the_homologs->get_num_homologs()/2+1;



	for(i=0; i<the_homologs->get_num_homologs(); i++)
		order[i]=the_trackings[0]->get_nth_tracking_list_pos(i);

	update_tracking();

	for(i=0; i<the_homologs->get_num_homologs(); i++) 
				delete[] new_scores[i];
	
	delete[] new_scores;
	delete[] full_break;

	for(i=0; i<the_genomes->get_num_genomes(); i++)
		the_trackings[i]->reset_internal_track_list();

	
}


int WGD_Tracks::section_insert(int start_pos, int end_pos, int rev, int loc, BOOL do_link, BOOL real_nums)
{
	int i, loc_val, retval=0;

	if (real_nums == FALSE) {
		if (loc > start_pos) 
			loc_val=move_locations[loc+(end_pos-start_pos)]-(move_locations[end_pos]-move_locations[start_pos]);
		else
			loc_val=move_locations[loc];

		for(i=0; i<the_genomes->get_num_genomes(); i++) 
			retval+=the_trackings[i]->do_section_insert(move_locations[start_pos], move_locations[end_pos], 
			loc_val, rev, do_link);
	}
	else {
		for(i=0; i<the_genomes->get_num_genomes(); i++) 
			retval+=the_trackings[i]->do_section_insert(start_pos, end_pos, loc, rev, do_link);

	}

	if (do_link == TRUE) {
		update_tracking();
		check_for_fully_connected_positions();
	}
	return(retval);
}




WGD_Tracks::~WGD_Tracks()
{
	delete[] the_trackings;
	delete[] move_locations;
	delete[] move_section_starts;
	delete[] move_section_ends;
	delete[] real_to_move_locs;
}



void WGD_Tracks::set_array_for_pos(int pos, int val, int size, int *array)
{
	int i;

	for(i=size-1; i>pos; i--) 
		array[i]=array[i-1];

	array[pos]=val;
}



void WGD_Tracks::check_for_fully_connected_positions()
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
			if (the_trackings[j]->has_tracked_double_break(i) == FALSE)
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
