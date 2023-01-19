#include "anneal_template.cpp"
#include <fstream>
#include <map>
#include "WGX_scaffold_anneal.h"


//#include <sys/time.h>
//#include <sys/resource.h>

//#define __SHOW_STEPS__

//unsigned long RunTime();  
using namespace::std;


extern Homolog_Set *global_homologs;



BOOL are_neighbors(Scaffold_Gene* gene1, Scaffold_Gene* gene2)
{
    BOOL retval=FALSE;
    
    if ((gene1->get_neighbor(0) == gene2) || (gene1->get_neighbor(1) == gene2)) retval=TRUE;
    return(retval);
}



Scaffold_Gene::Scaffold_Gene()
{
    int i;
    
    gene_num=0;
    name="";
    neighbors=new Scaffold_Gene*[2];
    for(i=0; i<2; i++) neighbors[i]=0;
    num_neighbors=0;
    keep=TRUE;
    used=TRUE;
    pillar_num=-1;
    
    tandem=FALSE;
    tandem_names=0;
    num_tandem=0;
    num_homologs=0;
}


Scaffold_Gene::Scaffold_Gene (int num, string new_name, int cnum) 
{
    int i;
    
    used=TRUE;
    contig_num=cnum;
    
    gene_num=num;
    name=new_name;
    neighbors=new Scaffold_Gene* [2];
    for(i=0; i<2; i++) neighbors[i]=0;
    num_neighbors=0;
    keep = TRUE;
    
    tandem=FALSE;
    tandem_names=0;
    num_tandem=0;
    num_homologs=0;
    pillar_num=-1;
}


Scaffold_Gene& Scaffold_Gene::operator=(Scaffold_Gene &assign_from)
{
    int i;
    
    gene_num=assign_from.get_gene_num();
    name=assign_from.get_name();
    used=assign_from.gene_used();
    contig_num=assign_from.get_contig_num();
    pillar_num=assign_from.pillar_num;
    
    num_homologs=assign_from.get_num_homologs();
    
    if (tandem_names !=0) delete[] tandem_names;
    
    num_tandem=assign_from.get_num_tandems();
    if (num_tandem >0) {
        tandem_names=new string[num_tandem];
        
        tandem=TRUE;
        for(i=0; i<num_tandem; i++) tandem_names[i]=assign_from.get_nth_tandem(i);
    }
    else {tandem=FALSE;}
    
    return(*this);
}



Scaffold_Gene * Scaffold_Gene::get_neighbor(int index)
{
    if (index == 0)
        return(neighbors[0]);
    else
        return(neighbors[1]);
}


void Scaffold_Gene::set_neighbor(Scaffold_Gene * new_neighbor, int index)
{
    Scaffold_Gene *old_neigh;
    
    if (index == 0) {
        old_neigh=neighbors[0];
        neighbors[0]=new_neighbor;
    }
    else{
        old_neigh=neighbors[1];
        neighbors[1]=new_neighbor;
    }
    if ((old_neigh==0) && (new_neighbor !=0))
        num_neighbors++;
    else {
        if ((old_neigh!=0) && (new_neighbor ==0))
            num_neighbors--;
    }
}

void Scaffold_Gene::add_tandem(string new_tandem)
{
    int i;
    string *new_names;
    
    if (num_tandem ==0) tandem=TRUE;
    
    new_names=new string[num_tandem+1];
    for(i=0; i<num_tandem; i++)
        new_names[i]=tandem_names[i];
    num_tandem++;
    
    new_names[num_tandem-1]=new_tandem;
    delete[] tandem_names;
    
    tandem_names=new_names;
}



Scaffold_Gene::~Scaffold_Gene()
{
    delete[] neighbors;
}

Scaffold_Genome::Scaffold_Genome()
{
    null_gene=new Scaffold_Gene();
    gene_list=0;
    num_genes=0;
}


Scaffold_Genome::Scaffold_Genome(string gname, int ngenes)
{
    int i;
    
    null_gene=new Scaffold_Gene();
    num_genes=ngenes;
    gene_list=new Scaffold_Gene*[num_genes];
    
    genome_name=gname;
    
    for(i=0; i<num_genes; i++) {
        gene_list[i]=new Scaffold_Gene(i, "Null",0);
    }
}


Scaffold_Genome::Scaffold_Genome(string gname, int ngenes, Scaffold_Gene **in_genelist)
{
    int i;
    
    
    null_gene=new Scaffold_Gene();
    num_genes=ngenes;
    gene_list=new Scaffold_Gene*[num_genes];
    
    genome_name=gname;
    
    for(i=0; i<num_genes; i++) {
        gene_list[i]=new Scaffold_Gene();
        (*gene_list[i])=(*in_genelist[i]);
    }
    
    for(i=0; i<num_genes; i++) {
        if (i!= (num_genes-1)) {
            if (gene_list[i]->get_contig_num() == gene_list[i+1]->get_contig_num()) {
                gene_list[i]->set_neighbor(gene_list[i+1], 1);
            }
        }
        if (i!=0) {
            if (gene_list[i]->get_contig_num() == gene_list[i-1]->get_contig_num()) {
                gene_list[i]->set_neighbor(gene_list[i-1], 0);
            }
        }
    }
    
}

Scaffold_Gene& Scaffold_Genome::operator[] (int index_num)
{
    if ((0<=index_num) && (index_num<num_genes))
        return(*gene_list[index_num]);
    else
        return(*null_gene);
}

Scaffold_Genome& Scaffold_Genome::operator= (Scaffold_Genome &assign_from)
{
    int i, neighbor_index;


    genome_name=assign_from.get_genome_name();
    
    if (gene_list !=0) {
        for(i=0; i<num_genes; i++) delete gene_list[i];
        delete[] gene_list;
    }
    
    num_genes=assign_from.get_num_genes();
    gene_list=new Scaffold_Gene * [num_genes];
    
    for(i=0; i<num_genes; i++) {
        gene_list[i]=new Scaffold_Gene();
        (*gene_list[i])=assign_from[i];
    }
    
    for(i=0; i<num_genes; i++)  {
        if (assign_from[i].gene_used() == FALSE) gene_list[i]->set_unused();
        
        if (assign_from[i].get_neighbor(0) !=0) {
            neighbor_index=assign_from[i].get_neighbor(0)->get_gene_num();
            gene_list[i]->set_neighbor(gene_list[assign_from[i].get_neighbor(0)->get_gene_num()],0);
        }
        if (assign_from[i].get_neighbor(1) !=0) {
            neighbor_index=assign_from[i].get_neighbor(1)->get_gene_num();
            gene_list[i]->set_neighbor(gene_list[assign_from[i].get_neighbor(1)->get_gene_num()],1);
        }
    }
    return(*this);
}

#if 0
void Scaffold_Genome::use_gene(int gene_num)
{
    int left_loc, right_loc;
    Scaffold_Gene *curr_gene, *left_neighbor, *right_neighbor;
    if (gene_list[gene_num]->gene_used() == FALSE) {
        curr_gene=gene_list[gene_num];
        
        curr_gene->set_used(-1);
        left_loc=gene_num-1;
        
        if (gene_num!=0) {
            left_neighbor=gene_list[left_loc];
            while((left_neighbor->gene_used()==FALSE) && (left_loc>0) && (left_neighbor->get_contig_num() == curr_gene->get_contig_num())) {
                left_loc--;
                left_neighbor=gene_list[left_loc];
            }
                  if (curr_gene->get_contig_num() != left_neighbor->get_contig_num()) left_neighbor=0;
            if ((left_loc ==0) &&(left_neighbor->gene_used() ==FALSE)) left_neighbor=0;
        }
        else
            left_neighbor=0;
        
        
        
        curr_gene->set_neighbor(left_neighbor, 0);
        if (left_neighbor!=0) left_neighbor->set_neighbor(curr_gene, 1);
        
        right_loc=gene_num+1;
        
        if (gene_num<(num_genes-1)) {
            right_neighbor=gene_list[right_loc];
            
            while((right_neighbor->gene_used()==FALSE) && (right_loc<(num_genes-1)) && (right_neighbor->get_contig_num() == curr_gene->get_contig_num())) {
                right_loc++;
                right_neighbor=gene_list[right_loc];
            }
            
            if (curr_gene->get_contig_num() != right_neighbor->get_contig_num()) right_neighbor=0;
            if ((right_loc ==(num_genes-1)) &&(right_neighbor->gene_used() ==FALSE)) right_neighbor=0;
            
        }
        else
            right_neighbor=0;
        
        
        
        curr_gene->set_neighbor(right_neighbor, 1);
        if (right_neighbor!=0) right_neighbor->set_neighbor(curr_gene, 0);
    }
}

#endif

void Scaffold_Genome::print_tandems (string filename)
{
    int my_gene, tandem;
    ofstream fout;
    Scaffold_Gene *curr_gene;
    
    fout.open(filename.c_str());
    
    if (fout.fail()) {
        cerr<<"ERROR: Cannot open tandem output file "<<filename<<endl;
        return;
    }
    else {
        for(my_gene=0; my_gene<num_genes; my_gene++) {
            curr_gene=gene_list[my_gene];
            if (curr_gene->get_num_tandems() > 0) {
                fout<<curr_gene->get_name()<<"\t"<<curr_gene->get_num_tandems();
                for(tandem =0; tandem< curr_gene->get_num_tandems(); tandem++) {
                    fout<<"\t"<<curr_gene->get_nth_tandem(tandem);
                }
                fout<<endl;
            }
        }
        fout.close();
    }
}

void Scaffold_Genome::omit_gene(int gene_num)
{
    Scaffold_Gene *curr_gene, *old_left, *old_right;
    
    if (gene_list[gene_num]->gene_used() == TRUE) {
        curr_gene=gene_list[gene_num];
        old_left=curr_gene->get_neighbor(0);
        old_right=curr_gene->get_neighbor(1);
        
        curr_gene->set_unused();
        
        if ((old_left !=0) && (old_right !=0)) {
            old_left->set_neighbor(old_right, 1);
            old_right->set_neighbor(old_left,0);
        }
        else {
            old_left->set_neighbor(0,1);
            old_right->set_neighbor(0,0);
        }
    }
}


void Scaffold_Genome::reset_neighbors()
{
    int my_gene, next_gene;
    
    for(my_gene=0; my_gene<num_genes; my_gene++) {
        gene_list[my_gene]->set_neighbor(0,0);
        gene_list[my_gene]->set_neighbor(0,1);
    }
    
    for(my_gene=0; my_gene<num_genes; my_gene++) {
        
        next_gene=my_gene-1;
        while((next_gene>=0) && (gene_list[next_gene]->gene_used() == FALSE)) next_gene--;
        
        if ((next_gene>=0) && (gene_list[next_gene]->gene_used() == TRUE)) {
            if (gene_list[next_gene]->get_contig_num() == gene_list[my_gene]->get_contig_num())
                gene_list[my_gene]->set_neighbor(gene_list[next_gene], 0);
        }
        
        next_gene=my_gene+1;
        while((next_gene<num_genes) && (gene_list[next_gene]->gene_used() == FALSE)) next_gene++;
        
        if ((next_gene<num_genes) && (gene_list[next_gene]->gene_used() == TRUE)) {
            if (gene_list[next_gene]->get_contig_num() == gene_list[my_gene]->get_contig_num())
                gene_list[my_gene]->set_neighbor(gene_list[next_gene], 1);
        }
        
    }
    
}


Scaffold_Genome::~Scaffold_Genome()
{
    int i;
    
    if (gene_list !=0) {
        for (i=0; i<num_genes; i++) delete gene_list[i];
        delete[] gene_list;
    }
    
    delete null_gene;
}


Scaffold_Genome* read_genome(string genome_file)
{
    int i, genecnt, num_genes=0, contig;
    string genome_name, gene_name;
    Scaffold_Gene **genelist;
    Scaffold_Genome *new_genome;
    ifstream fin;
    
    fin.open(genome_file.c_str());
    if (! fin.fail()) {
        fin>>genome_name;
        while(!(fin.eof())) {
            fin>>contig>>gene_name;
            if (!(fin.eof()))
                num_genes++;
        }
        fin.close();
        fin.clear();
        
        genelist=new Scaffold_Gene*[num_genes];
        fin.open(genome_file.c_str());
        fin>>genome_name;
        genecnt=0;
        while(!(fin.eof())) {
            fin>>contig>>gene_name;
            if (!(fin.eof())) {
                genelist[genecnt] = new Scaffold_Gene(genecnt, gene_name, contig);
                genecnt++;
            }
        }
        fin.close();
        
        new_genome=new Scaffold_Genome(genome_name, num_genes, genelist);
        
        for(i=0; i<num_genes; i++)
            delete genelist[i];
        delete[] genelist;
        
        return(new_genome);
    }
    else {
        cerr<<"ERROR: Invalid genome file "<<genome_file<<endl;
        return(0);
    }
}


Scaffold_Genome* collapse_tandems (string tandem_file, Scaffold_Genome *orig_genome)
{
    int i, j, new_num_genes=0, gene1, gene2, temp, cntnew, check_total=0, cnt=0;
    string name1, name2, tempname;
    BOOL *in_new_genome, has, already_merged;
    map<string, int> name_map, tandem_map;
    ifstream fin;
    Scaffold_Gene **new_genelist;
    Scaffold_Genome *new_genome;
    
    fin.open(tandem_file.c_str());
    
    if (!(fin.fail())) {
    
        in_new_genome=new BOOL[orig_genome->get_num_genes()];

        for(i=0; i<orig_genome->get_num_genes(); i++) {
            name_map[(*orig_genome)[i].get_name_string()]=i;
            tandem_map[(*orig_genome)[i].get_name_string()]=i;
            in_new_genome[i]=TRUE;
        }
        
        while(!(fin.eof())) {
            fin>>name1>>name2;
            cnt++;
            if (! (fin.eof())) {
                gene1=name_map[name1];
                gene2=name_map[name2];
                
                if (((*orig_genome)[gene1].get_neighbor(0) == &(*orig_genome)[gene2]) ||
                    ((*orig_genome)[gene1].get_neighbor(1) == &(*orig_genome)[gene2])) {
                    //cout<<cnt<<": Can merge "<<name1<<" and "<<name2<<" Raw ids are "<<gene1<<", "<<gene2;
                    gene1=tandem_map[name1];
                    gene2=tandem_map[name2];
                    //cout<<" Mapped ids are "<<gene1<<", "<<gene2<<endl;
                    
                    if (gene1 != gene2) {
                        if (gene1 > gene2) {
                            temp=gene1;
                            gene1=gene2;
                            gene2=temp;
                            tempname=name1;
                            name1=name2;
                            name2=tempname;
                        }
                        
                        already_merged=FALSE;
                        for(j=0; j<(*orig_genome)[gene1].get_num_tandems(); j++) {
                            if (name2 == (*orig_genome)[gene1].get_nth_tandem(j))
                                already_merged=TRUE;
                        }
                        
                        if (already_merged == FALSE) {
                            //cout<<"MErging "<<name1<<" and "<<name2<<". Ids are "<<gene1<<" and "<<gene2<<endl;
                            (*orig_genome)[gene1].add_tandem(name2);
                            for(i=0; i<(*orig_genome)[gene2].get_num_tandems(); i++) {
                                has=FALSE;
                                for(j=0; j<(*orig_genome)[gene1].get_num_tandems(); j++) {
                                    if ((*orig_genome)[gene2].get_nth_tandem(i) == (*orig_genome)[gene1].get_nth_tandem(j))
                                        has=TRUE;
                                }
                                
                                if (has == FALSE) {
                                    (*orig_genome)[gene1].add_tandem((*orig_genome)[gene2].get_nth_tandem(i));
                                }
                            }
                            tandem_map[name2]=gene1;
                            in_new_genome[gene2]=FALSE;
                            for(j=0; j<(*orig_genome)[gene1].get_num_tandems(); j++) {
                                //cout<<(*orig_genome)[gene1].get_name()<<" has tandem "<<(*orig_genome)[gene1].get_nth_tandem(j)<<endl;
                                tandem_map[(*orig_genome)[gene1].get_nth_tandem(j)]=gene1;
                            }
                        }
                        
                    }
                }
            }
        }
        fin.close();
        for(i=0; i<orig_genome->get_num_genes(); i++) {
            if (in_new_genome[i] ==TRUE) new_num_genes++;
        }
        
        new_genelist=new Scaffold_Gene*[new_num_genes];
        cntnew=0;
        for(i=0; i<orig_genome->get_num_genes(); i++) {
            if (in_new_genome[i] ==TRUE) {
                new_genelist[cntnew]=new Scaffold_Gene(cntnew, (*orig_genome)[i].get_name_string(),
                                                       (*orig_genome)[i].get_contig_num());
                cntnew++;
            }
        }
        
        
        
        
        new_genome=new Scaffold_Genome(orig_genome->get_genome_name(), new_num_genes, new_genelist);
        new_genome->reset_neighbors();
        
        cntnew=0;
        for(i=0; i<orig_genome->get_num_genes(); i++) {
            if (in_new_genome[i] == TRUE) {
                
                for(j=0; j< (*orig_genome)[i].get_num_tandems(); j++)
                    (*new_genome)[cntnew].add_tandem((*orig_genome)[i].get_nth_tandem(j));
                cntnew++;
            }
        }
        
        for(i=0; i<new_num_genes; i++) delete new_genelist[i];
        delete[] new_genelist;
        
        return(new_genome);
    }
    else {
        cerr<<"ERROR: Invalid tandem file "<<tandem_file<<endl;
        return(orig_genome);
    }
}

Gene_homologs::Gene_homologs()
{
    cerr<<"Invalid call to default constructor of Gene_homologs\n";
    me=0;
    num_homologs=0;
    the_homologs=0;
    the_vals=0;
}

Gene_homologs::Gene_homologs(Scaffold_Gene *to_me)
{
    me=to_me;
    num_homologs=0;
    the_homologs=0;
    the_vals=0;
}


Gene_homologs::Gene_homologs(int nhomologs, Scaffold_Gene* to_me, Scaffold_Gene** homologs, double *vals)
{
    int homolog;
    
    me=to_me;
    num_homologs=nhomologs;
    
    the_homologs=new Scaffold_Gene* [num_homologs];
    the_vals=new double [num_homologs];
    
    for(homolog=0; homolog<num_homologs; homolog++) {
        the_homologs[homolog]=homologs[homolog];
        the_vals[homolog]=vals[homolog];
    }
}


void Gene_homologs::add_homolog(Scaffold_Gene *new_homolog, double val)
{
    int i;
    double *new_vals;
    BOOL found=FALSE;
    Scaffold_Gene **new_list;
    
    for(i=0; i<num_homologs; i++) {
        if (the_homologs[i] == new_homolog)found =TRUE;
    }
    
    if (found == FALSE) {
        num_homologs++;
        new_vals=new double[num_homologs];
        new_list=new Scaffold_Gene*[num_homologs];
        
        for(i=0; i<num_homologs-1; i++) {
            new_vals[i]=the_vals[i];
            new_list[i]=the_homologs[i];
        }
        
        new_vals[num_homologs-1]=val;
        new_list[num_homologs-1]=new_homolog;
        
        delete[] the_vals;
        delete[] the_homologs;
        
        the_vals=new_vals;
        the_homologs=new_list;
    }
}

Scaffold_Gene* Gene_homologs::get_closest_homolog()
{
    int homolog=0, found_index=0;
    double min;
    
    min=the_vals[homolog];
    
    for(homolog=1; homolog<num_homologs; homolog++) {
        if (the_vals[homolog] <min) {
            min=the_vals[homolog];
            found_index=homolog;
        }
    }
    return(the_homologs[found_index]);
}



Gene_homologs::~Gene_homologs()
{
    if (the_homologs!=0) delete[] the_homologs;
    if (the_vals!=0) delete[] the_vals;
    
}


Homolog_Set::Homolog_Set()
{
    ances_genome=new Scaffold_Genome("NULL", 1);
    dupl_genome=new Scaffold_Genome("NULL", 1);
    
    the_homologs=0;
    null_homolog=new Gene_homologs(0);
    genomes_fixed=FALSE;
    has_order=FALSE;
    assigned_order=0;
}


Homolog_Set::Homolog_Set(Scaffold_Genome *agenome, Scaffold_Genome *dgenome)
{
    int i;
    
    ances_genome=new Scaffold_Genome("NULL", 1);
    dupl_genome=new Scaffold_Genome("NULL", 1);
    (*ances_genome)=(*agenome);
    (*dupl_genome)=(*dgenome);
    
    genomes_fixed=FALSE;

    null_homolog=new Gene_homologs(0);
    
    the_homologs =new Gene_homologs*[ances_genome->get_num_genes()];
    
    for(i=0; i<ances_genome->get_num_genes(); i++)
        the_homologs[i]=new Gene_homologs(&(*ances_genome)[i]);
    
    has_order=FALSE;
    assigned_order=0;
}



Gene_homologs& Homolog_Set::operator[] (int index_num)
{
    if ((the_homologs != 0) && (index_num<ances_genome->get_num_genes()) && (index_num>=0))
        return((*the_homologs[index_num]));
    else
        return((*null_homolog));
}

Homolog_Set& Homolog_Set::operator= (Homolog_Set &assign_from)
{
    int i, j;
    double *new_vals;
    Scaffold_Gene **new_homologs;
    
    if (genomes_fixed==FALSE) {
        if (the_homologs != 0) {
            for(i=0; i<ances_genome->get_num_genes(); i++) delete the_homologs[i];
            delete[] the_homologs;
        }
        
        (*ances_genome)=(*assign_from.get_ancestral_genome());
        (*dupl_genome)=(*assign_from.get_dupl_genome());
    }
    else {
        for(i=0; i<ances_genome->get_num_genes(); i++) {
                if ((*assign_from.get_ancestral_genome())[i].gene_used() ==TRUE)
                    (*ances_genome)[i].set_used((*assign_from.get_ancestral_genome())[i].get_pillar());
                else
                    (*ances_genome)[i].set_unused();
        }
        for(i=0; i<dupl_genome->get_num_genes(); i++) {
            if ((*assign_from.get_dupl_genome())[i].gene_used() ==TRUE)
                (*dupl_genome)[i].set_used((*assign_from.get_dupl_genome())[i].get_pillar());
            else
                (*dupl_genome)[i].set_unused();
        }
        ances_genome->reset_neighbors();
        dupl_genome->reset_neighbors();
    }
    
    if (genomes_fixed==FALSE) {
        the_homologs=new  Gene_homologs*[ances_genome->get_num_genes()];
        for(i=0; i<ances_genome->get_num_genes(); i++) {
            new_vals=new double[assign_from[i].get_num_homologs()];
            new_homologs=new Scaffold_Gene* [assign_from[i].get_num_homologs()];
            for(j=0; j<assign_from[i].get_num_homologs(); j++) {
                new_homologs[j]=&(*dupl_genome)[assign_from[i].get_nth_homolog(j)->get_gene_num()];
                new_vals[j]=assign_from[i].get_nth_val(j);
            }
            
            the_homologs[i]=new Gene_homologs(assign_from[i].get_num_homologs(), &(*ances_genome)[i], new_homologs, new_vals);
            delete[] new_vals;
            delete[] new_homologs;
        }
    }
    
    if (assign_from.stored_order()== TRUE) {
        has_order=TRUE;
        if (assigned_order ==0)
            assigned_order=new int [ances_genome->get_num_genes()];
        
        for(i=0; i<ances_genome->get_num_genes(); i++) assigned_order[i]=assign_from.get_order_pos(i);
    }
    
    return(*this);
}

void Homolog_Set::get_new_order(string order_file)
{
    int i, cnt;
    string ances_name, dupl_name, info, info2;
    ifstream fin;
    map<string, int> name_map;
    
    fin.open(order_file.c_str());
    if (fin.fail()) {
        cerr<<"ERROR: Could not open order file "<<order_file<<endl;
    }
    else {
        assigned_order =new int[ances_genome->get_num_genes()];
        has_order=TRUE;
        
        for(i=0; i<ances_genome->get_num_genes(); i++)
            name_map[(*ances_genome)[i].get_name()]=(*ances_genome)[i].get_gene_num();
        
        std::getline(fin, info);
        std::getline(fin, info);
        
        cnt=0;
        
        while(!fin.eof()) {
            fin>>ances_name;
            fin>>info;
            fin>>dupl_name>>info>>info2;
            fin>>dupl_name>>info>>info2;
            if (! fin.eof()) {
                assigned_order[name_map[ances_name]]=cnt;
                cout<<"Gene "<<ances_name<<" that was originally in position "<<name_map[ances_name]<<" is now in position "<<cnt<<endl;
                cnt++;
            }
        }
        fin.close();
    }
}
int Homolog_Set::get_order_pos(int orig_pos)
{
    if ((has_order==FALSE) || (abs(orig_pos) >= ances_genome->get_num_genes())) {
        cerr<<"ERROR in request for Homolog_Set order information\n";
        return (-1);
    }
    else return(assigned_order[abs(orig_pos)]);
}




Homolog_Set::~Homolog_Set()
{
    int i;
    
    if (the_homologs !=0) {
        for(i=0; i<ances_genome->get_num_genes(); i++) delete the_homologs[i];
        delete[] the_homologs;
    }

    
    delete ances_genome;
    delete dupl_genome;
    
    delete null_homolog;
    
}


Homolog_Set* read_homolog_set (string homolog_file, Scaffold_Genome *ances_genome, Scaffold_Genome *dupl_genome, double cutoff, int count_cutoff)
{
    int i, j, k, num_homo, num_new, num_new_dupl, cnt_homo, ances_id;
    double val;
    string ances_name, dupl_name, line;
    BOOL *has_homo;
    Scaffold_Gene **new_gene_list;
    Scaffold_Genome *new_ances_genome, *new_dupl_genome;
    map<string, int> ances_map, dupl_map, new_dupl_map, new_ances_map;
    ifstream fin;
    Homolog_Set *new_homologs, *final_homologs;
    
    fin.open(homolog_file.c_str());
    
    if (!(fin.fail())) {
        
    
        new_homologs=new Homolog_Set(ances_genome, dupl_genome);
        
        for(i=0; i<ances_genome->get_num_genes(); i++) {
            ances_map[(*ances_genome)[i].get_name_string()]=i;
            for(j=0; j<(*ances_genome)[i].get_num_tandems(); j++)
                ances_map[(*ances_genome)[i].get_nth_tandem(j)]=i;
        }
        
        for(i=0; i<dupl_genome->get_num_genes(); i++) {
            dupl_map[(*dupl_genome)[i].get_name_string()]=i;
            for(j=0; j<(*dupl_genome)[i].get_num_tandems(); j++)
                dupl_map[(*dupl_genome)[i].get_nth_tandem(j)]=i;
        }
        
        
        while((!fin.eof())) {
            fin>>ances_name>>num_homo;
            if (! fin.eof()){
                for(j=0; j<num_homo; j++) {
                    fin>>dupl_name>>val;
                    if (val <=cutoff) {
                        //cout<<"Adding "<<dupl_name<<" to "<<(*new_homologs->get_ancestral_genome())[ances_map[ances_name]].get_name()<<"("<<(*new_homologs)[ances_map[ances_name]].get_num_homologs()<<")"<<endl;
                        if ((ances_map.find(ances_name) == ances_map.end())  ||
                            (dupl_map.find(dupl_name) == dupl_map.end()))
                            cerr<<"ERROR: Could not find either "<<ances_name<<" or "<<dupl_name<<" in corresponding genome\n";
                        else
                            (*new_homologs)[ances_map[ances_name]].add_homolog(&(*new_homologs->get_dupl_genome())[dupl_map[dupl_name]], val);
                    }
                }
            }
            
        }
        fin.close();
        
        num_new=0;
        has_homo=new BOOL[new_homologs->get_ancestral_genome()->get_num_genes()];
        for(i=0; i<new_homologs->get_ancestral_genome()->get_num_genes(); i++) {
            if (((*new_homologs)[i].get_num_homologs() == 0) || ((count_cutoff!=-1) && ((*new_homologs)[i].get_num_homologs() > count_cutoff))) has_homo[i]=FALSE;
            else {
                has_homo[i]=TRUE;
                new_ances_map[(*new_homologs->get_ancestral_genome())[i].get_name()]=num_new;
                num_new++;
            }
        }
        
        for(i=0; i<new_homologs->get_dupl_genome()->get_num_genes(); i++) {
            cnt_homo=0;
            for(j=0; j<new_homologs->get_ancestral_genome()->get_num_genes(); j++) {
                if (has_homo[j] == TRUE) {
                    for(k=0; k<(*new_homologs)[j].get_num_homologs(); k++) {
                        if ((*new_homologs)[j].get_nth_homolog(k) == &(*new_homologs->get_dupl_genome())[i]) cnt_homo++;
                    }
                }
            }
            (*new_homologs->get_dupl_genome())[i].set_num_homologs(cnt_homo);
        }
        
        
        
        
        num_new_dupl=0;
        for(i=0; i<new_homologs->get_dupl_genome()->get_num_genes(); i++) {
             if ((*new_homologs->get_dupl_genome())[i].get_num_homologs() >0) {
                 new_dupl_map[(*new_homologs->get_dupl_genome())[i].get_name()]=num_new_dupl;
                 num_new_dupl++;
             }
         }
        
        cnt_homo=0;
        new_gene_list=new Scaffold_Gene*[num_new];
        for(i=0; i<new_homologs->get_ancestral_genome()->get_num_genes(); i++) {
            if (has_homo[i] == TRUE) {
                new_gene_list[cnt_homo] = new Scaffold_Gene();
                (*new_gene_list[cnt_homo])=(*new_homologs->get_ancestral_genome())[i];
                new_gene_list[cnt_homo]->set_gene_num(cnt_homo);
                cnt_homo++;
            }
        }
        new_ances_genome=new Scaffold_Genome(new_homologs->get_ancestral_genome()->get_genome_name(), num_new, new_gene_list);
        delete[] new_gene_list;

        new_gene_list=new Scaffold_Gene*[num_new_dupl];
        cnt_homo=0;
        for(i=0; i<new_homologs->get_dupl_genome()->get_num_genes(); i++) {
            if (new_dupl_map.find((*new_homologs->get_dupl_genome())[i].get_name())!=new_dupl_map.end()) {
                new_gene_list[cnt_homo] = new Scaffold_Gene();
                (*new_gene_list[cnt_homo])=(*new_homologs->get_dupl_genome())[i];
                new_gene_list[cnt_homo]->set_gene_num(new_dupl_map[(*new_homologs->get_dupl_genome())[i].get_name()]);
                cnt_homo++;
            }
        }
        new_dupl_genome=new Scaffold_Genome(new_homologs->get_dupl_genome()->get_genome_name(), num_new_dupl, new_gene_list);
        delete[] new_gene_list;
        
        
        cnt_homo=0;
        final_homologs= new Homolog_Set(new_ances_genome, new_dupl_genome);
        for(i=0; i<new_homologs->get_ancestral_genome()->get_num_genes(); i++) {
            if (has_homo[i] == TRUE) {
                for(j=0; j<(*new_homologs)[i].get_num_homologs(); j++) {
                    (*final_homologs)[cnt_homo].add_homolog(&(*final_homologs->get_dupl_genome())[new_dupl_map[(*new_homologs)[i].get_nth_homolog(j)->get_name()]], (*new_homologs)[i].get_nth_val(j));
                }
                cnt_homo++;
            }
        }
        
        
        delete new_ances_genome;
        delete new_dupl_genome;
        delete new_homologs;
        new_homologs=final_homologs;
        
        for(i=0; i<new_homologs->get_dupl_genome()->get_num_genes(); i++) {
            cnt_homo=0;
            for(j=0; j<new_homologs->get_ancestral_genome()->get_num_genes(); j++) {
                for(k=0; k<(*new_homologs)[j].get_num_homologs(); k++) {
                    if ((*new_homologs)[j].get_nth_homolog(k) == &(*new_homologs->get_dupl_genome())[i]) cnt_homo++;
                }
            }
            (*new_homologs->get_dupl_genome())[i].set_num_homologs(cnt_homo);
        }
        
        return(new_homologs);
    }
    else {
        cerr<<"ERROR: Invalid homolog file "<<homolog_file<<endl;
        return(0);
    }
}


Pillar::Pillar()
{
    cerr<<"Invalid call to base constructor of Pillar\n";
    out_neighbors=0;
    WGX_genes=0;
    WGX_neighbors=0;
    my_homologs=0;
    pillar_id=0;
}


Pillar::Pillar(Exchange *cexchange, Scaffold_Gene *out, Gene_homologs *homos, int id)
{
    int i;
    
    pillar_id=id;
    
    depth=cexchange->get_WGX_depth();
    curr_exchange=cexchange;
    
    outgroup=out;
    out_neighbors=new Scaffold_Gene* [2];
    
    out_neighbors[0]=out_neighbors[1]=0;
    
    WGX_genes=new Scaffold_Gene*[depth];
    WGX_neighbors=new Scaffold_Gene**[depth];
    
    for(i=0; i<depth; i++) {
        WGX_genes[i]=0;
        WGX_neighbors[i]=new Scaffold_Gene*[2];
        WGX_neighbors[i][0]=0;
        WGX_neighbors[i][1]=0;
    }
    my_homologs=homos;
}

void Pillar::set_WGX_gene(Scaffold_Gene* the_gene, int index_num)
{
    WGX_genes[abs(index_num)%depth]=the_gene;
}


void Pillar::set_outgoup_neighbor(Scaffold_Gene* new_neighbor, int index_num)
{
    if (new_neighbor == 0 )
        out_neighbors[abs(index_num)%2]=new_neighbor;
    else {
        if ((outgroup->get_neighbor(0) == new_neighbor) || (outgroup->get_neighbor(1) == new_neighbor))
        {
            out_neighbors[abs(index_num)%2]=new_neighbor;
        }
        else {
            cerr<<"ERROR: gene "<<new_neighbor->get_name_string()<<" is not the neighbor of "<<outgroup->get_name_string()<<endl;
        }
    }
}


void Pillar::set_WGX_neighbor(Scaffold_Gene* new_neighbor, int WGX_index, int neighbor_index)
{
    if (new_neighbor == 0) {
        WGX_neighbors[abs(WGX_index)%depth][abs(neighbor_index)%2]=new_neighbor;
    }
    else {
        if (WGX_genes[abs(WGX_index)%depth] !=0) {
            if ((WGX_genes[abs(WGX_index)%depth]->get_neighbor(0)==new_neighbor) || (WGX_genes[abs(WGX_index)%depth]->get_neighbor(1)==new_neighbor) ) {
                WGX_neighbors[abs(WGX_index)%depth][abs(neighbor_index)%2]=new_neighbor;
            }
            else {
                cerr<<"ERROR: gene "<<new_neighbor->get_name_string()<<" is not the neighbor of "<<WGX_genes[abs(WGX_index)%depth]<<endl;
            }
        }
        else {
            cerr<<"ERROR: TRying to set neighbor of null gene\n";
        }
    }
}

Pillar::~Pillar()
{
    int i;
    
    if (out_neighbors!=0) delete[] out_neighbors;
    if (WGX_genes!=0) delete[] WGX_genes;
    if (WGX_neighbors!=0) {
        for(i=0; i<depth; i++) delete[] WGX_neighbors[i];
        delete[] WGX_neighbors;
    }
}


WGX_scaffold::WGX_scaffold()
{
    cerr<<"Invalid call to base constructor of WGX_scaffold\n";
    null_pillar=new Pillar(curr_exchange, (Scaffold_Gene*)0, 0, 0);
    num_pillars=0;
    pillar_order=0;
    the_pillars=0;
    have_init_order=FALSE;
    
    the_homologs=new Homolog_Set();
    pillar_gene_fixed=0;
    num_blocks=0;
    block_set=0;
    num_empty=0;
    empty_index=0;
    curr_block_data=0;
    order_backref=0;
}




WGX_scaffold::WGX_scaffold(Homolog_Set *new_homologs, Exchange *cexchange)
{
    int i, j, k, cnt_open, old_score, new_score;
    BOOL all_one, changed;
    Gene_homologs *curr_homo;
    Scaffold_Gene *test_gene, *my_gene;
    
    curr_exchange=cexchange;
    
    num_empty=0;
    the_homologs=new Homolog_Set();
    (*the_homologs)=(*new_homologs);
    curr_exchange->set_num_sites(the_homologs->get_ancestral_genome()->get_num_genes());
    num_pillars=curr_exchange->get_num_sites();
    
    curr_block_data =new Gene_block[curr_exchange->get_WGX_depth()];
    
    if (the_homologs->stored_order() == FALSE)
        have_init_order=FALSE;
    else {
        have_init_order=TRUE;
        cout<<"Using stored initial pillar order\n";
    }
    
    
    //Avoid free/alloc of block list by setting it to max possible size.
    block_set=new Gene_block[num_pillars];
    
    //Avoid free/alloc of index by setting it to max possible size.
    empty_index=new int[num_pillars];
    
    null_pillar=new Pillar(curr_exchange, 0, 0, 0);
    
    the_pillars = new Pillar * [num_pillars];
    pillar_order = new int [num_pillars];
    order_backref = new int [num_pillars];
    
    for(i=0; i<num_pillars; i++) {
        if (the_homologs->stored_order() == FALSE) {
            pillar_order[i]=i;
            order_backref[i]=i;
        }
        else {
            pillar_order[i]=the_homologs->get_order_pos(i);
            order_backref[the_homologs->get_order_pos(i)]=i;
        }
        
        the_pillars[i]=new Pillar(curr_exchange, &(*the_homologs->get_ancestral_genome())[i], &(*the_homologs)[i], i);
    }
    
    pillar_gene_fixed=new BOOL [the_homologs->get_ancestral_genome()->get_num_genes()];
    
    num_open_pillars=0;
    
    for(i=0; i<the_homologs->get_dupl_genome()->get_num_genes(); i++) (*the_homologs->get_dupl_genome())[i].set_unused();
    
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        pillar_gene_fixed[i]=FALSE;
        
        if ((*the_homologs)[i].get_num_homologs() == 1) {
            //cout<<"Gene "<<(*the_homologs->get_ancestral_genome())[i].get_name()<<" has "<<(*the_homologs)[i].get_num_homologs()<<endl;
            if ((*the_homologs)[i].get_nth_homolog(0)->get_num_homologs() == 1) {
                if ((*the_homologs)[i].get_nth_homolog(0)->gene_used() ==FALSE) {
                    pillar_gene_fixed[i]=TRUE;
                    the_pillars[i]->set_WGX_gene((*the_homologs)[i].get_nth_homolog(0), 0);
                    (*the_homologs)[i].get_nth_homolog(0)->set_used(i);
                }
            }
        }
        if (pillar_gene_fixed[i] == FALSE) num_open_pillars++;
    }
    
    //check_pillars();
    
    cout<<num_open_pillars<<" of the ancestral genome can vary in homolog assignment\n";
    
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        if (the_pillars[i]->get_WGX_gene(0) == 0) {
            curr_homo=&(*the_homologs)[i];
            test_gene=(*the_homologs)[i].get_closest_homolog();
            
            if (test_gene->gene_used() ==FALSE) {
                the_pillars[i]->set_WGX_gene(test_gene, 0);
                test_gene->set_used(i);
            }
            else {
                j=0;
                test_gene=0;
                while((j<(*the_homologs)[i].get_num_homologs()) && (test_gene==0)) {
                    if ((*the_homologs)[i].get_nth_homolog(j)->gene_used() == FALSE) test_gene=(*the_homologs)[i].get_nth_homolog(j);
                    j++;
                }
                if (test_gene !=0) {
                    the_pillars[i]->set_WGX_gene(test_gene, 0);
                    test_gene->set_used(i);
                }
               
            }
                
        }
    }
    
   // check_pillars();
    
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        if (the_pillars[i]->get_WGX_gene(0) == 0) {
            curr_homo=&(*the_homologs)[i];
            if (curr_homo->get_num_homologs() == the_pillars[i]->get_WGX_depth()) {
                for(k=1; k<the_pillars[i]->get_WGX_depth(); k++) {
                    all_one=TRUE;
                    for(j=0; j<curr_homo->get_num_homologs(); j++) {
                        if (curr_homo->get_nth_homolog(j)->get_num_homologs() != 1) all_one =FALSE;
                    }
                    if (all_one == TRUE) {
                        test_gene=curr_homo->get_nth_homolog(0);
                        j=0;
                        while((test_gene->gene_used() ==TRUE) && (j<curr_homo->get_nth_homolog(j)->get_num_homologs())) {
                            j++;
                            test_gene=curr_homo->get_nth_homolog(j);
                        }
                        if ((j<curr_homo->get_num_homologs()) && (test_gene->gene_used() == FALSE)) {
                            the_pillars[i]->set_WGX_gene(test_gene, k);
                            test_gene->set_used(i);
                        }
                    }
                }
            }
        }
    }
    
    
    open_pillar_index=new Pillar * [num_open_pillars];
    cnt_open=0;
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        if (pillar_gene_fixed[i] == FALSE) {
            open_pillar_index[cnt_open]=the_pillars[i];
            cnt_open++;
        }
    }
    
    update_neighbors();
   
    //check_pillars();
    cout<<"Initial total joins: "<<get_neighbor_count()<<endl;
    
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        if (the_pillars[i]->get_WGX_gene(0) !=0) {
            my_gene=the_pillars[i]->get_WGX_gene(0);
            update_neighbors();
            old_score=get_neighbor_count();
            changed=FALSE;
            
            
            for(j=1; j<the_pillars[i]->get_WGX_depth(); j++) {
                if (changed == FALSE) {
                    the_pillars[i]->set_WGX_gene(0,0);
                    the_pillars[i]->set_WGX_gene(my_gene,j);
                    update_neighbors();
                    new_score=get_neighbor_count();
                    
                    if (new_score < old_score) {
                        the_pillars[i]->set_WGX_gene(0,j);
                        the_pillars[i]->set_WGX_gene(my_gene,0);
                    }
                    else changed=TRUE;
                  //  if (new_score > old_score)
                    //    cout<<"Swap improves from "<<old_score<<" to "<<new_score<<" at "<<i<<endl;
                //}
                }
            }
        }
    }
    
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        for(j=0; j<the_pillars[i]->get_WGX_depth(); j++) {
            if ((the_pillars[i]->get_WGX_gene(j)!=0) && (the_pillars[i]->get_WGX_gene(j)->get_pillar() !=i))
                cerr<<"ERROR in initial pillar back refs at "<<i<<" for "<<the_pillars[i]->get_WGX_gene(j)->get_name()<<" = "<<the_pillars[i]->get_WGX_gene(j)->get_pillar()<<endl;
        }
    }
    cout<<"Pillar assigned checked pre neighbor\n";
    //check_pillars();
    update_neighbors();
    
    for(i=0; i<the_homologs->get_ancestral_genome()->get_num_genes(); i++) {
        for(j=0; j<the_pillars[i]->get_WGX_depth(); j++) {
            if ((the_pillars[i]->get_WGX_gene(j)!=0) && (the_pillars[i]->get_WGX_gene(j)->get_pillar() !=i))
                cerr<<"ERROR in initial pillar back refs at "<<i<<" for "<<the_pillars[i]->get_WGX_gene(j)->get_name()<<" = "<<the_pillars[i]->get_WGX_gene(j)->get_pillar()<<endl;
        }
    }
    cout<<"Pillar assigned checked post-neighbor\n";
   // check_pillars();
    optimize_assigns();
   // check_pillars();
    update_neighbors();
    set_blocks();
    cout<<"Initial order has "<<get_neighbor_count()<<" neighbors and "<<get_num_blocks()<<" blocks\n";
}


Pillar& WGX_scaffold::operator[] (int element)
{
    
    if ((0<=element) &&(element<num_pillars))
        return(*the_pillars[element]);
    else
        return(*null_pillar);
}



WGX_scaffold& WGX_scaffold::operator= (WGX_scaffold &assign_from)
{
    int i, j;
    
    (*the_homologs)=(*assign_from.get_homolog_set());
    if (assign_from.get_homolog_set()->stored_order()==TRUE) have_init_order=TRUE;
    
    for(i=0; i<num_pillars; i++) pillar_order[i]=assign_from.get_order(i);
    
    for(i=0; i<num_pillars; i++) {
        for(j=0; j<assign_from[i].get_WGX_depth(); j++) {
            if (assign_from[i].get_WGX_gene(j) != 0)
                the_pillars[i]->set_WGX_gene(&(*the_homologs->get_dupl_genome())[assign_from[i].get_WGX_gene(j)->get_gene_num()],j);
            else
                the_pillars[i]->set_WGX_gene(0,j);
        }
    }
    update_neighbors();
    set_blocks();
    return(*this);
}


void WGX_scaffold::change_pillar_order (int *new_order)
{
    int i;
    
    for(i=0; i<num_pillars; i++) pillar_order[i]=new_order[i];
}

void WGX_scaffold::update_neighbors()
{
    int pillar, WGX_level, pillar_index, cnt_joins=0, i;
    Pillar *my_pillar, *next_pillar, *last_pillar;
    BOOL empty;
    
    
    for(pillar=0; pillar<num_pillars; pillar++)
        order_backref[pillar_order[pillar]]=pillar;
    
    //if (empty_index!=0) delete[] empty_index;
    num_empty=0;
    
    for(pillar=0; pillar<num_pillars; pillar++) {
        empty=TRUE;
        for(WGX_level=0; WGX_level<the_pillars[pillar]->get_WGX_depth(); WGX_level++) {
            if (the_pillars[pillar]->get_WGX_gene(WGX_level) !=0) empty=FALSE;
        }
        
        if (empty == TRUE) {
            //i=0;
            //while(pillar_order[i] != pillar) i++;
            empty_index[num_empty]=pillar;
            num_empty++;
        }
    }
    
   // if (num_empty >0)
   //     empty_index=new int[num_empty];
    //else
    //    empty_index=0;
    
    
    
    for(pillar=0; pillar<num_pillars; pillar++) {
        my_pillar=the_pillars[pillar_order[pillar]];
        my_pillar->set_outgoup_neighbor(0,0);
        my_pillar->set_outgoup_neighbor(0,1);
        
        for(WGX_level=0; WGX_level<my_pillar->get_WGX_depth(); WGX_level++) {
            if (my_pillar->get_WGX_gene(WGX_level) != 0) {
                my_pillar->set_WGX_neighbor(0, WGX_level,0);
                my_pillar->set_WGX_neighbor(0, WGX_level,1);
            }
        }
    }

    for(pillar=0; pillar<num_pillars; pillar++) {
        my_pillar=the_pillars[pillar_order[pillar]];
        
        if (pillar < (num_pillars-1)) {
            next_pillar=the_pillars[pillar_order[pillar+1]];
            if (are_neighbors(my_pillar->get_outgroup_gene(), next_pillar->get_outgroup_gene()) ==TRUE) {
                my_pillar->set_outgoup_neighbor(next_pillar->get_outgroup_gene(), 1);
                cnt_joins++;
            }
            
            for(WGX_level=0; WGX_level<my_pillar->get_WGX_depth(); WGX_level++) {
                if (my_pillar->get_WGX_gene(WGX_level) != 0) {
                    pillar_index=pillar+1;
                    next_pillar=the_pillars[pillar_order[pillar_index]];
                    while((next_pillar->get_WGX_gene(WGX_level)==0) && (pillar_index <(num_pillars-1))) {
                        pillar_index++;
                        next_pillar=the_pillars[pillar_order[pillar_index]];
                    }
                    
                    if (next_pillar->get_WGX_gene(WGX_level)!= 0) {
                        if (are_neighbors(my_pillar->get_WGX_gene(WGX_level), next_pillar->get_WGX_gene(WGX_level)) == TRUE) {
                            my_pillar->set_WGX_neighbor(next_pillar->get_WGX_gene(WGX_level), WGX_level, 1);
                            cnt_joins++;
                        }
                    }
                }
            }
            
            
        }
        if (pillar>0) {
            last_pillar=the_pillars[pillar_order[pillar-1]];
            if (are_neighbors(my_pillar->get_outgroup_gene(), last_pillar->get_outgroup_gene()) ==TRUE)
                my_pillar->set_outgoup_neighbor(last_pillar->get_outgroup_gene(), 0);
            
            for(WGX_level=0; WGX_level<my_pillar->get_WGX_depth(); WGX_level++) {
                if (my_pillar->get_WGX_gene(WGX_level) != 0) {
                    pillar_index=pillar-1;
                    last_pillar=the_pillars[pillar_order[pillar_index]];
                    
                    while((last_pillar->get_WGX_gene(WGX_level)==0) && (pillar_index >0)) {
                        pillar_index--;
                        last_pillar=the_pillars[pillar_order[pillar_index]];
                    }
                    
                    if (last_pillar->get_WGX_gene(WGX_level) !=0) {
                        if (are_neighbors(my_pillar->get_WGX_gene(WGX_level), last_pillar->get_WGX_gene(WGX_level)) == TRUE) {
                            my_pillar->set_WGX_neighbor(last_pillar->get_WGX_gene(WGX_level), WGX_level, 0);
                            cnt_joins++;
                        }
                    }

                }
            }
        }
    }
    //cout<<"Made "<<cnt_joins<<" while setting neighbors\n";
}
    
    

double WGX_scaffold::get_neighbor_count()
{
    int pillar, WGX_level;
    double retval=0;
    
    for(pillar=0; pillar<num_pillars; pillar++) {
        if (the_pillars[pillar_order[pillar]]->get_outgroup_neighbor(0) !=0) retval+=0.25;
        if (the_pillars[pillar_order[pillar]]->get_outgroup_neighbor(1) !=0) retval+=0.25;
        for(WGX_level=0; WGX_level<the_pillars[pillar_order[pillar]]->get_WGX_depth(); WGX_level++) {
            if (the_pillars[pillar_order[pillar]]->get_WGX_gene(WGX_level) != 0) {
                if (the_pillars[pillar_order[pillar]]->get_WGX_neighbor(WGX_level,0) !=0) retval+=1.0;
                if (the_pillars[pillar_order[pillar]]->get_WGX_neighbor(WGX_level,1) !=0) retval+=1.0;
            }
        }
        
    }
    return(retval);
}


void WGX_scaffold::optimize_assigns()
{
    int pillar, WGX_level, homolog, l_pillar, n_pillar, old_level;
    Scaffold_Gene *curr_gene, *orig_gene;
    Pillar *curr_pillar, *last_pillar, *next_pillar, *old_pillar;
    BOOL drop, changed;
    
    
    for(pillar=0; pillar<num_pillars; pillar++) {
        curr_pillar=the_pillars[pillar_order[pillar]];
        
        
        for(WGX_level=0; WGX_level<curr_pillar->get_WGX_depth(); WGX_level++) {
            drop=FALSE;
            
            l_pillar=pillar-1;
            
            while((l_pillar >=0) &&(the_pillars[pillar_order[l_pillar]]->get_WGX_gene(WGX_level)==0)) l_pillar--;
            
            if ((l_pillar >=0) && (the_pillars[pillar_order[l_pillar]]->get_WGX_gene(WGX_level)!=0))
                last_pillar=the_pillars[pillar_order[l_pillar]];
            else last_pillar=0;
            
            n_pillar=pillar+1;
            
            while((n_pillar < num_pillars) && (the_pillars[pillar_order[n_pillar]]->get_WGX_gene(WGX_level)==0)) n_pillar++;
            
            if ((n_pillar< num_pillars) && (the_pillars[pillar_order[n_pillar]]->get_WGX_gene(WGX_level)!=0))
                next_pillar =the_pillars[pillar_order[n_pillar]];
            else next_pillar=0;
            
        
            if ((curr_pillar->get_WGX_gene(WGX_level)!=0) && (curr_pillar->get_WGX_neighbor(WGX_level, 0) ==0)
                && (curr_pillar->get_WGX_neighbor(WGX_level,1)==0)) {
                //We have a gene, but no links
                orig_gene=curr_pillar->get_WGX_gene(WGX_level);
                
                if ((last_pillar != 0) && (next_pillar != 0)) {
                    if (are_neighbors(last_pillar->get_WGX_gene(WGX_level), next_pillar->get_WGX_gene(WGX_level)) == TRUE) {
                        orig_gene->set_unused();
                        curr_pillar->set_WGX_gene(0, WGX_level);
                      // cout<<"Clearing gene "<<orig_gene->get_name()<<" from pillar "<<curr_pillar->get_pillar_id()<<" level "<<WGX_level<<" to join last to next"<<endl;
                    }
                    else {
                        changed=FALSE;
                        homolog=0;
                        while((changed == FALSE) && (homolog<curr_pillar->get_my_homologs()->get_num_homologs())) {
                            //for(homolog=0; homolog < curr_pillar->get_my_homologs()->get_num_homologs(); homolog++) {
                            curr_gene=curr_pillar->get_my_homologs()->get_nth_homolog(homolog);
                            if (curr_gene != orig_gene) {
                                if ((are_neighbors(curr_gene, last_pillar->get_WGX_gene(WGX_level))==TRUE) ||
                                    (are_neighbors(curr_gene, next_pillar->get_WGX_gene(WGX_level))==TRUE)) {
                                    changed=TRUE;
                                    if (curr_gene->gene_used()== TRUE) {
                                        
                                        old_pillar=the_pillars[curr_gene->get_pillar()];
                                        old_level=0;
                                        while((old_level < old_pillar->get_WGX_depth()) && (old_pillar->get_WGX_gene(old_level) != curr_gene)) old_level++;
                                        
                                      //  cout<<"Clearing gene "<<curr_gene->get_name()<<" from old pillar "<<curr_gene->get_pillar()<<" level "<<old_level<<" to use at "<<pillar<<" level: "<<WGX_level<<endl;
                                        
                                        if (old_level >= old_pillar->get_WGX_depth()) cerr<<"ERROR: gene "<<curr_gene->get_name()<<" points to invalid pillar loc: "<<curr_gene->get_pillar()<<endl;
                                        
                                        old_pillar->set_WGX_gene(0, old_level);
                                    }
                                    
                                    drop=TRUE;
                                    curr_gene->set_used(curr_pillar->get_pillar_id());
                                    //cout<<pillar<<" Swap: Assigned gene "<<curr_gene->get_name()<<" to pillar "<<curr_pillar->get_pillar_id()<<" at "<<WGX_level<<endl;
                                    curr_pillar->set_WGX_gene(curr_gene, WGX_level);
                                }
                            }
                            homolog++;
                        }
                        if (drop == TRUE) orig_gene->set_unused();
                    }
                }
                
            }
            else {
                if (curr_pillar->get_WGX_gene(WGX_level)==0) {
                    if ((last_pillar != 0) && (next_pillar != 0)) {
                        if (are_neighbors(last_pillar->get_WGX_gene(WGX_level), next_pillar->get_WGX_gene(WGX_level)) == FALSE) {
                            changed=FALSE;
                            homolog=0;
                            while((changed == FALSE) && (homolog<curr_pillar->get_my_homologs()->get_num_homologs())) {
                            //for(homolog=0; homolog < curr_pillar->get_my_homologs()->get_num_homologs(); homolog++) {
                                curr_gene=curr_pillar->get_my_homologs()->get_nth_homolog(homolog);
                                if ((are_neighbors(curr_gene, last_pillar->get_WGX_gene(WGX_level))==TRUE) ||
                                    (are_neighbors(curr_gene, next_pillar->get_WGX_gene(WGX_level))==TRUE)) {
                                    changed=TRUE;
                                        if (curr_gene->gene_used()== TRUE) {
                                            old_pillar=the_pillars[curr_gene->get_pillar()];
                                            old_level=0;
                                            while((old_level < old_pillar->get_WGX_depth()) && (old_pillar->get_WGX_gene(old_level) != curr_gene)) old_level++;
                                            
                                            //cout<<"Clearing gene "<<curr_gene->get_name()<<" from old pillar "<<curr_gene->get_pillar()<<" level "<<old_level<<" to use at "<<pillar<<" level: "<<WGX_level<<endl;
                                            
                                            if (old_level >= old_pillar->get_WGX_depth()) cerr<<"ERROR: gene "<<curr_gene->get_name()<<" points to invalid pillar loc: "<<curr_gene->get_pillar()<<endl;
                                            
                                            old_pillar->set_WGX_gene(0, old_level);
                                        }
                                    
                                        curr_gene->set_used(curr_pillar->get_pillar_id());
                                       //cout<<pillar<<" level "<<WGX_level<< " is empty: Assigned gene "<<curr_gene->get_name()<<" to pillar "<<curr_pillar->get_pillar_id()<<" at "<<WGX_level<<endl;
                                        curr_pillar->set_WGX_gene(curr_gene, WGX_level);
                                }
                                homolog++;

                            }
                        }
                    }
                }
            }
        }
    }
}


int WGX_scaffold::get_block_start(int bl)
{
    return(block_set[bl].start);
}

int WGX_scaffold::get_block_end(int bl)
{
    return(block_set[bl].end);
}


void WGX_scaffold::set_blocks()
{
    int curr_pillar, next_pillar, WGX_level, curr_block_end, min_end;
    BOOL has_neighbors;
    Scaffold_Gene *neighbor;
    
    curr_pillar=0;
    
    num_blocks=0;
    
    while(curr_pillar < num_pillars) {
        for(WGX_level=0; WGX_level<the_pillars[pillar_order[curr_pillar]]->get_WGX_depth(); WGX_level++ ) {
            next_pillar=curr_pillar;
            while((next_pillar< num_pillars) && (the_pillars[pillar_order[next_pillar]]->get_WGX_gene(WGX_level) == 0)) next_pillar++;
            
            curr_block_data[WGX_level].start=next_pillar;
            
            if (next_pillar < num_pillars) {
                neighbor=the_pillars[pillar_order[next_pillar]]->get_WGX_neighbor(WGX_level, 1);
                

                if (neighbor !=0) {
                    if (neighbor->gene_used()==FALSE) cerr<<"ERROR: Neighboring gene not in use but assigned\n";
                    
                    while(the_pillars[neighbor->get_pillar()]->get_WGX_neighbor(WGX_level, 1) != 0) {
                        neighbor=the_pillars[neighbor->get_pillar()]->get_WGX_neighbor(WGX_level, 1);
                    }
                
                    curr_block_data[WGX_level].end=order_backref[neighbor->get_pillar()];
                }
                else {
                    curr_block_data[WGX_level].end=next_pillar;
                }
            }
            else curr_block_data[WGX_level].end=num_pillars-1;
        }
        
        min_end=curr_block_data[0].end;
        for(WGX_level=1; WGX_level<the_pillars[pillar_order[curr_pillar]]->get_WGX_depth(); WGX_level++ ) {
            if (curr_block_data[WGX_level].end < min_end) min_end=curr_block_data[WGX_level].end;
        }
        
        curr_block_end=min_end;
        
#if 0
        next_pillar=curr_pillar;
        while(next_pillar < curr_block_end) {
            empty=TRUE;
            for(WGX_level=0; WGX_level<the_pillars[pillar_order[next_pillar]]->get_WGX_depth(); WGX_level++ ) {
                if (the_pillars[pillar_order[next_pillar]]->get_WGX_gene(WGX_level) != 0) empty=FALSE;
            }
            
            if (empty == TRUE) curr_block_end=next_pillar;
            else next_pillar++;
        }
#endif
        block_set[num_blocks].start=curr_pillar;
        block_set[num_blocks].end=curr_block_end;
        if (curr_block_end > num_pillars) cerr<<"Error: block "<<num_blocks<<" ends after pillars: "<<curr_block_end<<endl;
        curr_pillar=curr_block_end+1;
        num_blocks++;
    }
    
    
    
       //cout<<"Setting "<<num_blocks<<" blocks. End is "<<block_set[num_blocks-1].end<<endl;
    
 
}


WGX_scaffold::~WGX_scaffold()
{
    int i;
    
    if (the_pillars != 0) {
        for(i=0; i<num_pillars; i++) delete[] the_pillars[i];
        delete[] the_pillars;
    }
    if (pillar_order!=0) delete[] pillar_order;
    if (order_backref!=0) delete[] order_backref;
    delete null_pillar;
    delete the_homologs;
    if (block_set !=0) delete[] block_set;
    //if (empty_index!=0) delete[] empty_index;
}


void WGX_scaffold::check_pillars()
{
    int pillar, level, level2;
    
    for(pillar=0; pillar<num_pillars; pillar++) {
        for(level=0; level<the_pillars[pillar_order[pillar]]->get_WGX_depth(); level++) {
            if (the_pillars[pillar_order[pillar]]->get_WGX_gene(level)!=0) {
                for(level2=level+1; level2<the_pillars[pillar_order[pillar]]->get_WGX_depth(); level2++) {
                    if (the_pillars[pillar_order[pillar]]->get_WGX_gene(level2)== the_pillars[pillar_order[pillar]]->get_WGX_gene(level)) {
                        cerr<<"Error at pillar "<<pillar<<" and levels "<<level<<", "<<level2<<" : Gene "<<the_pillars[pillar_order[pillar]]->get_WGX_gene(level)->get_name()<<" is present at both"<<endl;
                    }
                }
                if (the_pillars[pillar_order[pillar]]->get_WGX_gene(level)->gene_used()==FALSE)
                    cerr<<"Error: gene "<<the_pillars[pillar_order[pillar]]->get_WGX_gene(level)->get_name()<<" at pillar "<<pillar<<" and level "<<level<<" is set as unused\n";
            }
        }
    }
}

BOOL Scaffold_point::check_order()
{
    int i;
    BOOL retval=TRUE;
    
    for(i=0; i<my_scaffold->get_num_pillars(); i++) new_order[i]=0;
    for(i=0; i<my_scaffold->get_num_pillars(); i++) {
        new_order[pillar_order[i]]++;
        
        if (new_order[pillar_order[i]] >1) {
            retval=FALSE;
            cerr<<"ERROR: Found two entries for pillar "<<i<<" Current is "<<pillar_order[i]<<endl;
        }
    }
    
    for(i=0; i<my_scaffold->get_num_pillars(); i++) {
        if (new_order[i]==0)  {
            cerr<<"ERROR: Never found an entry for pillar "<<i<<endl;
            retval=FALSE;
        }
    }
    return(retval);
    
}


Scaffold_point::Scaffold_point(Exchange *cexchange)
{

    int i, j;
    Scaffold_Gene *new_homolog;
  
    curr_exchange=new Exchange();
    
    
    //Copies the exchange object into the new exchange
    //(Prevents conflicts between the exchanges of different walkers)
    (*curr_exchange)=(*cexchange);
  
    
    my_scaffold=new WGX_scaffold(global_homologs, curr_exchange);
    
    pillar_order=new int [curr_exchange->get_num_sites()];
    new_order=new int [curr_exchange->get_num_sites()];
    
    for(i=0; i<curr_exchange->get_num_sites(); i++)
	  pillar_order[i]=i;
    
    for(i=0; i<my_scaffold->get_homolog_set()->get_ancestral_genome()->get_num_genes(); i++) {
        j=0;
        new_homolog=0;
        while((j<(*my_scaffold->get_homolog_set())[i].get_num_homologs()) && (new_homolog==0)) {
            if ((*my_scaffold->get_homolog_set())[i].get_nth_homolog(j)->gene_used() == FALSE) {
                new_homolog=(*my_scaffold->get_homolog_set())[i].get_nth_homolog(j);
                
            }
            j++;
        }
        if (new_homolog !=0) {
            j=0;
            while((j<(*my_scaffold)[i].get_WGX_depth()) && ((*my_scaffold)[i].get_WGX_gene(j) !=0)) j++;
            
            if (j<(*my_scaffold)[i].get_WGX_depth()) {
                (*my_scaffold)[i].set_WGX_gene(new_homolog, j);
                new_homolog->set_used(i);
            }
        }
    }
}


Scaffold_point& Scaffold_point::operator=(Scaffold_point & assign_from)
{
    int i;
    
    (*my_scaffold)=(*assign_from.my_scaffold);
    
    for(i=0; i<my_scaffold->get_homolog_set()->get_num_ances_genes(); i++)
       pillar_order[i]=assign_from.pillar_order[i];
    
    return(*this);
}



Scaffold_point::~Scaffold_point ()
{
    delete[] pillar_order;
    delete[] new_order;
    delete my_scaffold;
    delete curr_exchange;
}







void Scaffold_anneal::move(int walker)
{
    int i, j,  dscore, old,  start, length, new_loc, second, empty_loc, order_loc, depth,
        rand1, rand2, num_full, rand_dupl, rand_slot, old_slot, old_pillar, block_start, block_length, new_block_loc;
    double val, rand_val;
    BOOL prev, do_reverse=FALSE, empty, goto_rand_loc;
    Scaffold_point *the_point;
    Pillar *change_pillar;
    Scaffold_Gene *gene1, *gene2, *curr_gene, *new_gene, *neighbor;


    pre_move->point->my_scaffold->get_homolog_set()->fix_genomes();
    *pre_move=*(*current_walkers)[walker];
    (*current_walkers)[walker]->point->my_scaffold->get_homolog_set()->fix_genomes();

    the_point=(*current_walkers)[walker]->point;

    rand_val=ranf();
    if (rand_val<0.3) {
        //start =ignuin(0, the_point->my_scaffold->get_homolog_set()->get_num_ances_genes()-1);
        block_start=ignuin(0, the_point->my_scaffold->get_num_blocks()-1);
        
        if (the_point->my_scaffold->order_preset() == FALSE) {
            if (ranf() <0.5)
                block_length = (int)gennor(3, 2);
            else
                block_length=1;
        }
        else {
            if (ranf() <0.6)
                block_length = (int)gennor(9, 7);
            else
                block_length=1;
        }
        //while ((length < 1) || (start+length > the_point->my_scaffold->get_homolog_set()->get_num_ances_genes())) {
        //    start =ignuin(0, the_point->my_scaffold->get_homolog_set()->get_num_ances_genes()-1);
         //   length = (int)gennor(12, 4);
        //}
        
        while ((block_length < 1) || ((block_start+block_length) > the_point->my_scaffold->get_num_blocks())) {
            block_start=ignuin(0, the_point->my_scaffold->get_num_blocks()-1);
            if (the_point->my_scaffold->order_preset() == FALSE)
                block_length = (int)gennor(3, 2);
            else
                block_length = (int)gennor(9, 7);
        }
        
        start=the_point->my_scaffold->get_block_start(block_start);
        length = the_point->my_scaffold->get_block_end(block_start+(block_length-1))-the_point->my_scaffold->get_block_start(block_start)+1;
        
        
        if (ranf()<0.25) {
          // cout<<"Moving far\n";
            //new_loc=ignuin(0, the_point->my_scaffold->get_homolog_set()->get_num_ances_genes()-(1+length));
            new_block_loc=ignuin(0,the_point->my_scaffold->get_num_blocks()-(1+block_length));
            
            if (new_block_loc >= block_start) new_block_loc+=block_length;
            
            new_loc=the_point->my_scaffold->get_block_start(new_block_loc);
            
           // cout<<" Start "<<start<<" len: "<<length<<" Loc: "<<new_loc<<"( BL: "<<new_block_loc<<" t: "<<the_point->my_scaffold->get_num_blocks()<<" -> "<<the_point->my_scaffold->get_block_start(new_block_loc)<<")"<<endl;
            
        }
        else {
           // cout<<"Moving near";
            //new_loc=gennor(0,18)+start;
            new_block_loc=gennor(0,18)+block_start;
        
            
            while(((new_block_loc >=block_start) && (new_block_loc < (block_start+block_length))) || (new_block_loc <0) || (new_block_loc >= the_point->my_scaffold->get_num_blocks()))
                new_block_loc=gennor(0,18)+block_start;
            
            new_loc=the_point->my_scaffold->get_block_start(new_block_loc);
            
            //if (new_loc >= start) new_loc+=length;
            
         // cout<<" Start "<<start<<" len: "<<length<<" Loc: "<<new_loc<<"( BL: "<<new_block_loc<<" t: "<<the_point->my_scaffold->get_num_blocks()<<" -> "<<the_point->my_scaffold->get_block_start(new_block_loc)<<")"<<endl;
            
        }
    
        if (ranf() > 0.3) do_reverse=TRUE;
        
        //if (do_reverse==FALSE) cout<<"Moving forward "<<start<<" for "<<length<<" to "<<new_loc<<"\t";
        //else cout<<"Moving reverse "<<start<<" for "<<length<<" to "<<new_loc<<"\t";
        
        //new_order=new int[the_point->my_scaffold->get_homolog_set()->get_num_ances_genes()];
        
        if (start > new_loc) {
            for(i=0; i<new_loc; i++) the_point->new_order[i]=the_point->pillar_order[i];
            if (do_reverse==FALSE) {
                for(i=new_loc; i<new_loc+length; i++) the_point->new_order[i] = the_point->pillar_order[start+(i-new_loc)];
            }
            else {
                for(i=new_loc; i<new_loc+length; i++) the_point->new_order[i] = the_point->pillar_order[(start+length-1)-(i-new_loc)];
            }
            for(i=new_loc+length; i<start+length; i++) the_point->new_order[i] = the_point->pillar_order[i-length];
            for(i=start+length; i<the_point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++)
                the_point->new_order[i]=the_point->pillar_order[i];
        }
        else {
            for(i=0; i<start; i++) the_point->new_order[i]=the_point->pillar_order[i];
            for(i=start+length; i<new_loc; i++) the_point->new_order[i-length]=the_point->pillar_order[i];
            if (do_reverse==FALSE) {
                for(i=new_loc; i<new_loc+length; i++) the_point->new_order[i-length]=the_point->pillar_order[start+(i-new_loc)];
            }
            else {
                for(i=new_loc; i<new_loc+length; i++) the_point->new_order[i-length]=the_point->pillar_order[(start+length-1)-(i-new_loc)];
            }
            for(i=new_loc; i<the_point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++)
                    the_point->new_order[i]=the_point->pillar_order[i];
            
        }
        for(i=0; i<the_point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++) the_point->pillar_order[i]=the_point->new_order[i];
        the_point->check_order();
        
       // delete[] new_order;
    }
    else {
        if (rand_val <0.8) {
           // cout<<"Swapping\t";
            //Swap two tracks across a block
            new_block_loc=ignuin(0, the_point->my_scaffold->get_num_blocks()-1);
            //new_loc=ignuin(0,the_point->my_scaffold->get_homolog_set()->get_num_ances_genes()-1);
            new_loc=the_point->my_scaffold->get_block_start(new_block_loc);
            
          //  cout<<" using block "<<new_block_loc<<" actual order coords: "<<new_loc<<endl;
            
            change_pillar=the_point->my_scaffold->get_nth_ordered_pillar(new_loc);
            rand1=ignuin(0,change_pillar->get_WGX_depth()-1);
            rand2=ignuin(0,change_pillar->get_WGX_depth()-1);
            while(rand1 == rand2) rand2=ignuin(0,change_pillar->get_WGX_depth()-1);

            if (rand1 > rand2) {
                old=rand1;
                rand1=rand2;
                rand2=old;
            }
            
            for(new_loc=the_point->my_scaffold->get_block_start(new_block_loc); new_loc<=the_point->my_scaffold->get_block_end(new_block_loc); new_loc++) {
                change_pillar=the_point->my_scaffold->get_nth_ordered_pillar(new_loc);
                gene1=change_pillar->get_WGX_gene(rand1);
                gene2=change_pillar->get_WGX_gene(rand2);
                change_pillar->set_WGX_gene(gene2, rand1);
                change_pillar->set_WGX_gene(gene1, rand2);
            }
        }
        else {
            if (ranf() < 0.3) {
              // cout<<"Filling empty\n";
            //Fill an empty pillar
                empty_loc=ignuin(0, the_point->my_scaffold->get_num_empty_pillars()-1);
                order_loc=the_point->my_scaffold->get_nth_empty_pillar_index(empty_loc);
                change_pillar=&(*the_point->my_scaffold)[order_loc];
                
                new_loc=0;
                while(the_point->pillar_order[new_loc] != change_pillar->get_pillar_id()) new_loc++;
                
                //cout<<"Empty pillar is "<<order_loc<<", which is "<<new_loc<<" in the order list\n";
                
                
                rand_slot=ignuin(0,change_pillar->get_WGX_depth()-1);
                rand_dupl=ignuin(0, change_pillar->get_my_homologs()->get_num_homologs()-1);
                
                
                old_pillar=-1;
                new_gene=change_pillar->get_my_homologs()->get_nth_homolog(rand_dupl);
                
                goto_rand_loc  = TRUE;
                
                if (new_gene->gene_used()==FALSE) {
                   // cout<<"Using unsued gene from "<<new_loc<<endl;
                    neighbor=0;
                    if (new_gene->get_neighbor(0) !=0) {
                        depth=0;
                        while((depth< change_pillar->get_WGX_depth()) && (neighbor != new_gene->get_neighbor(0))) {
                            neighbor =(*the_point->my_scaffold)[new_gene->get_neighbor(0)->get_pillar()].get_WGX_gene(depth);
                            if (neighbor != new_gene->get_neighbor(0)) depth++;
                        }
                        if (depth <change_pillar->get_WGX_depth() ){
                            if ((neighbor != new_gene->get_neighbor(0)) ||
                                (((*the_point->my_scaffold)[neighbor->get_pillar()].get_WGX_neighbor(depth, 0)!=0) &&
                                 ((*the_point->my_scaffold)[neighbor->get_pillar()].get_WGX_neighbor(depth, 1)!=0)))
                                neighbor=0;
                        }
                       // cout<<"For gene "<<new_gene->get_name()<<" neighbor 0 val is "<<neighbor<<" Pillar is "<<new_gene->get_neighbor(0)->get_pillar()<<" and depth is "<<depth<<endl;
                    }
                    if ((neighbor == 0) && (new_gene->get_neighbor(1) !=0)) {
                        depth=0;
                        while((depth< change_pillar->get_WGX_depth()) && (neighbor != new_gene->get_neighbor(1))) {
                            neighbor =(*the_point->my_scaffold)[new_gene->get_neighbor(1)->get_pillar()].get_WGX_gene(depth);
                            if (neighbor != new_gene->get_neighbor(1)) depth++;
                        }
                        if (depth <change_pillar->get_WGX_depth() ){
                            if ((neighbor != new_gene->get_neighbor(1)) ||
                                (((*the_point->my_scaffold)[neighbor->get_pillar()].get_WGX_neighbor(depth, 0)!=0) &&
                                 ((*the_point->my_scaffold)[neighbor->get_pillar()].get_WGX_neighbor(depth, 1)!=0)))
                                neighbor=0;
                        }
                     //   cout<<"For gene "<<new_gene->get_name()<<" neighbor 1 val is "<<neighbor<<" Pillar is "<<new_gene->get_neighbor(1)->get_pillar()<<" and depth is "<<depth<<endl;
                    }
                    if (neighbor !=0) goto_rand_loc=FALSE;
                }
             
                if (goto_rand_loc==FALSE) {
                    start= 0;
                    while(the_point->pillar_order[start] != neighbor->get_pillar()) start++;
                    
                   // cout<<"Trying targeted join: "<<new_loc<<" to "<<start<<endl;
                }
                else {
                    
                    block_start=ignuin(0, the_point->my_scaffold->get_num_blocks()-1);
                    
                    start=the_point->my_scaffold->get_block_start(block_start);
                    
                    while(start == new_loc) {
                        block_start=ignuin(0, the_point->my_scaffold->get_num_blocks()-1);
                        start=the_point->my_scaffold->get_block_start(block_start);
                    }
                }
                
                old=new_loc;
                new_loc=start;
                start=old;
                   // cout<<"Moving to "<<new_loc<<endl;
                    
                    //new_order=new int[the_point->my_scaffold->get_homolog_set()->get_num_ances_genes()];
                    
                if (start > new_loc) {
                    for(i=0; i<new_loc; i++) the_point->new_order[i]=the_point->pillar_order[i];
                    the_point->new_order[new_loc]=the_point->pillar_order[start];
                    
                    for(i=new_loc+1; i<start+1; i++) the_point->new_order[i] = the_point->pillar_order[i-1];
                    for(i=start+1; i<the_point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++)
                        the_point->new_order[i]=the_point->pillar_order[i];
                }
                else {
                    for(i=0; i<start; i++) the_point->new_order[i]=the_point->pillar_order[i];
                    for(i=start+1; i<new_loc; i++) the_point->new_order[i-1]=the_point->pillar_order[i];
                    the_point->new_order[new_loc-1]=the_point->pillar_order[start];
                    
                    for(i=new_loc; i<the_point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++)
                        the_point->new_order[i]=the_point->pillar_order[i];
                    
                }
                for(i=0; i<the_point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++) the_point->pillar_order[i]=the_point->new_order[i];
                
                the_point->check_order();
               // delete[] new_order;
                
              
                
               // cout<<"Assigninng homolog #"<<rand_dupl<<" which is "<<new_gene->get_name_string()<<" to slot "<<rand_slot<<endl;
                
                
                if (new_gene->gene_used() == TRUE) {
                    old_slot=0;
                    old_pillar=new_gene->get_pillar();
                    while((old_slot < (*the_point->my_scaffold)[new_gene->get_pillar()].get_WGX_depth()) &&
                          ((*the_point->my_scaffold)[new_gene->get_pillar()].get_WGX_gene(old_slot) != new_gene)) old_slot++;
                    
                    if (old_slot == (*the_point->my_scaffold)[new_gene->get_pillar()].get_WGX_depth())
                        cerr<<"ERROR in pillar mapping: gene "<<new_gene->get_name()<<" has wrong pillar address\n";
                    else {
                        (*the_point->my_scaffold)[new_gene->get_pillar()].set_WGX_gene(0, old_slot);
                    }
                }
                
                change_pillar->set_WGX_gene(new_gene,rand_slot);
                new_gene->set_used(change_pillar->get_pillar_id());
                
                
            }
            else {
               // cout<<"Exchanging\n";
                new_loc=ignuin(0,the_point->my_scaffold->get_num_open_pillars()-1);
                
                while (the_point->my_scaffold->get_nth_open_pillar(new_loc)->get_my_homologs()->get_num_homologs() ==1)
                     new_loc=ignuin(0,the_point->my_scaffold->get_num_open_pillars()-1);
                    
                change_pillar=the_point->my_scaffold->get_nth_open_pillar(new_loc);
                
                num_full=0;
                for(i=0; i<change_pillar->get_WGX_depth(); i++) {
                    if (change_pillar->get_WGX_gene(i) != 0) num_full++;
                }
                if ((num_full>1) && (rand_val>0.98)) {
                    
                    //Remove a duplicate
                    rand_slot=ignuin(0,change_pillar->get_WGX_depth()-1);
                    while(change_pillar->get_WGX_gene(rand_slot)==0)
                        rand_slot=ignuin(0,change_pillar->get_WGX_depth()-1);
                
                    curr_gene=change_pillar->get_WGX_gene(rand_slot);
                    change_pillar->set_WGX_gene(0,rand_slot);
                    curr_gene->set_unused();
                    empty=TRUE;
                    for(i=0; i<change_pillar->get_WGX_depth(); i++) {
                        if (change_pillar->get_WGX_gene(i) !=0) empty=FALSE;
                    }
                    if (empty==TRUE) change_pillar->get_outgroup_gene()->set_unused();
                }
                else {
                    //Add a homolog
                    rand_slot=ignuin(0,change_pillar->get_WGX_depth()-1);
                    rand_dupl=ignuin(0, change_pillar->get_my_homologs()->get_num_homologs()-1);
                    
                    while (change_pillar->get_WGX_gene(rand_slot) == change_pillar->get_my_homologs()->get_nth_homolog(rand_dupl))
                        rand_dupl=ignuin(0, change_pillar->get_my_homologs()->get_num_homologs()-1);
                    
                    curr_gene=change_pillar->get_WGX_gene(rand_slot);
                    if (curr_gene !=0)
                        curr_gene->set_unused();
                   
                    old_pillar=-1;
                    new_gene=change_pillar->get_my_homologs()->get_nth_homolog(rand_dupl);
                    if (new_gene->gene_used() == TRUE) {
                        old_slot=0;
                        old_pillar=new_gene->get_pillar();
                        while((old_slot < (*the_point->my_scaffold)[new_gene->get_pillar()].get_WGX_depth()) &&
                              ((*the_point->my_scaffold)[new_gene->get_pillar()].get_WGX_gene(old_slot) != new_gene)) old_slot++;
                        
                        if (old_slot == (*the_point->my_scaffold)[new_gene->get_pillar()].get_WGX_depth())
                            cerr<<"ERROR in pillar mapping: gene "<<new_gene->get_name()<<" has wrong pillar address\n";
                        else {
                            (*the_point->my_scaffold)[new_gene->get_pillar()].set_WGX_gene(0, old_slot);
                        }
                    }
                    
                    change_pillar->set_WGX_gene(new_gene,rand_slot);
                    new_gene->set_used(change_pillar->get_pillar_id());
                    
                    empty=TRUE;
                    for(i=0; i<change_pillar->get_WGX_depth(); i++) {
                        if (change_pillar->get_WGX_gene(i) !=0) empty=FALSE;
                    }
                    if (empty==TRUE) change_pillar->get_outgroup_gene()->set_unused();
                    
                    if ((old_pillar != -1)  && (&(*the_point->my_scaffold)[old_pillar] != change_pillar)) {
                        change_pillar=&(*the_point->my_scaffold)[old_pillar];
                        empty=TRUE;
                        for(i=0; i<change_pillar->get_WGX_depth(); i++) {
                            if (change_pillar->get_WGX_gene(i) !=0) empty=FALSE;
                        }
                        if (empty==TRUE) change_pillar->get_outgroup_gene()->set_unused();
                    }

                }

            }
        }
    }
    
    //for(i=0; i<the_point->the_homologs->get_num_homologs(); i++) cout<<new_order[i]<<" ";
    //cout<<endl;
    
    the_point->my_scaffold->change_pillar_order(the_point->pillar_order);
    the_point->my_scaffold->update_neighbors();
   // cout<<"Preopt neighbors: "<<the_point->my_scaffold->get_neighbor_count()<<endl;
    the_point->my_scaffold->optimize_assigns();
    the_point->my_scaffold->update_neighbors();
    the_point->my_scaffold->set_blocks();
    
    (*current_walkers)[walker]->new_score= -the_point->my_scaffold->get_neighbor_count();;


	//cout<<"Walker: "<<walker<<endl;
	//cout<<"Old score: "<<(*current_walkers)[walker]->score<<"\t";
	//cout<<"New score: "<<(*current_walkers)[walker]->new_score<<"\n";

}  //End Groups_anneal::move







void Scaffold_anneal::initialize_params(Space_point<Scaffold_point> *curr_condition)
{
    
	cout<<"Initial score: "<<-curr_condition->point->my_scaffold->get_neighbor_count()<<endl;
	
    curr_condition->new_score=-curr_condition->point->my_scaffold->get_neighbor_count();
    
    curr_condition->point->my_scaffold->get_homolog_set()->fix_genomes();
}


void Scaffold_anneal::after_move(int walker_num)
{
	

}


void Scaffold_anneal::output ()
{
    int i, j;
    ofstream fout;
    Pillar *my_pillar;
    
    cout<<"Final Score "<<best->score<<endl;
    fout.open(outfile.c_str());
    
    if (!(fout.fail())) {
        fout<<"Final Score "<<best->score<<endl;
        fout<<"OUTGROUPGENE\tOUTSYN";
        for(j=0; j<(*best->point->my_scaffold)[0].get_WGX_depth(); j++) {
            fout<<"\tDUPLGENE"<<j<<"\tDUPL"<<j<<"SYN";
        }
        fout<<endl;
        
        for(i=0; i<best->point->my_scaffold->get_homolog_set()->get_num_ances_genes(); i++) {
            my_pillar=&(*best->point->my_scaffold)[best->point->my_scaffold->get_order(i)];
            fout<<my_pillar->get_outgroup_gene()->get_name()<<"\t";
            if ((my_pillar->get_outgroup_neighbor(0) !=0)  &&
                (my_pillar->get_outgroup_neighbor(1) !=0) )
                fout<<"B";
            else {
                if (my_pillar->get_outgroup_neighbor(0) !=0) fout<<"L";
                else {
                    if (my_pillar->get_outgroup_neighbor(1) !=0) fout<<"R";
                    else fout<<"N";
                }
            }
            
            for(j=0; j<my_pillar->get_WGX_depth(); j++) {
                if (my_pillar->get_WGX_gene(j) != 0) {
                    fout<<"\t"<<my_pillar->get_WGX_gene(j)->get_name()<<"\t";
                    if ((my_pillar->get_WGX_neighbor(j,0) !=0) &&
                        (my_pillar->get_WGX_neighbor(j,1) !=0)) fout<<"B";
                    else {
                        if (my_pillar->get_WGX_neighbor(j,0) !=0) fout<<"L";
                        else {
                            if (my_pillar->get_WGX_neighbor(j,1) !=0) fout<<"R";
                            else fout<<"N";
                        }
                    }
                    if (my_pillar->get_WGX_gene(j)->get_num_homologs() ==1) fout<<"\tS";
                    else fout<<"\tM";
                    
                }
                else {
                    fout<<"\tNONE\t*\t*";
                }
            }
            fout<<endl;
            
        }
        fout.close();
    }
    else {
        cerr<<"ERROR: Cannot open output file "<<outfile<<endl;
    }
  

}
