#include "phylo_model_matrix.h"
#include <string>
#include <algorithm>
#include <cctype>
#include <sstream>


TransProb_Matrix::TransProb_Matrix (int n_states)
{
    int i, j;
    
    num_rates=1;
    num_states=n_states;
    
    transprobs =new double ** [num_rates];
    condprobs= new long double [num_states];
    for(i=0; i<num_rates; i++) {
        transprobs[i]=new double * [num_states];
        for (j=0; j<num_states; j++)
            transprobs[i][j]=new double [num_states];
    }
}


TransProb_Matrix::TransProb_Matrix (int n_states, int n_rates)
{
    int i, j;
    
    num_rates=n_rates;
    num_states=n_states;
    
    transprobs =new double ** [num_rates];
    condprobs=new long double [num_states];
    for(i=0; i<num_rates; i++) {
        transprobs[i]=new double * [num_states];
        for (j=0; j<num_states; j++)
            transprobs[i][j]=new double [num_states];
    }

}

TransProb_Matrix& TransProb_Matrix::operator=(TransProb_Matrix & assign_from)
{
    int i, j, k;
    
    for(i=0; i<num_rates; i++) {
        for(j=0; j<num_states; j++) {
            for (k=0; k<num_states; k++)
                transprobs[i][j][k]=assign_from.get_transprob(j,k,i);
        }
    }
    return(*this);
}


double TransProb_Matrix::get_transprob (int start, int end)
{
    return(transprobs[0][start][end]);
}


double TransProb_Matrix::get_transprob (int start, int end, int rate)
{
    return(transprobs[rate][start][end]);
}


void TransProb_Matrix::set_transprob(int start, int end, double val)
{
    transprobs[0][start][end]=val;
}


void TransProb_Matrix::set_transprob(int start, int end, int rate, double val)
{
    transprobs[rate][start][end]=val;
}


TransProb_Matrix::~TransProb_Matrix()
{
    int i, j;
    
    if (transprobs !=0 ) {
        for (i=0; i<num_rates; i++) {
            for (j=0; j<num_states; j++) delete[] transprobs[i][j];
            delete[] transprobs[i];
        }
        delete[] transprobs;
    }
    
    if (condprobs != 0) delete[] condprobs;
}

Model_State::Model_State()
{
    state_num=-1;
    redund_assigns=0;
    name="NONE";
    std::cerr<<"Error: call to default constructor of Model_State\n";
}


Model_State::Model_State(int s_id,  string new_name)
{
    state_num=s_id;
    name=new_name;
    num_redund=0;
    redund_assigns=0;
    num_redund_from=0;
    observed_state=TRUE;
    root_state=FALSE;
    state_level=0;
    positions_present=0;
}


Model_State::Model_State(int s_id,  string new_name, int s_level, int b_mask, int tlevels)
{
    int i, *pow2;
    
    state_num=s_id;
    name=new_name;
    num_redund=0;
    redund_assigns=0;
    num_redund_from=0;
    observed_state=TRUE;
    root_state=FALSE;
    state_level=s_level;
    binary_mask=b_mask;
    total_levels=tlevels;
    
    pow2=new int[total_levels];
    pow2[0]=1;
    for(i=1; i<total_levels; i++) pow2[i]=pow2[i-1]*2;
    
    positions_present=new BOOL [total_levels];
    
    for(i=0; i<total_levels; i++) {
        if (pow2[i] & binary_mask)
            positions_present[i]=TRUE;
        else
            positions_present[i]=FALSE;
    }
    delete[] pow2;
}



int Model_State::get_state_id()
{
    return(state_num);
}


string Model_State::get_state_name()
{
    return(name);
}


int Model_State::num_state_redunds()
{
    return(num_redund);
}


int Model_State::get_ith_redund_state(int redund_num)
{
    return(redund_assigns[redund_num]);
}


void Model_State::assign_redunds(std::map<int,int> redund_hash)
{
    int i;
    std::map<int, int>::iterator it;
    
    for (it=redund_hash.begin(); it!=redund_hash.end(); ++it) {
        if (it->second == state_num) num_redund++;
        
        if (it->first == state_num) num_redund_from++;
    }
    
    if (num_redund_from > 0) observed_state = FALSE;
    
    //Accounts for the state itself
    num_redund++;
    redund_assigns=new int[num_redund];

    redund_assigns[0]=state_num;
    
    i=1;
    for (it=redund_hash.begin(); it!=redund_hash.end(); ++it) {
        if (it->second == state_num) {
            redund_assigns[i]=it->first;
            i++;
        }
    }
}

void Model_State::initialize_cross_ref(int nc)
{
    int i;
    
    num_cross_refs=nc;
    cross_refs=new Model_State * [num_cross_refs];
    
    for(i=0; i<num_cross_refs; i++)
        cross_refs[i]=this;
}


void Model_State::set_cross_ref(int crnum, Model_State *other_state)
{
    cross_refs[crnum]=other_state;
}


Model_State::~Model_State()
{
    if (redund_assigns !=0) delete[] redund_assigns;
    if (positions_present!=0) delete[] positions_present;
}


Param_Set::Param_Set()
{
    vals=new double [1];
    vals[0]=0.0;
    real_set_size=1;
    is_global=TRUE;
    global_name="NULLPARAM";
    param_type=MINUS_INF_TO_INF;
}


Param_Set::Param_Set(string name, Numerical_Param_Type type, double init_val, int set_size, BOOL g)
{
    int i;
    
    is_global=g;
    
    param_type=type;
    global_name=name;
    
    if (is_global==TRUE) {
        vals=new double[1];
        real_set_size=1;
        vals[0]=init_val;
    }
    else {
        real_set_size=set_size;
        vals=new double[real_set_size];
        
        for(i=0; i<real_set_size; i++)
            vals[i]=init_val;
    }
}

Param_Set::Param_Set(string name, Numerical_Param_Type type, double* init_vals, int set_size, BOOL g)
{
    int i;
    
    is_global=g;
    
    param_type=type;
    global_name=name;
    
    if (is_global==TRUE) {
        vals=new double[1];
        real_set_size=1;
        vals[0]=init_vals[0];
    }
    else {
        real_set_size=set_size;
        vals=new double[real_set_size];
        
        for(i=0; i<real_set_size; i++)
            vals[i]=init_vals[i];
    }
}


string Param_Set::get_global_param_name()
{
    return(global_name);
}



string Param_Set::get_param_set_name(int set_id)
{
    char temp_num[30];
   string retval;
    
    if (is_global==TRUE) retval=global_name;
    else {
        retval =global_name + "_";
        int_to_string(temp_num, 29, set_id);
        retval=retval+temp_num;
    }
    
    return(retval);
}




double Param_Set::get_parameter(int set_num)
{
    if (is_global==TRUE) return(vals[0]);
    else return(vals[set_num]);
}


double Param_Set::get_scaled_parameter(int my_set)
{
    double my_val;
    
    if (is_global==TRUE) my_val=vals[0];
    else my_val=vals[my_set];
    
    switch (param_type) {
        case MINUS_INF_TO_INF:
            return(my_val);
        case ZERO_TO_INF:
            return(fabs(my_val));
        case ONE_TO_INF:
        case ONE_PLUS_TO_INF:
            return(fabs(my_val)-1);
        case ZERO_TO_ONE:
            return(log(my_val));
    }
    
    return(0);
}


void Param_Set::set_parameter(double scaled_val, int set_id)
{
    switch (param_type) {
        case MINUS_INF_TO_INF:
            if (is_global == TRUE)
                vals[0]=scaled_val;
            else
                vals[set_id]=scaled_val;
            break;
        case ZERO_TO_INF:
            if (is_global == TRUE)
                vals[0]=fabs(scaled_val);
            else
                vals[set_id]=fabs(scaled_val);
            break;
        case ONE_TO_INF:
        case ONE_PLUS_TO_INF:
            if (is_global == TRUE)
                vals[0]=fabs(scaled_val)+1.0;
            else
                vals[set_id]=fabs(scaled_val)+1.0;
            break;
        case ZERO_TO_ONE:
            if (is_global == TRUE)
                vals[0]=exp(-1.0*fabs(scaled_val));
            else
                vals[set_id]=exp(-1.0*fabs(scaled_val));
            break;
    }
}


Param_Set::~Param_Set()
{
    
    if (vals !=0) delete[] vals;
}

Phylo_Matrix::Phylo_Matrix()
{
    std::cerr<<"Error: Call to default constructor of Phylo_Matrix\n";
    num_states=0;
    num_parameters=0;
    matrix_descript=0;
    matrix_params=0;
    rate_matrix=0;
    receive_matrix=0;
    the_parameters=0;
    the_states=0;
    curr_exchange=0;
    the_matrices=0;
    curr_lin_alg=0;
    num_root_states=0;
    level_states=0;
    brnlen_one=FALSE;
    num_param_sets=1;
    num_states_by_level=0;
    total_states_by_level=0;
    is_symm=FALSE;
    
}



Phylo_Matrix::Phylo_Matrix(Exchange *cexchange, string model_file)
{
    int i, j, k, l, state_level, cnt, pow2[32], *bin_string=0, binary_mask, setcnt, state_cnt,
        c1_state, c2_state, con1_state, con2_state, dupl_state;
    double val, *vals=0;
    BOOL swap=FALSE, global, use_set, same, foundcopy;
    string range_name, dummy, temp_name, param_name, search_param, descript, redund_master, redund_child, state_name, limit_type, model_type, param_branch_range, param_line;
    Numerical_Param_Type temp_type;
    Param_Set *temp_set;
    ifstream fin;
    std::map<string, int> state_ids;
    std::map<int, int> redund_hash;
    std::size_t found;
    
    curr_exchange=cexchange;
    brnlen_one=FALSE;
    is_symm=FALSE;
    
    num_param_sets=1;
    
    pow2[0]=1;
    for(i=1; i<32; i++) pow2[i]=pow2[i-1]*2;
    
    
    fin.open(model_file.c_str());
    
    if (!fin.fail()) {
        rate_matrix_invalid=TRUE;
        num_root_states=0;
        fin>>model_name;
        fin>>model_type;
        std::cout<<"Using model "<<model_name<<" and type "<<model_type<<std::endl;
        std::transform(model_type.begin(), model_type.end(), model_type.begin(), ::tolower);

        if (model_type == "hierarchical") {
            hierarchical_model =TRUE;
            fin>>num_levels;
            bin_string=new int[num_levels];
            for(i=0; i<num_levels; i++) bin_string[i]=0;
            std::cout<<"Model has "<<num_levels<<" levels"<<std::endl;
        }
        else hierarchical_model =FALSE;
        
        fin>>descript;
        
        std::transform(descript.begin(), descript.end(), descript.begin(), ::tolower);
        //cout<<"Next description is "<<descript<<endl;
        if(descript == "#numberofparametersets") {
            fin>>num_param_sets;
            cout<<"Model is branch-specific with: "<<num_param_sets<<" sets of model catagories\n";
            
            vals=new double[num_param_sets];
            
            fin>>descript;
        }
        
        fin>>num_states;
        std::cout<<"Model has "<<num_states<<" states"<<std::flush<<std::endl;
        fin>>descript;
        
        fin>>num_params_per_set;
        std::cout<<"Model has "<<num_params_per_set<<" parameters per parameter set and "<<num_param_sets<<" sets"<<std::flush<<std::endl;
        the_states=new Model_State*[num_states];
        the_parameters=new Param_Set*[num_params_per_set];
        
        
       
        matrix_descript=new Entry_Type * [num_states];
        matrix_params=new Param_Set ***[num_states];
        num_params_per_entry=new int *[num_states];
        rate_matrix=new double* [num_states];
        receive_matrix = new double* [num_states];
        for(i=0; i<num_states; i++) {
            matrix_descript[i]=new Entry_Type[num_states];
            matrix_params[i]=new Param_Set **[num_states];
            num_params_per_entry[i]=new int [num_states];
            rate_matrix[i] = new double [num_states];
            receive_matrix[i]=new double[num_states];
        }
        
        the_matrices =new TransProb_Matrix * [curr_exchange->get_num_branches()];
        
        for(i=0; i<curr_exchange->get_num_branches(); i++)
            the_matrices[i]=new TransProb_Matrix(num_states,curr_exchange->get_num_rates());
        
        curr_lin_alg=new  Generic_Linear_Algebra(num_states);
        
        fin>>descript;
        for (i=0; i<num_states; i++) {
            if (hierarchical_model == FALSE) {
                fin>>state_name;
                the_states[i]=new Model_State(i, state_name);
            }
            else {
                fin>>state_name;
                for(j=0; j<num_levels; j++) fin>>bin_string[j];
                binary_mask=0;
                state_level=0;
                //std::cout<<state_name<<"\t";
                for (j=0; j<num_levels; j++) {
                    //std::cout<<j<<"&"<<bin_string[j]<<"\t";
                    if (bin_string[j] ==1) {
                        binary_mask += pow2[j];
                        state_level++;
                    }
                }
                //std::cout<<"\nRead state "<<state_name<<" with mask "<<binary_mask<<" and level "<<state_level<<std::endl;
                the_states[i]=new Model_State(i, state_name, state_level, binary_mask, num_levels);
            }
            
            state_ids[state_name]=i;
        }
        fin>>descript;
        //cout<<"Now reading "<<descript<<endl;
        for (i=0; i<num_params_per_set; i++) {
            use_set=FALSE;
            
            param_branch_range="";
            if (num_param_sets ==1) {
                fin>>param_name>>range_name>>val;
                global=TRUE;
            }
            else {
                if (i ==0) std::getline(fin, param_line);
                
                std::getline(fin, param_line);
                //cout<<"Line is "<<param_line<<endl;
                //fin>>param_name>>range_name>>val>>param_branch_range;
            
                std::stringstream ss(param_line);
                ss >>param_name>>range_name>>val>>param_branch_range;
                
                std::transform(param_branch_range.begin(), param_branch_range.end(), param_branch_range.begin(), ::tolower);
                
               //std::cout<<"Initial: read parameter "<<i<<": |"<<param_name<<"| range: "<<range_name<<" Value: "<<val<<" Global "<<global<<endl;
                
                if (param_branch_range == "global") global=TRUE;
                else {
                    global=FALSE;
                    for(setcnt=0; setcnt<num_param_sets; setcnt++) vals[setcnt]=val;
                    setcnt=1;
                    while(ss>>val) {
                        std::cout<<"Read extra value "<<setcnt<<" as "<<val<<endl;
                        vals[setcnt]=val;
                        setcnt++;
                    }
                    
                    use_set=TRUE;
                    
                    
                    
                }
            }
            
            std::cout<<"Read parameter "<<i<<": |"<<param_name<<"| range: "<<range_name<<" Value: "<<val<<" Global "<<global<<endl;
            temp_type=get_param_type(range_name);

            if (use_set ==FALSE) the_parameters[i]=new Param_Set(param_name, temp_type, val, num_param_sets, global);

            else the_parameters[i]=new Param_Set(param_name, temp_type, vals, num_param_sets, global);
            
        }
        
        //Hack to make sure switch prob is last for WGD/WGX models
        for (i=0; i<num_params_per_set-1; i++) {
            if (the_parameters[i]->get_global_param_name() == "SwitchProb") {
                swap=TRUE;
                temp_set=the_parameters[i];
                the_parameters[i]=the_parameters[num_params_per_set-1];
            }
        }
        if (swap == TRUE) {
            the_parameters[num_params_per_set-1]=temp_set;
        }
        
        if (vals !=0) delete[] vals;
        
        initialize_linear_param_list();
        
        fin>>descript;
        //Read the matrix by names
        for (i=0; i<num_states; i++) {
            for (j=0; j<num_states; j++) {
                if (i==j) {
                    fin>>dummy;
                    num_params_per_entry[i][j]=1;
                    matrix_params[i][j]=new Param_Set* [1];
                    matrix_params[i][j][0]=0;
                    matrix_descript[i][j]=ONE;
                }
                else {
                    fin>>param_name;
                    if(param_name != "Default") {
                        if (param_name != "Zero") {
                            matrix_descript[i][j]=PARAM;
                            search_param=param_name;
                            found=search_param.find("*");

                            if (found ==std::string::npos) {
                                num_params_per_entry[i][j]=1;
                                matrix_params[i][j]=new Param_Set * [1];
                                matrix_params[i][j][0]=get_param(param_name);
                               // cout<<"Entry "<<i<<", "<<j<<" has a single parameter: "<<param_name<<" == "<<matrix_params[i][j][0]->get_global_param_name()<<endl;
                            }
                            else {
                                num_params_per_entry[i][j]=0;
                               // cout<<"Entry "<<i<<", "<<j<<" has several parameters: "<<param_name<<endl;
                                search_param=param_name;
                                while((found=search_param.find("*")) !=std::string::npos)
                                {
                                    num_params_per_entry[i][j]++;
                                    search_param=search_param.substr(found+1, string::npos);
                                }
                                num_params_per_entry[i][j]++;
                                
                               // cout<<"Found "<<num_params_per_entry[i][j]<<" parameters in entry "<<param_name<<endl;
                                matrix_params[i][j]=new Param_Set * [num_params_per_entry[i][j]];
                                cnt=0;
                                search_param=param_name;
                                while((found=search_param.find("*")) !=std::string::npos)
                                {
                                    param_name=search_param.substr(0, found);
                                    matrix_params[i][j][cnt]=get_param(param_name);
                                    
                                    search_param=search_param.substr(found+1, string::npos);
                                   // cout<<"Setting parameter "<<cnt<<" to "<<param_name<<" (remaining: "<<search_param<<")"<<endl;
                                    cnt++;
                                }
                                matrix_params[i][j][cnt]=get_param(search_param);
                               // cout<<"Setting last parameter "<<cnt<<" to "<<search_param<<endl;
                                
                            }
                        }
                        else {
                            matrix_descript[i][j]=ZERO;
                            num_params_per_entry[i][j]=1;
                            matrix_params[i][j]=new Param_Set * [1];
                            matrix_params[i][j][0]=0;
                        }
                    }
                    else {
                        matrix_descript[i][j]=ONE;
                        num_params_per_entry[i][j]=1;
                        matrix_params[i][j]=new Param_Set * [1];
                        matrix_params[i][j][0] = 0;
                    }
                }
            }
        }
        
        if (num_levels == 2) {
            same=TRUE;
            foundcopy=TRUE;
            
            std::map<string,int>::iterator it;
            
            it = state_ids.find("Copy1");
            if (it != state_ids.end())
                c1_state=state_ids["Copy1"];
            else
                foundcopy=FALSE;
            
            it = state_ids.find("Copy2");
            if (it != state_ids.end())
                c2_state=state_ids["Copy2"];
            else
                foundcopy=FALSE;
            
            it = state_ids.find("Dupl");
            if (it != state_ids.end())
                dupl_state=state_ids["Dupl"];
            else
                foundcopy=FALSE;
            
            con1_state=-1;
            con2_state=-1;
            it = state_ids.find("DuplC1");
            if (it != state_ids.end())
                con1_state=state_ids["DuplC1"];
            
            
            it = state_ids.find("DuplC2");
            if (it != state_ids.end())
                con2_state=state_ids["DuplC2"];
            
            
            if (foundcopy==TRUE) {
                if (num_params_per_entry[dupl_state][c1_state] != num_params_per_entry[dupl_state][c2_state])
                    same=FALSE;
                else {
                    for(l=0; l<num_params_per_entry[dupl_state][c1_state]; l++) {
                        if (matrix_params[dupl_state][c1_state][l] != matrix_params[dupl_state][c2_state][l]) same=FALSE;
                    }
                }
                if (con1_state != -1) {
                    if (num_params_per_entry[con1_state][c1_state] != num_params_per_entry[con2_state][c2_state])
                        same=FALSE;
                    else {
                        for(l=0; l<num_params_per_entry[dupl_state][c1_state]; l++) {
                            if (matrix_params[con1_state][c1_state][l] != matrix_params[con2_state][c2_state][l]) same=FALSE;
                        }
                    }
                }
                
            }
            
            if ((same == TRUE) &&(foundcopy == TRUE)) is_symm=TRUE;
         
        }
        
       
        
        fin>>descript;
        redund_hash.clear();

        while(! (fin.eof())) {
            fin>>limit_type;
            std::transform(limit_type.begin(), limit_type.end(), limit_type.begin(), ::tolower);
            if (limit_type == "redundancy") {
                fin>>redund_master>>redund_child;
                 if (!(fin.eof())) {
                     redund_hash[state_ids[redund_child]]=state_ids[redund_master];
                 }
            }
            if (limit_type == "rootstate") {
                fin>>state_name;
                if (!(fin.eof())) {
                    cout<<"Setting root state to "<<state_name<<" =="<<state_ids[state_name]<<endl;
                    the_states[state_ids[state_name]]->make_root_state();
                }
            }
        }
        fin.close();
              
        for (i=0; i<num_states; i++) the_states[i]->assign_redunds(redund_hash);
       
        for (i=0; i<num_states; i++)
            if (the_states[i]->is_root_state() == TRUE) num_root_states++;
            
        
        
        if (num_root_states ==0) {
            for (i=0; i<num_states; i++) {
                the_states[i]->make_root_state();
                num_root_states++;
            }
        }
        
        
        root_states=new Model_State *[num_root_states];
        state_cnt=0;
        for(i=0; i<num_states; i++) {
            if (the_states[i]->is_root_state() == TRUE) {
                root_states[state_cnt]=the_states[i];
                state_cnt++;
            }
        }
        
        
        i=0;
        while (the_states[i]->is_root_state() == FALSE) i++;
        first_root_state=the_states[i];
        
        if (hierarchical_model == TRUE) {
            level_states=new Model_State ** [num_levels+1];
            full_level_states=new Model_State ** [num_levels+1];
            
            num_states_by_level=new int [num_levels+1];
            total_states_by_level = new int [num_levels+1];
            for (i=0; i<=num_levels; i++) {
                cnt=0;
                for(j=0; j<num_states; j++) {
                    if ((the_states[j]->get_state_level() == i) &&(the_states[j]->is_observed_state() )) cnt++;
                }
                num_states_by_level[i]=cnt;
                
                cnt=0;
                for(j=0; j<num_states; j++) {
                    if (the_states[j]->get_state_level() == i) cnt++;
                }
                total_states_by_level[i]=cnt;
            
                
                level_states[i]=new Model_State*[num_states_by_level[i]];
                cnt=0;
                for(j=0; j<num_states; j++) {
                    if ((the_states[j]->get_state_level() == i) &&(the_states[j]->is_observed_state() ))  {
                        level_states[i][cnt]=the_states[j];
                        the_states[j]->set_level_state_id(cnt);
                        cnt++;
                    }
                }
                
                full_level_states[i] = new Model_State*[total_states_by_level[i]];
                cnt=0;
                for(j=0; j<num_states; j++) {
                    if (the_states[j]->get_state_level() == i) {
                        full_level_states[i][cnt]=the_states[j];
                        cnt++;
                    }
                }
            }
            //for (i=0; i<=num_levels; i++) {
            //    for(j=0; j<num_states_by_level[i]; j++) cout<<"Level "<<i<<" state number "<<j<<" is "<<level_states[i][j]->get_state_name()<<endl;
            //}
            
            setup_masked_states();
            delete[] bin_string;
        }
        
    }
    else {
        std::cerr<<"Error: Model description file "<<model_file<<" does not exist\n";
        num_states=0;
        matrix_descript=0;
        matrix_params=0;
        rate_matrix=0;
        receive_matrix=0;
        num_parameters=0;
        num_params_per_set=0;
        the_states=0;
        the_parameters=0;
        curr_exchange=0;
        the_matrices=0;
        curr_lin_alg=0;
        receive_matrix=0;
        level_states=0;
    }
    

}

TransProb_Matrix* Phylo_Matrix::get_tp_matrix_num(int matrix_num)
{
    return(the_matrices[matrix_num]);
}

Model_State * Phylo_Matrix::get_nth_level_ith_state(int level, int position)
{
    return(level_states[level][position]);
}

Model_State * Phylo_Matrix::get_ith_full_state_level(int level, int position)
{
    return(full_level_states[level][position]);
}

    
Model_State * Phylo_Matrix::get_masked_state(int mask)
{
    return(masked_states[mask]);
}

double Phylo_Matrix::get_param(int param_num)
{
    return(the_parameters[param_set_lookup[param_num]]->get_parameter(param_set_id[param_num]));
}


double Phylo_Matrix::get_scaled_param(int param_num)
{
    return(the_parameters[param_set_lookup[param_num]]->get_scaled_parameter(param_set_id[param_num]));
}

void Phylo_Matrix::set_param(int param_num, double scaled_val)
{
    the_parameters[param_set_lookup[param_num]]->set_parameter(scaled_val, param_set_id[param_num]);
    rate_matrix_invalid=TRUE;
}


string Phylo_Matrix::get_param_name(int param_num)
{
    return(the_parameters[param_set_lookup[param_num]]->get_param_set_name(param_set_id[param_num]));
}


int Phylo_Matrix::get_num_level_states (int level)
{
    return(num_states_by_level[(level%(num_levels+1))]);
}

int Phylo_Matrix::get_num_full_level_states (int level)
{
    return (total_states_by_level[(level%(num_levels+1))]);
}


Entry_Type Phylo_Matrix::get_trans_type(int state1, int state2)
{
    return(matrix_descript[abs(state1)%num_states][abs(state2)%num_states]);
}

int Phylo_Matrix::num_trans_params(int state1, int state2)
{
    return(num_params_per_entry[abs(state1)%num_states][abs(state2)%num_states]);
}

Param_Set* Phylo_Matrix::get_nth_trans_param(int state1, int state2, int param_num)
{
    return(matrix_params[abs(state1)%num_states][abs(state2)%num_states][abs(param_num)%num_params_per_entry[abs(state1)%num_states][abs(state2)%num_states]]);
}

void Phylo_Matrix::initialize_rate_matrix(int set_id)
{
    int i, j, k;
    double sum;
    
    for(i=0; i<num_states; i++) {
        sum=0;
        for(j=0; j<num_states; j++) {
            if (i != j) {
                switch (matrix_descript[i][j]) {
                    case PARAM:
                        rate_matrix[i][j]=1.0;
                        for (k=0; k<num_params_per_entry[i][j]; k++)
                            rate_matrix[i][j] *= matrix_params[i][j][k]->get_parameter(set_id);
                        break;
                    case ZERO:
                        rate_matrix[i][j]=0.0;
                        break;
                    case ONE:
                        rate_matrix[i][j]=1.0;
                        break;
                }
                
                
                sum+=rate_matrix[i][j];
            }
        }
        rate_matrix[i][i]=-1.0*sum;
    }
    
    curr_lin_alg->find_eigen_matrix(rate_matrix);
    rate_matrix_invalid=FALSE;
}



void Phylo_Matrix::calc_transprobs(Branch *my_branch)
{
    int i, j, k;
    
    //std::cout<<"Calculating tps for "<<my_branch->get_name()<<" using length "<<my_branch->get_brnlen()<<endl;
    
    //if (rate_matrix_invalid == TRUE) std::cout<<"Re-initializing rate matrix\n";
    if (brnlen_one==TRUE) {
        if ((rate_matrix_invalid == TRUE) || (num_param_sets > 1)) initialize_rate_matrix((my_branch->get_param_set()%num_param_sets));
        
        for(i=0; i<curr_exchange->get_num_rates(); i++) {
            curr_lin_alg->find_scaled_matrix(receive_matrix, 1.0, 0);
        }
        
        for(i=0; i<curr_exchange->get_num_rates(); i++) {
            for(j=0; j<num_states; j++) {
                for(k=0; k<num_states; k++) {
                    the_matrices[my_branch->get_brn_num()]->set_transprob(j,k,i, receive_matrix[j][k]);
                    //cout<<the_matrices[my_branch->get_brn_num()]->get_transprob(j,k)<<"\t";
                }
                // cout<<endl;
            }
        }

    }
    else {
        if (my_branch->get_brnlen()  > MIN_BRLEN) {
            if ((rate_matrix_invalid == TRUE) || (num_param_sets > 1)) initialize_rate_matrix((my_branch->get_param_set()%num_param_sets));
            
            for(i=0; i<curr_exchange->get_num_rates(); i++) {
                curr_lin_alg->find_scaled_matrix(receive_matrix, my_branch->get_brnlen(), 0);
            }
            
            for(i=0; i<curr_exchange->get_num_rates(); i++) {
                for(j=0; j<num_states; j++) {
                    for(k=0; k<num_states; k++) {
                        
                        the_matrices[my_branch->get_brn_num()]->set_transprob(j,k,i, receive_matrix[j][k]);
                        //cout<<the_matrices[my_branch->get_brn_num()]->get_transprob(j,k)<<"\t";
                    }
                    //cout<<endl;
                }
            }
            
        }
        else {
            my_branch->set_brnlen(MIN_BRLEN);
            for(i=0; i<curr_exchange->get_num_rates(); i++) {
                for(j=0; j<num_states; j++) {
                    for(k=0; k<num_states; k++) {
                        if (j==k) the_matrices[my_branch->get_brn_num()]->set_transprob(j,k,i, 1.0);
                        else the_matrices[my_branch->get_brn_num()]->set_transprob(j,k,i, 0.0);
                    }
                }
            }
        }
    }
}




void Phylo_Matrix::calc_all_transprobs(Tree *curr_tree)
{
    int i;
    
    for(i=0; i<curr_exchange->get_num_branches(); i++) calc_transprobs((*curr_tree)[i]);
}


Phylo_Matrix::~Phylo_Matrix()
{
    int i, j;
    
    if (the_parameters != 0) {
        for(i=0; i<num_param_sets; i++) delete the_parameters[i];
        delete[] the_parameters;
    }
    if (the_states !=0) {
        for(i=0; i<num_states; i++)
            delete the_states[i];
        delete[] the_states;
    }
    
    if (matrix_params != 0) {
        for(i=0; i<num_states; i++)  {
            for( (j=0); j<num_states; j++)
                delete matrix_params[i][j];
            delete[] num_params_per_entry[i];
            delete[] matrix_params[i];
        }
        delete[] matrix_params;
    }
 
    if (matrix_descript !=0) {
        for(i=0; i<num_states; i++) delete[] matrix_descript[i];
        delete[] matrix_descript;
    }
    if (rate_matrix !=0) {
        for(i=0; i<num_states; i++) delete[] rate_matrix[i];
        delete[] rate_matrix;
    }
    
    if (receive_matrix !=0) {
        for(i=0; i<num_states; i++) delete[] receive_matrix[i];
        delete[] receive_matrix;
    }
    
    if (the_matrices != 0) delete[] the_matrices;
    if (curr_lin_alg !=0) delete curr_lin_alg;
    
    if (level_states !=0) {
        for (i=0; i<=num_levels; i++) delete[] level_states[i];
        delete[] level_states;
    }
    
    if (num_states_by_level !=0) delete[] num_states_by_level;
    if (total_states_by_level !=0) delete[] total_states_by_level;
}



Numerical_Param_Type Phylo_Matrix::get_param_type(string type_name)
{
    Numerical_Param_Type retval=MINUS_INF_TO_INF;
    
    if (type_name == "MINUS_INF_TO_INF") retval=MINUS_INF_TO_INF;
    
    if (type_name == "ZERO_TO_INF") retval=ZERO_TO_INF;
    
    if (type_name == "ONE_TO_INF") retval=ONE_TO_INF;
    
    if (type_name == "ZERO_TO_ONE") retval =ZERO_TO_ONE;
    
    return(retval);
}



Param_Set * Phylo_Matrix::get_param(string param_name)
{
    int i=0;
    
    //cout<<"Looking for parameter: "<<param_name<<endl;
    while((i<num_params_per_set) && (param_name != the_parameters[i]->get_global_param_name())) i++;
    
    //cout<<"Found position "<<i<<": "<<the_parameters[i]->get_global_param_name()<<endl;
    if (i == num_params_per_set) return(0);
    else return(the_parameters[i]);
}


void Phylo_Matrix::setup_masked_states()
{
    int i,  pow2[32], max_mask;
    
    pow2[0]=1;
    for (i=1; i<32; i++) pow2[i]=2*pow2[i-1];
    
    max_mask=0;
    for(i=0; i<num_states; i++) if (max_mask < the_states[i]->get_binary_rep()) max_mask=the_states[i]->get_binary_rep();
    
    i=0;
    while(pow2[i] < max_mask) i++;
    mask_array_size=pow2[i];
    
    masked_states=new Model_State* [mask_array_size];
    
    for(i=0; i<mask_array_size; i++) masked_states[i]=0;
    
    for(i=0; i<num_states; i++) {
            if (the_states[i]->is_observed_state() == TRUE)
                masked_states[the_states[i]->get_binary_rep()]= the_states[i];
    }
    
}


void Phylo_Matrix::initialize_all_cross_refs(int level, int nc)
{
    int i;
    
    for(i=0; i<num_states; i++) {
        if (the_states[i]->get_state_level() == level) {
            the_states[i]->initialize_cross_ref(nc);
        }
    }
}



void Phylo_Matrix::initialize_linear_param_list()
{
    int i, j, para_id;
    
    num_parameters=0;
    
    for(i=0; i<num_params_per_set; i++) {
        if (the_parameters[i]->is_global_parameter() == TRUE) num_parameters++;
        else num_parameters+=num_param_sets;
    }
    param_set_lookup= new int[num_parameters];
    param_set_id = new int [num_parameters];
    
    cout<<"Linear parameter list has "<<num_parameters<<" elements\n";
    
    para_id=0;
    for(i=0; i<num_params_per_set; i++) {
        if (the_parameters[i]->is_global_parameter() == TRUE) {
            param_set_lookup[para_id]=i;
            param_set_id[para_id]=0;
            para_id++;
        }
        else {
            for(j=0; j<num_param_sets; j++) {
                param_set_lookup[para_id]=i;
                param_set_id[para_id]=j;
                para_id++;
            }
        }
    }
}


Branch_Ex::Branch_Ex() : Branch()
{
    std::cerr<<"Error: call to default construction of Branch_Ex class\n";
    my_transprobs=0;
}


Branch_Ex::Branch_Ex(Exchange *cexchange, TransProb_Matrix *new_mat, int b_num) : Branch(cexchange, b_num)
{
    my_transprobs=new_mat;
    extern_trpb=TRUE;
}

Branch_Ex  & Branch_Ex::operator=(Branch_Ex & assign_from)
{
    int i, j, k;
    
    set_name(assign_from.get_name());
    
    brlen=assign_from.get_brnlen();
    expect_numsubs_site=assign_from.expect_subs_site();
    has_name=assign_from.branch_has_name();
    
    uninitialized=assign_from.is_uninitialized();
    param_set_id=assign_from.get_param_set();
    
    tip=assign_from.is_tip();
    if (tip == TRUE)
        children[0]=children[1]=0;
    
    taxa_id=assign_from.get_taxa_id();
    pns_num=assign_from.get_p_nonsyn_num();
    pitg_num=assign_from.get_p_intergroup_num();
    for(i=0; i<curr_exchange->get_num_aa_props()+1; i++)
        aa_prop_num[i]=assign_from.get_aa_prop_num(i);
    
    *(my_transprobs) = *(assign_from.get_tb_matrix());
    
    return(*this);
}


double Branch_Ex::get_trpb(int rate, int start, int end)
{
    return(my_transprobs->get_transprob(start, end, rate));
}


void Branch_Ex::set_trpb(int rate, int start, int end, double trpb)
{
    std::cerr<<"Error: call to invalid set_trpb from Branch_Ex object\n";
}


long double Branch_Ex::get_cond_prob (int state)
{
    return(my_transprobs->get_cond_prob(state));
}



void Branch_Ex::set_cond_prob (int state, long double prob)
{
    my_transprobs->set_cond_prob(state, prob);
}


Tree_Ex::Tree_Ex() : Tree()
{
    tree_Ex =0;
    std::cerr<<"Error: Call to default constructor of Tree_Ex object\n";
}


Tree_Ex::Tree_Ex(Exchange *cexchange, Phylo_Matrix *my_mat, BOOL rooted) : Tree()
{
    int i;
    curr_exchange=cexchange;
    
    //We flag that this tree contains its own local data
    tree_is_local_mem=TRUE;
    is_rooted_tree=rooted;
    if(is_rooted_tree==TRUE)
        null_root_brlns=FALSE;
    else
        null_root_brlns=TRUE;
    
    start_extra_brn=new Branch_list;
    extra_brn=start_extra_brn;
    start_tips=new Branch_list;
    tips=start_tips;
    
    start_tips->next=0;
    start_extra_brn->next=0;
    
    
    
    if (curr_exchange->get_num_taxa() <= 3)
        three_taxa_tree=TRUE;
    else
        three_taxa_tree=FALSE;
    
    prune_root=0;
    constrain_brn=0;
    
    //Declares a new array of branches large enough for the tree
    tree = new Branch * [curr_exchange->get_num_branches()];
    tree_Ex =new Branch_Ex *[curr_exchange->get_num_branches()];
    for (i=0; i<curr_exchange->get_num_branches(); i++) {
        tree_Ex[i]=new Branch_Ex(curr_exchange, my_mat->get_tp_matrix_num(i), i);
        tree[i]=tree_Ex[i];
    }
}


Tree_Ex & Tree_Ex::operator= (Tree_Ex & assign_from)
{
    int i, j;
    
    for (i=0; i<curr_exchange->get_num_branches(); i++) {
        *(tree_Ex[i])=*(assign_from.get_nth_branch(i));
        tree[i]=tree_Ex[i];
        
        if (assign_from[i]!=assign_from.find_root())   {
            //Added a bunch of checks so that we can copy trees that are partially assembled:
            //i.e. where some branches do not have parents/siblings/children yet
            j=0;
            
            if(assign_from[i]->get_sibling() != 0 ){
                while (assign_from[j]!=assign_from[i]->get_sibling()) j++;
                set_as_siblings(tree[i], tree[j]);
            }
            
            if (tree[i]->is_tip()==FALSE) {
                j=0;
                
                if(assign_from[i]->get_child(0) != 0) {
                    while (assign_from[j]!=assign_from[i]->get_child(0)) j++;
                    set_as_parent_child(tree[i], tree[j], 0);
                }
                
                j=0;
                if(assign_from[i]->get_child(1) != 0) {
                    while (assign_from[j]!=assign_from[i]->get_child(1)) j++;
                    set_as_parent_child(tree[i], tree[j], 1);
                }
            }
            j=0;
            
            if (assign_from[i]->get_parent() != 0) {
                while (assign_from[j]!=assign_from[i]->get_parent()) j++;
                if (assign_from[i]->get_parent()->get_child(0)==assign_from[i])
                    set_as_parent_child(tree[j], tree[i], 0);
                else
                    set_as_parent_child(tree[j], tree[i], 1);
            }
        }
        
        else {
            j=0;
            
            if (assign_from[i]->get_child(0) != 0) {
                while (assign_from[j]!=assign_from[i]->get_child(0)) j++;
                set_as_parent_child(tree[i], tree[j], 0);
            }
            j=0;
            
            if (assign_from[i]->get_child(1) != 0) {
                while (assign_from[j]!=assign_from[i]->get_child(1)) j++;
                set_as_parent_child(tree[i], tree[j], 1);
            }
            tree[i]->null_parent();
            tree[i]->null_sibling();
        }
    }
    set_root(tree[0]);
    
    return(*this);
}

Branch_Ex * Tree_Ex::get_nth_branch(int n)
{
    return(tree_Ex[n]);
}

Branch* Read_PAUP_Tree_Ex::initialize_branch ()
//This function is used to create new internal branches.  It finds the first
//"unused" branch in the current tree array, gives it the next unused id number
//and returns a pointer to it.

{
    
    return(curr_tree->initialize_branch());
}  //End Read_PAUP_Tree::initialize_branch



Tree_Ex * Read_PAUP_Tree_Ex::create_tree_from_file (Exchange *cexchange, Sequence_dataset *curr_data, Phylo_Matrix *curr_matrix, BOOL rooted)
{
    Branch  *curbranch;
    ifstream treein;
    char dump;
    string line, read_name, val_string;
    int i,j=0, curid;
 
    
    curr_exchange=cexchange;
    
    curr_tree=new Tree_Ex(curr_exchange, curr_matrix, rooted);
    
    treein.open(curr_exchange->get_treefile());
    if (treein.fail())
    {
        std::cerr<<"Cannot find file "<<curr_exchange->get_treefile()<<endl;
        return(0);
    }
    else
    {
        
        treein>>line;
        std::transform(line.begin(), line.end(), line.begin(), ::tolower);
       // std::cout<<"Read: "<<line<<endl;
        while (line.find("translate") == std::string::npos) {
            treein>>line;
            std::transform(line.begin(), line.end(), line.begin(), ::tolower);
            //std::cout<<"Read: "<<line<<endl;
           // std::cout<<"Find returns "<<line.find("translate")<<endl;
        }
        
        
        //Uses Paup's translate table to name the tip branches
        for (i=0; i<curr_exchange->get_num_taxa(); i++)
        {
            
            treein>>curid>>read_name;
            if (read_name.find(",") != std::string::npos)
                read_name=read_name.substr(0, read_name.length()-1);
            if (read_name[0] == '\'') read_name=read_name.substr(1, read_name.length());
            if (read_name.find("\'") != std::string::npos)
                read_name=read_name.substr(0, read_name.length()-1);
            //std::cout<<"Read "<<read_name<<" for taxa "<<i<<endl;

            
            (*curr_tree)[i]->set_name(read_name.c_str());
            (*curr_tree)[i]->initialized();
            
            if (curr_exchange->have_data()==TRUE)
            {
                j=0;
                
                while ((strcmp((*curr_data)[j].Sequence_name(), read_name.c_str())!=0) && (j<curr_exchange->get_num_taxa()))
                    j++;
                
                if (j==curr_exchange->get_num_taxa())
                    std::cerr<<"Number of taxa exceeded and sequence "<<read_name<<" still not found in sequence file\n";
                
                (*curr_tree)[i]->set_taxa_id(j);
            }
            else
                (*curr_tree)[i]->set_taxa_id(i);
            
        }
        
        treein.get(dump);
        treein>>line;
        //std::cout<<"Read: "<<line<<endl;
        treein>>line;
        //std::cout<<"Read: "<<line<<endl;
        treein>>line;
        //std::cout<<"Read: "<<line<<endl;
        treein>>line;
        //std::cout<<"Read: "<<line<<endl;
        treein>>line;
        //std::cout<<"Read: "<<line<<endl;
        
        if (line == "[&R]")
             curr_tree->set_rooted_tree(TRUE);
        else
            curr_tree->set_rooted_tree(FALSE);
           
        

        if (curr_tree->rooted_tree()==TRUE) std::cout<<"ROOTED TREE\n";
        else std::cout<<"UNROOTED TREE\n";

        
        nest_level=0;
        //****************ROOTED*****************
        //Handles rooted trees
        if (curr_tree->rooted_tree()==TRUE)
        { 	
            treein>>rel;
            //cout<<"First tree char: "<<rel<<endl;
            curbranch=make_subtree(treein);
            

            //std::cout<<curbranch->get_name()<<" Length: "<<curbranch->expect_subs_site()<<endl;

            
            
        } //End Rooted Section
        
        
        //****************UNROOTED*****************
        //Handles unrooted trees by finding the three siblings at the base and joining them
        else
        {
            treein>>rel;
            
            curbranch=make_base_unrooted(treein);
            
            
        } //End Unrooted Section
        
        
        //curr_tree->set_root((*curr_tree)[0]);
        curr_tree->set_root(curbranch);
        return(curr_tree);
    }
} //End Read_PAUP_Tree::create_tree_from_file


Branch* Read_PAUP_Tree_Ex::get_tip(ifstream &treein)
{
    int curid;
    double blen;
    Branch *new_tip;
    
    curid=rel-48;
    
    //Handles Taxa IDs greater than 9
    treein>>rel1;
    
    while(rel1>=48 && rel1<=57)
    {
        curid=curid*10+(rel1-48);
        treein>>rel1;
    }
    
    new_tip=(*curr_tree)[curid-1];
    new_tip->set_tip(TRUE);
    treein>>blen;
    new_tip->set_expect_subs_site(blen);
    
    return(new_tip);
    
}  //End Read_PAUP_Tree::get_tip



Branch* Read_PAUP_Tree_Ex::get_interior(ifstream &treein)
{
    double blen;
    Branch *curbranch;
    
    treein>>rel1;
    treein>>rel1;
    
    curbranch=initialize_branch();
    curbranch->set_tip(FALSE);
    treein>>blen;
    curbranch->set_expect_subs_site(blen);
    
    return(curbranch);
}  //End Read_PAUP_Tree::get_interior




Branch* Read_PAUP_Tree_Ex::make_subtree(ifstream &treein)
{
    int i;
    double brn;
    Branch *sibling1, *sibling2, *parent;
    
    if (rel=='(')
        nest_level++;
    
    treein>>rel;
    
    if (rel=='(')
        sibling1=make_subtree(treein);                                     //Recursion
    else
        sibling1=get_tip(treein);
    
    treein>>rel;
    treein>>rel;
    
    if (rel=='(')
        sibling2=make_subtree(treein);                                 //Recursion
    else
        sibling2=get_tip(treein);
    
#if defined (DEBUG)
    cout<<" Sibling 1: "<<sibling1->get_name()<<" Sibling 2: "<<sibling2->get_name()<<" Nest level: "<<nest_level<<" \n"<<flush;
#endif
    
    curr_tree->set_as_siblings(sibling1, sibling2);
    
    if (nest_level>1)
        parent=get_interior(treein);
    else
    {
        parent=initialize_branch();
        parent->set_tip(FALSE);
        if (curr_tree->rooted_tree() == FALSE) {
            parent->set_brnlen(0.0);
            parent->set_expect_subs_site(0.0);
        }
        else {
            treein>>rel;
            treein>>rel;
            treein>>brn;
            parent->set_expect_subs_site(brn);
        }
        parent->set_parent(0);
        parent->set_name("Root");
        
    }
    curr_tree->set_as_parent_child(parent, sibling1, 0);
    curr_tree->set_as_parent_child(parent, sibling2, 1);
    
    nest_level--;
    return (parent);
}  //End Read_PAUP_Tree::make_subtree

Branch* Read_PAUP_Tree_Ex::make_base_unrooted(ifstream &treein)
{
    int i;
    Branch *sibling1, *sibling2, *sibling3, *parent;
    
    if (rel=='(')
        nest_level++;
    
    treein>>rel;
    
    if (rel=='(')
        sibling1=make_subtree(treein);                                     //Recursion
    else
        sibling1=get_tip(treein);
    
    treein>>rel;
    treein>>rel;
    
    if (rel=='(')
        sibling2=make_subtree(treein);                                 //Recursion
    else
        sibling2=get_tip(treein);
    
    treein>>rel;
    treein>>rel;
    
    if (rel=='(')
        sibling3=make_subtree(treein);                                 //Recursion
    else
        sibling3=get_tip(treein);
    
#if defined (DEBUG)
    cout<<"Unrooted base: Sibling 1: "<<sibling1->get_name()<<" Sibling 2: "<<sibling2->get_name()
    <<" Sibling 3: "<<sibling3->get_name()<<endl;
#endif
    
    curr_tree->set_as_siblings(sibling1, sibling2);
    
    parent=initialize_branch();
    parent->set_tip(FALSE);
    parent->set_brnlen(0.0);
    parent->set_expect_subs_site(0.0);
    
    curr_tree->set_as_parent_child(parent, sibling1, 0);
    curr_tree->set_as_parent_child(parent, sibling2, 1);
    
    sibling1=sibling3;
    sibling2=parent;
    
    curr_tree->set_as_siblings(sibling1, sibling2);
    
    parent=initialize_branch();
    parent->set_tip(FALSE);
    parent->set_brnlen(0.0);
    parent->set_expect_subs_site(0.0);
    parent->set_name("Root");
    
    curr_tree->set_as_parent_child(parent, sibling1, 1);
    curr_tree->set_as_parent_child(parent, sibling2, 0);
    
#if defined (DEBUG)
    cout<<" Parent: "<<parent->get_name()<<")";
#endif
    
    nest_level--;
    return (parent);
}  //End Read_PAUP_Tree_Ex::make_base_unrooted


void Write_Tree_Arb_Model::write_tree(string output_file, string prog_name, Phylo_Matrix *cmatrix, Phylo_Matrix *rmatrix, Tree *ctree, Exchange *cexchange)
{
    
    curr_exchange=cexchange;
    tree_object=ctree;
    the_matrix=cmatrix;
    root_matrix=rmatrix;
    
    treefile = new ofstream(output_file.c_str());
    
    //strcpy(filename,output_file.c_str());
    
    strcpy(program_name, prog_name.c_str());
    strcpy(matrix_file, "\0");
    write_descript_string();
    write_tree_string();
    delete treefile;
}


void Write_Tree_Arb_Model::write_descript_string()
{
    int i, j, pos;
    char temp_num[42];
    
    if (root_matrix ==0)
        dlines = the_matrix->get_num_params()+3;
    else
        dlines = the_matrix->get_num_params() + root_matrix->get_num_params()+5;
    
    descript=new char*[dlines];
    for(i=0; i<dlines; i++)
        descript[i]=new char[120];
    
    cout<<"Writing tree for model: "<<the_matrix->get_model_name()<<endl;
    strcpy(descript[0], program_name);
    strcpy(descript[1], the_matrix->get_model_name().c_str());
    
    pos=2;
    for(i=0; i<the_matrix->get_num_params(); i++) {
        cout<<"Storing parameter "<<the_matrix->get_param_name(i)<<endl;
        strcpy(descript[pos], the_matrix->get_param_name(i).c_str());
        strcat(descript[pos], ": ");
        double_to_string (temp_num, 39, 4,the_matrix->get_param(i));
        strcat(descript[pos], temp_num);
        pos++;
    }
    
    if (root_matrix !=0) {
        strcpy(descript[pos], root_matrix->get_model_name().c_str());
        pos++;
        
        for(i=0; i<root_matrix->get_num_params(); i++) {
            strcpy(descript[pos], root_matrix->get_param_name(i).c_str());
            strcat(descript[pos], ": ");
            double_to_string (temp_num, 39, 4,root_matrix->get_param(i));
            strcat(descript[pos], temp_num);
            pos++;
        }
        
        strcpy(descript[pos],"");
        for(i=0; i<root_matrix->get_num_states(); i++) {
            if (root_matrix->get_nth_state(i)->is_root_state()==TRUE) {
                for (j=0; j<root_matrix->get_num_states(); j++) {
                    if (root_matrix->get_tp_matrix_num(tree_object->find_root()->get_brn_num())->get_transprob(i, j, 0) !=0) {
                        strcat(descript[pos], root_matrix->get_nth_state(i)->get_state_name().c_str());
                        strcat(descript[pos], "->");
                        strcat(descript[pos], root_matrix->get_nth_state(j)->get_state_name().c_str());
                        strcat(descript[pos], "=");
                        double_to_string (temp_num, 39, 4,root_matrix->get_tp_matrix_num(tree_object->find_root()->get_brn_num())->get_transprob(i, j, 0));
                        strcat(descript[pos], temp_num);
                        strcat(descript[pos], ", ");
                    }
                }
            }
        }
        pos++;
    }
    
    strcpy(descript[pos], "Final ln likelihood: ");
    double_to_string (temp_num, 39, 4, curr_exchange->get_saved_lnL());
    strcat(descript[pos], temp_num);
        
    
}




Branch* Read_PAUP_Settings_Tree_Ex::make_subtree(ifstream &treein)
{
    int i, model_id;
    char dump;
    double brn;
    Branch *sibling1, *sibling2, *parent;
    
    if (rel=='(')
        nest_level++;
    
    treein>>rel;
    
    if (rel=='(')
        sibling1=make_subtree(treein);                                     //Recursion
    else
        sibling1=get_tip(treein);
    
    treein>>rel;
    treein>>rel;
    
    if (rel=='(')
        sibling2=make_subtree(treein);                                 //Recursion
    else
        sibling2=get_tip(treein);
    
//#if defined (DEBUG)
    cout<<" Sibling 1: "<<sibling1->get_name()<<" Sibling 2: "<<sibling2->get_name()<<" Nest level: "<<nest_level<<" \n"<<flush;
//#endif
    
    curr_tree->set_as_siblings(sibling1, sibling2);
    
    if (nest_level>1)
        parent=get_interior(treein);
    else
    {
        parent=initialize_branch();
        parent->set_tip(FALSE);
        if (curr_tree->rooted_tree() == FALSE) {
            parent->set_brnlen(0.0);
            parent->set_expect_subs_site(0.0);
        }
        else {
            treein>>rel;
            treein>>rel;
            treein>>brn>>dump>>model_id;
            parent->set_expect_subs_site(brn);
            parent->set_param_set(model_id-1);

        }
        parent->set_parent(0);
        parent->set_name("Root");
        
    }
    curr_tree->set_as_parent_child(parent, sibling1, 0);
    curr_tree->set_as_parent_child(parent, sibling2, 1);
    
    nest_level--;
    return (parent);
}  //End Read_PAUP_Tree::make_subtree

Branch* Read_PAUP_Settings_Tree_Ex::get_tip(ifstream &treein)
{
    int i, curid, model_id;
    char dump;
    double blen;
    Branch *new_tip;
    
    curid=rel-48;
    
    //Handles Taxa IDs greater than 9
    treein>>rel1;
    
    while(rel1>=48 && rel1<=57)
    {
        curid=curid*10+(rel1-48);
        treein>>rel1;
    }
    
    new_tip=(*curr_tree)[curid-1];
    new_tip->set_tip(TRUE);
    treein>>blen>>dump>>model_id;
    
    cout<<"Creating tip: "<<curid<<" len "<<blen<<" model: "<<model_id<<endl;
    
    new_tip->set_expect_subs_site(blen);
    new_tip->set_param_set(model_id-1);
    
    return(new_tip);
    
}  //End Read_PAUP_Tree::get_tip



Branch* Read_PAUP_Settings_Tree_Ex::get_interior(ifstream &treein)
{
    int i, model_id;
    char dump;
    double blen;
    Branch *curbranch;
    
    treein>>rel1;
    treein>>rel1;
    
    curbranch=initialize_branch();
    curbranch->set_tip(FALSE);
    treein>>blen>>dump>>model_id;
    curbranch->set_expect_subs_site(blen);
    curbranch->set_param_set(model_id-1);
    
    cout<<"Creating interior: len "<<blen<<" model: "<<model_id<<endl;
    
    
    return(curbranch);
}  //End Read_PAUP_Tree::get_interior

