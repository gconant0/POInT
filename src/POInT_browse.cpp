#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

#include "cgicc/Cgicc.h"
#include "cgicc/HTTPHTMLHeader.h"
#include "cgicc/HTMLClasses.h"
#include <cgicc/HTTPContentHeader.h>
//#include "base64.hpp"
#include "Base64.h"

using namespace std;
using namespace cgicc;




#define MAX_NAME_LEN 100
#define BUFFER_SIZE 2048
#define JUMP_FRAMES 100

enum MSG_TYPE {MESSAGE, DATABLOCK, SIZE, PING, ENDCOMM, PHYLOTREE, BATCHBLOCK, MODELDIAG, TREEDIAG};
enum IMAGE_SIZE {SMALL, MEDIUM, LARGE};
//Receive message types from POInT
MSG_TYPE get_message_type (char* buffer);

void get_image_data(string gene_name, string window_size, string my_socketname, char imageformat, bool use_position, bool &error, string &error_msg, std::stringstream *plot_ss, int &pillar, int &start_pillar, int &end_pillar, int &last_pillar, IMAGE_SIZE mysize, int *map_coords, int tree_pillar, string &treestring, string tree_type_string, bool is_random);

int get_pdf_image(std::stringstream& image_instream, std::stringstream& image_outstream);
int clean_image(std::stringstream& image_instream, std::stringstream& image_outstream);

bool ping_image_server(string my_socketname, string &error_msg, int &num_genomes, int &num_pillars, int &num_genes, string &version);
void get_modeltree_data(string my_socketname, std::stringstream *model_ss, bool is_model);



bool does_have_events(std::map<std::string, std::string> socket_hash, std::map<std::string, bool> &valid_map, int &num_events,
                      std::map<std::string, int> &genome_cnt, std::map<std::string, int> &pillar_cnt, std::map<std::string, int> &gene_cnt, string &version)
{
    int genomes, pillars, genes;
    string error_msg;
 
    
    num_events=0;
    for (std::map<string,string>::iterator it=socket_hash.begin(); it!=socket_hash.end(); ++it) {
        //cout<<"Testing "<<it->first<<": "<<it->second<<std::endl;
        if (ping_image_server(it->second, error_msg,  genomes, pillars, genes, version)==true) {
            //if (!(fin.fail())) {
             //cout<<" Success"<<endl;
            valid_map[it->first]=true;
            genome_cnt[it->first]=genomes;
            pillar_cnt[it->first]=pillars;
            gene_cnt[it->first]=genes;
            num_events++;
        }
        
    }
    if (num_events==0)
        return(false);
    
    else return(true);
}

// Print the form for this CGI
void printForm(const Cgicc& cgi, std::map<std::string, std::string> socket_hash, std::map<std::string, bool> valid_map, int num_events, std::string server_path, bool is_random, int rand_pillar)
{
    int frame_size;
    bool have_event=false, have_size=false;
    string error_msg, default_gene, default_event, default_size, last_pillar, href_path, space_string;
    ifstream fin;
    IMAGE_SIZE my_image_size=MEDIUM;
    
    last_pillar="-1";
    
    default_gene="";
    
    cout << "<form name=\"drawform\" id=\"drawform\" method=\"post\" action=\""
    << cgi.getEnvironment().getScriptName()
        << "\" enctype=\"multipart/form-data\">" << endl;
   
    const_form_iterator name = cgi.getElement("genename");
    if(name != cgi.getElements().end()) {
        //cout << "Your name: " << **name << endl;
        default_gene=**name;
    }
    
    const_form_iterator event  = cgi.getElement("event");
    if (event != cgi.getElements().end()) {
        default_event=**event;
        have_event=true;
        //cout<<"Using socket "<<socket_name<<endl;
    }
    
    const_form_iterator lastpill = cgi.getElement("pillarid");
    if (lastpill != cgi.getElements().end()) {
        last_pillar=**lastpill;
    }
    
    const_form_iterator imgsize = cgi.getElement("imagesize");
    if (imgsize != cgi.getElements().end()) {
        if (**imgsize == "Large") my_image_size=LARGE;
        else {
            if (**imgsize == "Medium") my_image_size=MEDIUM;
            else {
                if (**imgsize == "Small") my_image_size=SMALL;
            }
        }
    }
    
    const_form_iterator size  = cgi.getElement("wind_size");
    if (size != cgi.getElements().end()) {
        default_size=**size;
        have_size=true;
        
        try {
            frame_size=std::stoi(default_size,nullptr);
        }
        catch (std::invalid_argument const &e) {have_size=false;}
        catch (std::out_of_range const &e) {have_size=false;}
    }
    
    switch (my_image_size) {
        case LARGE:
            cout << "<table width=1500>" << endl;
            space_string="&nbsp;&nbsp;";
            break;
        case MEDIUM:
            cout << "<table width=1340>" << endl;
            space_string="&nbsp;";
            break;
        case SMALL:
            cout << "<table width=950>" << endl;
            space_string="&nbsp;";
            break;
    }
    
    //cout << "<table>" << endl;
    href_path="http://"+server_path+"/help.html#genesearch_ref";
    //cout << "<tr><td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Gene name:</a></b></td>"
        cout<< "<td class=\"form\">"<<"<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Gene name:</a></b>"<<space_string
    << "<input type=\"text\" name=\"genename\" accept=\"text/plain\" value=\""<<default_gene<<"\" />"<<space_string;
    
    href_path="http://"+server_path+"/help.html#randompillar_ref";
    cout<<" or <b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Random</a></b>: "<<"<input type=\"checkbox\" name=\"random_pillar\" value=\"Random\">"
            << "</td>";
   
    

    if (num_events>0) {
        href_path="http://"+server_path+"/help.html#event_ref";
        //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Event:</a></b></td>"<<"
        cout<<"<td class=\"form\">"<<"<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Event:</a></b>"<<space_string
            <<"<select id=\"event\" name = \"event\">";
        for (std::map<string,string>::iterator it=socket_hash.begin(); it!=socket_hash.end(); ++it) {
            if (valid_map[it->first]==true) {
                if (have_event==true) {
                    if (it->first==default_event)
                        cout<<"<option selected=\"selected\" value=\""<<it->first<<"\">"<<it->first<<"</option> ";
                    else
                        cout<<"<option  value=\""<<it->first<<"\">"<<it->first<<"</option> ";
                    
                }
                else {
                    cout<<"<option value=\""<<it->first<<"\">"<<it->first<<"</option> ";
                }
            }
        }
        cout<< "</td>";
        

    }
    
    href_path="http://"+server_path+"/help.html#window_size_ref";
        //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Window size:</a></b></td>"
        cout<< "<td class=\"form\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Window size:</a></b>"<<space_string
    //<< "<input type=\"text\" name=\"wind_size\" accept=\"text/plain\" />"
            << "<select id=\"wind_size\" name = \"wind_size\">";
    for( int i=10; i<26; i++) {
        if (have_size==true) {
            if (frame_size==i)
                cout<<"<option selected=\"selected\"  value=\""<<i<<"\">"<<i<<"</option> ";
            else
                cout<<"<option value=\""<<i<<"\">"<<i<<"</option> ";
        }
        else  cout<<"<option value=\""<<i<<"\">"<<i<<"</option> ";
    }
    cout<< "</td>";
    
    if (my_image_size==SMALL) cout<<"<tr>";
    href_path="http://"+server_path+"/help.html#viz_size_ref";
    //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Image Size:</a></b></td>"
    cout<<"<td class=\"form\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Image Size:</a></b>";
    cout<<" <input type=\"radio\" name=\"imagesize\" value=\"Small\" ";
    if (my_image_size == SMALL) cout<<" checked=\"checked\" />Small";
    else cout<<"/>Small";
    cout<<" <input type=\"radio\" name=\"imagesize\" value=\"Medium\" ";
    if (my_image_size == MEDIUM) cout<<" checked=\"checked\" />Medium";
    else cout<<"/>Medium";
    cout<<" <input type=\"radio\" name=\"imagesize\" value=\"Large\" ";
    if (my_image_size == LARGE) cout<<" checked=\"checked\" />Large";
    else cout<<"/>Large";
    
    //cout<<"<br />"<<endl;
    
    cout<<"</td>" <<endl;
    
    
    href_path="http://"+server_path+"/help.html#treeformat_ref";
    //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Tree Format:</a></b></td>"
    cout<<"<td class=\"form\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Tree Format:</a></b>";
    cout<<" <input type=\"radio\" name=\"treeformat\" value=\"PAUP\" checked=\"checked\">NEXUS";
    cout<<" <input type=\"radio\" name=\"treeformat\" value=\"Newick\" >Phylip</td>";
    if (my_image_size==SMALL) cout<<"</tr>";
    //cout<<"<td class=\"form\"><input type=\"hidden\" id=\"pillarid\" name=\"pillarid\" value=\"-1\"></td>"<<endl;
    if (is_random ==false)
        cout<<"<td class=\"form\"><input type=\"hidden\" id=\"pillarid\" name=\"pillarid\" value=\""<<last_pillar<<"\"></td>"<<endl;
    else
        cout<<"<td class=\"form\"><input type=\"hidden\" id=\"pillarid\" name=\"pillarid\" value=\""<<rand_pillar<<"\"></td>"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"treepillar\" name=\"treepillar\" value=\"-1\"></td>"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"cdspillar\" name=\"cdspillar\" value=\"-1\"></td>"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"imageformat\" name=\"imageformat\" value=\"Y\"></td>"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"drawmodel\" name=\"drawmodel\" value=\"N\"></td>"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"drawtreediag\" name=\"drawtreediag\" value=\"N\"></td>"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"printstats\" name=\"printstats\" value=\"N\"></td>"<<endl;
    cout<<"</tr>" << endl;
  
    
    cout << "</table>" << endl;
    
    cout << "<div class=\"center\"><p>"
        << "<input type=\"submit\" name=\"drawbutton\" id=\"drawbutton\" value=\"Draw\" onclick=\"reset_tree_pillar()\" />"<< "</p></div></form>" << endl;
}


int main(int argc, char **argv)
{
    int i, frame_size=12, pos, pillar=0, start_pillar=0, end_pillar, last_pillar, new_pillar, my_pillar=-1, map_coords[7], cent_pillar, endx,
    tree_pillar=-1, num_events, cds_pillar, x1, x2, button_size, total, arrow_x, arrow_y, arrow_xs[6], arrow_ys[6], total_genomes, total_pillars, total_genes, line_x_start, line_y_start, line_len, line_height;
    double jframe_size, startx;
    char myimageformat='Y', printchar, draw_model='N', draw_treediag='N', print_stats='N';
    bool size_valid=true, load_data=false, load_error, have_name=false, have_events, have_pillar=false, get_model=false, is_random=false;
    string error_msg, gene_name, window_size, image_string, image_string2, encoded_string2, encoded_string, default_gene, tree_pillar_string, link_string2, image_string3, link_string3, encoded_string3, blank_model_string, link_string4, encoded_string4, image_string4, stat_image,
    link_string, new_name, socket_name, polyploidy_files, polyploidy_name, polyploidy_socket, dump, help_string, cite_string, server_path, image_path, href_path, modelstring, treediag_string, random_string, stats_string, version,
        new_pillar_string, chosen_pillar, tree_string, tree_type_string, socket_id, cds_file, cds_pillar_string, cds_num, line, formatstring, click_string, tree_pict_string;
    std::stringstream plot_ss, plot_ss_cleaned, plot_ss_pdf, model_ss, model_ss_cleaned,
        treediag_ss, treediag_ss_cleaned;
    std::ifstream pp_in, cds_in;
    std::ofstream fout;
    std::map<std::string, bool> valid_map;
    std::map<std::string, std::string> socket_hash, tree_pict_hash;
    std::map<std::string, int> genome_cnts, pillar_cnts, gene_cnts;
    IMAGE_SIZE my_image_size=MEDIUM;
    
    default_gene="";
    
    window_size="12";
    
    tree_type_string="NEXUS";
    tree_string="";
    cds_file="";
    formatstring="YES";
    
    if (const char* server_p=std::getenv("HTTP_HOST"))
        server_path=server_p;
    else server_path="UNKNOWN";
    
    if (argc <3) {
        //if (argc == 2) {
            //polyploidy_files=argv[1];
     
        polyploidy_files="/data/POInT/POInT_browser_events_cleaned_prefix_commonname.txt";
        pp_in.open(polyploidy_files.c_str());
            
        if (pp_in.fail()) {
            socket_hash["GrassRho"]="/data/POInT/POInTGrass";
        }
        else {
            while (! pp_in.eof()) {
                polyploidy_name="";
                pp_in>>polyploidy_name>>polyploidy_socket;
                std::getline(pp_in, line);
                
                std::istringstream ssp(line);
#if 0
                while (! ssp.eof())
                    ssp>>tree_pict_string;
                
                tree_pict_hash[polyploidy_name]=tree_pict_string;
                //cout<<"Got event "<<polyploidy_name<<" at "<<polyploidy_socket<<endl;
#endif
                if(polyploidy_name!= "") {
                    socket_hash[polyploidy_name]=polyploidy_socket;
                }
            }
            pp_in.close();
        }
            
        
    
        
       try {
           Cgicc cgi;
           
           //cout << HTTPHTMLHeader() << endl;
           
           // Set up the HTML document
           //cout << html() << head(title("POInT browse")) << endl;
           //cout << body() << endl;
           
           form_iterator imgsize = cgi.getElement("imagesize");
           if (imgsize != cgi.getElements().end()) {
               if (**imgsize == "Large")
                   my_image_size=LARGE;
               else {
                   if (**imgsize == "Medium")
                       my_image_size=MEDIUM;
                   else {
                       if (**imgsize == "Small")
                           my_image_size=SMALL;
                   }
               }
           }
           
           form_iterator name = cgi.getElement("genename");
           if(name != cgi.getElements().end()) {
               //cout << "Your name: " << **name << endl;
               gene_name=**name;
               
               new_name="";
               //pos=0;
               for (string::iterator it = gene_name.begin(); it < gene_name.end(); it++){
                   if ((((int)*it)>=32) &&(((int)*it)<=126)) {
                       new_name+=*it;
                       //pos++;
                   }
               }
               
               gene_name=new_name;
               
               if (gene_name.length() >MAX_NAME_LEN) {
                   gene_name=gene_name.substr(0,MAX_NAME_LEN);
               }
               have_name=true;
              // cout << "Your name: |" << **name << "| gene name = |"<<gene_name<<"| len = "<<gene_name.length()<<"|<br>"<<endl;
               
               if (gene_name.length() == 0) {
                   have_name=false;
               }
           }
           
           form_iterator size  = cgi.getElement("wind_size");
           if (size != cgi.getElements().end()) {
               window_size=**size;
               
              // cout<<"Current window size value "<<window_size<<"<br>"<<endl;
               
               try {
                   frame_size=std::stoi(window_size,nullptr);
               }
               catch (std::invalid_argument const &e) {size_valid=false;}
               catch (std::out_of_range const &e) {size_valid=false;}
           }
           
           form_iterator hidpillar  = cgi.getElement("pillarid");
           if (hidpillar != cgi.getElements().end()) {
               chosen_pillar=**hidpillar;
               
               my_pillar=std::stoi(chosen_pillar, nullptr);
               if (my_pillar >= 0) {
                   have_pillar=true;
                   //cout<<"Plotting from pillar "<<my_pillar<<endl;
               }
           }
           
           form_iterator randpillar  = cgi.getElement("random_pillar");
           if (randpillar != cgi.getElements().end()) {
               if (**randpillar=="Random")
                   is_random=true;
                   //cout<<"Plotting from pillar "<<my_pillar<<endl;
               
           }
           
           form_iterator treepill = cgi.getElement("treepillar");
           if (treepill != cgi.getElements().end()) {
               tree_pillar_string=**treepill;
               tree_pillar=std::stoi(tree_pillar_string, nullptr);
               
               //cout<<"REquesing geen tree for "<<tree_pillar<<std::endl<<std::flush;
           }

           form_iterator imageformat = cgi.getElement("imageformat");
           if (imageformat != cgi.getElements().end()) {
               formatstring=**imageformat;
               if (formatstring[0] == 'N') myimageformat='N';
           }

           form_iterator drawmodel = cgi.getElement("drawmodel");
           if (drawmodel != cgi.getElements().end()) {
               modelstring=**drawmodel;
               draw_model=modelstring[0];
           }
           
           form_iterator drawtreediag = cgi.getElement("drawtreediag");
           if (drawtreediag != cgi.getElements().end()) {
               treediag_string=**drawtreediag;
               draw_treediag=treediag_string[0];
           }
           
           form_iterator printstats = cgi.getElement("printstats");
           if (printstats != cgi.getElements().end()) {
               stats_string=**printstats;
               print_stats=stats_string[0];
           }
           
           form_iterator treetype = cgi.getElement("treeformat");
           if (treetype != cgi.getElements().end()) {
               tree_type_string=**treetype;
           }
           
           form_iterator event  = cgi.getElement("event");
           if (event != cgi.getElements().end()) {
               socket_id=**event;
               socket_name=socket_hash[**event];
               //cout<<"Using socket "<<socket_name<<endl;
           }
           
           
           form_iterator cdspill = cgi.getElement("cdspillar");
           if (cdspill != cgi.getElements().end()) {
               cds_pillar_string=**cdspill;
               cds_pillar=std::stoi(cds_pillar_string, nullptr);
               if (cds_pillar != -1) {
                   //if ((cds_pillar<0) || (cds_pillar>=last_pillar)) {
                   //    load_error=true;
                   //}
                   //else {
                       stringstream ss5;
                       ss5 << cds_pillar;
                       cds_num = ss5.str();
                       
                       cds_file = "/var/www/html/Downloads/" + socket_id + "/Pillars/Pillar" + cds_num + "_CDS.fas";
                       
                       cds_in.open(cds_file.c_str());
                       
                       if (cds_in.fail()) {
                           load_error=true;
                           error_msg="Requested CDS file does not exist";
                       }
                        else load_data=true;
                       cds_in.close();
                       cds_in.clear();
                       
                   //}
               }
               
               //cout<<"REquesing geen tree for "<<tree_pillar<<std::endl<<std::flush;
           }
           
           load_error=true;
           have_events=does_have_events(socket_hash, valid_map, num_events, genome_cnts, pillar_cnts, gene_cnts, version);
           total_genomes=0;
           total_pillars=0;
           total_genes=0;
           for (std::map<string,int>::iterator it=genome_cnts.begin(); it!=genome_cnts.end(); ++it) {
               total_genomes +=it->second;
           }
           for (std::map<string,int>::iterator it=pillar_cnts.begin(); it!=pillar_cnts.end(); ++it) {
               total_pillars +=it->second;
           }
           for (std::map<string,int>::iterator it=gene_cnts.begin(); it!=gene_cnts.end(); ++it) {
               total_genes +=it->second;
           }
           
           
           //cout<<"MAp for "<<socket_name<<" is "<<valid_map[socket_id]<<" hn hp he le"<<have_name<<have_pillar<<have_events<<load_error<<endl;
           
           if (have_events == true) {
               load_error=false;
               
               if (size_valid == false) {
                   error_msg= "Cannot convert ";
                   error_msg += **size;
                   error_msg +=" to a valid window size<br>";
                   //cout<<"ERROR: Cannot convert "<< **size <<" to a valid window size<br>"<<endl;
                   load_error=true;
                   load_data=false;
               }
               else {
                 if ((load_error == false) && (valid_map[socket_id]==true) && (cds_file == "") && ((have_name == true) || (have_pillar==true) || (is_random==true))) {
                       window_size=**size;
                       plot_ss.str(std::string());
                     
                     if (draw_model=='Y') {
                         get_modeltree_data(socket_name, &model_ss, true);
                         clean_image(model_ss, model_ss_cleaned);
                       
                     }
                     
                     if (draw_treediag=='Y') {
                         get_modeltree_data(socket_name, &treediag_ss, false);
                         clean_image(treediag_ss, treediag_ss_cleaned);
                     }
                       
                       //cout<<"Sending data to plotter: "<<gene_name<<" size: |"<<window_size<<"|"<<" socket: "<<socket_name<<" tree pillar "<<tree_pillar<<endl<<flush;
                       tree_string="";
                       get_image_data(gene_name, window_size, socket_name, myimageformat, have_pillar, load_error, error_msg, &plot_ss, my_pillar, start_pillar, end_pillar, last_pillar, my_image_size, map_coords, tree_pillar, tree_string, tree_type_string, is_random);
                       pillar=my_pillar;
                       load_data=true;
                       default_gene=gene_name;
                   }
               }
           }
           
           //cout << HTTPHTMLHeader() << endl;
           //cout << html() << head(title("POInT browse")) << endl;
           //cout << body() << endl;
           //cout<<"Load choice: le ld: "<<load_error<<load_data<<" TS: "<<tree_string<<endl<<flush;
           //cout << body() << html();
           
           if ((cds_file != "") && (load_error == false)  && (load_data==true) ) {
               cout<<"Content-Type:application/x-download\n";
               cout<<"Content-Disposition:attachment;filename=Pillar"<<cds_pillar<<"_CDS.fas\n\n";
               
               cds_in.open(cds_file.c_str());
               
               while (! cds_in.eof()) {
                   std::getline(cds_in, line);
                   std::cout<<line<<std::endl;
               }
               cds_in.close();
               //cout<<cds_file<<std::endl;
               
           }
           else {
               if ((tree_string != "") && (load_error == false)  && (load_data==true) ) {
                   //cout << HTTPHTMLHeader() << endl;
                   
                   // Set up the HTML document
                  
                   
                   //cout<<"Dummy tree download\n";
                   cout<<"Content-Type:application/x-download\n";
                   cout<<"Content-Disposition:attachment;filename=Pillar"<<tree_pillar<<"_genetree.tre\n\n";
                   cout<<tree_string<<std::endl;
                   
               }
               else {
                   if (formatstring[0] == 'N') {
                       cout<<"Content-Type:application/x-download\n";
                       //cout<<"Content-Disposition:attachment;filename=POInTSnap.ps\n\n";
                       cout<<"Content-Disposition:attachment;filename=POInTSnap.pdf\n\n";
                       
                      
                      
                       get_pdf_image(plot_ss, plot_ss_pdf);
                    //plot_ss.seekg(0, std::ios::beg);
                       
                       
                       //std::cout<<"Size of ps stream: "<<plot_ss.str().length()<<std::endl;
                       //std::cout<<"Size of pdf stream: "<<plot_ss_pdf.str().length()<<std::endl;
                       

#if 1
                       while(!(plot_ss_pdf.eof())) {
                           printchar=(char)plot_ss_pdf.get();
                           if (!(plot_ss_pdf.eof()))
                               cout<<printchar;
                       }
#else
                       while(!(plot_ss.eof())) {
                           printchar=(char)plot_ss.get();
                           if (!(plot_ss.eof()) && (printchar>0))
                               cout<<printchar;
                       }
                       cout<<std::endl;
#endif
                       
                       //cout<<plot_ss.str()<<std::endl;
                   }
                   else {
                       // Send HTTP header
                       cout << HTTPHTMLHeader() << endl;
                       
                       // Set up the HTML document
                       cout << html() << head(title("POInT browse")) << endl;
                       cout << body() << endl;
                       
                       // if (cds_file != "") {
                       //     cout<<"Should have downloaded "<<cds_file<<" | "<<cds_pillar_string<<" | "<<cds_pillar<<endl;
                       // }
                       //if (cds_file == "")
                       //     cds_file="/var/www/html/Downloads/GrassRho/Pillar0_CDS.fas";
                       // cds_in.open(cds_file.c_str());
                       
                       //while (! cds_in.eof()) {
                       //  std::getline(cds_in, line);
                       //  std::cout<<line<<std::endl;
                       //}
                       //cds_in.close();
                       //  }
                       //  if (tree_string != "") {
                       //      cout<<"Content-Type:application/x-download\n";
                       //      cout<<"Content-Disposition:attachment;filename=Pillar"<<tree_pillar<<"_genetree.tre\n\n";
                       //      cout<<tree_string<<std::endl;
                       //  }
                       
                       
                       
                       switch (my_image_size) {
                           case LARGE:
                               image_path ="http://"+server_path + "/POInT_browse_header.png";
                               
                               break;
                           case MEDIUM:
                               image_path ="http://"+server_path + "/POInT_browse_header_M.png";
                               break;
                           case SMALL:
                               image_path ="http://"+server_path + "/POInT_browse_header_S.png";
                               break;
                       }
                       cout << img().set("SRC", image_path)
                           .set("ALT", "Bad Header image") << endl;
                       
                       printForm(cgi, socket_hash, valid_map, num_events, server_path, is_random, pillar);
                       
                       //cout<<"Host name is "<<server_path<<endl;
                       //cout<<"LD is "<<load_data<<" and le is "<<load_error<<endl<<flush;
                       
                       if ((load_data==false) || (load_error==true) ) {
                           //cout<<"LD is "<<load_data<<" and le is "<<load_error<<endl<<flush;
                           switch (my_image_size) {
                               case LARGE:
                                   image_path="http://" +server_path+"/POInT_browse_default_screen.png";
                                   break;
                               case MEDIUM:
                                   image_path="http://" +server_path+"/POInT_browse_default_screen_M.png";
                                   break;
                               case SMALL:
                                   image_path="http://" +server_path+"/POInT_browse_default_screen_S.png";
                                   break;
                           }
                           cout << img().set("SRC", image_path)
                               .set("ALT", "Bad Header image") << endl;
                           cout <<"<br>"<<endl;
                           //if (load_error==true) {
                           //    cout<<"<h3>Messages</h3><br>"<<endl;
                           //    cout<<error_msg<<"<br>"<<endl;
                               
                           //}
                       }
                       else {
                           
                           
                           cout << "<table>" << endl;
                           cout << "<tr><td>";
                           
                           
                           switch (my_image_size) {
                               case LARGE:
                                   arrow_x=78;
                                   arrow_y=907;
                                   arrow_xs[0]=2;
                                   arrow_xs[1]=71;
                                   arrow_xs[2]=11;
                                   arrow_xs[3]=73;
                                   arrow_xs[4]=2;
                                   arrow_xs[5]=71;
                                   
                                   arrow_ys[0]=26;
                                   arrow_ys[1]=71;
                                   arrow_ys[2]=319;
                                   arrow_ys[3]=585;
                                   arrow_ys[4]=108;
                                   arrow_ys[5]=160;
                                   
                                   line_x_start=15;
                                   line_y_start=896;
                                   line_len=1770;
                                   line_height=12;
                                   
                                   break;
                               case MEDIUM:
                                   arrow_x=50;
                                   arrow_y=603;
                                   arrow_xs[0]=1;
                                   arrow_xs[1]=46;
                                   arrow_xs[2]=6;
                                   arrow_xs[3]=46;
                                   arrow_xs[4]=1;
                                   arrow_xs[5]=46;
                                  
                                   
                                   arrow_ys[0]=16;
                                   arrow_ys[1]=54;
                                   arrow_ys[2]=214;
                                   arrow_ys[3]=388;
                                   arrow_ys[4]=62;
                                   arrow_ys[5]=94;
                                   
                                   line_x_start=7;
                                   line_y_start=598;
                                   line_len=1179;
                                   line_height=9;
                                   break;
                               case SMALL:
                                   arrow_x=38;
                                   arrow_y=452;
                                   arrow_xs[0]=2;
                                   arrow_xs[1]=35;
                                   arrow_xs[2]=3;
                                   arrow_xs[3]=38;
                                   arrow_xs[4]=2;
                                   arrow_xs[5]=35;
                                   
                                   arrow_ys[0]=17;
                                   arrow_ys[1]=39;
                                   arrow_ys[2]=158;
                                   arrow_ys[3]=290;
                                   arrow_ys[4]=60;
                                   arrow_ys[5]=85;
                                   
                                   line_x_start=4;
                                   line_y_start=447;
                                   line_len=886;
                                   line_height=6;
                                   break;
                           }
                           cout<<"<map NAME=\"arrowmap2\">";
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<arrow_xs[0]<<", "<<arrow_ys[0]<<", "<<arrow_xs[1]<<", "<<arrow_ys[1]
                           <<"\" HREF=\"#\"  onclick=\"draw_treediag()\">";
                           
                           
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<arrow_xs[4]<<", "<<arrow_ys[4]<<", "<<arrow_xs[5]<<", "<<arrow_ys[5]
                           <<"\" HREF=\"#\"  onclick=\"draw_model()\">";
                           
                           
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<arrow_xs[2]<<", "<<arrow_ys[2]<<", "<<arrow_xs[3]<<", "<<arrow_ys[3]
                           <<"\" HREF=\"#\"  onclick=\"move_left()\">";
                           
                           
                           cout<<"</map>\n";
                           switch (my_image_size) {
                               case LARGE:
                                   image_path="http://"+server_path+"/LeftArrow_L.png";
                                   break;
                               case MEDIUM:
                                   image_path="http://"+server_path+"/LeftArrow_M.png";
                                   break;
                               case SMALL:
                                   image_path="http://"+server_path+"/LeftArrow_S.png";
                                   break;
                           }
                           cout<<"<img usemap=\"#arrowmap2\" src=\""<<image_path<<"\" alt=\"Bad broswer image\"/>\n";
                           
                           cout <<"</td><td>";
                           
                           cout<<"<map NAME=\"genemap\">";
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<map_coords[0]<<", "<<map_coords[3]<<", "<<map_coords[0]+map_coords[1]<<", "<<map_coords[4]
                           <<"\" HREF=\"#\"  onclick=\"move_to_pillar("<<0<<")\" >";
                           
                           
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<map_coords[0]<<", "<<map_coords[5]<<", "<<map_coords[0]+(map_coords[1]/2)<<", "<<map_coords[6]<<"\" HREF=\"#\"  onclick=\"draw_pillar_tree("<<0<<")\" >";
                           
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<map_coords[0]+(map_coords[1]/2)<<", "<<map_coords[5]<<", "<<map_coords[0]+map_coords[1]<<", "<<map_coords[6]<<"\" HREF=\"#\"  onclick=\"download_cds("<<0<<")\" >";
                           cout<<std::endl;
                           
                           startx=(double)line_x_start;
                           
                           jframe_size=(double)line_len/(double)JUMP_FRAMES;
                           for (i=0; i<JUMP_FRAMES; i++) {
                               endx=(int)(startx+jframe_size);
                               cent_pillar=(int)(((i+0.5)/JUMP_FRAMES)*last_pillar);
                               cout<<"<area SHAPE=\"RECT\" COORDS=\""<<(int)startx<<", "<<line_y_start<<", "<<endx<<", "<<line_y_start-line_height<<"\" HREF=\"#\"  onclick=\"move_to_abs_pillar("<<cent_pillar<<")\" >";
                               cout<<std::endl;
                               startx+=jframe_size;
                           }
                           
                           x1=map_coords[0];
                           x2=map_coords[0]+map_coords[1];
                           
                           for(i=1; i<frame_size; i++) {
                               x1 = x1 + map_coords[1]+map_coords[2];
                               if (i%2==1) x1++;
                               
                               x2 = x2 + map_coords[1]+map_coords[2];
                               if (i%2 ==1) x2++;
                               
                               cout<<"<area SHAPE=\"RECT\" COORDS=\""<<x1<<", "<<map_coords[3]<<", "<<x2<<", "<<map_coords[4]<<"\" HREF=\"#\"  onclick=\"move_to_pillar("<<i<<")\" >";
                               
                               cout<<"<area SHAPE=\"RECT\" COORDS=\""<<x1<<", "<<map_coords[5]<<", "<<x1+(map_coords[1]/2)<<", "<<map_coords[6]<<"\" HREF=\"#\"  onclick=\"draw_pillar_tree("<<i<<")\" >";
                               
                               cout<<"<area SHAPE=\"RECT\" COORDS=\""<<x1+(map_coords[1]/2)<<", "<<map_coords[5]<<", "<<x2<<", "<<map_coords[6]<<"\" HREF=\"#\"  onclick=\"download_cds("<<i<<")\" >";
                               cout<<std::endl;
                           }
                           cout<<"</map>\n";
                           
                           
                           
                           clean_image(plot_ss, plot_ss_cleaned);
                           
                           
                           macaron::Base64 encoder2;
                           image_string2=plot_ss_cleaned.str();
                           encoded_string2=encoder2.Encode(image_string2);
                           link_string2="data:image/png;base64,";
                           link_string2+=encoded_string2;
                           //cout<<"Encoded: "<<encoded_string<<endl;
                           //auto encoded = base64::encode(data);
                           
                           
                           
                           //cout<<"<img usemap=\"#genemap\"  src=\""<<link_string<<"\" alt=\"Bad broswer image\"/>\n";
                           cout<<"<img usemap=\"#genemap\" src=\""<<link_string2<<"\" alt=\"Bad broswer image\"/>\n";
                           //cout<<plot_ss.str()<<endl;
                           //cout << html();
                           //cout << img().set("SRC", plot_ss.str())
                           //.set("ALT", "Bad draw request") << endl;
                           
                           cout <<"</td><td>";
                           
                           
                           switch (my_image_size) {
                               case LARGE:
                                   arrow_x=76;
                                   arrow_y=907;
                                   arrow_xs[0]=6;
                                   arrow_xs[1]=67;
                                   arrow_xs[2]=14;
                                   arrow_xs[3]=78;
                                   arrow_xs[4]=8;
                                   arrow_xs[5]=73;
                                   
                                   arrow_ys[0]=28;
                                   arrow_ys[1]=68;
                                   arrow_ys[2]=323;
                                   arrow_ys[3]=586;
                                   arrow_ys[4]=110;
                                   arrow_ys[5]=154;
                                   break;
                               case MEDIUM:
                                   arrow_x=56;
                                   arrow_y=605;
                                   arrow_xs[0]=6;
                                   arrow_xs[1]=49;
                                   arrow_xs[2]=9;
                                   arrow_xs[3]=48;
                                   arrow_xs[4]=2;
                                   arrow_xs[5]=53;
                                   
                                   arrow_ys[0]=19;
                                   arrow_ys[1]=54;
                                   arrow_ys[2]=213;
                                   arrow_ys[3]=393;
                                   arrow_ys[4]=69;
                                   arrow_ys[5]=106;
                                   break;
                               case SMALL:
                                   arrow_x=43;
                                   arrow_y=453;
                                   arrow_xs[0]=0;
                                   arrow_xs[1]=39;
                                   arrow_xs[2]=5;
                                   arrow_xs[3]=38;
                                   arrow_xs[4]=2;
                                   arrow_xs[5]=40;
                                   
                                   arrow_ys[0]=11;
                                   arrow_ys[1]=33;
                                   arrow_ys[2]=158;
                                   arrow_ys[3]=294;
                                   arrow_ys[4]=56;
                                   arrow_ys[5]=83;
                                   break;
                           }
                           cout<<"<map NAME=\"arrowmap\">";
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<arrow_xs[0]<<", "<<arrow_ys[0]<<", "<<arrow_xs[1]<<", "<<arrow_ys[1];
                           href_path="http://"+server_path+"/Key.png";
                           cout<<"\" HREF=\""<<href_path<<"\"  onclick=\"return open_help(this)\">";
                           
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<arrow_xs[4]<<", "<<arrow_ys[4]<<", "<<arrow_xs[5]<<", "<<arrow_ys[5]
                           <<"\" HREF=\"#\"  onclick=\"do_print_stats()\">";
                           
                           cout<<"<area SHAPE=\"RECT\" COORDS=\""<<arrow_xs[2]<<", "<<arrow_ys[2]<<", "<<arrow_xs[3]<<", "<<arrow_ys[3]
                           <<"\" HREF=\"#\"  onclick=\"move_right()\">";
                           
                           
                           
                           
                           cout<<"</map>\n";
                           switch (my_image_size) {
                               case LARGE:
                                   image_path="http://"+server_path+"/RightArrow_L.png";
                                   break;
                               case MEDIUM:
                                   image_path="http://"+server_path+"/RightArrow_M.png";
                                   break;
                               case SMALL:
                                   image_path="http://"+server_path+"/RightArrow_S.png";
                                   break;
                           }
                           cout<<"<img usemap=\"#arrowmap\" src=\""<<image_path<<"\" alt=\"Bad browser image\"/>\n";
                           
                           
                           cout <<"</td></tr>";
                           
                           
                           cout << "</table>" << endl;
                           
                           
                           
                       }
                       
                       
                       switch (my_image_size) {
                           case LARGE:
                               button_size=220;
                               total = 1950;
                               break;
                           case MEDIUM:
                               button_size=180;
                               total = 1290;
                               break;
                           case SMALL:
                               button_size=140;
                               total = 965;
                               break;
                       }
                       
                       help_string="http://"+server_path+"/help.html";
                       cite_string="http://"+server_path+"/cite.html";
                       blank_model_string="http://"+server_path+"/Blankmodel.png";
                       cout << "<table class=\"center\">" << endl;
                       cout << "<tr><td style=\"width:"<<button_size<<"px\"><a href=\""<<help_string<<"\" onclick=\"return open_help(this)\">";
                       image_path="http://"+server_path+"/Help.png";
                       cout << img().set("SRC", image_path)
                           .set("ALT", "Bad image")  <<"</a></td>"<< endl;
                       //.set("ALT", "Bad image") .set("HREF", help_string) .set("onclick", click_string) << endl;
                       //.set("ALT", "Bad image") .set("onclick", "return open_help(this)") <<"</a>"<< endl;
                       cout<<"<td style=\"width:"<<button_size<<"px\"><a href=\"mailto:gconant@ncsu.edu\">";
                       image_path="http://"+server_path+"/Contact.png";
                       cout << img().set("SRC", image_path)
                           .set("ALT", "Bad image") <<"</a>"<<endl;
                       //.set("ALT", "Bad image") .set("HREF", "mailto:gconant@ncsu.edu") << endl;
                       cout<<"</td>";
                       cout<<"<td style=\"width:"<<button_size<<"px\">";
                       image_path="http://"+server_path+"/Cite.png";
                       cout << img().set("SRC", image_path)
                           .set("ALT", "Bad image") .set("onclick", "return open_cite()") << endl;
                       cout<<"</td>";
                       href_path="http://"+server_path+"/Downloads";
                       cout<<"<td style=\"width:"<<button_size<<"px\"><a href=\""<<href_path<<"\">";
                       image_path="http://"+server_path+"/Data.png";
                       cout << img().set("SRC", image_path)
                           .set("ALT", "Bad image") <<"</a>"<<endl;
                       cout<<"</td>";
                       
                       href_path="http://"+server_path+"/cgi-bin/POInT_download";
                       cout << "<td style=\"width:"<<button_size<<"px\"><a href=\""<<href_path<<"\" onclick=\"return open_batch(this)\">";
                       image_path="http://"+server_path+"/Batch.png";
                       cout << img().set("SRC", image_path)
                           .set("ALT", "Bad image")  <<"</a></td>"<< endl;
                       
                       cout<<"<td style=\"float:right;width:"<<total-(5*button_size)<<"px\">";
                       image_path="http://"+server_path+"/Download.png";
                       cout << img().set("SRC", image_path).set("onclick", "return download_ps_image()").set("ALIGN", "right")
                           .set("ALT", "Bad image") <<endl;
                       cout<<"</td></tr></table>"<<endl;
                       
#if 0
                       if (draw_model=='Y') {
                           macaron::Base64 encoder3;
                           image_string3=model_ss_cleaned.str();
                           encoded_string3=encoder3.Encode(image_string3);
                           link_string3="data:image/png;base64,";
                           link_string3+=encoded_string3;
                           //cout<<"Window-target: new_window\n";
                           //cout<<"Content-type: text/html\n\n";
                           cout<<"<img  src=\""<<link_string3<<"\" alt=\"Bad broswer image\"/>\n";
                       }
#endif
                       cout<<"<script type=\"text/javascript\">\n var newpillar = document.getElementById(\"pillarid\")\n var myform = document.getElementsByName(\"DrawForm\")\n";
                       cout<<" var treepillar = document.getElementById(\"treepillar\")\n";
                       cout<<" var cdspillar = document.getElementById(\"cdspillar\")\n";
                       cout<<" var model_draw = document.getElementById(\"drawmodel\")\n";
                       cout<<" var tree_draw = document.getElementById(\"drawtreediag\")\n";
                       cout<<" var print_stats = document.getElementById(\"printstats\")\n";
                       cout<<" var imageformat = document.getElementById(\"imageformat\")\n";
                       cout<<"function move_left() {\n";
                       cout<<"\tnewpillar.setAttribute(\'value\', \'";
                       new_pillar = pillar - (frame_size/2);
                       if ((new_pillar-(frame_size/2)) <0) new_pillar=(frame_size/2);
                       if ((new_pillar-(frame_size/2)+frame_size) >= last_pillar) new_pillar = ((last_pillar-1)-frame_size)+(frame_size/2);
                       
                       stringstream ss3;
                       ss3 << new_pillar;
                       new_pillar_string = ss3.str();
                       
                       
                       cout<<new_pillar_string<<"\');\n";
                       //cout<<"\tgenename.setAttribute(\'value\', \'\');\n";
                       //cout<<"\tdocument.getElementById('DrawForm').submit();\n";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n";
                       //cout<<"\tmyform[0].submit();\n";
                       
                       cout<<"}\nfunction move_right() {\n";
                       cout<<"\tnewpillar.setAttribute(\'value\', \'";
                       new_pillar = pillar + (frame_size/2);
                       if ((new_pillar-(frame_size/2)) <0) new_pillar=(frame_size/2);
                       if ((new_pillar-(frame_size/2)+frame_size) >= last_pillar) new_pillar = ((last_pillar-1)-frame_size)+(frame_size/2);
                       
                       stringstream ss4;
                       ss4 << new_pillar;
                       new_pillar_string = ss4.str();
                       
                       
                       cout<<new_pillar_string<<"\');\n";
                       //cout<<"\tgenename.setAttribute(\'value\', \'\');\n";
                       //cout<<"\tdocument.getElementById('DrawForm').submit();\n";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n";
                       //cout<<"\tmyform[0].submit();\n";
                       
                       
                       cout<<"}\n";
                       cout<<"function move_to_pillar(pos) {\n";
                       cout<<"new_val=pos+"<<start_pillar<<";\n";
                       cout<<"\tnewpillar.setAttribute(\'value\', new_val);";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n";
                       cout<<"}\n";
                       
                       cout<<"function move_to_abs_pillar(pos) {\n";
                       cout<<"new_val=pos;\n";
                       cout<<"\tnewpillar.setAttribute(\'value\', new_val);";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n";
                       cout<<"}\n";
                       
                       cout<<"function draw_pillar_tree(pos) {\n";
                       cout<<"new_tval=pos+"<<start_pillar<<";\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n}\n";
                       
                       cout<<"function download_cds(pos) {\n";
                       cout<<"new_tval=pos+"<<start_pillar<<";\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n}\n";
                       
                       
                       cout<<"function download_ps_image() {\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'N\');";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n}\n";
                       
                       cout<<"function draw_model() {\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');\n";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'Y\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       //cout<<"\tvar modelwin_var= window.open(\""<<blank_model_string<<"\",\'modelwin\', \'left=20,top=20,width=500,height=500,toolbar=0,scrollbars=1,resizable=0\');\n";
                       cout<<"\tvar modelwin_var= window.open(\"\",\'modelwin\', \'left=20,top=20,width=550,height=550,toolbar=0,scrollbars=1,resizable=0\');\n";
                       //cout<<"\treturn false;\n}\n";
                       
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n}\n";
                       
                       cout<<"function draw_treediag() {\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');\n";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'Y\');\n";
                       //cout<<"\tvar treewin_var= window.open(\""<<blank_model_string<<"\",\'treewin\', \'left=20,top=20,width=1000,height=500,toolbar=0,scrollbars=1,resizable=0\');\n";
                       cout<<"\tvar treewin_var= window.open(\" \",\'treewin\', \'left=20,top=20,width=1050,height=550,toolbar=0,scrollbars=1,resizable=0\');\n";
                       //cout<<"\treturn false;\n}\n";
                       
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n}\n";
                       
                       cout<<"function do_print_stats() {\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');\n";
                       cout<<"\tnew_tval=-1\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'Y\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       //cout<<"\tvar modelwin_var= window.open(\""<<blank_model_string<<"\",\'modelwin\', \'left=20,top=20,width=500,height=500,toolbar=0,scrollbars=1,resizable=0\');\n";
                       cout<<"\tvar statwin_var= window.open(\"\",\'statwin\', \'left=20,top=20,width=510,height=350,toolbar=0,scrollbars=1,resizable=0\');\n";
                       //cout<<"\treturn false;\n}\n";
                       
                       cout<<"\tdocument.forms[\"drawform\"].submit();\n}\n";
                       
                       
                       cout<<"function reset_tree_pillar() {\n";
                       cout<<"\timageformat.setAttribute(\'value\', \'Y\');\n";
                       cout<<"\tnew_tval=-1;\n";
                       cout<<"\ttreepillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tnewpillar.setAttribute(\'value\', new_tval);\n";
                       cout<<"\tcdspillar.setAttribute(\'value\', new_tval);";
                       cout<<"\tmodel_draw.setAttribute(\'value\',\'N\');";
                       cout<<"\tprint_stats.setAttribute(\'value\',\'N\');\n";
                       cout<<"\ttree_draw.setAttribute(\'value\',\'N\');\n";
                       cout<<"\t}\n";
                       
                       //cout<<"function open_help(helphttp) {\n";
                       //cout<<"\tif (! window.focus) return true;\n";
                       //cout<<"\tvar href;\n";
                       //cout<<"\tif (typeof(helphttp) == 'string') href=helphttp;\n";
                       //cout<<"\telse href=helphttp.href;\n";
                       //cout<<"\twindow.open(href,\'helpwin\', \'left=20,top=20,width=500,height=500,toolbar=0,scrollbars=1,resizable=0\');\n";
                       //cout<<"\treturn false;\n}\n";
                       cout<<"function open_help(helphttp) {\n";
                       cout<<"\tif (! window.focus) return true;\n";
                       //cout<<"\t var href;\n";
                       //cout<<"\t href=\""<<help_string<<"\";\n";
                       cout<<"var href;\n";
                       cout<<"if (typeof(helphttp) == \'string\') href=helphttp;\n";
                       cout<<"else href=helphttp.href;\n";
                       cout<<"\twindow.open(href,\'helpwin\', \'left=20,top=20,width=1350,height=800,toolbar=0,scrollbars=1,resizable=0\');\n";
                       cout<<"\treturn false;\n}\n";
                       
                       cout<<"function open_batch(batchhttp) {\n";
                       cout<<"\tif (! window.focus) return true;\n";
                       //cout<<"\t var href;\n";
                       //cout<<"\t href=\""<<help_string<<"\";\n";
                       cout<<"var href;\n";
                       cout<<"if (typeof(batchhttp) == \'string\') href=batchhttp;\n";
                       cout<<"else href=batchhttp.href;\n";
                       cout<<"\twindow.open(href,\'batchwin\', \'left=20,top=20,width=1000,height=300,toolbar=0,scrollbars=1,resizable=0\');\n";
                       cout<<"\treturn false;\n}\n";
                       
                       
                       cout<<"function open_cite() {\n";
                       cout<<"\tif (! window.focus) return true;\n";
                       cout<<"\t var href;\n";
                       cout<<"\t href=\""<<cite_string<<"\";\n";
                       cout<<"\twindow.open(href,\'citewin\', \'left=20,top=20,width=1400,height=500,toolbar=1,scrollbars=1,resizable=1\');\n";
                       cout<<"\treturn false;\n}\n";
                       
                       cout<<"function printerror() {\n";
                       cout<<"alert(\""<<error_msg<<"\");}\n";
                       
#if 1
                       if (load_error==true) {
                           cout<<"document.onload = printerror()\n";
                       }
                       if (draw_model=='Y') {
                           macaron::Base64 encoder3;
                           image_string3=model_ss_cleaned.str();
                           encoded_string3=encoder3.Encode(image_string3);
                           link_string3="data:image/png;base64,";
                           link_string3+=encoded_string3;
                           //cout<<"Window-target: new_window\n";
                           //cout<<"Content-type: text/html\n\n";
                           //cout<<"<img  src=\""<<link_string3<<"\" alt=\"Bad broswer image\"/>\n";
                           
                           //cout<<"window.addEventListener(\"load\", show_model());\n";
                           cout<<"document.onload = show_model()\n";
                           
                           cout<<"function show_model() {\n";
                           cout<<"var modelwin_var=window.open(\"\",\'modelwin\');\n";
                           cout<<"modelwin_var.document.write (\'<img src=\""<<link_string3
                            <<"\" alt=\"Bad broswer image\"/>\');\n";
                           cout<<"modelwin_var.document.close();\n";
                           cout<<"}\n";
                       }
                       if (draw_treediag == 'Y') {
                           macaron::Base64 encoder4;
                           image_string4=treediag_ss_cleaned.str();
                           encoded_string4=encoder4.Encode(image_string4);
                           link_string4="data:image/png;base64,";
                           link_string4+=encoded_string4;
                           //cout<<"Window-target: new_window\n";
                           //cout<<"Content-type: text/html\n\n";
                           //cout<<"<img  src=\""<<link_string3<<"\" alt=\"Bad broswer image\"/>\n";
                           
                           //cout<<"window.addEventListener(\"load\", show_model());\n";
                           cout<<"document.onload = show_treediag()\n";
                           
                           cout<<"function show_treediag() {\n";
                           cout<<"var treewin_var=window.open(\"\",\'treewin\');\n";
                           cout<<"treewin_var.document.write (\'<img src=\""<<link_string4
                            <<"\" alt=\"Bad broswer image\"/>\');\n";
                           cout<<"treewin_var.document.close();\n";
                           cout<<"}\n";
                       }
                       if (print_stats == 'Y') {
                           
                           stat_image="http://"+server_path + "/Statistics.png";
                           
                           cout<<"document.onload = show_stats()\n";
                           
                           cout<<"function show_stats() {\n";
                           cout<<"var statwin_var=window.open(\"\",\'statwin\');\n";
                           cout<<"statwin_var.document.write (\'<img src=\""<<stat_image
                            <<"\" alt=\"Bad browser image\"/>\');\n";
                           cout<<"statwin_var.document.write (\'<b>POInT version "<<version<<"</b><br><br>\')\n";
                           cout<<"statwin_var.document.write (\'<b>Total events: </b> "<<num_events<<"<br><br>\')\n";
                           cout<<"statwin_var.document.write (\'<table> <tr><td class=\"title\"><b>Dataset</b></td>\');\n";
                           cout<<"statwin_var.document.write (\'<td class=\"title\"><b>Total Genomes</b></td>\');\n";
                           cout<<"statwin_var.document.write (\'<td class=\"title\"><b>Total Pillars</b></td>\');\n";
                           cout<<"statwin_var.document.write (\'<td class=\"title\"><b>Total Genes</b></td></tr>\');\n";
                           cout<<"statwin_var.document.write (\'<tr><td>"<<socket_id<<"</td>\');\n";
                           cout<<"statwin_var.document.write (\'<td>"<<genome_cnts[socket_id]<<"</td>\');\n";
                           cout<<"statwin_var.document.write (\'<td>"<<pillar_cnts[socket_id]<<"</td>\');\n";
                           cout<<"statwin_var.document.write (\'<td>"<<gene_cnts[socket_id]<<"</td></tr>\');\n";
                           cout<<"statwin_var.document.write (\'<tr><td>All events</td>\');\n";
                           cout<<"statwin_var.document.write (\'<td>"<<total_genomes<<"</td>\');\n";
                           cout<<"statwin_var.document.write (\'<td>"<<total_pillars<<"</td>\');\n";
                           cout<<"statwin_var.document.write (\'<td>"<<total_genes<<"</td></tr></table>\');\n";
                           cout<<"statwin_var.document.close();\n";
                           cout<<"}\n";
                       }
#endif
                       
                       cout<<"</script>\n";
                       
                       
                       
                      // Close the HTML document
                      cout << body() << html();
                   }
               }
           }
       }
       catch(exception& e) {
          // handle any errors - omitted for brevity
       }
    }
#if 0
    else {
        gene_name=argv[1];
        socket_name=argv[2];
        
        if (ping_image_server(socket_name, error_msg) == false) cout<<error_msg<<endl<<flush;
        else {
        
            get_image_data(gene_name, "12", socket_name, load_error, error_msg, &plot_ss);
         
            if (load_error==false) {
                fout.open("test.png");
                while(!(plot_ss.eof()))
                fout<<(char)plot_ss.get();
            //fout<<ss.str();
                fout.close();
            }
            else {
                cout<<"GOT error: "<<error_msg<<endl;
            }
        }
    }
#endif
}


void get_image_data(string gene_name, string window_size, string my_socketname, char imageformat, bool use_position, bool &error, string &error_msg, std::stringstream *plot_ss, int &pillar, int &start_pillar, int &end_pillar, int &last_pillar, IMAGE_SIZE mysize, int *map_coords, int tree_pillar, string &tree_string, string tree_type_string, bool is_random)
{
    int i, sock, rval, stop, end_data, data_size, data_cnt, tree_pos=0, tree_size;
    char dummy;
    struct sockaddr_un server;
    std::string my_message, pillar_string, tree_pillar_string;
    char buffer[BUFFER_SIZE], recv_buffer[BUFFER_SIZE], t_buf[BUFFER_SIZE];
    bool have_size=false;
    MSG_TYPE msg_type;
    std::stringstream ss;
    //std::ofstream fout;
    
    
    if (use_position==false)
        pillar=-1;
    
    stringstream ss2;
    ss2 << pillar;
    pillar_string = ss2.str();
    
    stringstream ss3;
    ss3 <<tree_pillar;
    tree_pillar_string=ss3.str();
    
    error=false;
    sock = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock < 0) {
        perror("opening stream socket");
        return;
    }
    
    server.sun_family = AF_UNIX;
    strcpy(server.sun_path, my_socketname.c_str());
    
    //cout<<"Trying to connect to "<<my_socketname<<endl<<flush;
    if (connect(sock, (struct sockaddr *) &server, sizeof(struct sockaddr_un)) < 0) {
        close(sock);
        //perror("connecting stream socket");
        //cout<<"FAiled connection"<<endl<<flush;
        error=true;
        error_msg="Error: failed to connect to the image server for event "+my_socketname;
        return;
    }
    
    
    strcpy(buffer, "D\t");
    if (use_position==false) {
        if (is_random==false)
            strcat(buffer, gene_name.c_str());
        else
            strcat(buffer, "NONE");
    }
    else
        strcat(buffer, "NONE");
    strcat(buffer, "\t");
    strcat(buffer, window_size.c_str());
    if (imageformat == 'N')
        strcat(buffer, "\tN");
    else
        strcat(buffer, "\tY");
    strcat(buffer, "\t");
    strcat(buffer, pillar_string.c_str());
    
    
    switch (mysize) {
        case LARGE:
            strcat(buffer, "\tX");
            break;
        case MEDIUM:
            strcat(buffer, "\tB");
            break;
        case SMALL:
            strcat(buffer, "\tL");
            break;
    }
        
    strcat(buffer, "\t");
    strcat(buffer, "R\t");
    strcat(buffer, tree_pillar_string.c_str());
  
    
    if (tree_type_string == "Phylip")
        strcat(buffer, "\tP");
    else
        strcat(buffer, "\tN");
    
    if (is_random==true)
        strcat(buffer, "\tY");
    else
        strcat(buffer, "\tN");
    
    
    if (write(sock, buffer, sizeof(buffer)) < 0)
        perror("writing on stream socket");
    //cout<<"Sent gene name "<<buffer<<endl<<flush;
    
    bzero(recv_buffer, sizeof(recv_buffer));
    
    rval = read(sock, recv_buffer, sizeof(recv_buffer));
    
    //cout<<"TRyed to read back: "<<rval<<endl<<flush;
    
    if (rval < 0)
        perror("reading stream message");
    else
        msg_type=get_message_type(recv_buffer);
    
    while (msg_type != ENDCOMM){
        
        switch (msg_type){
            case SIZE:
                //std::stringstream size_ss (recv_buffer);
                //cout<<"Got rev buffer: "<<recv_buffer<<std::endl<<std::flush;
                ss<<recv_buffer;
                ss>>dummy>>data_size>>pillar>>start_pillar>>end_pillar>>last_pillar;
                for(i=0; i<7; i++) ss>>map_coords[i];
                ss>>tree_size;
                ss.clear();
                have_size=true;
                //cout<<"Expecting "<<data_size<<" in data stream"<<" c0 is "<<map_coords[0]<<" Tree size "<<tree_size<<std::endl<<std::flush;
                data_cnt=0;
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                
                break;
                
                break;
            case PHYLOTREE:
                i=1;
                while((i<BUFFER_SIZE) && (tree_pos<tree_size)){
                        //tree_string[tree_pos]=recv_buffer[i];
                        t_buf[tree_pos]=recv_buffer[i];
                        tree_pos++;
                        i++;
                }
                t_buf[tree_pos]='\0';
                tree_string+=t_buf;
                
                //cout<<"Got part of tree: "<<" P: "<<tree_pos<<" i: "<<i<<" RB: "<<recv_buffer<<std::endl<<std::flush<<tree_string<<std::endl<<std::flush;
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                
                break;
            case DATABLOCK:
                //std::cout<<"Data message received"<<std::endl<<std::flush;
                //std::stringstream ss3 (recv_buffer);
                //std::cout<<"Received message of size "<<sizeof(recv_buffer)<<std::endl<<std::flush;
                
                if (have_size==false) {
                    end_data=BUFFER_SIZE-1;
                    
                    while(recv_buffer[end_data] == '\0') end_data--;
                    end_data++;
                    
                    
                    //std::cout<<"Reading until "<<end_data<<" of "<<BUFFER_SIZE<<" in mssg"<<std::endl<<std::flush;
                    for(i=1; i<end_data; i++) {
                        //if (recv_buffer[i] != '\0')
                            *plot_ss<<recv_buffer[i];
                    }
                    //std::cout<<"Waiting for next message"<<std::endl<<std::flush;
                }
                else {
                    i=1;
                    while((i<BUFFER_SIZE) && (data_cnt < data_size)) {
                        //if ((int)recv_buffer[i] > 0)
                            *plot_ss<<recv_buffer[i];
                        data_cnt++;
                        i++;
                    }
                    
                }
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                
                break;
                
            case MESSAGE:
                //std::cout<<"Message received: "<<std::endl<<std::flush;
                //std::stringstream ss (recv_buffer);
                ss<<recv_buffer;
                error_msg=ss.str();
                //ss>>error_msg;
                error_msg=error_msg.substr(1, error_msg.length()-1);
            
                //std::cout<<"Message= "<<error_msg<<std::endl<<std::flush;
                ss.clear();
                error=true;
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                break;
                
            case BATCHBLOCK:
                //std::cout<<"Message received: "<<std::endl<<std::flush;
                //std::stringstream ss (recv_buffer);
                
                error_msg="ERROR: Browser received batch download data";
                //ss>>error_msg;
                //error_msg=error_msg.substr(1, error_msg.length()-1);
                
                //std::cout<<"Message= "<<error_msg<<std::endl<<std::flush;
                ss.clear();
                error=true;
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                break;
        }
    }
    
    
}

bool ping_image_server(string my_socketname, string &error_msg, int &num_genomes, int &num_pillars, int &num_genes, string &version)
{
    int i, sock, rval;
    bool retval=false;
    struct sockaddr_un server;
    std::string buf_data;
    char dummy, buffer[BUFFER_SIZE], recv_buffer[BUFFER_SIZE];
    MSG_TYPE msg_type;
    std::stringstream ss;
    
    
    sock = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock < 0) {
        //perror("opening stream socket");
        error_msg="Error openning stream socket";
        return(retval);
    }
    
    server.sun_family = AF_UNIX;
    strcpy(server.sun_path, my_socketname.c_str());
    
    //cout<<"Trying to PING "<<my_socketname<<endl<<flush;
    if (connect(sock, (struct sockaddr *) &server, sizeof(struct sockaddr_un)) < 0) {
        close(sock);
        //perror("connecting stream socket");
        //cout<<"FAiled connection"<<endl<<flush;
        error_msg="Error connecting to stream socket";
        return(retval);
    }
    
    strcpy(buffer, "PPINGSEND");
    
    if (write(sock, buffer, sizeof(buffer)) < 0)
        //perror("writing on stream socket");
        error_msg="Error writing on stream socket";
        //cout<<"Sent ping--getting response"<<endl<<flush;
    
    bzero(recv_buffer, sizeof(recv_buffer));
    
    rval = read(sock, recv_buffer, sizeof(recv_buffer));
    
    //cout<<"TRyed to read back: "<<rval<<endl<<flush;
    
    if (rval < 0)
        //perror("reading stream message");
        error_msg="Error reading stream message";
    else {
        msg_type=get_message_type(recv_buffer);
     
        if (msg_type ==PING) {
            retval=true;
            ss<<recv_buffer;
            ss>>dummy>>buf_data>>num_genomes>>num_pillars>>num_genes>>version;
        }
        bzero(recv_buffer, sizeof(recv_buffer));
        
        rval = read(sock, recv_buffer, sizeof(recv_buffer));
        
        //cout<<"Got PING, waiting for END: "<<rval<<endl<<flush;
        
        if (rval < 0) {
            //perror("reading stream message");
            error_msg="Error reading stream message";
            retval=false;
        }
        else {
            
            msg_type=get_message_type(recv_buffer);
            if (msg_type !=ENDCOMM) retval=false;
            else close(sock);
        }
    }
    
    
    
    return(retval);
    
}
    
MSG_TYPE get_message_type (char* buffer)
{
    MSG_TYPE msg_type=ENDCOMM;
    
    if (buffer[0] == 'M') msg_type=MESSAGE;
    if (buffer[0] == 'D') msg_type=DATABLOCK;
    if (buffer[0] == 'S') msg_type=SIZE;
    if (buffer[0] == 'E') msg_type=ENDCOMM;
    if (buffer[0] == 'P') msg_type=PING;
    if (buffer[0] == 'T') msg_type=PHYLOTREE;
    if (buffer[0] == 'B') msg_type=BATCHBLOCK;
    if (buffer[0] == 'I') msg_type=MODELDIAG;
    if (buffer[0] == 'V') msg_type=TREEDIAG;
    
    //std::cout<<"MSG KEY: "<<buffer[0]<<std::endl<<std::flush;
    
    switch (msg_type) {
        case DATABLOCK:
            //std::cout<<"Data message received"<<std::endl<<std::flush;
            break;
        case MESSAGE:
            //std::cout<<"Message received"<<std::endl<<std::flush;
            break;
        case ENDCOMM:
            //std::cout<<"End connection message received"<<std::endl<<std::flush;
            break;
        case PHYLOTREE:
            //std::cout<<"Message is tree data"<<std::endl<<std::flush;
            break;
    }
    
    return(msg_type);
}


#define MAXLINE 1000
int get_pdf_image(std::stringstream& image_instream, std::stringstream& image_outstream)
{
    pid_t pid;
    int i, rv, readval, cnt=0;
    int    commpipe1[2], commpipe2[2];
    char dummy[1], buf[MAXLINE];
    std::string retval="ERROR", passdata;
    //std::stringstream image_outstream;
    
    
    image_outstream.clear();
    //std::cout<<"Initial Length of cleaned stream is "<<image_outstream.str().length()<<std::endl<<std::flush;
    
    
    /* Setup communication pipeline first */
    if ( (pipe(commpipe1) < 0) || (pipe(commpipe2) < 0) )
    {
        std::cerr << "PIPE ERROR" << std::endl;
        return(-1);
    }
    
    
    /* Attempt to fork and check for errors */
    if( (pid=fork()) == -1){
        std::cerr<<"Fork error.\n";  /* something went wrong */
        return(-1);
    }
    
    if(pid){
        /* A positive (non-negative) PID indicates the parent process */
        close(commpipe1[0]);
        close(commpipe2[1]);
        //std::cout<<"Intermed. Length of cleaned stream is "<<image_outstream.str().length()<<std::endl;
        
#if 0
        while (!(image_instream.eof())) {
            i=0;
            strcpy(buf, "");
            while((!(image_instream.eof())) && (i<MAXLINE)) {
                dummy[0]=(char)image_instream.get();
                strcat(buf, dummy);
                i++;
            }
                
           // if (!(image_instream.eof()))
           if (i<0)
               write(commpipe1[1],buf,i);
        }
#else
        while(!(image_instream.eof())) {
            dummy[0]=(char)image_instream.get();
            if (!(image_instream.eof()))
                write(commpipe1[1],dummy,1);
        }
#endif
     
        close(commpipe1[1]);
        //wait();
        readval = read(commpipe2[0], dummy, 1);
       
        if ( readval  < 0 )
        {
            std::cerr << "READ ERROR FROM PIPE" << std::endl;
            //return(-1);
        }
        else if (readval == 0)
        {
            std::cerr << "Child Closed Pipe" << std::endl;
            return(-1);
        }
       
        while(readval >0) {
            image_outstream<<dummy[0];
            
            cnt++;
            readval=read(commpipe2[0], dummy, 1);
        }
        
        
        
        close(commpipe2[0]);
        
        //std::cout << "OUTPUT of PROGRAM B is: " << line;
        //std::cout<<"Put "<<cnt<<" items in stream\n"<<std::flush;
        //std::cout<<"New Length of cleaned stream is "<<image_outstream.str().length()<<std::endl;
        
        
        return(0);
        
    }
    else{
        close(commpipe1[1]);
        close(commpipe2[0]);
        
        if (commpipe1[0] != STDIN_FILENO)
        {
            if (dup2(commpipe1[0], STDIN_FILENO) != STDIN_FILENO)
            {
                std::cerr << "dup2 error to stdin" << std::endl;
            }
            close(commpipe1[0]);
            //return(-1);
        }
        
        if (commpipe2[1] != STDOUT_FILENO)
        {
            if (dup2(commpipe2[1], STDOUT_FILENO) != STDOUT_FILENO)
            {
                std::cerr << "dup2 error to stdout" << std::endl;
            }
            close(commpipe2[1]);
            //return(-1);
        }
        
        if ( execl("/usr/bin/ps2pdf", "ps2pdf", "-dEPSCrop", "-", "-",  (char *)0) < 0 )
            //if ( execl("/usr/local/bin/magick", "magick", (char *)0) < 0 )
        {
            //printf("execl returned! errno is [%d]\n",errno);
            //perror("The error message is :");
            
            std::cerr << "system error" << std::endl;
            return(-1);
        }
        
        return(-1);
    }
    return(-1);
    
    
    
    
}


int clean_image(std::stringstream& image_instream, std::stringstream& image_outstream)
{
    pid_t pid;
    int i, rv, readval, cnt=0;
    int    commpipe1[2], commpipe2[2];
    char dummy[1];
    std::string retval="ERROR", passdata;
    //std::stringstream image_outstream;
    
    
    image_outstream.clear();
    //std::cout<<"Initial Length of cleaned stream is "<<image_outstream.str().length()<<std::endl<<std::flush;
    
    
    /* Setup communication pipeline first */
    if ( (pipe(commpipe1) < 0) || (pipe(commpipe2) < 0) )
    {
        std::cerr << "PIPE ERROR" << std::endl;
        return(-1);
    }
    
    
    /* Attempt to fork and check for errors */
    if( (pid=fork()) == -1){
        std::cerr<<"Fork error.\n";  /* something went wrong */
        return(-1);
    }
    
    if(pid){
        /* A positive (non-negative) PID indicates the parent process */
        close(commpipe1[0]);
        close(commpipe2[1]);
        //std::cout<<"Intermed. Length of cleaned stream is "<<image_outstream.str().length()<<std::endl;
        
        
        while (!(image_instream.eof())) {
            dummy[0]=(char)image_instream.get();
            if (!(image_instream.eof()))
                write(commpipe1[1],dummy,1);
        }
       
        close(commpipe1[1]);
        //wait();
        readval = read(commpipe2[0], dummy, 1);
       
        if ( readval  < 0 )    {
            std::cerr << "READ ERROR FROM PIPE" << std::endl;
            //return(-1);
        }
        else if (readval == 0)  {
            std::cerr << "Child Closed Pipe" << std::endl;
            return(-1);
        }
        //image_outstream<<dummy;
        while(readval >0) {
            image_outstream<<dummy[0];
            cnt++;
            readval=read(commpipe2[0], dummy, 1);
            //readval = read(commpipe2[0], buf, MAXLINE);
        }
        
        
        
        close(commpipe2[0]);
        
        //std::cout << "OUTPUT of PROGRAM B is: " << line;
        //std::cout<<"Put "<<cnt<<" items in stream\n"<<std::flush;
        //std::cout<<"New Length of cleaned stream is "<<image_outstream.str().length()<<std::endl;
        
        
        return(0);
        
    }
    else{
        close(commpipe1[1]);
        close(commpipe2[0]);
        
        if (commpipe1[0] != STDIN_FILENO)
        {
            if (dup2(commpipe1[0], STDIN_FILENO) != STDIN_FILENO)
            {
                std::cerr << "dup2 error to stdin" << std::endl;
            }
            close(commpipe1[0]);
            //return(-1);
        }
        
        if (commpipe2[1] != STDOUT_FILENO)
        {
            if (dup2(commpipe2[1], STDOUT_FILENO) != STDOUT_FILENO)
            {
                std::cerr << "dup2 error to stdout" << std::endl;
            }
            close(commpipe2[1]);
            //return(-1);
        }
        
        if ( execl("/usr/bin/convert", "magick", "png:-", "-resize", "50%", "-antialias", "png:-", (char *)0) < 0 ) {
            std::cerr << "system error" << std::endl;
            return(-1);
        }
        
        return(-1);
    }
    return(-1);
    
    
    
}


void get_modeltree_data(string my_socketname, std::stringstream *model_ss, bool is_model)
{
    int i, sock, rval, stop, end_data, data_size, data_cnt;
    char dummy;
    struct sockaddr_un server;
    std::string my_message, error_msg;
    char buffer[BUFFER_SIZE], recv_buffer[BUFFER_SIZE], t_buf[BUFFER_SIZE];
    bool have_size=false;
    MSG_TYPE msg_type;
    std::stringstream ss;
    
    
    //error=false;
    sock = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock < 0) {
        perror("opening stream socket");
        return;
    }
    
    server.sun_family = AF_UNIX;
    strcpy(server.sun_path, my_socketname.c_str());
    
    //cout<<"Trying to connect to "<<my_socketname<<endl<<flush;
    if (connect(sock, (struct sockaddr *) &server, sizeof(struct sockaddr_un)) < 0) {
        close(sock);
        //perror("connecting stream socket");
        //cout<<"FAiled connection"<<endl<<flush;
        //error=true;
        error_msg="Error: failed to connect to the image server for event "+my_socketname;
        return;
    }
    
    if (is_model == true)
        strcpy(buffer, "I\n");
    else
        strcpy(buffer, "V\n");
    
    if (write(sock, buffer, sizeof(buffer)) < 0)
        perror("writing on stream socket");
    
    bzero(recv_buffer, sizeof(recv_buffer));
    rval = read(sock, recv_buffer, sizeof(recv_buffer));
    if (rval < 0)
        perror("reading stream message");
    else
        msg_type=get_message_type(recv_buffer);
    
    while (msg_type != ENDCOMM){
        
        switch (msg_type){
            case SIZE:
                //std::stringstream size_ss (recv_buffer);
                //cout<<"Got rev buffer: "<<recv_buffer<<std::endl<<std::flush;
                ss<<recv_buffer;
                ss>>dummy>>data_size;
                ss.clear();
                have_size=true;
                //cout<<"Expecting "<<data_size<<" in data stream"<<" c0 is "<<map_coords[0]<<" Tree size "<<tree_size<<std::endl<<std::flush;
                data_cnt=0;
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                
                break;
                
            case MODELDIAG:
            case TREEDIAG:
                //std::cout<<"Data message received"<<std::endl<<std::flush;
                //std::stringstream ss3 (recv_buffer);
                //std::cout<<"Received message of size "<<sizeof(recv_buffer)<<std::endl<<std::flush;
                
                if (have_size==false) {
                    end_data=BUFFER_SIZE-1;
                    
                    while(recv_buffer[end_data] == '\0') end_data--;
                    end_data++;
                    
                    
                    //std::cout<<"Reading until "<<end_data<<" of "<<BUFFER_SIZE<<" in mssg"<<std::endl<<std::flush;
                    for(i=1; i<end_data; i++) {
                        //if (recv_buffer[i] != '\0')
                        *model_ss<<recv_buffer[i];
                    }
                    //std::cout<<"Waiting for next message"<<std::endl<<std::flush;
                }
                else {
                    i=1;
                    while((i<BUFFER_SIZE) && (data_cnt < data_size)) {
                        //if ((int)recv_buffer[i] > 0)
                            *model_ss<<recv_buffer[i];
                        //std::cout<<"RB: "<<i<<" - "<<(int)recv_buffer[i]<<std::endl;
                        //if ((int)recv_buffer[i]<=0) cnt_neg++;
                        data_cnt++;
                        i++;
                    }
                    
                }
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                
                break;
                
            case MESSAGE:
                //std::cout<<"Message received: "<<std::endl<<std::flush;
                //std::stringstream ss (recv_buffer);
                ss<<recv_buffer;
                ss>>error_msg;
                //std::cout<<"Message= "<<error_msg<<std::endl<<std::flush;
                ss.clear();
                //error=true;
                
                bzero(recv_buffer, sizeof(recv_buffer));
                rval = read(sock, recv_buffer, sizeof(recv_buffer));
                if (rval < 0)
                    perror("reading stream message");
                else msg_type=get_message_type(recv_buffer);
                break;
        }
    }
    
    close(sock);
    //printf("Ending connection\n");
    //std::cout<<"Data count was "<<data_cnt<<" neg was : "<<cnt_neg<<std::endl;
    
    
}
