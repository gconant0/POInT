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

#include <boost/filesystem.hpp>

using namespace std;
using namespace cgicc;




#define MAX_NAME_LEN 100
#define BUFFER_SIZE 2048

enum MSG_TYPE {MESSAGE, DATABLOCK, SIZE, PING, ENDCOMM, PHYLOTREE, BATCHBLOCK, MODELDIAG, TREEDIAG};
enum EVENT_TYPE {WGD, WGT};

//Receive message types from POInT
MSG_TYPE get_message_type (char* buffer);

void get_data (std::string tree_header, std::string CDS_header, double ortho_cut, int num_miss, std::string my_socketname, std::string tree_type_string,std::string subgenome, bool orthos_only, bool &error, std::string &error_msg, std::stringstream *tar_ss);
int tar_files(int num_pillars, std::string *treefiles, std::string *CDS_files, std::string tree_header, std::string CDS_header, std::stringstream *tar_outstream);


bool ping_image_server(string my_socketname, string &error_msg);

bool does_have_events(std::map<std::string, std::string> socket_hash, std::map<std::string, bool> &valid_map, int &num_events)
{
    string error_msg;
 
    
    num_events=0;
    for (std::map<string,string>::iterator it=socket_hash.begin(); it!=socket_hash.end(); ++it) {
        //cout<<"Testing "<<it->first<<": "<<it->second<<std::endl;
        if (ping_image_server(it->second, error_msg)==true) {
            //if (!(fin.fail())) {
             //cout<<" Success"<<endl;
            valid_map[it->first]=true;
            num_events++;
        }
        
    }
    if (num_events==0)
        return(false);
    
    else return(true);
}

// Print the form for this CGI
void printForm(const Cgicc& cgi, std::map<std::string, std::string> socket_hash, std::map<std::string, bool> valid_map, std::map<std::string, bool> is_degen, std::map<std::string, int> event_sizes, int num_events, std::string server_path, std::map<std::string, EVENT_TYPE> type_hash)
{
    int i;
    bool have_event=false, orthos_only=false, real_submit=false, picked_event=false;
    string error_msg, default_event, my_event, default_thresh, default_missing_genes, submit_val, href_path;
    ifstream fin;
    
   
    
    default_thresh="0.9";
    
    cout << "<form name=\"downloadform\" id=\"downloadform\" method=\"post\" action=\""
    << cgi.getEnvironment().getScriptName()
        << "\" enctype=\"multipart/form-data\">" << endl;
   
    const_form_iterator event  = cgi.getElement("event");
    if (event != cgi.getElements().end()) {
        default_event=**event;
        have_event=true;
        //cout<<"Using socket "<<socket_name<<endl;
    }
   
    
    const_form_iterator orthosonly = cgi.getElement("orthos_only");
    if (orthosonly != cgi.getElements().end()) {
        if (**orthosonly == "Orthos") orthos_only=true;
    }
    
    const_form_iterator realsub = cgi.getElement("eventsubmit");
    if (realsub != cgi.getElements().end()) {
        if (**realsub == "Y")  real_submit=true;
    }
    
    cout << "<table width=900>" << endl;
    cout << "<tr>";
    if (num_events>0) {
        href_path="http://"+server_path+"/help.html#event_ref";
        //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Event:</a></b></td>\n"
        cout<<"<td class=\"form\" style=\"text-align:left\">"<<"<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Event:</a></b>&nbsp;&nbsp;<select id=\"event\" name = \"event\" onchange=\"redraw_form()\">";
        
        
        for (std::map<string,string>::iterator it=socket_hash.begin(); it!=socket_hash.end(); ++it) {
            if (valid_map[it->first]==true) {
                if (have_event==true) {
                    if (it->first==default_event) {
                        cout<<"<option selected=\"selected\" value=\""<<it->first<<"\">"<<it->first<<"</option> \n";
                        picked_event=true;
                        my_event=default_event;
                    }
                    else
                        cout<<"<option  value=\""<<it->first<<"\">"<<it->first<<"</option> \n";
                    
                }
                else {
                    cout<<"<option value=\""<<it->first<<"\">"<<it->first<<"</option> \n";
                }
            }
        }
        
        if (picked_event==false) {
            for (std::map<string,string>::iterator it=socket_hash.begin(); it!=socket_hash.end(); ++it) {
                if (valid_map[it->first]==true) {
                    if (picked_event==false) {
                        picked_event=true;
                        my_event=it->first;
                    }
                }
            }
        }
        
        cout<< "</select></td>\n";
        
        
    }
    
    href_path="http://"+server_path+"/POInTDownload_help.html#singleortho_ref";
    //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Single-copy orthologs only?:</a></b></td>\n";
    cout<<"<td class=\"form\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Single-copy orthologs only?</a></b>&nbsp;&nbsp;<input type=\"checkbox\" name=\"orthos_only\" value=\"Orthos\"";
    if (orthos_only == true) cout<<"checked=\"checked\"";
    cout<<" onchange=\"redraw_form()\"></td>\n";
    
    if (((have_event==true) || (picked_event==true)) && (orthos_only ==false)) {
        href_path="http://"+server_path+"/POInTDownload_help.html#missingcnt_ref";
       // cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Maximum number of missing genes:</a></b></td>\n"
        cout<<"<td class=\"form\">"<< "<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Maximum number of missing genes:</a></b>&nbsp;&nbsp;<select id=\"missing_genes\" name = \"missing_genes\" ongotfocus=\"redraw_form()\">";
        if (type_hash[my_event] == WGD) {
            for(i=0; i<=event_sizes[my_event]; i++) {
                cout<<"<option ";
                if (i ==0) cout<<" selected=\"selected\" ";
                cout<<" value=\""<<i<<"\">"<<i<<"</option> ";
            }
        }
        else {
            for(i=0; i<=2*event_sizes[my_event]; i++) {
                cout<<"<option ";
                if (i ==0) cout<<" selected=\"selected\" ";
                cout<<" value=\""<<i<<"\">"<<i<<"</option> ";
            }
        }
        cout<< "</td>\n";
    }
    else {
        if ((picked_event==true) && (orthos_only == true)) {
            href_path="http://"+server_path+"/POInTDownload_help.html#subgenomes_ref";
            //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Subgenomes:</a></b></td>\n"
            cout<<"<td class=\"form\">"<<"<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Subgenomes:</a></b>&nbsp;&nbsp;<select id=\"subgenome\" name = \"subgenome\" ongotfocus=\"redraw_form()\">";
                cout<<"<option "<<" selected=\"selected\" value=\"All\">All</option>\n";
            
            if (type_hash[my_event] == WGD) {
                if (is_degen[my_event] == false) {
                    cout<<"<option "<<"  value=\"Least\">Least fractionated</option>\n";
                    cout<<"<option "<<"  value=\"Most\">Most fractionated</option>\n";
                }
            }
            else {
                cout<<"<option "<<"  value=\"Least\">Least fractionated</option>\n";
                cout<<"<option "<<"  value=\"Intermed\">Intermed. fractionated</option>\n";
                cout<<"<option "<<"  value=\"Most\">Most fractionated</option>\n";
            }
            cout<< "</td>\n";
        }
    }
    
   
    
    cout<<"</tr><tr>";
    
    href_path="http://"+server_path+"/POInTDownload_help.html#cutoff_ref";
    //cout <<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Orthology Confidence cutoff:</a></b></td>\n"
    cout<< "<td class=\"form\">"<< "<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Orthology Confidence cutoff:</a></b>&nbsp;&nbsp;<select id=\"ortho_cut\" name = \"ortho_cut\">";
    //<< "<input type=\"text\" name=\"wind_size\" accept=\"text/plain\" />"
   // if ( (have_event== false) || (is_degen[my_event] == false)) {
        cout<<"<option value=\"0.5\">"<<0.5<<"</option> \n";
        cout<<"<option value=\"0.6\">"<<0.6<<"</option> \n";
        cout<<"<option value=\"0.7\">"<<0.7<<"</option> \n";
        cout<<"<option value=\"0.8\">"<<0.8<<"</option> \n";
        cout<<"<option selected=\"selected\" value=\"0.9\">"<<0.9<<"</option> \n";
        cout<<"<option value=\"0.95\">"<<0.95<<"</option> \n";
        cout<<"<option value=\"0.99\">"<<0.99<<"</option> \n";
        //cout<<"<option value=\"1.1\">"<<1.1<<"</option> \n";
   // }
#if 0
    else {
        cout<<"<option value=\"0.25\">"<<0.25<<"</option> \n";
        cout<<"<option value=\"0.3\">"<<0.3<<"</option> \n";
        cout<<"<option value=\"0.35\">"<<0.35<<"</option> \n";
        cout<<"<option value=\"0.4\">"<<0.4<<"</option> \n";
        cout<<"<option selected=\"selected\" value=\"0.45\">"<<0.45<<"</option> \n";
        cout<<"<option value=\"0.475\">"<<0.475<<"</option> \n";
        cout<<"<option value=\"0.495\">"<<0.495<<"</option> \n";
        //cout<<"<option value=\"0.495\">"<<0.555<<"</option> \n";
    }
#endif
    cout<< "</select></td>\n";
    href_path="http://"+server_path+"/help.html#treeformat_ref";
    //cout<<"<td class=\"title\"><b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Tree Format:</a></b></td>\n";
    cout<<"<td class=\"form\">";
    cout<<"<b><a href=\""<<href_path<<"\" onClick=\"return open_help(this)\">Tree Format:</a></b>&nbsp;&nbsp;";
    cout<<" <input type=\"radio\" name=\"treeformat\" value=\"NEXUS\" checked=\"checked\">NEXUS";
    cout<<" <input type=\"radio\" name=\"treeformat\" value=\"Phylip\" >Newick</td>\n";
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"eventsubmit\" name=\"eventsubmit\" value=";
    if (real_submit==true) cout<<"\"Y\"";
    else cout<<"\"N\"";
    cout<<"></td>\n"<<endl;
    cout<<"<td class=\"form\"><input type=\"hidden\" id=\"orthossubmit\" name=\"orthossubmit\" value=\"-1\"></td>"<<endl;
    cout<<"</tr>" << endl;
  
    
    cout << "</table>\n" << endl;
    
    cout << "<div class=\"center\"><p>"
        << "<input type=\"submit\" name=\"drawbutton\" id=\"drawbutton\" value=\"Download\" onclick=\"submit_form()\" >"<< "</p></div></form>" << endl;
}


int main(int argc, char **argv)
{
    int i, num_miss=-1, num_events, cds_pillar, x1, x2, button_size, total, event_size;
    char myimageformat='Y', printchar, dummy;
    double ortho_cut=-1.0;
    bool  load_data=false, load_error, have_events, orthos_only=false, real_submit=false;
    string error_msg, new_name, socket_name, polyploidy_files, polyploidy_name, polyploidy_socket, dump, help_string, cite_string,
        socket_id, ortho_cut_string, cds_num, line, formatstring, click_string, treefile, modelfile, homologfile, num_miss_string, tree_type_string, tree_header, *tree_files, *full_tree_files, CDS_header, *CDS_strings, pid_string, server_path, href_path, image_path, subgenome_string;
    std::stringstream tar_ss;
    std::ifstream pp_in;
    std::ofstream fout;
    std::map<std::string, bool> valid_map, is_degen;
    std::map<std::string, std::string> socket_hash;
    std::map<std::string, int> event_size_hash;
    std::map<std::string, EVENT_TYPE> type_hash;
    pid_t my_pid;
    //EVENT_TYPE my_type=WGD;
 
    
     
    polyploidy_files="/data/POInT/POInT_browser_events_cleaned_prefix_commonname.txt";
    pp_in.open(polyploidy_files.c_str());
    
    if (const char* server_p=std::getenv("HTTP_HOST"))
        server_path=server_p;
    else server_path="UNKNOWN";
    
    if (pp_in.fail()) {
        socket_hash["GrassRho"]="/data/POInT/POInTGrass";
    }
    else {
        while (! pp_in.eof()) {
            polyploidy_name="";
            pp_in>>polyploidy_name>>polyploidy_socket>>treefile>>modelfile>>homologfile;
            event_size=0;
            
           
            
            std::getline(pp_in, line);
            std::istringstream ssp(line);
            
            while (! ssp.eof()) {
                ssp>>dump;
                if (dump[0] != '-') {
                    event_size++;
                }
            }
            if (polyploidy_name != "YeastWGD")
            event_size = event_size/2;
          
            
            //cout<<"Got event "<<polyploidy_name<<" at "<<polyploidy_socket<<endl;
            
            if(polyploidy_name!= "") {
                socket_hash[polyploidy_name]=polyploidy_socket;
                event_size_hash[polyploidy_name]=event_size;
                
                std::size_t found3 = modelfile.find("WGT");
                if (found3 != std::string::npos)
                    type_hash[polyploidy_name]=WGT;
                else
                    type_hash[polyploidy_name]=WGD;
                
                std::size_t found = modelfile.find("bias");
                std::size_t found2 = modelfile.find("WGD");
                if (found2 == std::string::npos)
                    is_degen[polyploidy_name]=false;
                else {
                    if (found!=std::string::npos) is_degen[polyploidy_name]=false;
                    else is_degen[polyploidy_name]=true;
                }
                
            }
        }
        pp_in.close();
    }
    
   try {
       Cgicc cgi;
       
       
       form_iterator event  = cgi.getElement("event");
       if (event != cgi.getElements().end()) {
           socket_id=**event;
           socket_name=socket_hash[**event];
       }
       
       const_form_iterator orthosonly = cgi.getElement("orthos_only");
       if (orthosonly != cgi.getElements().end()) {
           if (**orthosonly == "Orthos") orthos_only=true;
       }
       
       form_iterator orthocut  = cgi.getElement("ortho_cut");
       if (orthocut != cgi.getElements().end()) {
           ortho_cut_string=**orthocut;
           ortho_cut=std::stod(ortho_cut_string);
       }
       
       form_iterator nummiss = cgi.getElement("missing_genes");
       if (nummiss != cgi.getElements().end()) {
           num_miss_string=**nummiss;
           num_miss=std::stoi(num_miss_string, nullptr);
       }
       
       form_iterator treetype = cgi.getElement("treeformat");
       if (treetype != cgi.getElements().end()) {
           tree_type_string=**treetype;
       }
       
       form_iterator subgen = cgi.getElement("subgenome");
       if (subgen != cgi.getElements().end()) {
           subgenome_string=**subgen;
       }
       
       const_form_iterator realsub = cgi.getElement("eventsubmit");
       if (realsub != cgi.getElements().end()) {
           if (**realsub == "Y")  real_submit=true;
       }
       
       load_error=true;
       have_events=does_have_events(socket_hash, valid_map, num_events);
      
       if (have_events == true) {
           load_error=false;
           
           if ((ortho_cut<0.0) ||(ortho_cut>1.0)) {
               error_msg= "Invalid orthology confidence cutoff selected\n";
               load_error=true;
               load_data=false;
           }
           
           if ((num_miss ==-1) && (orthos_only==false)) {
               error_msg= "Invalid cutoff for the number of duplicates selected\n";
               load_error=true;
               load_data=false;
           }
           if ((load_error == false) && (valid_map[socket_id]==true) && (real_submit==true)) {
                tar_ss.str(std::string());
               
            my_pid=getpid();

            std::stringstream ss_d;
            ss_d << my_pid;
            pid_string = ss_d.str();

            //cout<<"PID string :"<<pid_string<<endl;
            //pid_string ="fred";
               
               tree_header="/media/ramdisk/POInT/" + pid_string + "/";
               CDS_header="/var/www/html/Downloads/" + socket_id + "/Pillars/";
               //CDS_header="/var/www/html/Downloads/" +
               
               
               get_data (tree_header, CDS_header, ortho_cut, num_miss, socket_name, tree_type_string, subgenome_string, orthos_only,  load_error, error_msg, &tar_ss);
               load_data=true;
               
               //load_data=false;
               
           }
           
       }

       
       if ((load_error == false)  && (load_data==true) &&(real_submit==true)  ) {
           cout<<"Content-Type:application/x-download\n";
           cout<<"Content-Disposition:attachment;filename=POInT_data.tar\n\n";
           //cout << HTTPHTMLHeader() << endl;
           
           // Set up the HTML document
#if 0
           cout << html() << head(title("POInT download")) << endl;
           cout <<"TRyinf to download tar\n";
           cout<<"Tree dir: "<<tree_header<<endl;
           cout<<" SI: "<<socket_id<<endl;
           cout<<"LE: "<<load_error<<endl;
           cout<<"realsub:"<<real_submit<<endl;
           cout<<"CDS dir; " <<CDS_header<<endl;
           
           cout << body() << endl;
#else
           while(!(tar_ss.eof())) {
               dummy=(char)tar_ss.get();
               
               if (!(tar_ss.eof())) {
                   std::cout<<dummy;
               }
           }
           
           boost::filesystem::remove_all(tree_header);
#endif
          //cout << body() << html();
       }
       else {
          
        // Send HTTP header
        cout << HTTPHTMLHeader() << endl;

        // Set up the HTML document
        cout << html() << head(title("POInT download")) << endl;
           
           //cout<<"Tree dir: "<<tree_header<<endl;
           //cout<<"CDS dir; " <<CDS_header<<endl;
           //cout<<" SI: "<<socket_id<<endl;
           //cout<<"LE: "<<load_error<<endl;
           //cout<<"realsub:"<<real_submit<<endl;
           //cout<<"PID: "<<pid_string<<endl;
        cout << body() << endl;
           
           href_path="http://"+server_path+"/POInTDownload_help.html";
           image_path="http://"+server_path+"/Download_Header.png";
           
           cout<<"<map NAME=\"headmap\">";
           cout<<"<area SHAPE=\"RECT\" COORDS=\"783,68,963,108\" HREF=\""<<href_path<<"\"  onClick=\"return open_help(this)\">";
           cout<<"</map>\n";
           cout<<"<img usemap=\"#headmap\" src=\""<<image_path<<"\" alt=\"Bad broswer image\"/>\n";
           
           
           //cout << img().set("SRC", image_path)
           //.set("ALT", "Bad Header image") << endl;
        printForm(cgi, socket_hash, valid_map, is_degen, event_size_hash, num_events, server_path, type_hash);
           
           
        if ((real_submit == true) &&  (load_error==true) ) {
               //cout <<"<br>"<<endl;
            
            cout<<"<b>Error: </b>";
            cout<<error_msg<<"<br>"<<endl;
            
        }
           
           
           
           
       cout<<"<script type=\"text/javascript\">\n";
           cout<<"var realsubval = document.getElementById(\"eventsubmit\")\n";
       cout<<"function redraw_form() {\n";
       cout<<"\trealsubval.setAttribute(\'value\', \'N\');\n";
       cout<<"\tdocument.forms[\"downloadform\"].submit();\n";
       cout<<"}\n";
       cout<<"function submit_form(val) {\n";
       cout<<"\trealsubval.setAttribute(\'value\', \'Y\');\n";
       cout<<"\tdocument.forms[\"downloadform\"].submit();\n";
       cout<<"}\n";
       
       cout<<"function open_help(helphttp) {\n";
       cout<<"\tif (! window.focus) return true;\n";
       //cout<<"\t var href;\n";
       //cout<<"\t href=\""<<help_string<<"\";\n";
       cout<<"var href;\n";
       cout<<"if (typeof(helphttp) == \'string\') href=helphttp;\n";
       cout<<"else href=helphttp.href;\n";
       cout<<"\twindow.open(href,\'helpwin\', \'left=20,top=20,width=1050,height=700,toolbar=1,scrollbars=1,resizable=1\');\n";
       cout<<"\treturn false;\n}\n";
       
      
       cout<<"</script>\n";
           
           
           
       // Close the HTML document
       cout << body() << html();
        
       }
   }
   catch(exception& e) {
      // handle any errors - omitted for brevity
   }

}


void get_data (std::string tree_header, std::string CDS_header, double ortho_cut, int num_miss, std::string my_socketname, std::string tree_type_string, std::string subgenome,  bool orthos_only,   bool &error, std::string &error_msg, std::stringstream *tar_ss)
{
    int i, sock, rval, stop, end_data, data_size, pillar, num_pillars,  new_pillar, pillar_cnt, data_cnt, tree_size, tree_pos, last_pillar, val[1], cnt_neg=0;
    char dummy;
    char dummy2[1];
    
    struct sockaddr_un server;
    std::string socketname, my_message, tree_string, tree_file, CDS_file, cutoff_string, missing_string, *tree_strings, *full_tree_strings, *CDS_strings;
    char buffer[BUFFER_SIZE], recv_buffer[BUFFER_SIZE], t_buf[BUFFER_SIZE];
    bool have_size;
    MSG_TYPE msg_type;
    std::stringstream ss, data_ss;
    std::ofstream fout;
    int errorp;
    
    error=false;
    error_msg="";
    
    sock = socket(AF_UNIX, SOCK_STREAM, 0);
    if (sock < 0) {
        perror("opening stream socket");
        return;
    }
    
    socketname=my_socketname;
    server.sun_family = AF_UNIX;
    strcpy(server.sun_path, socketname.c_str());
    
    
    if (connect(sock, (struct sockaddr *) &server, sizeof(struct sockaddr_un)) < 0) {
        close(sock);
        error=true;
        error_msg="Error: failed to connect to the image server for event "+my_socketname;
        return;
    }
    
    
    stringstream ss2;
    ss2 << ortho_cut;
    cutoff_string = ss2.str();
    
    stringstream ss3;
    ss3 << num_miss;
    missing_string = ss3.str();
    
    strcpy(buffer, "B\t");
    strcat(buffer, cutoff_string.c_str());
    strcat(buffer, "\t");
    strcat(buffer, missing_string.c_str());
    strcat(buffer, "\t");
    if (orthos_only==true) strcat(buffer, "Y");
    else strcat(buffer, "N");
    
    if (tree_type_string == "Phylip")
        strcat(buffer, "\tP");
    else
        strcat(buffer, "\tN");
    
    if (subgenome == "All")
        strcat(buffer, "\tA");
    else {
        if (subgenome == "Least")
            strcat(buffer, "\tL");
        else {
            if (subgenome == "Intermed")
                strcat(buffer, "\tI");
            else
                strcat(buffer, "\tM");
        }
    }
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
                ss>>dummy>>data_size>>num_pillars;
                ss.clear();
                have_size=true;
                //std::cout<<"Expecting "<<data_size<<" in data stream. There are "<<num_pillars<<" pillars"<<std::endl<<std::flush;
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
            case BATCHBLOCK:
                //  std::cout<<" More batch data: "<<data_cnt<<std::endl;
                i=1;
                error=false;
                while((i<BUFFER_SIZE) && (data_cnt < data_size)) {
                    //if ((int)recv_buffer[i] > 0)
                    data_ss<<recv_buffer[i];
                    //std::cout<<"RB: "<<i<<" - "<<(int)recv_buffer[i]<<std::endl;
                    if ((int)recv_buffer[i]<=0) cnt_neg++;
                    data_cnt++;
                    i++;
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
                if (error_msg.length()>8)
                    error_msg=error_msg.substr(8,error_msg.length()-8);
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
    
    //std::cout<<"Out of message receive block"<<std::endl;
    if (error == false) {
            boost::filesystem::path p(tree_header);
            
            boost::filesystem::create_directories(p);
        
    
        pillar_cnt=0;
        tree_string="";
        
        //tree_header="/media/ramdisk/";
        //CDS_header="/data/Polyploidy_analyses/Legume_WGD/Pillars/AllPillars/";
        
        tree_strings=new std::string[num_pillars];
        full_tree_strings=new std::string[num_pillars];
        CDS_strings=new std::string[num_pillars];
        
        data_ss>>pillar;
        dummy=(char)data_ss.get();
        //std::cout<<"First pillar is "<<pillar<<std::endl;
        while(!(data_ss.eof())) {
            dummy=(char)data_ss.get();
            //std::cout<<"Dummy: "<<dummy<<" = "<<dummy+0<<std::endl;
            if (!(data_ss.eof())) {
                if (dummy == '\0') {
                    data_ss>>new_pillar;
                    dummy=(char)data_ss.get();
                    std::ostringstream t;
                    t<<tree_header<<"Pillar"<<pillar<<".tre";
                    tree_file =t.str();
                    
                    //std::cout<<"Writing "<<tree_string<<" to "<<tree_file<<std::endl;
                    full_tree_strings[pillar_cnt]=tree_file;
                    fout.open(tree_file.c_str());
                    fout<<tree_string;
                    fout.close();
                    std::ostringstream v;
                    v<<"Pillar"<<pillar<<".tre";
                    tree_strings[pillar_cnt]=v.str();
                    std::ostringstream u;
                    u<<"Pillar"<<pillar<<"_CDS.fas";
                    CDS_file = u.str();
                    CDS_strings[pillar_cnt]=CDS_file;
                    tree_string = "";
                    pillar_cnt++;
                    pillar=new_pillar;
                }
                else {
                    tree_string=tree_string+dummy;
                    //std::cout<<tree_string<<std::endl;
                }
            }
        }
        
        

        tar_files(num_pillars, tree_strings, CDS_strings, tree_header, CDS_header, tar_ss);
        
   
        
        
        for(i=0; i<num_pillars; i++) remove(full_tree_strings[i].c_str());
        
        
        delete[] full_tree_strings;
        delete[] tree_strings;
        delete[] CDS_strings;
    }
    close(sock);
    return;
}



bool ping_image_server(string my_socketname, string &error_msg)
{
    int i, sock, rval;
    bool retval=false;
    struct sockaddr_un server;
    char buffer[BUFFER_SIZE], recv_buffer[BUFFER_SIZE];
    MSG_TYPE msg_type;
    
    
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
     
        if (msg_type ==PING) retval=true;
        
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

int tar_files(int num_pillars, std::string *treefiles, std::string *CDS_files,  std::string tree_header, std::string CDS_header, std::stringstream* tar_outstream)
{
    pid_t pid;
    int i, rv, readval, cnt=0, max_arg=0, tree_end;
    int    commpipe1[2], commpipe2[2];
    char dummy[1],  **para_list;
    std::string retval="ERROR", passdata;
    //std::stringstream image_outstream;
    
    max_arg=tree_header.length();
    if (CDS_header.length()>max_arg) max_arg=CDS_header.length();
    
    for(i=0; i<num_pillars; i++) {
        if (treefiles[i].length() > max_arg) max_arg=treefiles[i].length();
        if (CDS_files[i].length()> max_arg) max_arg=CDS_files[i].length();
    }
    
    
    //std::cout<<"Length of max is "<<max_arg<<std::endl;
    
    para_list=new char*[2*num_pillars+8];
    
    para_list[0]=new char[max_arg+1];
    para_list[1]=new char[max_arg+1];
    para_list[2]=new char[max_arg+1];
    para_list[3]=new char[max_arg+1];
    para_list[4]=new char[max_arg+1];
    
    strcpy(para_list[0], "tar");
    strcpy(para_list[1], "-cvf");
    strcpy(para_list[2], "-");
    strcpy(para_list[3], "-C");
    strcpy(para_list[4], tree_header.c_str());
    
    for(i=0; i<num_pillars; i++) {
        para_list[i+5]=new char[max_arg+1];
        //para_list[2*i+5]=new char[max_arg+1];
        strcpy(para_list[i+5], treefiles[i].c_str());
        //strcpy(para_list[2*i+5], CDS_files[i].c_str());
        
        //std::cout<<"For i: "<<i<<" at "<<2*i+3<<" added "<<para_list[2*i+3]<<", "<<para_list[2*i+4]<<std::endl;
    }
    tree_end=5+num_pillars;
    
    para_list[tree_end]=new char[max_arg+1];
    para_list[tree_end+1]=new char[max_arg+1];
    
    strcpy(para_list[tree_end], "-C");
    strcpy(para_list[tree_end+1], CDS_header.c_str());
    
    
    
    
    for(i=0; i<num_pillars; i++) {
        para_list[i+tree_end+2]=new char[max_arg+1];
        //para_list[2*i+5]=new char[max_arg+1];
        strcpy(para_list[i+tree_end+2], CDS_files[i].c_str());
        //strcpy(para_list[2*i+5], CDS_files[i].c_str());
        
        //std::cout<<"For i: "<<i<<" at "<<2*i+3<<" added "<<para_list[2*i+3]<<", "<<para_list[2*i+4]<<std::endl;
    }
    
    
    para_list[2*num_pillars+7]=new char;
    
    para_list[2*num_pillars+7]=NULL;
    
    //for(i=0; i<2*num_pillars+8; i++) std::cout<<" "<<para_list[i];
    //std::cout<<std::endl;
    tar_outstream->clear();
    //std::cout<<"Initial Length of cleaned stream is "<<tar_outstream.str().length()<<std::endl<<std::flush;
    
    
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
        
        
        
        close(commpipe1[1]);
        //wait();
        readval = read(commpipe2[0], dummy, 1);
        //readval = read(commpipe2[0], buf, MAXLINE);
        //cnt++;
#
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
        //image_outstream<<dummy;
        while(readval >0) {
            *tar_outstream<<dummy[0];
            //image_outstream<<buf;
            //std::cout<<"D "<<cnt<<": "<<dummy<<std::endl;
            cnt++;
            readval=read(commpipe2[0], dummy, 1);
            //readval = read(commpipe2[0], buf, MAXLINE);
        }
        
        
        
        close(commpipe2[0]);
        
        //std::cout << "OUTPUT of PROGRAM B is: " << line;
        //std::cout<<"Put "<<cnt<<" items in stream\n"<<std::flush;
        //std::cout<<"New Length of cleaned stream is "<<tar_outstream.str().length()<<std::endl;
        
        //for(i=0; i<2*num_pillars+3; i++) delete[] para_list[i];
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
        
        if ( execv("/bin/tar", para_list) < 0 ) {
            std::cerr<<"execv returned! errno is "<<errno<<std::endl;
            perror("The error message is :");
            std::cerr << "system error" << std::endl;
            return(-1);
        }
        
        return(-1);
    }
    return(-1);
    
    
    
}



