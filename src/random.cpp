#include "random.h"


Random_Gen::Random_Gen()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    
    seedin.open("seed.txt");
    if (seedin.fail()) {
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        
        strftime(buffer,sizeof(buffer),"%A, %B %d, %Y %H:%M:%S",timeinfo);
        
        std::cout<<"No seed file found: seeding from time phrase "<<buffer<<std::endl;
        
        phrtsd(buffer, &seed1, &seed2);
    }
    else {
        seedin>>seed1>>seed2;
        seedin.close();
    }
    setall(seed1, seed2);
    
}

Random_Gen::~Random_Gen()
{
    seedout.open("./seed.txt");
    getsd(&seed1, &seed2);
    seedout<<seed1<<"\t"<<seed2<<std::endl;
    seedout.close();
}
