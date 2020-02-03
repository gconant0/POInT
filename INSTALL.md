# POInT installation instructions:

### Basic Instructions:
The following 4 commands should generate POInT  from POInT.tar-- <br>
the % indicates your prompt--don't type it :)


### Type
`% tar xvf POInT.tar`

Now, you will need to make the package.  

### \*\*\*\*\*BUILDING POInT\*\*\*\*\*
The previous command will create a new directory called POInT. Enter that directory<br>
`% cd POInT`

Now the configure script will make some simple checks of your OS and build the files needed for compiling POInT<br>
`%./configure.pl`

You can now make all of the needed programs by typing<br>
`% make`



### Options from the configure script:
- If you have install the Gnu Plotutils package and would like to use it for illustrations with POInT, use<br>
`%./configure.pl -p:<PATH TO libplot>`

- If you want to compile the OpenMP parallel version of POInT, use<br>
`%./configure.pl -omp`


Please see the included [POInT documentation](http://conantlab.org/POInT/POInT.html) for instructions and tips on using POInT<br>


### Comments on basic installation process:

To build and use the package, you will need a c and a c++ compiler (preferably the GNu gcc and g++). The make command will build liblapack (a partial lapack library), libf2c (because lapack was originally built in Fortran) and POInT. 
  
The configure.pl script will check the standard locations to see if your system already has the lapack and blas linear algebra libraries installed (/usr/lib/, /usr/local/lib and /lib).  If the libraries are detected, they will be used for the installation: otherwise the script will invoke the compilation of the included copy of the necessary parts of lapack and the f2c routines needed to use it.


### Other platforms:
This package was developed and tested using Mac OX 10.13 and Ubuntu Linux with the clang and gcc compilers.  It should work on 
other platforms, but has not been tested on them.  

