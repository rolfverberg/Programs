# Programs

# ConvertFVtoTecplot:
A routine to convert binary Fieldview data into binary or ascii Tecplot format. A linux makefile is provided in the Compile directory.
Build the executable in the Compile directory and execute convertfvtotecplot.
It will first ask whether you want to convert a single Fieldview file or a two seperate files, one with the grid plus one with the solution. Next it will ask you to enter the filename(s). Finally it will ask you to choose binary or ascii Tecplot output.
Alternatively you can provide one or two command line arguments with the filename(s). In the latter case the Tecplot output is binary by default.

# Source:
General source code files:

clComplex.cpp/h: A class to work with complex numbers and arrays and 2D matrices of complex numbers.

clFieldview.cpp/h: A class to work with binary Fieldview data files.

clMatrix.h: A template class to work with 2D matrices of an arbitrary type.

clMatrix3D.h: A template class to work with 3D matrices of an arbitrary type.

clMatrix4D.h: A template class to work with 4D matrices of an arbitrary type.

clMpi.cpp/h: A class with routines to interface between your code and MPI.

clString.cpp/h: A class to work with strings and arrays and 2D matrices of strings.

clTecplot.cpp/h: A class to work with Tecplot data files (binary or ascii).

clUtils.cpp/h: A collection of general utilities.

clVector.h: A template class to work with vectors 2D matrices of an arbitrary type.

# TestTecplot:
Test routine for clTecplot, a class with routines to read or write data in Tecplot ascii or binary format. It also contains routines for some basic operations on the data variables and zones. Build the executable in the Compile directory, copy the executable testtecplot to the base TestTecplot directory and execute it. It will run abunch of tests on files in the Input directory. Output will be written to the Output directory.
