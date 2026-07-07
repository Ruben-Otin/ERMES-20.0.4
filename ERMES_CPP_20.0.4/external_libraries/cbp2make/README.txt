************************************************************************************
* ERMES 20.0 - Linux compilation instructions
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - July 2024
************************************************************************************

To compile ERMES 20.0 for Linux, follow these steps:

1) Export the Visual Studio solution "ERMES.sln" to a CodeBlocks project "ERMES.cbp" using CodeBlocks GUI.  

2) Say NO to default options and then select "GCC compiler" and "Release_x64" options.

3) Copy&Paste "ERMES.cbp" and "ERMES.layout" files from "/windows" folder to "/unix" folder. 

4) Execute in Windows the batch file "cbp2makefile.bat" located in "/external_libraries/cbp2make" with a double click. 

5) Check "Make_Template" and correct CodeBlocks generated "Makefile" compiler options.

6) It may be necessary to load/unload some modules before compiling on Linux, for instance:

   >> module unload pgi
   >> module load gcc/7.3.0
   >> module load openmpi/3.1.1
   >> module load pgi/18.1

   or:

   >> module unload ifort 
   >> module unload pgi
   >> module load gcc/7.3.0
   >> module load openmpi
   >> module load pgi 
    
7) Compile ERMES on Linux using the "make" command inside the "/unix" folder.

8) Rename the executable "ERMES" with the current version name "ERMES_20.0_c7.4". 

9) Copy&Paste the compiled Linux executable "ERMES_20.0_c7.4" into the GiD problemtype folder "/ERMES_20.0/ERMES.gid/execs/unix".

10) Edit "ERMES.unix.bat" in the problemtype folder "/ERMES_20.0/ERMES.gid" with the proper executable name and path:

    >> "$3/execs/unix/ERMES_20.0_c7.4" "$1" > "$2/$1.info"
    
11) After compilation, it may be necessary to change Linux permissions on the executable files:

    >> chmod a+rwx ERMES.unix.bat
    >> chmod a+rwx ERMES
   
12) If using GiD, remember to load the required Linux modules before executing "./gid &".