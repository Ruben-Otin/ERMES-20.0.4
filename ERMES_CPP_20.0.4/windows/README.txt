************************************************************************************
* ERMES 20.0 - Windows compilation instructions
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - July 2024
************************************************************************************

To compile ERMES 20.0 for Windows, follow these steps: 

1) Open the file "ERMES.sln" located in "\ERMES_CPP_20.0\windows" with Microsoft Visual Studio 2022.

2) Click the "Build" button on Microsoft Visual Studio 2022.

3) The compiled executable "ERMES.exe" will be located in the folder: "\windows\x64\Release".

4) Rename "ERMES.exe" with the current version name "ERMES_20.0.exe". 

5) Copy&Paste "ERMES_20.0.exe" on the GiD problemtype folder "\ERMES_20.0\ERMES.gid\execs\win".

6) Edit "ERMES.win.bat" on "\ERMES_20.0\ERMES.gid" to call "ERMES_20.0.exe" from GiD:

   >> %3\execs\win\ERMES_20.0.exe %1 > %2\%1.info

7) To save disk space, delete after compilation the following folders:

   - "\windows\x64"
   - "\windows\Release"
   - "\windows\.vs\ERMES\FileContentIndex"
   - "\windows\.vs\ERMES\v17\ipch"
  
8) To compile ERMES for Linux, read the "README.txt" file located in the "\unix" folder.