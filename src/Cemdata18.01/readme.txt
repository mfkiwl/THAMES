To use the Cemdata2018 data base in GEMS-PSI package, please download it to your hard disk and perform the following steps:

    
Unzip the downloaded zip file (contains a directory named "cemdata18") into a temporary directory, e.g. as /Tempfiles/cemdata
 
    

Find where you have GEMS installed (on Windows, usually under C:\GEMS335) the Gems3-app\Resources\DB.default directory.  Under Linux, this may need a root password.
 
    
close GEMS. Copy all files from /Tempfiles/cemdata into the ems3-app\Resources\DB.default diractory.

    

Start GEMS and create a new project. In the "Selection of Independent
    Components..." dialog, keep the "psi-nagra" thermodynamic database
and turn on
 "3rd party". Tick on "cemdata"; all the data within cemdata are selected , 
depending on whether you calculate 1) Portland and blended cements or 2) alkali activated cements select a different subset

keep "pc", ".", "ss" and "ss-fe" (and deselect "aam") if you want to calculate Portland and blended cements; 
choose within "pc":"csh" one of the four alternative CSH model "cshq", "cshkn", "csh3t" or "csh2o") 

for alkali activated systems, select "aam" and deselect "pc". 


This will add the cemdata database
 to the default Nagra-PSI database. 
   
 
Select Independent Components to form the system and click "Ok" to proceed
 as usual.

    


version 1.1: 8-1-2019
definition and data for INFCNA corrected to (CaO)1(SiO2)1.1875(Al2O3)0.15625(Na2O)0.34375(H2O)1.3125 as given in Table 4