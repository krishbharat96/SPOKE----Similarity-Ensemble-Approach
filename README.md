# SEA Data within SPOKE

In an effort to integrate the latest assay information from ChEMBL in order to utilize Dr. Michael Keiser's latest pipeline for making predictions on compound-protein binding data, we developed a way to 
In order to extract assay data from SPOKE, we used parts of ChEMBL db to obtain compound to protein assay data. These data include different types of assays such as IC50, EC50, logEC50 etc. along with filtration parameters to obtain the best curated data possible. There are two main files to look at:
  1. parse_sea.py: This is a module that utilizes the latest copy of the sqlite ChEMBL db file (which can be found at: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/) to extract the appropriate SEA assay information
  2. Update_SEA_Assays.py: Utilizes info from the previous module to extract relevant SEA information and update the latest SPOKE db (you have to input SPOKE address etc.)
  3. The extracted information is then run using the SEA algorithm which can be extracted from Dr. Michael Keiser's Lab.
