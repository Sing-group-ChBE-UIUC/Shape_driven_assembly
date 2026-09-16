## File names
- Each files named by shape contains raw data files for structure analysis (g(r), S(q) and conformation) and simulation set-up (LAMMPS).
- Names are not organized (Too lazy to edit all files), but file name contains "RDF*" is g(r) file with specific pair indexes and "Sq*" is structure factor files. The index in the name of the structure factor files indicate the nth trials to get the structure factor.
- "input" is the lammps input file, "PMF.txt" or "pot0.table" is lammps pair potential file.
- "D_cfg*" or "D_list.txt" files are the list of diameter of beads in a chain where the third column is Kuhn length conversion. (didn't use it at all though)
- Structure factors are generated first iterating three q vectors, then sorted into absolute q bins. "Sq_raw.txt" or "Sq_*" files indicate the raw data, and "Sq_by_bin.txt" or "S_q*_bin" indicates sorted data.


This data is related to the paper "Local structures of Shape-Defined Bottlebrush Homopolymers in Melts" submitted to Soft Matter (2026.9 current). Hope it is published. I'll update the DOI later. 

  If you find any issues or need further details, contact through my personal email (khs12137@gmail.com)
   

