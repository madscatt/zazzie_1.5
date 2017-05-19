<META http-equiv="Content-Style-Type" content="text/css">
<LINK href="../sassie_style.css" rel="stylesheet" type="text/css">

####[Return to Main Documents Page](../../sassie_docs.html)  


##Docking

Generates structures by performing molecular docking of input structure using the open-source molecular docking package [RosettaDock](http://https://www.rosettacommons.org/docs/latest/application_documentation/docking/docking-protocol).

---

###Accessibility

The Docking module is accessible from the [Alpha](../../alpha/alpha.html) section of the main menu.

The Docking module is accessible from the [Simulate](../simulate.html) section of the main menu.


---

###Basic Usage
The purpose of the module is to perform molecular docking for a protein complex, with  one protein partner fixed and a second mobile. The Docking stage generates a set of decoy structures which are used to generate SAS profiles.  The SAS profiles are filtered using a interpolated data file to generate reduced Chi-Square values which along with docking scores can be used to select a subset of decoy structures.

---

###Notes

* The starting structure must be a complete structure without missing residues.  Atom and residue naming must be compatable with those defined in the CHARMM force field See [Notes on Starting Structures and Force Fields](../../structures_and_force_fields/notes.html) and [PDB Scan](../../build/pdbscan/pdbscan.html) for further details.


* The output file format is DCD since in most cases many structures are generated.  There is no option to save the output files in PDB format.  One can use [Extract Utilities](../../tools/extract_utilities/extract_utilities.html) to convert DCD files to multi-frame PDB files.
* Either Docking or Partition may be selected (not both).

* The number of q-points, range of q, and the spacing of the q-points used to create interpolated data files MUST match the input settings that you use in this module to calculate SAS profiles for comparison of theoretical and experimental data (Chi-Square Filter).
* Docking proceeds in cycles, each cycle attempts to generate an additional 10,000 decoy structures.  After each cycle a test for spacial convergence is performed. When the average difference between successive spacial voxels is less than **spatial voxel convergence tolerance**, no additional docking cycles are performed.  Side chains are then restored and SAS profiles are generated.  The default setting for the **maximum  number of docking cycles** is 20.  If the input value for this variable is changed, then the program will NOT check for spatial convergence, and, instead, will perform the specified number of docking cycles, generating total number of LR decoy structures equal to approximately_ 10,000 x **maximum  number of docking cycles** (usually approximate since 10,000/**number of processors** is truncated to integer)

* The structures in the output DCD and the corresponding reference PDB file produced in the docking stage, are to be used as input DCD and reference PDB files for SAS-Partition.  In addition,  file named "x2_scores.txt", produced in SAS-Docking stage, is required for the subsequent SAS-Partition.



---

#SAS-Docking

###Screen Shots and Description of Input Fields

---

This example generates a set of docking decoy structures using the RosettaDock Full-Protocol.  This involves a global search in space using a low resolution stage of docking scoring in which one protein partner is randomly re-positioned around the other partner. The re-positioning involves a random reorientation of both partners, and the output structure is then transformed to the reference frame of the fixed partner.  Side chains are removed and replaced with single atoms for efficiency.  This is followed by a high resolution stage with side-chains restored.  The selected protein is a tissue factor/FAB complex (PDB ID: 1AHW), in which protein chains A and B are fixed, and chain C is globally re-postioned.

->![inline image](images/docking_LR_input.png)<-


---



* **run name:** user defined name of folder that will contain the results.  

* **reference pdb:** PDB file with naming information and coordinates of the starting structure.  

* **number of processors:** Number of processors to assign to Docking.
* **maximum number of docking cycles:** Number of docking cycles to perform, each cycle generates approximately 10,000 docking decoys.  If set to default (20) program will check for spacial convergence after each cycle, and not do additional cycles if convergence criteria is met.  For any other setting, it will do the specified number of cycles, without check for convergence.
* **partner1_partner2:**  Protein chain designation(s) for partner 1 (fixed), and mobile partner 2, _seperated by the symbol_ _.
* **SAS option:** Scattering choice for use with SASCALC called internally. Must be *either* _neutron_ or _xray_.
* **spatial voxel convergence tolerance:** Convergence criteria for spacial voxels corresponding to each accumulated cycle of docking structures.  If average difference between successive voxels over last 2000 voxels is less than this value, program will exit docking cycles, unless "maximum number of docking cycles" is not set to 20.
* **interpolated data file:** Name of input file with interpolated data with at least three columns: q, I(q), and error in I(q)
* **number of q values:** The number of individual q-points for which I(q) is calculated, including q=0.
* **maximum q value:** The maximum value to calculate I(q).
* **I(0):** Experimentally determined value of scattering intensity at q = 0.
   



---

###Example Output

-->![inline image](images/docking_LR_output.png)<-

The output will show a plot of Chi-Square values versus docking score, for use in selecting cutoffs for the SAS-Partition stage.

Results are written to a new directory within the given "run name" as noted in the output.  

Several files are generated and saved to the "run name" docking directory.  The output DCD file containing decoy structures in "centroid" form,  a reference PDB file in "centroid" form, the file with Chi-Square values and LR docking scores, and a file with cumulative spatial voxels for the LR decoy set.




*./run_ 0/docking/run_ 0_dock.dcd*

*./run_ 0/docking/run_ 0_ dock.pdb*

*./run_ 0/docking/run_ 0_ x2_ vs_ scores.txt*

*./run_ 0//run_ 0_ occupied_ voxels.npy*
  



---

####Visualization

None

---
###Files Used and Created in Example

* input files

	[1ahw.pdb](files/1ahw.pdb)  
	[1ahw_iq.dat](files/1ahw_iq.dat)

* output files

	** caution: DCD file is > 1 GB ** 


	[run_ 0_ dock.dcd](files/run_0_dock.dcd)        
	[run_ 0_ dock.pdb](files/run_0_dock.pdb)  
	[run_ 0_ x2_ vs_ scores_ vs_ x2.txt](files/run_0_scores_vs_x2.txt)  
	[run_ 0_ occupied_ voxels.npy](files/run_0_occupied_voxels.npy)  

---


##SAS-Partition

###Screen Shots and Description of Input Fields

---

 The program will partition the original input set of structures into a subset based favorable Chi-Square and docking score from an input file and using provided input cutoff values. Required input are structures in a DCD file, generated from a SAS-Docking run and a associated reference pdb file. 

---

-->![inline image](images/docking_hr_input.png)<-


---
* **run name:** user defined name of folder that will contain the results.  

* **reference PDB:** PDB file with naming information and coordinates of the full-atom structure

* **trajectory file filename:** DCD file with structures in full-atom mode. 

* **number of processors:** Number of processors to assign to HR Docking.
* **partner1_partner2:**  Protein chain designation(s) for partner 1 (fixed), and mobile partner 2, _seperated by the symbol_ _.
 
* **scores chi-square file:**  file with chi-square values and docking scores in columns one and two, respectively.
* **scores cutoff:** minimum value of absolute value docking score to use to partition the input set of LR decoy structures into smaller subset.
* **chi-sqaure cutoff:** maximum value of chi-square to use to partition the input set of decoy structures in smaller  subset.

---
Several files are generated and saved to the "run name" partition directory.  The output DCD file containing partitioned set of docking structures, and a file with corresponding docking scores and chi-square values



*./run_ 0/partition/run_ 0.dcd*

*./run_ 0/partition/x2_scores _partition.txt*

---

####Visualization

None

---
###Files Used and Created in Example

* input files

	[run_ 0_dock.pdb](files/run_0_dock.pdb)  
	[run_ 0_dock.dcd](files/run_0_dock.dcd)   
	[run_ 0_ scores_ vs _x2.txt](files/run_0_scores_vs_x2.txt)


* output files

	** caution: DCD file is > 6 GB ** 


	[run_ 0.dcd](files/run_0.dcd)        
	[x2_ scores_partition.txt](files/x2_scores_partition.txt)  

---
###Reference(s) and Citations


1. [Application of Small-Angle Scattering Data to Improve Protein-Protein Docking Success Rate Using Rosetta](http://scripts.iucr.org/cgi-bin/paper?a12999) J.A. Snyder, P. Ghatty,  A. McAuley, J. E. Curtis [BIBTeX](references/bibtex_snyder_2016.html), [EndNote](references/endnote_snyder_2016.html),  [Plain Text](references/plain_text_snyder_2016.html)  


5. [SASSIE: A program to study intrinsically disordered biological molecules and macromolecular ensembles using experimental scattering restraints](http://www.sciencedirect.com/science/article/pii/S0010465511003286)  J. E. Curtis, S. Raghunandan, H. Nanda, S. Krueger, Comp. Phys. Comm. 183, 382-389 (2012). [BIBTeX](../../references/bibtex_sassie.html), [EndNote](../../references/endnote_sassie.html),  [Plain Text](../../references/plain_text_sassie.html)  

---

####[Return to Simulate](../simulate.html)
####[Return to Main Documents Page](../../sassie_docs.html)  


<a href=#>Go to top</a>

<footer>
  <ul>
  Supported via CCP-SAS a joint EPSRC (EP/K039121/1) and NSF (CHE-1265821) grant
  </ul>
</footer> 
