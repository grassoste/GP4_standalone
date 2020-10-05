# GP4 standalone
## Gram-Positive Protein Prediction Pipeline

Standalone version of the prediction tool presented in the publication: 


### Requisites

To run GP4 locally you first need to download and install the following programs on your local machine:
	
	* SignalP 4.1 (https://services.healthtech.dtu.dk/software.php)
	* SignalP 5.0 (https://services.healthtech.dtu.dk/software.php)
	* TatP (https://services.healthtech.dtu.dk/software.php)
	* LipoP (https://services.healthtech.dtu.dk/software.php)
	* Phobius (http://phobius.sbc.su.se/data.html)
	* TMHMM (https://services.healthtech.dtu.dk/software.php)
	* InterproScan 5.33 (InterproDB 72.0): ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.33-72.0/

Additionally you need to have a working version of Python (>=3.5) and Biopython (>=1.74). 
Please also make sure you have the following packages: pandas, numpy, decouple.

### Installation instructions

First, clone this repository in a folder of your choice:

	git clone https://github.com/grassoste/GP4_standalone.git

Then modify the file conf.env adding the absolute path to the required programs.

To see the help with the possible options just type:

	./GP4.py -h

We suggest to first test the installation by running a fasta file with a single protein, simply doing:

	./GP4.py -i Fasta_file


### Cite us!

If you use our tool in your paper please cite: 

### Troubleshooting

If you run in any issue specific with GP4, please report it. We'll try to help you! 
In the meanwhile you can run your predictions at: http://gp4.hpc.rug.nl/
