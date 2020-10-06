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


### Examples of usage

```
usage: GP4 [-h] -i INPUT [-o OUTPUT] [-n NAME] [-t THREADS]
           [-r {data,pred,all}] [-v]

Script to perform a proteomic based consensus prediction of protein
localization in Gram+. Written by Stefano Grasso (c).

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file. If you chose -r "pred" then indicate the
                        folder with the generated data.
  -o OUTPUT, --output OUTPUT
                        Output folder. Default: Results
  -n NAME, --name NAME  Output file name. Default: same as input name. If
                        already existing a short string will be attached to
                        it.
  -t THREADS, --threads THREADS
                        Number of threads to be used. Default: all.
  -r {data,pred,all}, --run {data,pred,all}
                        'data' to generate data; 'pred' requires generated
                        data as input; 'all' (default).
  -v, --version         show program's version number and exit
```

The most basic command is: 

	./GP4.py -i Fasta_file

Results will be saved by default in 'Results' within the 'Fasta_file' folder. If already existing a string will be attached in order to avoid overwriting. By default both data are generated and predictions performed.

To specify an output folder and a new name run:

	./GP4.py -i Fasta_file -o New-results-folder -n New-name

Results will be saved in New-results-folder/New-name.

You can also decide to run only the module to generate the data (for instance if you have a big dataset):

	./GP4.py -i Fasta_file -o New-results-folder -n New-name -r data

Then when you are ready to generate the final predictions you can run:

	./GP4.py -i New-results-folder/New-name -r pred

*Note*: now the input is not the Fasta_file anymore, but the folder containing the data.
By default the mode is 'all' so data and prediction will be generated in a single run.

You can also select the number of threads to generate the data in parallel (this affect only the 'data' module):

	./GP4.py -i Fasta_file -o New-results-folder -n New-name -t 8

*Note*: by deafult all available threads will be used. 


### Cite us!

If you use our tool in your paper please cite: 

### Troubleshooting

If you run in any issue specific with GP4, please report it. We'll try to help you! 
In the meanwhile you can run your predictions at: http://gp4.hpc.rug.nl/
