# SPOT-RNA2
Improved RNA Secondary Secondary Structure Prediction using Evolutionary Profile, Mutational Coupling and Two-dimensional Transfer Learning.

SYSTEM REQUIREMENTS
====
Hardware Requirments:
----
SPOT-RNA2 predictor requires only a standard computer with around 32 GB RAM to support the in-memory operations for RNAs sequence length less than 500.

Software Requirments:
----
* [Python3.6](https://docs.python-guide.org/starting/install3/linux/)
* [Perl-5.4 or later](https://www.perl.org/get.html)
* [virtualenv](https://virtualenv.pypa.io/en/latest/installation/) or [Anaconda](https://anaconda.org/anaconda/virtualenv)
* [CUDA 10.0](https://developer.nvidia.com/cuda-10.0-download-archive) (Optional If using GPU)
* [cuDNN (>= 7.4.1)](https://developer.nvidia.com/cudnn) (Optional If using GPU)

SPOT-RNA2 has been tested on Ubuntu 14.04, 16.04, and 18.04 operating systems.

USAGE
====

Installation:
----

To install SPOT-RNA2 and it's dependencies following commands can be used in terminal:

1. `git clone https://github.com/jaswindersingh2/SPOT-RNA2.git`
2. `cd SPOT-RNA2`

Either follow **virtualenv** column steps or **conda** column steps to create virtual environment and to install SPOT-RNA2 python dependencies given in table below:<br />

|  | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; virtualenv | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; conda |
| :- | :-------- | :--- |
| 3. | `virtualenv -p python3.6 venv` | `conda create -n venv python=3.6` |
| 4. | `source ./venv/bin/activate` | `conda activate venv` | 
| 5. | `pip install -r requirements.txt && deactivate` | `while read p; do conda install --yes $p; done < requirements.txt && conda deactivate` | 

If Infernal tool is already installed in the system, please add path to the folder contains binary files in line no. 12 of `run_spotrna2.sh` file. In case, Infernal tool is not installed in the system, please use follwing two commands to download and extract it. In case of any problem regarding Infernal download, please refer to [Infernal webpage](http://eddylab.org/infernal/) as following commands only tested on Ubuntu 18.04, 64 bit system.

6. `wget 'eddylab.org/infernal/infernal-1.1.3-linux-intel-gcc.tar.gz'`
7. `tar -xvzf infernal-*.tar.gz && rm infernal-*.tar.gz`

If BLASTN tool is already installed in the system, please add path to the folder contains binary files in line no. 10 of `run_spotrna2.sh` file. In case, BLASTN tool is not installed in the system, please use follwing two commands to download and extract it. In case of any problem regarding BLASTN download, please refer to [BLASTN webpage](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlasstDocs&DOC_TYPE=Download) as following commands only tested on Ubuntu 18.04, 64 bit system.

8. `wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*+-x64-linux.tar.gz'`
9. `tar -xvzf ncbi-blast-*+-x64-linux.tar.gz && rm ncbi-blast-*+-x64-linux.tar.gz`

To install **SPOT-RNA** predictor to obtain consensus secondary structure for the MSA-1 from BLASTN, the following three command can be used.<br />

10. `git clone https://github.com/jaswindersingh2/SPOT-RNA.git && cd SPOT-RNA`

11. `wget 'https://www.dropbox.com/s/dsrcf460nbjqpxa/SPOT-RNA-models.tar.gz' || wget -O SPOT-RNA-models.tar.gz 'https://app.nihaocloud.com/f/fbf3315a91d542c0bdc2/?dl=1'`

12. `tar -xvzf SPOT-RNA-models.tar.gz && rm SPOT-RNA-models.tar.gz && cd ../`

If NCBI's nt database already available in your system, please add path to database in line no. 11 and line 13 of `run_spotrna.sh` file.  Otherwise, download the reference database ([NCBI's nt database](ftp://ftp.ncbi.nlm.nih.gov/blast/db/)) for BLASTN and INFERNAL. The following command can used for NCBI's nt database download. Make sure there is enough space on the system as database is of size around 270 GB after extraction and it can take couple of hours to download depending on the internet bandwidth. In case of any problem, please refer to [NCBI's database website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

13. `wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" -O ./nt_database/nt.gz && gunzip ./nt_database/nt.gz`

The database need to formated to use with BLASTN tool. To format it, the following command can be used. Please make sure system have enough space as formated database is of size around 120 GB in addition to appox. 270 GB from previous step and it can few hours for it.

14. `./ncbi-blast-2.10.0+/bin/makeblastdb -in ./nt_database/nt -dbtype nucl`

To install the DCA predictor, the following 2 command can be used:<br />

15. `git clone "https://github.com/sokrypton/GREMLIN_CPP"`

16. `cd GREMLIN_CPP && g++ -O3 -std=c++0x -o gremlin_cpp gremlin_cpp.cpp -fopenmp && cd ../`

To install the LinearPartition, the following 2 command can be used:<br />

17. `git clone 'https://github.com/LinearFold/LinearPartition.git'`

18. `cd LinearPartition/ && make && cd ../`

To run the SPOT-RNA2
-----

```
./run_spotrna.sh sample_run/6ufj.fasta 
```

The above command creates two folder `6ufj_features` and `6ufj_outputs` under the folder contains input file (`sample_run/` folder in this case). `6ufj_features/` contains all the alignments (MSA-1, MSA-2) and features (PSSM, DCA, bps probability) generated from SPOT-RNA2 pipeline. `6ufj_outputs/` contains predicted secondary structure in bpseq format (`6ufj.bpseq`), ct format (`6ufj.ct`), dbn format (`6ufj.st`) with secondary structure motifs, and base-pair probability (`6ufj.prob`). The verify the results, it should be as in `sample_seq_features` and `sample_seq_outputs` folder because both sequences (`sample_seq.fasta` and `6ufj.fasta`) are same.

References
====

**If you use SPOT-RNA2 for your research please cite the following papers:**

Singh, J., Paliwal, K., Zhang, T., Singh, J., Litfin, T., Zhou, Y., 2020. Improved RNA secondary structure prediction using evolutionary profile, mutational coupling and two-dimensional transfer learning. (submitting soon)

**If you use SPOT-RNA2 data sets and/or input feature pipeline, please consider citing the following papers:**

[1] Zhang, T., Singh, J., Litfin, T., Zhan, J., Paliwal, K. and Zhou, Y., 2020. RNAcmap: A Fully Automatic Method for Predicting Contact Maps of RNAs by Evolutionary Coupling Analysis. bioRxiv.

[2] Zhang, H., Zhang, L., Mathews, D.H. and Huang, L., 2020. LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities. Bioinformatics, 36(Supplement_1), pp.i258-i267.

[3] Singh, J., Hanson, J., Paliwal, K. and Zhou, Y., 2019. RNA secondary structure prediction using an ensemble of two-dimensional deep neural networks and transfer learning. Nature communications, 10(1), pp.1-13.

[4] Nawrocki, E.P. and Eddy, S.R., 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29(22), pp.2933-2935.

[5] H.M. Berman, J. Westbrook, Z. Feng, G. Gilliland, T.N. Bhat, H. Weissig, I.N. Shindyalov, P.E. Bourne. (2000) The Protein Data Bank Nucleic Acids Research, 28: 235-242.

[6] Padideh Danaee, Mason Rouches, Michelle Wiley, Dezhong Deng, Liang Huang, David Hendrix, bpRNA: large-scale automated annotation and analysis of RNA secondary structure, Nucleic Acids Research, Volume 46, Issue 11, 20 June 2018, Pages 5381–5394, https://doi.org/10.1093/nar/gky285

[7] Kamisetty, H., Ovchinnikov, S. and Baker, D., 2013. Assessing the utility of coevolution-based residue–residue contact predictions in a sequence-and structure-rich era. Proceedings of the National Academy of Sciences, 110(39), pp.15674-15679.

[8] Chiu, J.K.H. and Chen, Y.P.P., 2014. Efficient conversion of RNA pseudoknots to knot-free structures using a graphical model. IEEE Transactions on Biomedical Engineering, 62(5), pp.1265-1271.

Licence
====
Mozilla Public License 2.0


Contact
====
jaswinder.singh3@griffithuni.edu.au, yaoqi.zhou@griffith.edu.au
