#!/bin/bash

start=`date +%s`

input="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
input_dir=$(dirname $input)
seq_id=$(basename $(basename $input) | cut -d. -f1)
program_dir=$(dirname $(readlink -f $0))

path_blastn=$program_dir/ncbi-blast-2.10.*+/bin       				# set here path to the folder contains executable binary files of Blast package
path_blastn_database=$program_dir/nt_database/nt      				# set here path to the formatted NCBI's database file without extension 
path_infernal=$program_dir/infernal-1.1.3-linux-intel-gcc/binaries  # set here path to the folder contains executable binary files Infernal package
path_infernal_database=$program_dir/nt_database/nt					# set here path to the NCBI's database database file

mkdir -p $input_dir/${seq_id}_features && mkdir -p $input_dir/${seq_id}_outputs
echo ">"$seq_id > $input_dir/${seq_id}_features/$seq_id.fasta
awk -i inplace '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $input 
tail -n1 $input >> $input_dir/${seq_id}_features/$seq_id.fasta

feature_dir=$input_dir/${seq_id}_features
output_dir=$input_dir/${seq_id}_outputs

#exit 1

if [ ! -f $path_blastn_database ];  then
    echo ""
    echo "========================================================================================"
    echo "            Looks like nt database doesn't exists in the path $path_blastn_database.    "
    echo "            If you want to download the database now, please make sure you have enough  "
    echo "            space in mounted directory and internet connection have enough bandwidth as "
    echo "            file is of size 270 GBs after unzip. It may take forever to download if     "
    echo "                                internet is slow!                                       "
    echo "========================================================================================"
    echo ""

    echo -n "Type 'y' for download or any other key to exit: "    
    read userinput

    if [[ $(echo $userinput | tr '[A-Z]' '[a-z]') == 'y' ]]; then

		echo ""
		echo "=============================================================================================="
		echo "       Downloading NCBI's database form ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz link. "
		echo "                                 May take few hours to download.                              "
		echo "=============================================================================================="
		echo ""
		wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" -O $program_dir/nt_database/nt.gz


		if [[ $? -eq 0 ]]; then 
	        echo ""
	        echo "======================================================================="
	        echo "            nt database is completed successfully.                     "
	        echo "======================================================================="
	        echo ""
		else
	        echo ""
	        echo "======================================================================="
	        echo "            Error! Unable to download database sucessfully.            "
	        echo "            Check wget command or internet connection.            "
	        echo "======================================================================="
	        echo ""
	        exit 1        
		fi

		echo ""
		echo "======================================================================"
		echo "            Unziping the downloaded nt database.                      "
		echo "       May take few hours as size of unzipped file is around 270 GBs. "
		echo "======================================================================"
		echo ""
		
	############ unzip the nt data base file ############
		gunzip $program_dir/nt_database/nt.gz

		if [[ $? -eq 0 ]]; then 
	        echo ""
	        echo "======================================================================="
	        echo "            nt database unzip completed successfully.                  "
	        echo "======================================================================="
	        echo ""
		else
	        echo ""
	        echo "======================================================================="
	        echo "            Error! unable to unzip database sucessfully.               "
	        echo "            Please check if gunzip program exists!                     "
	        echo "======================================================================="
	        echo ""
	        exit 1        
		fi

    else
		echo ""
		echo "==========================================================="
		echo "      Existing the program because nt database is missing! "
		echo "==========================================================="
		echo ""
        exit 1
    fi

fi


###### check if aligned homologous sequences file already exists ############
if [ -f $feature_dir/$seq_id.a2m ];	then
        echo ""
        echo "======================================================================"
        echo "    MSA file $feature_dir/$seq_id.a2m from Infernal Pipeline already  "
        echo "    exists for query sequence $feature_dir/$seq_id.fasta.             "
        echo "                                                                      "
        echo "    Delete existing $feature_dir/$seq_id.a2m if want to generate new  "
        echo "    alignment file                                                    "
        echo "======================================================================"
    	echo ""
else

   #### check if formatted nt database exists or not ##### 
    if [[ ! -f "$path_blastn_database.nal" ]]; then
        echo ""
        echo "====================================================================="
        echo "    Nucleotide database file $path_database/nt need to formated      "
        echo "    formated to use with 'makeblastdb' program in BLAST-N program.   "  
        echo ""          
		echo "    Formatting may take 2-3 hours as size of file is around 270 GBs. "
        echo "====================================================================="
        echo ""
        $path_blastn/makeblastdb -in $path_database/nt -dbtype nucl
        
        if [[ $? -eq 0 ]]; then
                echo ""
                echo "======================================================="
                echo "          nt database formatted successfully.          "
                echo "======================================================="
                echo ""
        else
                echo ""
                echo "=================================================================="
                echo "        Error occured while formatting the nt database.           "
                echo ""
                echo "  Check for '$path_blastn/makeblastdb' program in BLAST package   "
                echo "=================================================================="
                echo ""
                exit 1
        fi                      
    fi


    #################### check if blastn alignment file ready exists ######################
    if [ -f $feature_dir/$seq_id.bla ];       then
	    echo ""
	    echo "======================================================================="
	    echo "    MSA-1 file $feature_dir/$seq_id.bla from Infernal Pipeline already "
	    echo "    exists for query sequence $feature_dir/$seq_id.fasta.              "
	    echo "                                                                       "
	    echo "    Delete existing $feature_dir/$seq_id.a2m if want to generate new   "
	    echo "    alignment file.                                                    "
	    echo "======================================================================="
		echo ""
    else
        echo ""
        echo "==========================================================================================================================="
        echo "      Start Running BLASTN for first round of homologous sequence search for query sequence $feature_dir/$seq_id.fasta.    "
        echo "      May take 5 mins to few hours depending on sequence length and no. of homologous sequences in database.               "
        echo "==========================================================================================================================="
        echo ""
        $path_blastn/blastn -db $path_blastn_database -query $feature_dir/$seq_id.fasta -out $feature_dir/$seq_id.bla -evalue 0.001 -num_descriptions 1 -num_threads 8 -line_length 1000 -num_alignments 50000
    fi
			
	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      First round of MSA-1 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "=================================================================="
        echo "        Error occured while formatting the nt database.           "
        echo ""
        echo "  Check for '$path_blastn/makeblastdb' program in BLAST package   "
        echo "=================================================================="
        echo ""
        exit 1
    fi

	######## reformat the output ################
    echo ""
    echo "========================================================================================"
    echo "         Converting $feature_dir/$seq_id.bla from BLASTN to $feature_dir/$seq_id.sto.   "
    echo "========================================================================================"
    echo ""
	$program_dir/utils/parse_blastn_local.pl $feature_dir/$seq_id.bla $feature_dir/$seq_id.fasta $feature_dir/$seq_id.aln
	$program_dir/utils/reformat.pl fas sto $feature_dir/$seq_id.aln $feature_dir/$seq_id.sto


	if [ $? -eq 0 ]; then
	    echo ""
	    echo "=========================================="
        echo "      Converison completed successfully.  "
	    echo "=========================================="
	    echo ""
	else
        echo ""
        echo "============================================================================================="
        echo "   Error occured while Converting $feature_dir/$seq_id.bla to $feature_dir/$seq_id.sto       "
        echo " "
        echo "  Check for $program_dir/utils/parse_blastn_local.pl and $program_dir/utils/reformat.pl file."
        echo "============================================================================================="
        echo ""
        exit 1
    fi

	######## predict secondary structure from SPOT-RNA ################
    echo ""
    echo "==============================================================================================================================="
    echo "       Predicting Consensus Secondary Structure (CSS) of query sequence $feature_dir/$seq_id.fasta using SPOT-RNA predictor.   "
    echo "==============================================================================================================================="
    echo ""
	source $program_dir/venv/bin/activate || conda activate venv
	cd $program_dir/SPOT-RNA
	python3 SPOT-RNA.py --inputs $feature_dir/$seq_id.fasta --outputs $feature_dir	
	cd -

	export PERL5LIB=$program_dir/utils/FreeKnot
	perl $program_dir/utils/FreeKnot/remove_pseudoknot.pl -i bpseq -s bp $feature_dir/$seq_id.bpseq > $feature_dir/$seq_id.bpseq.unknotted
	python3 $program_dir/utils/bpseq2dbn.py --inputs $feature_dir --outputs $feature_dir --rna_id $seq_id
	tail -n +3 $feature_dir/$seq_id.dbn > $feature_dir/$seq_id.db

	deactivate || conda deactivate

	################ reformat ss with according to gaps in reference sequence of .sto file from blastn ################
	for i in `awk '{print $2}' $feature_dir/$seq_id.sto | head -n5 | tail -n1 | grep -b -o - | sed 's/..$//'`; do sed -i "s/./&-/$i" $feature_dir/$seq_id.db; done

	#########  add reformated ss from last step to .sto file of blastn ##############
	head -n -1 $feature_dir/$seq_id.sto > $feature_dir/temp.sto
	echo "#=GC SS_cons                     "`cat $feature_dir/$seq_id.db` > $feature_dir/temp.txt
	cat $feature_dir/temp.sto $feature_dir/temp.txt > $feature_dir/$seq_id.sto
	echo "//" >> $feature_dir/$seq_id.sto

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "=================================================================="
        echo "      Consensus Secondary Structure (CSS) generated successfully. "
	    echo "=================================================================="
	    echo ""
	else
        echo ""
        echo "=============================================================================="
        echo "             Error occured while generating structure from SPOT-RNA.          "
        echo " "
        echo "  Please raise issue at 'https://github.com/jaswindersingh2/SPOT-RNA2/issues'."
        echo "=============================================================================="
        echo ""
        exit 1
    fi

	######## run infernal ################
    echo ""
    echo "=============================================================================================================="
    echo "      Building Covariance Model from BLASTN alignment (with SS from SPOT-RNA) from $feature_dir/$seq_id.sto file.         "
    echo "=============================================================================================================="
    echo ""
	$path_infernal/cmbuild --hand -F $feature_dir/$seq_id.cm $feature_dir/$seq_id.sto

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "============================================================================"
        echo "    Covariance Model (CM) built successfully from $feature_dir/$seq_id.sto. "
	    echo "============================================================================"
	    echo ""
	else
        echo ""
        echo "==============================================================================================="
        echo "     Error occured while building Covariance Model (CM) from $path_infernal/cmbuild.           "
        echo " "
        echo "     Please check for $path_infernal/cmbuild program.      "
        echo "==============================================================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "===================================================================="
    echo "       Calibrating the Covariance Model $feature_dir/$seq_id.cm.    "
    echo "===================================================================="
    echo ""
	$path_infernal/cmcalibrate $feature_dir/$seq_id.cm

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "    CM calibrated $feature_dir/$seq_id.cm successfully.    "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "==============================================================="
        echo "     Error occured while calibrating $feature_dir/$seq_id.cm.  "
        echo " "
        echo "     Please check for $path_infernal/cmcalibrate program.      "
        echo "==============================================================="
        echo ""
        exit 1
    fi

    echo ""
    echo "======================================================================================================================"
    echo "        Second round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm.    "
    echo "                 May take 15 mins to few hours for this step.                                                         "
    echo "======================================================================================================================"
    echo ""
	$path_infernal/cmsearch -o $feature_dir/$seq_id.out -A $feature_dir/$seq_id.msa --cpu 24 --incE 10.0 $feature_dir/$seq_id.cm $path_infernal_database

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "      Second round of MSA-2 search completed successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "===================================================================================="
        echo "     Error occured during the second round search using CM $feature_dir/$seq_id.cm. "
        echo " "
        echo "     Please check for $path_infernal/cmsearch program.                              "
        echo "===================================================================================="
        echo ""
        exit 1
    fi

	######### reformat the alignment without gaps and dashes  ###############
    echo ""
    echo "======================================================================="
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa   "
    echo "          for PSSM and DCA features by removing the gaps and dashes.   "
    echo "======================================================================="
    echo ""
	$path_infernal/esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa > $feature_dir/temp.a2m

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==========================================================="
        echo "   Reformatted the $feature_dir/$seq_id.msa successfully.  "
	    echo "==========================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the refomatting the alignment file $feature_dir/$seq_id.msa.  "
        echo " "
        echo "     Please check for $path_infernal/esl-reformat program.                              "
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s $feature_dir/temp.a2m > $feature_dir/$seq_id.a2m

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================="
        echo "   Duplicate sequences removed successfully.   "
	    echo "==============================================="
	    echo ""
	else
        echo ""
        echo "========================================================================================"
        echo "     Error occured during the removel of duplicates from MSA-2.  "
        echo " "
        echo "     Please check for $program_dir/utils/seqkit program.                              "
        echo "========================================================================================"
        echo ""
        exit 1
    fi

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.a2m | sed '/^$/d' > $feature_dir/temp.a2m 
	############# add query sequence at the top of MSA file  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m > $feature_dir/$seq_id.a2m 

fi

############# check if pssm file already exists otherwise generate from alignment file #############
if [ -f $feature_dir/$seq_id.pssm ];	then
        echo ""
        echo "=============================================================================================================================================="
        echo "    PSSM feature file $feature_dir/$seq_id.pssm already exists for query sequence $feature_dir/$seq_id.fasta.  "
        echo "=============================================================================================================================================="
    	echo ""
else
	echo ""
	echo "======================================================================================"
	echo "          Extracting PSSM features from the alignment $feature_dir/$seq_id.a2m.       "
	echo "======================================================================================"
	echo ""
	$program_dir/utils/getpssm.pl $feature_dir/$seq_id.fasta $feature_dir/$seq_id.a2m $feature_dir/$seq_id.pssm

	if [ $? -eq 0 ]; then
	    echo ""
	    echo "==============================================================="
        echo "   PSSM extracted successfully from $feature_dir/$seq_id.a2m.  "
	    echo "==============================================================="
	    echo ""
	else
        echo ""
        echo "========================================================================="
        echo "     Error occured while extracting PSSM from $feature_dir/$seq_id.a2m.  "
        echo " "
        echo "     Please check for $program_dir/utils/getpssm.pl program.             "
        echo "========================================================================="
        echo ""
        exit 1
    fi
fi

######### run linearpartition RNA secondary structure base-pair probability predictor ###############
echo ""
echo "============================================================================"
echo "          Running LinearPartition-V for base-pair probabilty features.      "
echo "============================================================================"
echo ""
tail -n +2 $feature_dir/$seq_id.fasta | $program_dir/LinearPartition/linearpartition -V -r $feature_dir/$seq_id.prob

if [ $? -eq 0 ]; then
    echo ""
    echo "===================================================================="
    echo "   Base-pair probabilty successfully obtained from LinearPartition. "
    echo "===================================================================="
    echo ""
else
    echo ""
    echo "============================================================================="
    echo "                Error occured while running LinearPartition.  "
    echo " "
    echo "     Please check for $program_dir/LinearPartition/linearpartition program.  "
    echo "============================================================================="
    echo ""
    exit 1
fi

############# check if dca file already exists otherwise generate from alignment file #############
if [ -f $feature_dir/$seq_id.dca ];	then
        echo ""
        echo "==============================================================="
        echo "    GRELMLIN feature file $feature_dir/$seq_id.dca already     "
        echo "    exists for query sequence $feature_dir/$seq_id.fasta.      "
        echo " "
        echo "    Delete the existing file if want to generate new dca file. "
        echo "==============================================================="
    	echo ""
else
	echo ""
	echo "============================================================================"
	echo "          Running GREMLIN for DCA features.                                 "
	echo "============================================================================"
	echo ""
	$program_dir/GREMLIN_CPP/gremlin_cpp -alphabet rna -i $feature_dir/$seq_id.a2m -o $feature_dir/$seq_id.dca > $feature_dir/$seq_id.log_gremlin
	if [ $? -eq 0 ]; then
		echo ""
		echo "===================================================="
		echo "   DCA features successfully obtained from GREMLIN. "
		echo "===================================================="
		echo ""
	else
		echo ""
		echo "============================================================================="
		echo "                Error occured while running GREMLIN.  "
		echo " "
		echo "     Please check for $program_dir/GREMLIN_CPP/gremlin_cpp program.  "
		echo "============================================================================="
		echo ""
		exit 1
	fi
fi


echo ""
echo "============================================================================"
echo "          Running SPOT-RNA2 for RNA secondary structure prediction.                                 "
echo "============================================================================"
echo ""
source $program_dir/venv/bin/activate || conda activate venv
python3 $program_dir/utils/SPOT-RNA2.py --inputs $feature_dir/$seq_id.fasta --outputs $output_dir --motifs True
deactivate || conda deactivate

end=`date +%s`

runtime=$((end-start))

echo -e "\ncomputation time = "$runtime" seconds"

