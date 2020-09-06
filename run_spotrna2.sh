#!/bin/bash

start=`date +%s`

input="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
input_dir=$(dirname $input)
seq_id=$(basename $(basename $input) | cut -d. -f1)
program_dir=$(dirname $(readlink -f $0))

path_blastn=$program_dir/ncbi-blast-2.10.*+/bin
path_blastn_database=$program_dir/nt_database/nt
path_infernal=$program_dir/infernal-1.1.3-linux-intel-gcc/binaries
path_infernal_database=$program_dir/nt_database/nt

mkdir -p $input_dir/${seq_id}_features && mkdir -p $input_dir/${seq_id}_outputs
echo ">"$seq_id > $input_dir/${seq_id}_features/$seq_id.fasta
awk -i inplace '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $input 
tail -n1 $input >> $input_dir/${seq_id}_features/$seq_id.fasta

feature_dir=$input_dir/${seq_id}_features
output_dir=$input_dir/${seq_id}_outputs

#exit 1

if [ ! -f $path_blastn_database ];  then
        echo ""
        echo "=============================================================================================================================================="
        echo "            Looks like nt database doesn't exists in the directory mounted to docker container directory /mnt/ "
        echo "            If you want to download the database now, please make sure you have enough space in mounted directory and internet connection have"
        echo "            enough bandwidth as file is of size 270 GBs after unzip. It may take forever to download if internet is slow!"
        echo "=============================================================================================================================================="
        echo ""

        echo -n "Type 'y' for download or any other key to exit: "    
        read userinput

        if [[ $(echo $userinput | tr '[A-Z]' '[a-z]') == 'y' ]]; then

                echo ""
                echo "======================================================================================================================================"
                echo "       Start downloading nt database form ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz link. May take few hours to download. "
                echo "======================================================================================================================================"
                echo ""
                wget -c "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz" -O /nt_database/nt.gz

        
                if [[ $? -eq 0 ]]; then 
                        echo ""
                        echo "======================================================================="
                        echo "            Download of nt database is complete.  "
                        echo "======================================================================="
                        echo ""
                else
                        echo ""
                        echo "======================================================================="
                        echo "            Error! Unable to download sucessfully.  "
                        echo "========================================================================"
                        echo ""
                        exit 1        
                fi

                echo ""
                echo "================================================================================================================================"
                echo "            Start unziping the downloaded nt database. May take few hours as size of unzipped file is around 270 GBs. "
                echo "================================================================================================================================"
                echo ""
                
        ############ unzip the nt data base file ############
                gunzip /nt_database/nt.gz

        else
                exit 1
        fi

fi


if [ -f $feature_dir/$seq_id.a2m ];	then
        echo ""
        echo "=============================================================================================================================================="
        echo "    MSA file $feature_dir/$seq_id.a2m from Infernal Pipeline already exists for query sequence $feature_dir/$seq_id.fasta.  "
        echo "=============================================================================================================================================="
    	echo ""
else

    if [[ ! -f "$path_blastn_database.nal" ]]; then
            echo ""
            echo "=========================================================================================================================="
            echo "    Nucleotide database file $path_database/nt need to formated to use with 'makeblastdb' program in BLAST package                "
            echo "    Formatting the nt database may take 2-3 hours as size of file is around 270 GBs"
            echo "=========================================================================================================================="
            echo ""
            $path_blastn/makeblastdb -in $path_database/nt -dbtype nucl
            
            if [[ $? -eq 0 ]]; then
                    echo ""
                    echo "==============================================================================="
                    echo "                    nt database formatted successfully!"
                    echo "==============================================================================="
                    echo ""
            else
                    echo ""
                    echo "======================================================================================================================================"
                    echo "  Error occured while formatting the nt database in $path_database/ mounted directory                                                 "
                    echo "  Please, try fomatting the nt database on your system outside the docker container using 'makeblastdb' program in BLAST package      "
                    echo "======================================================================================================================================"
                    echo ""
                    exit
            fi                      
    fi


        #################### run blastn ######################
        if [ -f $feature_dir/$seq_id.bla ];       then
                echo ""
                echo "=========================================================================================="
                echo "        MSA file $feature_dir/$seq_id.bla from BLASTN already exists for query sequence.  "
                echo "=========================================================================================="
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
			

	######## reformat the output ################
    echo ""
    echo "========================================================================================"
    echo "         Converting $feature_dir/$seq_id.bla from BLASTN to $feature_dir/$seq_id.sto.   "
    echo "========================================================================================"
    echo ""
	$program_dir/utils/parse_blastn_local.pl $feature_dir/$seq_id.bla $feature_dir/$seq_id.fasta $feature_dir/$seq_id.aln
	$program_dir/utils/reformat.pl fas sto $feature_dir/$seq_id.aln $feature_dir/$seq_id.sto

	######## predict secondary ################
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

	######## run infernal ################
    echo ""
    echo "=============================================================================================================="
    echo "      Building Covariance Model from BLASTN alignment (with SS from SPOT-RNA) from $feature_dir/$seq_id.sto file.         "
    echo "=============================================================================================================="
    echo ""
	$path_infernal/cmbuild --hand -F $feature_dir/$seq_id.cm $feature_dir/$seq_id.sto

    echo ""
    echo "======================================================================="
    echo "          Calibrating the Covariance Model $feature_dir/$seq_id.cm.    "
    echo "======================================================================="
    echo ""
	$path_infernal/cmcalibrate $feature_dir/$seq_id.cm

    echo ""
    echo "======================================================================================================================"
    echo "        Second round of homologous sequences search using the calibrated covariance model $feature_dir/$seq_id.cm.    "
    echo "                 May take 15 mins to few hours for this step.                                                         "
    echo "======================================================================================================================"
    echo ""
	$path_infernal/cmsearch -o $feature_dir/$seq_id.out -A $feature_dir/$seq_id.msa --cpu 24 --incE 10.0 $feature_dir/$seq_id.cm $path_infernal_database

	######### reformat the output for dca input ###############
    echo ""
    echo "======================================================================="
    echo "          Reformatting the output alignment $feature_dir/$seq_id.msa   "
    echo "          for PSSM and DCA features by removing the gaps.   "
    echo "======================================================================="
    echo ""
	$path_infernal/esl-reformat --replace acgturyswkmbdhvn:................ a2m $feature_dir/$seq_id.msa > $feature_dir/temp.a2m

	######### remove duplicates sequences from the alignment ###############
    echo ""
    echo "======================================================================="
    echo "          Removing duplicates from the alignment.                      "
    echo "======================================================================="
    echo ""
	$program_dir/utils/seqkit rmdup -s $feature_dir/temp.a2m > $feature_dir/$seq_id.a2m

	############# multiline fasta to single line fasta file   #############
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $feature_dir/$seq_id.a2m | sed '/^$/d' > $feature_dir/temp.a2m 

	############# add query sequence at the top of MSA file  #############
    cat $feature_dir/$seq_id.fasta $feature_dir/temp.a2m > $feature_dir/$seq_id.a2m 

fi

############# getting PSSM from alignment  #############
echo ""
echo "======================================================================================"
echo "          Extracting PSSM features from the alignment $feature_dir/$seq_id.a2m.       "
echo "======================================================================================"
echo ""
$program_dir/utils/getpssm.pl $feature_dir/$seq_id.fasta $feature_dir/$seq_id.a2m $feature_dir/$seq_id.pssm

######### run linearpartition RNA secondary structure base-pair probability predictor ###############
echo ""
echo "============================================================================"
echo "          Running LinearPartition-V for base-pair probabilty features.      "
echo "============================================================================"
echo ""
tail -n +2 $feature_dir/$seq_id.fasta | $program_dir/LinearPartition/linearpartition -V -r $feature_dir/$seq_id.prob

######### run GREMLIN for DCA features ###############
echo ""
echo "============================================================================"
echo "          Running GREMLIN for DCA features.                                 "
echo "============================================================================"
echo ""
$program_dir/GREMLIN_CPP/gremlin_cpp -alphabet rna -i $feature_dir/$seq_id.a2m -o $feature_dir/$seq_id.dca > $feature_dir/$seq_id.log_gremlin

############ save output in ct, bpseq and base-pair matrix #############
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

