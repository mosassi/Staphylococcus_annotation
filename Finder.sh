#! /bin/bash

#Read input parameters
while [[ $# > 1 ]]
do
arg="$1"
case $arg in
	-i|--input)
	input="$2"
	shift
	;;
	-o|--output)
	output="$2"
	shift
	;;
	-db)
	dataBase="$2"
	shift
	;;
	-py|--script)
	pyScript="$2"
	shift
	;;
	-s|--solo)
	next=false
	;;
	-t|--type)
	finder_type="$2"
	;;
esac
shift
done


#Create output directory if it didn't already exist

if [ ! -d $output ]
then
  mkdir $output
fi

#Merge together the files of the database if it was not already done

if [ ! -f $dataBase"Finder_DB_cat.fasta" ]
then
  cat $dataBase*.fsa > $dataBase"Finder_DB_cat.fasta"
fi

#Select all assembled genomes in the input directory
if [ $finder_type = "ResFinder" ] || [ $finder_type = "VirFinder" ]
then
	genome_list=$(ls $input  | egrep "\.ffn$")
elif [ $finder_type = "sRNAFinder" ]
then
	genome_list=$(ls $input  | egrep "\.fasta$")
else
	genome_list=$(ls $input  | egrep "\.fsa$")
fi

#genere la bbd pour Blast "> /dev/null " permet de ne rien afficher en console
eval makeblastdb -in $dataBase"Finder_DB_cat.fasta" -dbtype nucl > /dev/null

#Iterate blast on each genome
for genome in $genome_list
do
  genome_name=${genome/.fsa/}
  echo $input$genome #Prints the name to record progress

  fichierSortie=$output"sortieBlast_"$genome_name".txt" #Name of the output file

	if [ $finder_type = "ProphageFinder" ] || [ $finder_type = "sRNAFinder" ]
	then
		eval blastn -query $input$genome -outfmt '"6 sseqid qseqid score evalue pident sstart send slen qstart qend qlen length qseq"' -out $fichierSortie -db $dataBase"Finder_DB_cat.fasta" -evalue 0.0001 -perc_identity 90 -qcov_hsp_perc 60
	else
		eval blastn -query $input$genome -outfmt '"6 sseqid qseqid score evalue pident sstart send slen qstart qend qlen length qseq"' -out $fichierSortie -db $dataBase"Finder_DB_cat.fasta" -evalue 0.0001 -perc_identity 90 -qcov_hsp_perc 60 -num_alignments 1
	fi

	echo "blastn -query" $input$genome "-out" $fichierSortie "-db" $dataBase"Finder_DB_cat.fasta"
	echo \n
done

# Cleanup blast for plamid and prophage finder

if [ $finder_type = "ProphageFinder" ] || [ $finder_type = "sRNAFinder" ]
then
	if [ -f  $pyScript ]
	then
		echo $pyScript "-i" $output "-o" $output "-b" "-db" $finder_type ">" $output"python.log 2>&1"
		echo \n
		python $pyScript -i $output -o $output -b -db $finder_type > $output/python.log 2>&1
	fi
fi


if [ $finder_type = "ProphageFinder" ] || [ $finder_type = "sRNAFinder" ]
then
	cat $output"Clean/"*.txt > $output"blast.txt"
else
	cat $output*.txt > $output"blast.txt"
fi

# Gene list preparation

if [ $finder_type = "ResFinder" ] || [ $finder_type = "VirFinder" ]
then
	gene_list=$(tr "_" "\\t" < $output"blast.txt" | cut -f 1)
	#gene_list=$(grep -G -o '^[a-z]\{3\}' $dataBase"notes.txt" | sort -u )
elif [ $finder_type = "ProphageFinder" ]
then
	gene_list=$(ls $dataBase | grep -G -o ^uid[0-9]\* | sort -u)
fi

if [ ! $finder_type = "sRNAFinder" ]
then

	# fasta files preparation

	if [ $finder_type = "ResFinder" ]
	then
		tr "_" "\\t" < $output"blast.txt" | cut -f 1,4,16 | awk '{print ">"$1"_"$2"\n"$3}' > $output"Genes_annotation.fasta"
	elif [ $finder_type = "VirFinder" ]
	then
		tr "_" "\\t" < $output"blast.txt" | tr ":" "\\t" | cut -f 1,4,16 | awk '{print ">"$1"_"$2"\n"$3}' > $output"Genes_annotation.fasta"
	elif [ $finder_type = "ProphageFinder" ]
	then
		tr "_" "\\t" < $output"blast.txt" | cut -f 1,3,15 | awk '{print ">"$1"_"$2"\n"$3}' > $output"Genes_annotation.fasta"
	fi

	# genomes fasta files preparation

	mkdir $output"genomes_fasta"
	for genome in $genome_list
	do
	 genome_name=${genome/.fsa/}
	 grep -A 1 $genome_name $output"Genes_annotation.fasta" | grep -v "\-\-" > $output"genomes_fasta/"$genome_name".fasta"
	done

	# genes fasta files preparation

	mkdir $output"genes_fasta"
	for gene in $gene_list
	do
		grep -A 1 $gene $output"Genes_annotation.fasta" | grep -v "\-\-" > $output"genes_fasta/"$gene".fasta"
	done

	# create matrix

	if [ $finder_type = "ResFinder" ] || [ $finder_type = "VirFinder" ]
	then
		echo "\\t\\t" $genome_list > $output"matrix.txt"

		for gene in $gene_list
		do
			declare -a array=()
			array+=( $(grep $gene $dataBase"notes.txt" | tr ":" "\\t" | cut -f 2) )
			array+=( $gene )
			for genome in $genome_list
			do
				################
				genome_gene_list=$(grep -o -G '^[a-z]\{3\}')
				for $genome_gene in $genome_gene_list
				do
					if [ "$(grep -o $genome_gene $output"genes_fasta/"$gene".fasta")" ]
					then
						array+=( "1" )
					else
						array+=( "0" )
					fi
				done
			done
			echo ${array[*]} >> $output"matrix.txt"
		done

	else
		echo "\\t" $genome_list > $output"matrix.txt"

		for gene in $gene_list
		do
			declare -a array=()
			array+=( $gene )
			for genome in $genome_list
			do
				genome_name=${genome/.fsa/}
				if [ "$(grep -o $genome_name $output"genes_fasta/"$gene".fasta")" ]
				then
					array+=( "1" )
				else
					array+=( "0" )
				fi
			done
			echo ${array[*]} >> $output"matrix.txt"
		done

	fi

fi
