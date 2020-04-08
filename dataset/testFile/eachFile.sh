#for file in /dir/*
#do
#  cmd [option] "$file" >> results.out
#done
#!/bin/bash

filename='/home/david/Documents/juliaProjet/optimiseSearchProtDB/dataset/test_dataset.list_family.complete'

filename2='/home/david/Documents/GitProjet/LU3IN013_Projet_Source/dataset/fastas'

declare -a file_array_test

declare -a file_array_fastas

cpt=0
while read -r line
do
  #Take the 7 first letter of each line
  id=$(cut -c-7 <<< "$line")
  makeblastdb -in /home/david/Documents/GitProjet/LU3IN013_Projet_Source/dataset/fastas/$id.fasta -out $id.db -dbtype prot
  #code for passing id to other script file as parameter
  cpt=$((cpt+1))
  echo "$cpt"
done < "$filename"




#for entry in "/home/david/Documents/GitProjet/LU3IN013_Projet_Source/dataset/fastas"/*
#do
#	filename=$(basename -- "$entry")
#	filename="${filename%.*}"
#  	#echo "coucou"
#        #basename "$filename"
#	file_array_fastas+=("${filename%.*}")
#	echo "ok"
#done

#for i in "${file_array_test[@]}"
#  do
#	for j in "${file_array_fastas[@]}"
#		do
#			if [ $i == $j ] then
#				echo "yes"
#			else
#				echo "no"
#			fi
#done

#for i in "${file_array_test[@]}"; do 
#	echo "$i";
#done



#for thelist in "${file_array_test[@]}" ; do
#	echo "enter"
#	for key in "${file_array_fastas[@]}" ; do
#        	echo "the key is: $key"
#    done
#done'
