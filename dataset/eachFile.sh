#!/bin/bash

filename='/home/david/Documents/juliaProjet/optimiseSearchProtDB/dataset/test_dataset.list_family.complete'


cpt=0
while read -r line
do
	#Take the 7 first letter of each line and put the result in id variable
	id=$(cut -c-7 <<< "$line")
	#Execute the makeblastdb command
	makeblastdb -in /home/david/Documents/GitProjet/LU3IN013_Projet_Source/dataset/fastas/$id.fasta -out d1a0aa_.db -dbtype prot
 	#code for passing id to other script file as parameter
  	cpt=$((cpt+1))
done < "$filename" #Passing the argument for the while loop




