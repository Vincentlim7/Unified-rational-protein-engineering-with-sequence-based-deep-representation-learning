#!/bin/bash

for filePath in ./dataset/new_psiblast_res/min_val_psiblast_res/* ; 
do 
	fileName=$(echo $filePath | rev | cut -d'/' -f1 | rev) ; 

	stringToSplit=`cat $filePath`
	IFS='	' read -ra my_array <<< "$stringToSplit"

	res=`grep ${my_array[1]} ./dataset/fullProtein.list | tr '\n' '\t'`
	res+=`grep ${my_array[0]} ./dataset/fullProtein.list | tr '\n' '\t'`
	res+=`echo ${my_array[10]}`
	echo $res > ./dataset/psiblast_res/min_files_processed/$fileName ; 
done





