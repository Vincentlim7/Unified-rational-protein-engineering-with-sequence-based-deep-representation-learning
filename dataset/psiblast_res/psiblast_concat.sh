#!/bin/bash

#allo to concatenate the results file of psiblast output
function concat_function(){
	#just an echo of your sequence
	echo "$1"
	#Execute ths psiblast command
	psiblast -db $1.db -in_pssm $1.pssm -outfmt 6 > $1.results
	#Put in a variable the name of the results file create by psiblast
	resultFile="${1}.results"
	#Cat the file then concatenate the global results file with the result of the cat
	cat "$resultFile">>"$2"
}

#Call the function with your arguments pass in command line
concat_function $1 $2

#Example:

#bash concat.sh d1a4fb_ concatResults.results 

#here $1 is d1a4fb_ and $2 is concatResults.results
