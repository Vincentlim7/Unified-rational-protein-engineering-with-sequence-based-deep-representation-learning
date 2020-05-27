for filePath in ./dataset/psiblast_res/full_files/* ; 
do 
	fileName=$(echo $filePath | rev | cut -d'/' -f1 | rev) ; 
	sort -k 11 -g $filePath  | head -n1 > ./dataset/psiblast_res/min_files/$fileName;
	
done


