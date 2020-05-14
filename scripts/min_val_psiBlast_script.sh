for filename in ../dataset/psiblast_res/* ; 
do 
	FILE=$(echo $filename | rev | cut -d'/' -f1 | rev) ; 
	sort -k 11 -g $filename  | head -n1 > ../dataset/min_val_psiblast_res/$FILE;
	
done


