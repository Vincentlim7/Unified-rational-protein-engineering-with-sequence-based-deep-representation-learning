filename='./test_dataset.list'

filename2='./fastas'

declare -a file_array_test

declare -a file_array_fastas

while read -r line
do
  #Take the 7 first letter of each line
  id=$(cut -c-7 <<< "$line")
  makeblastdb -in ./fastas/$id.fasta -out ./testFiles/$id.db -dbtype prot
done < "$filename"
