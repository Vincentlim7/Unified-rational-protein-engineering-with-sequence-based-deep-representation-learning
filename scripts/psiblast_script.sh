./ncbi-blast-2.10.0/ncbi-blast-2.10.0+/bin/psiblast -db ./dataset/queries/$1.db -in_pssm ./dataset/pssmModels/$2.pssm -outfmt 6 >> ./dataset/psiblast_res/$1.results

