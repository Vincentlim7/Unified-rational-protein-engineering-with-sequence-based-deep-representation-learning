#!/bin/bash

prot=$1
touch res/$1_res
for f in res/*results; do
	cat $f >> res/$1_res
done
