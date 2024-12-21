#!/usr/bin/env bash

model="blosum" #swap out for "default" and "blosum"

eval=(0.00001 0.00010 0.00100 0.01000 0.10000 1.00000 10.00000)
for i in "${eval[@]}"
do
  hmmsearch --tblout "${model}_test/e_val/${i}_results.tbl" -o "${model}_test/e_val/${i}_results.txt" -E $i --domE $i --seed 42 "models/${model}_model.hmm" complete_sequences.fasta
done

score=(0 10 50 100)
for y in "${score[@]}"
do
  hmmsearch --tblout "${model}_test/score/${y}_results.tbl" -o "${model}_test/score/${y}_results.txt" -T $y --domT $y --seed 42 "models/${model}_model.hmm" complete_sequences.fasta
done