#!/bin/bash

n=$1

mkdir "Dom${n}_${n}"

for i in {1..7} ; do
    octave-cli gauss.m && mv "PesosEPontosIntegrecaoFinal.txt" "Dom${n}_${n}/${i}_PesosEPontosIntegrecaoFinal.txt"
done
