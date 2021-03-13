#!/bin/bash

a=$1
b=$2

mkdir "Dom${a##*[+-]}_${b##*[+-]}"

for i in {1..7} ; do
    octave-cli gauss.m && mv "PesosEPontosIntegrecaoFinal.txt" "Dom${a##*[+-]}_${b##*[+-]}/${i}_PesosEPontosIntegrecaoFinal.txt"
done
