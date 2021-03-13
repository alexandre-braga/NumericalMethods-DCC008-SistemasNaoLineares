#!/bin/bash

for i in {1..7} ; do
    octave-cli gauss.m && mv "PesosEPontosIntegrecaoFinal.txt" "Dom2_2/${i}_PesosEPontosIntegrecaoFinal.txt"
done
