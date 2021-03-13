#!/snap/bin/octave.octave-cli

function [I] = pontoMedio(f, a, b, n)
    h = (b-a)/n;
    soma = 0;
    x = a;
    for i = 1:n,
        xnext = a + i*h;
        soma += h*f((x + xnext)/2);
        x = xnext;
    endfor
    I = soma;
    save a.txt a b n I;
endfunction
