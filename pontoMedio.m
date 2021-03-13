#!/snap/bin/octave.octave-cli

function [I] = pontoMedio(f, a, b, n)
    h = (b-a)/n;
    I = 0;
    x = a;
    for i = 1:n,
        xnext = a + i*h;
        I += h*f((x + xnext)/2);
        x = xnext;
    endfor
endfunction
