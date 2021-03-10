#!/snap/bin/octave.octave-cli

function [I] = simpson38(f, a, b, n)
    h = (b-a)/(3*n);
    soma = f(a);
    for i = 1:3*n-1,
        x = a + i*h;
        if(mod(i, 3) == 0)
            soma += 2*f(x);
        else
            soma += 3*f(x);
        endif
    endfor
    soma += f(b);
    I = 3/8 * h * soma;
endfunction
