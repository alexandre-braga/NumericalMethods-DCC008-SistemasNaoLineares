#!/snap/bin/octave.octave-cli -qf
1;
function gauss(f, a, b, n)
    w = zeros(n, 1);
    t = zeros(n, 1);
    for i = 1:n/2,
        w[i] = i*(b-a)/(2*n);
        t[i] = a + i*w[i]/2;
    endfor
endfunction

addpath(pwd);
args = argv();

f = @(x) sqrt(1 + x.^2);
gauss(@f, 0, 1, 10);

printf("Fim do programa\n");
