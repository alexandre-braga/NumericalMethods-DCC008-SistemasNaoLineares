#!/snap/bin/octave.octave-cli
1;
function integra(f, a, b, n)
    w = zeros(n, 1);
    t = zeros(n, 1);
    printf("%4s    %4s    %4s\n", "w", "t", "f(t)");
    for i = 1:n/2,
        w(i) = i*(b-a)/(2*n);
        t(i) = a + i*w(i)/2;
        printf("%.4f %.4f %.4f\n", w(i), t(i), f(t(i)));
    endfor
endfunction

addpath(pwd);
args = argv();

f = @(x) 1/x;
I = simpson38(f, 1, 2, 6);
printf("I = %f\n", I);

printf("Fim do programa\n");
