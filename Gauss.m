#!/home/kuroneko/snap/octave -qf
1;

function [w,t,g] = metodoDeIntegracao(a,b,N)
    w = zeros(N, 1);
    t = zeros(N, 1);
    for i = 1:N/2,
        w(i) = i*(b-a)/(2*n);
        w(N-i) = w(i);
        t(i) = a + i*w(i)/2;
        t(N-i) = t(i);
    endfor
    g = simpson38();
endfunction

function [f] = funcaoWTdw(w,t,g,i,j,N,E)
  f = 0;
  for k = 1:N
    if k == j
      f += (w(k) + E) * t(k)^i - g(k);
    else 
      f += w(k) * t(k)^i - g(k);
    endif
  endfor
endfunction

function [f] = funcaoWTdt(w,t,g,i,j,N,E)
  f = 0;
  for k = 1:N
    if k == j
      f += w(k) * (t(k) + E)^i - g(k);
    else 
      f += w(k) * t(k)^i - g(k);
    endif
  endfor
endfunction

function [df] = derivadaParcial(w,t,i,g,j,N,E)
  if j <= N
    f1 = funcaoWTdw(w,t,g,i,jN,E);
    f2 = funcaoWTdw(w,t,g,i,0,N,0);
    df = (f1-f2)/E;
  endif
  else
    f1 = funcaoWTdt(w,t,g,i,j,N,E);
    f2 = funcaoWTdt(w,t,g,i,0,N,0);
    df = (f1-f2)/E;
  endif
endfunction

function [J] = jacobiana(w,t,g,N)
    J = zeros(N, N);
    E = 10e-8;
    for j = 1:(2*N)
       for i = 1:(2*N)
          J(i,j) = derivadaParcial(w,t,g,i,j,E); 
       endfor
    endfor
endfunction


addpath(pwd);
args = argv();

printf("Iniciando o programa\n");
tol = 10e-8;
a = -1; #input('Insira o valor de a: ');
b =  1; #input('Insira o valor de b: ');
N =  7; input('Defina o número de pontos de integração: ');

[w,t,g] = metodoDeIntegracao(a,b,N);
[J] = jacobiana(w,t,g,N);


while norm(s(i), inf) > tol
    
endwhile

save PesosEPontosIntegrecao.text w t;
save 

printf("Fim programa\n");

