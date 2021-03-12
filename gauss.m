#!/snap/bin/octave.octave-cli
1;

function [w,t] = metodoDeIntegracao(a,b,N)
    w = zeros(N, 1);
    t = zeros(N, 1);
    for i = 1:ceil(N/2)
        w(i) = i*(b-a)/(2*N);
        w(N+1-i) = w(i);
        t(i) = a + i*w(i)/2;
        t(N+1-i) = t(i);
        if i == N+1-i
          t(i) = (a+b)/2;
        endif
    endfor
endfunction

function [F] = funcaoWT(w,t,g,N)
  F = zeros(2*N,1);
  for i = 1:(2*N),
    f = 0;
    for k = 1:N,
      f += w(k) * t(k).^i - g(k);
    endfor
    F(i) = f;
  endfor
endfunction

function [f] = funcaoWTdw(w,t,g,i,j,N,E)
  f = 0;
  for k = 1:N,
    if k == j
      f += (w(k) + E) * t(k).^(i-1);
    else 
      f += w(k) * t(k).^(i-1);
    endif
  endfor
  f -= g(i);
endfunction

function [f] = funcaoWTdt(w,t,g,i,j,N,E)
  f = 0;
  for k = 1:N,
    if k == j
      f += w(k) * (t(k) + E).^(i-1);
    else 
      f += w(k) * t(k).^(i-1);
    endif
  endfor
  f -= g(i);
endfunction

function [df] = derivadaParcial(w,t,g,i,j,N,E)
  if j <= N
    f1 = funcaoWTdw(w,t,g,i,j,N,E);
    f0 = funcaoWTdw(w,t,g,i,0,N,0);
    df = (f1-f0)/E;
  else
    f1 = funcaoWTdt(w,t,g,i,j,N,E);
    f0 = funcaoWTdt(w,t,g,i,0,N,0);
    df = (f1-f0)/E;
  endif
endfunction

function [J] = jacobiana(w,t,g,N,E)
    J = zeros(2*N,2*N);
    for i = 1:(2*N),
       for j = 1:(2*N),
          J(i,j) = derivadaParcial(w,t,g,i,j,N,E); 
       endfor
    endfor
endfunction


addpath(pwd);
args = argv();

printf("Iniciando o programa\n");
tol = 10e-8;
E = 10e-8;
a = -1;
b =  1;
N =  7;
f = @(x) x.^2 + log(x);

#tol = input('Insira a tolerância: ');
#E = input('Insira a perturbação: ');
#N = input('Defina o número de pontos de integração: ');
#a = input('Insira o valor de a: ');
#b = input('Insira o valor de b: ');
#f = str2func(input('Defina a função na forma f = @(x) expressão: ', 's'));
printf("f = %s\n", func2str(f));

[w,t] = metodoDeIntegracao(a,b,N);

g = zeros(2*N,1);
for i = 1:(2*N),
    y = @(x) x.^(i-1);
    g(i) = simpson38(y,a,b,2*N);
endfor

printf("g: ");
disp(g);

[w,t] = metodoDeIntegracao(a,b,N);
printf("W: ");
disp(w);
printf("T: ");
disp(t);
while true
    [F] = funcaoWT(w,t,g,N);
    [J] = jacobiana(w,t,g,N,E);
    s = J\(-F);
    for i = 1:(2*N)
      if i <= N
        wnext(i) = s(i) + w(i);
      else
        tnext(i-N) = s(i) + t(i-N);
      endif
    endfor
    w = wnext;
    t = tnext;
    if(norm(s, inf) <= tol)
        break;
    endif
endwhile

save PesosEPontosIntegrecao.txt a b tol E h N w t;

printf("Fim do programa\n");

