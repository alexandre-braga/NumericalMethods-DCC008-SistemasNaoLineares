#!/snap/bin/octave.octave-cli
1;
format long;

function [w,t] = metodoDeIntegracao(a,b,N)
    w = double(zeros(N, 1));
    t = double(zeros(N, 1));
    for i = 1:ceil(N/2)
        w(i) = i*(b-a)/(2*N);
        w(N+1-i) = w(i);
        if i == N+1-i
          t(i) = (a+b)/2;
        else
          t(i) = a + i*w(i)/2;
          t(N+1-i) = -t(i);
        endif
    endfor
endfunction

function [F] = funcaoWT(w,t,g,N)
  F = double(zeros(2*N,1));
  for i = 1:(2*N),
    for k = 1:N,
      F(i) += w(k) * t(k).^(i-1);
    endfor
    F(i) -= g(i);
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
    if k == j-N
      f += w(k) * (t(k) + E).^(i-1);
    else 
      f += w(k) * t(k).^(i-1);
    endif
  endfor
  f -= g(i);
endfunction

function [df] = derivadaParcial(w,t,g,i,j,N,E)
  if j <= N
    df = funcaoWTdw(w,t,g,i,j,N,E);
  else
    df = funcaoWTdt(w,t,g,i,j,N,E);
  endif
endfunction

function [J] = jacobiana(w,t,g,N,E,F)
    J = double(zeros(2*N,2*N));
    for i = 1:(2*N),
       for j = 1:(2*N),
          J(i,j) = (derivadaParcial(w,t,g,i,j,N,E)-derivadaParcial(w,t,g,i,j,N,0))/E;
       endfor
    endfor
endfunction


addpath(pwd);
args = argv();

printf("Iniciando o programa\n");
tol = 10e-8;
E = 10e-8;

N = input('Defina o número de pontos de integração: ');
a = input('Insira o valor de a: ');
b = input('Insira o valor de b: ');


[w,t] = metodoDeIntegracao(a,b,N);

g = zeros(2*N,1);
for i = 1:(2*N),
    y = @(x) x.^(i-1);
    g(i) = pontoMedio(y,a,b,1000);
endfor

save PesosIniciais.txt w t g a b N tol E;

cont = 0;
while true
    [F] = funcaoWT(w,t,g,N);
    [J] = jacobiana(w,t,g,N,E,F);
    
    s = J\(-F);

    for i = 1:(2*N)
      if i <= N
        wnext(i) = s(i) + w(i);
      else
        tnext(i-N) = s(i) + t(i-N);
      endif
    endfor
    
    erro = norm(s, inf);
    if(erro <= tol)
        break;
    endif

    nome = sprintf("log/PesosEPontosIntegrecao%d.txt", cont);
    save(nome, 'F', 'J', 's', 'w', 't', 'erro', 'tol', 'E', 'N');
    w = wnext;
    t = tnext;
    cont++;

endwhile

save PesosEPontosIntegrecaoFinal.txt F J s w t erro tol E N;

printf("Fim do programa\n");

