clear all
close all
syms s z r s0 s1 s2

dt = 0.1;


%processo
Hs = tf([2], [1 2 1]);
%processo em z
Hz = c2d(Hs, dt)

%definindo A e B
[Bz, Az] = tfdata(Hz, 'v') %tf em vetor

%% MÉTODO PAR DE POLOS DOMINANTES
Mp = 0.05; % Sobressinal máximo 0.1 -> 5%
amortecimento = sqrt( log(Mp)^2 / (pi^2 + log(Mp)^2) ) %calculo do amortecimento
ts = 2; % Tempo de acomodação 
wn = 4/(amortecimento*ts)

Acs = s^2 + 2*amortecimento*wn*s + wn^2
raizes = roots([1 2*amortecimento*wn wn^2])

polo_s1 = raizes(1) % Raíz 1
polo_s2 = raizes(2) % Raíz 2


%% Acz
A(z) = poly2sym(Az, z)
B(z) = poly2sym(Bz, z)

Ao(z) = z^2;

z1 = exp(dt*polo_s1)
z2 = exp(dt*polo_s2)

Ac(z) = expand((z-z1)*(z-z2))

%% RTS
R(z) = (z-1)*(z+r);
S(z) = s0*z^2 + s1*z + s2


polinomio_esq = expand(A(z)*R(z) + B(z)*S(z))
polinomio_dir = expand(Ac(z)*Ao(z))

C = coeffs(polinomio_esq - polinomio_dir, z); % Matriz de coeficientes resultantes

eq1 = C(1) == 0;
eq2 = C(2) == 0;
eq3 = C(3) == 0;
eq4 = C(4) == 0;

[r, s0, s1, s2] = solve([eq1, eq2, eq3, eq4], [r, s0, s1, s2]);

R(z) = (z-1)*(z+r);
%R(z) = z + r
S(z) = s0*z^2 + s1*z + s2;
%S(z) = s0*z + s1
to = Ac(1)/B(1);
T(z) = to*Ao(z);

coeff_R = sym2poly(R(z))
coeff_T = sym2poly(T(z))
coeff_S = sym2poly(S(z))

%% Controlador

% Plota o controlador melhorado
sim('DiscreteControlRTS'); % Simula no Simulink

tempo=reference(:,1)
referencia=reference(:,2);
saida=output(:,2);
controle=control(:,2);

IAE = round( sum(abs(referencia-saida)) , 2) % Calcula IAE - Integral Absolute Error
ITAE = round( sum(tempo.*abs(referencia-saida)) , 2) % Calcula ITAE - Integral Time Absolute Error
figure(1)
stairs(tempo, referencia,'r','linewidth',2); ; % Referência
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Referência')
title('Resposta Controlador RTS - Par de polos')
grid on;

figure(2)
stairs(tempo, controle,'b','linewidth',1);  % Sinal de controle
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Sinal de controle')
title('Resposta Controlador RTS - Par de polos')
grid on;


figure(3)
stairs(tempo, saida,'g','linewidth',1); ; % Processo com controlador projetado
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Saida')
title('Resposta Controlador RTS - Par de polos')
grid on;
disp('Sobressinal máximo (%):')
disp((max(saida)-1)*100)

t1 = 0; t2 = 0;
for i = 1:length(saida)
       if saida(i) > 0.0999 && t1 == 0
           t1 = tempo(i);
       end
       if saida(i) > 0.9
           t2 = tempo(i);
           break
       end
end
disp('Tempo de subida (10% até 90%) em segundos:')
disp(t2 - t1)


%%FINISH
%% Restart
clear all

syms s z r s0 s1 s2

dt = 0.1;


%processo
Hs = tf([2], [1 2 1]);
%processo em z
Hz = c2d(Hs, dt)

%definindo A e B
[Bz, Az] = tfdata(Hz, 'v') %tf em vetor

%% MÉTODO LUGAR DAS RAÍZES RECÍPROCO (LRR)

ro = 100;

%mudar implementação
num = poly2sym(Bz, z)*poly2sym(flip(Bz), z);
den = ro*poly2sym(Az, z)*poly2sym(flip(Az), z);
eq = num + den;


raizes = roots(sym2poly(eq))

polo1 = -raizes(1) % Raíz 1
polo2 = -raizes(2) % Raíz 2


%% Acz
A(z) = poly2sym(Az, z)
B(z) = poly2sym(Bz, z)

Ao(z) = z^2;

z1 = exp(dt*polo1)
z2 = exp(dt*polo2)



Ac(z) = expand((z-z1)*(z-z2))

%% RTS
R(z) = (z-1)*(z+r);
S(z) = s0*z^2 + s1*z + s2


polinomio_esq = expand(A(z)*R(z) + B(z)*S(z))
polinomio_dir = expand(Ac(z)*Ao(z))

C = coeffs(polinomio_esq - polinomio_dir, z); % Matriz de coeficientes resultantes
eq1 = C(1) == 0;
eq2 = C(2) == 0;
eq3 = C(3) == 0;
eq4 = C(4) == 0;


[r, s0, s1, s2] = solve([eq1, eq2, eq3, eq4], [r, s0, s1, s2]);


Rz = (z-1)*(z+r);
Sz = s0*z^2 + s1*z + s2;
to = Ac(1)/B(1);
Tz = to*Ao(z);

coeff_R = sym2poly(Rz)
coeff_T = sym2poly(Tz)
coeff_S = sym2poly(Sz)
%% Controlador

% Plota o controlador melhorado
sim('DiscreteControlRTS'); % Simula no Simulink

tempo=reference(:,1)
referencia=reference(:,2);
saida=output(:,2);
controle=control(:,2);

IAE = round( sum(abs(referencia-saida)) , 2) % Calcula IAE - Integral Absolute Error
ITAE = round( sum(tempo.*abs(referencia-saida)) , 2) % Calcula ITAE - Integral Time Absolute Error
figure(4)
stairs(tempo, referencia,'r','linewidth',2); ; % Referência
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Referência')
title('Resposta Controlador RTS - LRR')
grid on;

figure(5)
stairs(tempo, controle,'b','linewidth',1);  % Sinal de controle
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Sinal de controle')
title('Resposta Controlador RTS - LRR')
grid on;


figure(6)
stairs(tempo, saida,'g','linewidth',1); ; % Processo com controlador projetado
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Saida')
title('Resposta Controlador RTS - LRR')
grid on;
disp('Sobressinal máximo (%):')
disp((max(saida)-1)*100)

t1 = 0; t2 = 0;
for i = 1:length(saida)
       if saida(i) > 0.0999 && t1 == 0
           t1 = tempo(i);
       end
       if saida(i) > 0.9
           t2 = tempo(i);
           break
       end
end
disp('Tempo de subida (10% até 90%) em segundos:')
disp(t2 - t1)


