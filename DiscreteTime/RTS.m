clear all
syms s z r s0 s1 s2

dt = 0.1;
tempo = 0:dt:30;

%processo
Hs = tf([2], [1 2 1]);
%processo em z
Hz = c2d(Hs, dt)

%definindo A e B
[Bz, Az] = tfdata(Hz, 'v'); %tf em vetor

%% MÉTODO PAR DE POLOS DOMINANTES
Mp = 0.05; % Sobressinal máximo 0.1 -> 5%
amortecimento = sqrt( log(Mp)^2 / (pi^2 + log(Mp)^2) ) %calculo do amortecimento
ts = 2; % Tempo de acomodação 
wn = 4/(amortecimento*ts)

Acs = s^2 + 2*amortecimento*wn*s + wn^2
raizes = roots([1 2*amortecimento*wn wn^2])

polo_s1 = raizes(1) % Raíz 1
polo_s2 = raizes(2) % Raíz 2

%% MÉTODO LUGAR DAS RAÍZES RECÍPROCO (LRR)

ro = 100;

%H(z) = poly2sym(coeff_B, z)/poly2sym(coeff_A, z)
%eq(z) = 1 + 1/ro*H(z)*H(1/z)
%raizes = roots(sym2poly(eq(z)))

num = poly2sym(Bz, z)*poly2sym(flip(Bz), z);
den = ro*poly2sym(Az, z)*poly2sym(flip(Az), z);
num = num + den;
raizes = roots(sym2poly(num));
polo1 = -raizes(3) % Raíz 1
polo2 = -raizes(4) % Raíz 2

%% Acz
A(z) = poly2sym(Az, z)
B(z) = poly2sym(Bz, z)

%Ao(z) = z^2;
Ao(z) = (z + 1.5)^2; %verificar Ao(z)

z1 = exp(dt*polo_s1)
z2 = exp(dt*polo_s2)

Ac(z) = expand((z-z1)*(z-z2))

%% RTS
R(z) = (z-1)*(z+r);
%R(z) = z+r
S(z) = s0*z^2 + s1*z + s2
%S(z) = s0*z + s1

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




