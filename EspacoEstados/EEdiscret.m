%at2
%projeto por alocação dee polos: controlador 
%pag 71

%%Definindo valores
T = 0.1; % Intervalo de amostragem
Hs = tf([2],[1 2 1]);
EE_continuo = ss(Hs);
EE_discreto = c2d(EE_continuo,T,'zoh'); %tira os valores de A B C D

wn=1;
amortecimento=0.7;

%%Tirando raizes
raizes = roots([1 2*amortecimento*wn wn^2]);
s1 = raizes(1); % Raíz 1 -0.7000 + 0.7141i
s2 = raizes(2); % Raíz 2 -0.7000 - 0.7141i
z1 = exp(s1*T); % Raíz discreta 1 0.9300 + 0.0665i
z2 = exp(s2*T); % Raíz discreta 2 0.9300 - 0.0665i

syms z
Pk(z) = expand((z-z1)*(z-z2)); %z^2 - 1.8600 + 0.8694
vpa(Pk(z), 4); %arredonda pra 4 digitos

%% Matrizes de EE organizadas de acordo com o MATLAB
A = EE_discreto.A; B = EE_discreto.B; C = EE_discreto.C; D = EE_discreto.D;
%A = [0.9983 0.0565;-0.0565 0.8853]; B = [0.0035; 0.1130]; C = [1 0]; D = [0];

PkA = A^2 - 1.915979346*A + 0.9194312561*eye(2); %   [ 0.0141    0.0096]
                                                 %   [-0.0096   -0.0051]

%% Definição para o controlador
Wc = [B A*B];                                     %   [0.1810    0.1465]
                                                  %   [0.0094    0.0257]



K = [0 1]*Wc^(-1)*PkA % -0.5714   -0.3094

%% Dead-Beat
Wo = [C; C*A];                                     %     [0    1.0000     ]
                                                   %     [0.0905    0.9953]
PL(z)=z^2;
PLA=A^2;
L = PLA*Wo^-1*[0;1] %L = [7.2387 1.8097]


%% Ganho KC
Kc = 1/( C*(eye(2)-(A-B*K))^-1*B ) %0.1906


%% sistema MF

Amf = [ A-B*K , B*K; zeros(size(A)) , A-L*C];
Bmf = [B*Kc;[0;0]];
Cmf=  [C 0 0];
Lmf=  [7.2387;1.8097;0;0]
sys=ss(Amf,Bmf,Cmf,D,T)
figure(1)
step(tf([1],[1]),sys);
grid on

title('Controlador MF por EE')
xlabel('Tempo');
ylabel('Saída');
legend('Referência', 'Sistema em EE')

%% COMEÇO DA 3
%Definindo valores
Kmf = [-0.5714,-0.3094];


t=0:0.06:30-0.1;
N = length(t);

w=ones(size(t)); % referencia
x = [0; 0; 0; 0] % x
xz = [0; 0; 0; 0]; % x~
u = [0; 0]; % u
y = zeros(1,N);
v = [zeros(1,250),ones(1,N-250)];
%=zeros(1,N);

for k=1:N
    x(:,k+1) = Amf*x(:,k) + Bmf*(u(k) + v(k));
    u(:,k) = -Kmf*x(k) + Kc*w(k);  
    xz(:, k+1) = Amf*xz(:,k) + Bmf*u(k) + Lmf*( y(k)-Cmf*xz(:,k) ); 
    y(k) = Cmf*x(:,k) + D*u(k);
end
figure(2)
plot(y)
hold on
%plot(ones(size(y)))
hold off
%grid on
title('Disturbio na entrada do controlador EE')
xlabel('Iterações');
ylabel('Saída');
legend('Sistema em EE')



