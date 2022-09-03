%% Importando biblioteca para controle de processos.
addpath 'C:\TFS\Controladores Digitais\ProcessosIndustriais\src'

%% Definindo o processo e resposta em malha aberta
num = 2;
den = [1 2 1];

processo = tf(num, den);

%               2
%  H(s) = --------------
%         s^2 + 2 s + 1

% Resposta em malha aberta:
figure(1);
step(processo);

% Tempo de simulação:
tempo = 0:0.1:30;

%% Process Dynamics

% Nesta seção obtemos os parâmetros que definem
% a dinâmica do processo:

% theta: atraso de transporte
% tau: constante de tempo
% k: ganho estático

dynamics = ProcessDynamics(processo, tempo);
dynamics_parameters = dynamics.getDynamicsParameters();

%% Ziegler Nichols

zn = ZieglerNichols(dynamics_parameters);
controller_parameters = zn.getPIDParameters();

%% CC

cc = CCTunning(dynamics_parameters);
controller_parameters = cc.getPIDParameters();

%% Controlador ZN

% Aqui utilizaremos o modelo do simulink disponi-
% bilizado para simulações: BaseControl.slx

% Dentro do simulink adicionei um bloco de con-
% trolador PID. Sendo que cada ganho recebe o 
% valor de uma variável do workspace:

% P : PROPORTIONAL_GAIN
% I : INTEGRAL_GAIN
% D : DERIVATIVE_GAIN

% Lembrando que I = Kp / Ti e D = Kp * Td
% Onde: 
% Kp: Ganho proporcional
% Ti: Tempo integrador
% Td: Tempo derivador

% Como estamos projetando um controlador PI, logo:
PROPORTIONAL_GAIN = controller_parameters.Kp;
INTEGRAL_GAIN = controller_parameters.Kp / controller_parameters.Ti;
DERIVATIVE_GAIN = controller_parameters.Kp * controller_parameters.Td;

%% Rodando a Simulação

% É importante ressaltar que os dados obtidos com
% a simulação serão salvos no workspace

sim('CustomBaseControl');

%% Plotando as respostas no tempo

figure(2);
% A classe SimulationVisualizer possui métodos para plotar 
% os dados da simulação.
sgtitle('Controller Tunning')
%% Saída: y(t)

subplot(311);
SimulationVisualizer.plotOutput(Reference, OutputRead);

%% Sinal de Controle: u(t)

subplot(312);
SimulationVisualizer.plotControlSignal(Input);
fprintf('O valor máximo do sinal de controle é %f\n', max(Input(:,2)));


%% Error: e(t)

subplot(313);
SimulationVisualizer.plotError(Error);

%% Controle de Performance

Rt=Reference(:,2);
Yt=OutputRead(:,2);
Ut=Input(:,2);

IAE = sum(abs(Rt-Yt));

ITAE = 0;
tempo_simulacao = 0:0.01:20;
for i = 1:size(tempo_simulacao, 2)
    ITAE = ITAE + abs(Rt(i)-Yt(i)) * tempo_simulacao(i);
end

TV = sum(abs(diff(Ut)));

disp('Para o método, obteve-se:')
fprintf('IAE: %f   ITAE: %f  TV: %f\n', IAE, ITAE, TV);

