%% Importando biblioteca para controle de processos.
addpath '..\ProcessosIndustriais\src'
addpath '..\continous time'

%% Definindo referência simulação incial
reference_amplitude = 1;

%% Removendo distrubio para simulação inicial

input_noise_amplitude = 0;
input_noise_time = 0;
ruido_simulation_input = timeseries();
ruido_simulation_input.Time = [0];
ruido_simulation_input.Data = [0]';

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

%% Tunning

tunning_method = 'CHRSR';

switch tunning_method
    case 'ZN'
        tunning = ZieglerNichols(dynamics_parameters);
    case 'CC'
        tunning = CCTunning(dynamics_parameters);
    case 'ITAERT'
        tunning = ITAERTunning(dynamics_parameters);
    case 'CHRRT'
        tunning = CHRRTunning(dynamics_parameters);
    case 'CHR20'
        tunning = CHR20Tunning(dynamics_parameters);
    case 'IAER'
        tunning = IAERTunning(dynamics_parameters);
    case 'IAESR'
        tunning = IAESRTunning(dynamics_parameters);
    case 'ITAEST'
        tunning = ITAESTunning(dynamics_parameters);
    case 'CHRSR'
        tunning = CHRSRTunning(dynamics_parameters);
end

controller_parameters = tunning.getPIDParameters();

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

%% Frequency Analysis

% Para encontrar a frequência de corte e largura de banda
% precisamos fazer uma análise do diagrama de bode em malha aberta e
% fechada
controller = pid(PROPORTIONAL_GAIN, INTEGRAL_GAIN, DERIVATIVE_GAIN);
open_loop_tf = series(processo, controller);
closed_loop_tf = feedback(open_loop_tf, 1);

margins = allmargin(open_loop_tf);

figure(30)
bode(open_loop_tf);grid
figure(40)
bode(closed_loop_tf);grid

%% Definindo período de amostragem

% Pela frequência de corte: 0.064 < T < 0.213
% Pela largura de banda: T = 0.0784
% Considerando os valores, optou-se por utilizar T = 0.1s

sample_time = 0.1;

%% Rodando a Simulação Discreta

% É importante ressaltar que os dados obtidos com
% a simulação serão salvos no workspace

sim('DiscreteBaseControl');

%% Rodando a Simulação Continua

% É importante ressaltar que os dados obtidos com
% a simulação serão salvos no workspace

sim('CustomBaseControl');

%% Plotando as respostas no tempo

figure(2);
% A classe SimulationVisualizer possui métodos para plotar 
% os dados da simulação.
sgtitle('Controller Tunning')

%% Saída: y(t)

SimulationVisualizer.plotOutput(Reference, OutputRead, 'r');
SimulationVisualizer.plotDiscreteOutput(DiscreteReference, DiscreteOutputRead);

%% Disturbio constante na entrada

input_noise_amplitude = 1;
input_noise_time = 10;

%% Rodando a Simulação Discreta com disturbio na entrada

%sim('DiscreteBaseControl');
figure(3);
SimulationVisualizer.plotDiscreteOutput(DiscreteReference, DiscreteOutputRead);

%% Noise

reference_amplitude = 0;
input_noise_amplitude = 0;
input_noise_time = 0;

load dados_ruido.mat 
% tempo: contém o vetor tempo 
% ruido: contém o vetor das amplitudes do sinal de ruído 

ruido_simulation_input = timeseries();
ruido_simulation_input.Time = tempo;
ruido_simulation_input.Data = ruido';
%% Rodando simulação com referência nula e ruido na saída

sim('DiscreteBaseControl');

%% Sinal de controle e ruido

figure(5);
%subplot(211);
SimulationVisualizer.plotDiscreteControlSignal(DiscreteInput);
%%
subplot(212);
plot(tempo,ruido, 'r');
xlabel('Time [s]');
ylabel('Noise [-]');
grid on;

%% Sinal de ruido e saida

figure(6);
subplot(211);
plot(tempo,ruido, 'r');
xlabel('Time [s]');
ylabel('Noise [-]');
grid on;
subplot(212);
SimulationVisualizer.plotDiscreteOutput(DiscreteReference, DiscreteOutputRead);

%% Valor médio absoluto

noise_absolute_mean = mean(abs(ruido));

input_absolute_mean = mean(abs(DiscreteInput(:,2)));

disp("O valor absoluto da média do ruído e sinal de controle são\n");
fprintf("Ruído: %f       Sinal de Controle: %f\n", noise_absolute_mean, input_absolute_mean);

