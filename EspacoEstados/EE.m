% at 1
%Completa, referente a questão numero 1
T = 0.1;

Hs = tf([2],[1 2 1]);

EE_continuo = ss(Hs); %Espaço de estado continuo
EE_discreto = c2d(EE_continuo,T,'zoh') %continuo para discreto
step(EE_continuo, EE_discreto);
grid on
title('')
xlabel('Tempo');
ylabel('Saída');
legend('EE Contínuo','EE Discreto')
