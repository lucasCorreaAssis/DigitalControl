%% Processo em tempo contínuo e interavalo de amostragem

Hs = tf(2, conv([1,1],[1,1]));
dt = 0.1;

%% Controlador em tempo contínuo

Kp=2.8;
Ki=1;
Kd=0.42;

s = tf('s');
z = tf('z', dt);

Cs = Kp + Ki/s + (Kd*s / (Kd/20*s + 1));

%% Discretização - Aproximação forward e backward

sf = (z-1)/dt;
sb = (z-1)/z*dt;

Cz = Kp + Ki/sf + (Kd*sb / (Kd/20*sb + 1))

%% Discretização segurador de ordem 0

Hz = c2d(Hs, dt);

%% Equação a diferenças

t = 0:dt:30;

w = ones(size(t));

e = [0,0,0];
u = [0,0,0];
y = [0,0,0];

timeSize = size(t);
N = timeSize(1,2) - 1;

for k=3:1:N
    e(k+1) = w(k) - y(k);
    u(k+1) = u(k) - 0.002*u(k-1) + 2.848*e(k+1) - 2.796*e(k) + 0.0476*e(k-1);
    y(k+1) = 1.81*y(k) - 0.8187*y(k-1) + 0.009358*u(k) + 0.008754*u(k-1);
end

%% Saída - Eq. a diferenças

figure(3)
hold on;
stairs(t,y,'r');
grid on;
hold off;

%% Sinal de controle - Eq. a diferenças

figure(5)
hold on;
stairs(t,u,'r')
hold off;
grid on;






