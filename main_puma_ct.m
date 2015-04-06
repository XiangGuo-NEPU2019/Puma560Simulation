function main_puma_ct 
[t,x] = ode45(@puma_ct,0:0.001:20,[0 0 0 0 0 0]);

figure(1),
plot(t,x(:,1),'-',t,x(:,2),'--',t,x(:,3),'-.','linewidth',2.5);
legend('q1','q2','q3');
hold on
plot(t,pi/2+sin(t),'--');
xlabel('t/s'),ylabel('q/rad');
hold off 
figure(2),
plot(t,x(:,4),'-',t,x(:,5),'--',t,x(:,6),'-.','linewidth',2.5);
legend('w1','w2','w3');
hold on
plot(t,cos(t),'--');
xlabel('t/s'),ylabel('w/rad/s');
hold off
end

function xdot = puma_ct(t,q)
xdot = zeros(6,1);
qd = [pi/2+sin(t); pi/2+sin(t); pi/2+sin(t)];
qddot = [cos(t); cos(t); cos(t)];
qd2dot = -[sin(t); sin(t); sin(t)];

A = [500 0 0; 0 180 0; 0 0 50];
B = [350 0 0; 0 50 0; 0 0 20];

xdot(1) = q(4);
xdot(2) = q(5);
xdot(3) = q(6);

xdot(4:6) = qd2dot - A*(q(4:6)-qddot) -B*(q(1:3)-qd);
   
end