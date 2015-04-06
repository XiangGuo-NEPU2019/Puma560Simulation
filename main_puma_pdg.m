function main_puma_pdg 
[t,x] = ode45(@puma_pdg,0:0.001:20,[0 0 0 0 0 0]);

plot(t,x(:,1),'-',t,x(:,2),'--',t,x(:,3),'-.','linewidth',2);
legend('q1','q2','q3');
hold on
plot(t,pi*ones(1,length(t))/2,'--');
xlabel('t/s'),ylabel('q/rad');
hold off 
end

function xdot = puma_pdg(~,q)
xdot = zeros(6,1);
qd = [pi/2; pi/2; pi/2];

A = [500 0 0; 0 180 0; 0 0 50];
B = [350 0 0; 0 50 0; 0 0 20];
tau = -A*([q(1);q(2);q(3)]-[qd(1);qd(2);qd(3)])-B*[q(4);q(5);q(6)];

xdot(1) = q(4);
xdot(2) = q(5);
xdot(3) = q(6);

H = [22+0.9*cos(q(2))*cos(q(2)), 1.17+1.92*sin(q(2)), -0.3*cos(q(2)+q(3));
    1.17+1.92*sin(q(2)), 1.66, -0.29;
    -0.3*cos(q(2)+q(3)), -0.29, 0.11];
xdot(4:6) = H\(-...
    [-1.8*cos(q(2))*sin(q(2))*q(4)*q(5)+1.92*cos(q(2))*q(5)*q(5)+0.3*sin(q(2)+q(3))*(q(5)+q(6))*q(6);
    0.9*sin(q(2))*cos(q(2))*q(4)*q(4)-0.3*sin(q(2)+q(3))*q(4)*q(6);
    0.3*sin(q(2)+q(3))*q(4)*q(5)]+tau);
    
end