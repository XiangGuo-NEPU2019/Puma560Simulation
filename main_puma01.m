function main_puma01

tau1 = zeros(1,20001);
tau2 = zeros(1,20001);
tau3 = zeros(1,20001);

for i = 0:1:20000;
    
    t = i*0.001;
    q = sin(t);
    qdot = cos(t);
    qddot = -sin(t);

    cq = cos(q);
    sq = sin(q);
    c2q = cos(q+q);
    s2q = sin(q+q);

    R01 = [cq,-sq,0;sq,cq,0;0,0,1];
    R02 = [cq*cq,-cq*sq,-sq;sq*cq,-sq*sq,cq;-sq,-cq,0];
    R03 = [cq*c2q,-cq*s2q,-sq;sq*c2q,-sq*s2q,cq;-s2q,-c2q,0];

    z1 = [0;0;1];
    z2 = [-sq;cq;0];
    z3 = z2;

    s1 = 0.243*z2;
    s2 = R02*[0.432;0;-0.093];
    r2 = R02*[0.068;0.006;-0.016];
    r3 = R03*[0;-0.07;0.014];

    w1 = [0;0;qdot];
    w1dot = [0;0;qddot];
    w2 = qdot*(z1+z2);
    w2dot = qddot*(z1+z2)+qdot*qdot*[-cq;-sq;0];
    w3 = qdot*(z1+2*z2);
    w3dot = qddot*(z1+2*z2)+qdot*qdot*[-2*cq;-2*sq;0];

    v2 = cross(w1,s1);
    v3 = v2+cross(w2,s2);

    a2 = cross(w1dot,s1)+cross(w1,cross(w1,s1));
    a3 = a2+cross(w2dot,s2)+cross(w2,cross(w2,s2));

    ac2 = a2+cross(w2dot,r2)+cross(w2,cross(w2,r2));
    ac3 = a3+cross(w3dot,r3)+cross(w3,cross(w3,r3));

    m2 = 17.4;
    m3 = 4.8;
    I1 = 0.35;
    I2 = R02*diag([0.13,0.542,0.539])*R02';
    I3 = R03*diag([0.066,0.0125,0.086])*R03';
    g = 9.8;

    f3 = m3*ac3+[0;0;m3*g];
    n3 = I3*w3dot+cross(w3,I3*w3)+cross(r3,f3);
    t3 = z3'*n3;

    f2 = m2*ac2+f3+[0;0;m2*g];
    n2 = I2*w2dot+cross(w2,I2*w2)+cross(r2,f2)+n3+cross((s2-r2),f3);
    t2 = z2'*n2;

    n1 = I1*w1+n2+cross(s1,f2);
    t1 = z1'*n1;

    tau1(i+1) = t1;
    tau2(i+1) = t2;
    tau3(i+1) = t3;

end

t = 0:0.001:20;
plot(t,tau1,'-',t,tau2,'--',t,tau3,'-.');
xlabel('t/s'),ylabel('\tau /Nm');
legend('\tau 1','\tau 2','\tau 3');
