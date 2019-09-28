m = 0.100;
d = 0.025;
I = m*d^2;

A = [-bx/m 0;
    0 -bth/I];
B = [1 1;1 -1];
C = [1 1];
D = [0 0];
% state_space = ss(A,B,C,D);
[G1_num G1_den] = ss2tf(A,B,C,D,1)
G1 = tf(G1_num,G1_den)
[G2_num G2_den] = ss2tf(A,B,C,D,2)
G2 = tf(G2_num,G2_den)

figure
bode(G1)
hold on
bode(G2)
legend('Axial','Rotational')

figure
pzplot(G1)
figure
pzplot(G2)
