clc; clear;
load('Pvsrequire_long1.mat') 

A4=[mean(P_Sch4(1,:,1)) mean(P_Sch4(2,:,1)) mean(P_Sch4(3,:,1)) mean(P_Sch4(4,:,1)) mean(P_Sch4(5,:,1))];
A3=[mean(P_Sch3(1,:,end)) mean(P_Sch3(2,:,end)) mean(P_Sch3(3,:,end)) mean(P_Sch3(4,:,end)) mean(P_Sch3(5,:,end))];
A2=[mean(P_Sch2(1,:,end)) mean(P_Sch2(2,:,end)) mean(P_Sch2(3,:,end)) mean(P_Sch2(4,:,end)) mean(P_Sch2(5,:,end))];
A1=[mean(P_Sch1(1,:,end)) mean(P_Sch1(2,:,end)) mean(P_Sch1(3,:,end)) mean(P_Sch1(4,:,end)) mean(P_Sch1(5,:,end))];
A5=[mean(P_Sch5(1,:,end)) mean(P_Sch5(2,:,end)) mean(P_Sch5(3,:,end)) mean(P_Sch5(4,:,end)) mean(P_Sch5(5,:,end))];

plot(1:5,10*log10(A1), '-o', 'linewidth', 1.5);
hold on
plot(10*log10((A2)), '-s', 'linewidth', 1.5);
% plot(10*log10(A3), '-<', 'linewidth', 1.5);
plot(10*log10(A5), '-<', 'linewidth', 1.5);
plot(10*log10(A4), '-p', 'linewidth', 1.5);
grid on;
axis([1 5 -20 12]);
legend('Proposed Rot-Pol RA + DBF','Rotation-Only RA + DBF','Boresight-Only RA + DBF','Fixed UPA + DBF','FontSize',12)
xlabel('Rate requirement, ${\bar R}$','Interpreter','latex');
ylabel('Transmit power (dBm)');

%zheli p=2;