clear all; close all; clc

% s = tf('s');
% sys = (s+2)/((s*(s+10)*(s^2+2*s+2)));
% bode(sys);
% 
% title("Magnitude vs. Angular velocity");
% xlabel("Angular Velocity (rad/s)")
% ylabel("Phase (Theta)")
% grid on;



% w = 0.1:1e3;
% 
% H = (1i.*w+2)./((1i.*w).*(1i.*w+10).*(((1i.*w).^2)+2*1i.*w+2));
% % H = 10*(1i.*w+4)./((1i.*w).*(1i.*w+1).*(1i.*w+200));
% amplitude = abs(H);
% phase = angle(H)*(180/pi);
% 
% figure(1)
% loglog(w,amplitude,'LineWidth',1);
% grid on;
% title("Magnitude vs. Angular velocity");
% xlabel("Angular Velocity (rad/s)")
% ylabel("Magnitude")
% 
% figure(3)
% semilogx(w,phase,'LineWidth',1);
% title("Phase vs. Angular Velocity");
% xlabel("Angular Velocity (rad/s)")
% ylabel("Phase (Theta)")
% grid on;
