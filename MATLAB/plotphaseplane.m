function [ output_args ] = plotphaseplane( theta, thetadot, c1, c2, deadband, u )
%PLOTPHASEPLANE Summary of this function goes here
%   Detailed explanation goes here

figure;
hold on;
thrust = max(abs(u));
k1 = 1 / (2*thrust);
a3 = -deadband/2;
a6 = deadband/2;
a2 = a3 - 0.1;
a1 = -k1 * c2^2 + a2;
a4 = k1 * c1^2 + a3;
a5 = -k1 * c1^2 + a6;
a7 = a6 + 0.1;
a8 = k1 * c2^2 + a7;
plot([-20 a5], [c1 c1], 'k');
plot(-k1 .* (c1:-0.0001:0).^2 + a6, c1:-0.0001:0, 'k');
plot([a6 a7], [0 0], 'k');
plot(k1 .* (0:-0.0001:-c2).^2 + a7, 0:-0.0001:-c2, 'k');
plot([a8 20], [-c2 -c2], 'k');
plot([-20 a1], [c2 c2], 'k');
plot(-k1 .* (c2:-0.0001:0).^2 + a2, c2:-0.0001:0, 'k');
plot([a2 a3], [0 0], 'k');
plot(k1 .* (0:-0.0001:-c1).^2 + a3, 0:-0.0001:-c1, 'k');
plot([a4 20], [-c1 -c1], 'k');

ii = find(u);
plot(theta(ii), thetadot(ii), 'r.');
jj = find(~u);
plot(theta(jj), thetadot(jj), 'b.');

xl = xlim;
yl = ylim;
xlim([-1 1] * max(max(abs(xl)), 2*a8));
ylim([-1 1] * max(max(abs(yl)), 2*c1));
xl = xlim;
yl = ylim;
plot(xl, [0 0], 'k');
plot([0 0], yl, 'k');
xlabel('$\theta_e$', 'Interpreter', 'latex');
ylabel('$\dot{\theta}_e$', 'Interpreter', 'latex');
% set(gca, 'XAxisLocation', 'origin');
% set(gca, 'YAxisLocation', 'origin');
xticks([a3 0 a6]);
yticks([-c1 -c2 0 c2 c1]);

end

