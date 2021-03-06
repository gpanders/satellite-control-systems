% load mc.mat

Nt = length(tVec);

biasErrorAvg = zeros(Nt, 3);
biasErrorAvg(:, 1) = mean(squeeze(biasEstError(:, 1, :)), 2);
biasErrorAvg(:, 2) = mean(squeeze(biasEstError(:, 2, :)), 2);
biasErrorAvg(:, 3) = mean(squeeze(biasEstError(:, 3, :)), 2);

figure;
for k = 1:3
    subplot(3, 1, k), ...
    plot(tVec, rad2deg(biasErrorAvg(:, k))), ...
    legend(['\beta_' num2str(k)]), ...
    xlabel('seconds'), ...
    ylabel('deg / sec'), ...
    xlim([tVec(1) tVec(end)]);
end

attErrorAvg = zeros(Nt, 3);
attErrorAvg(:, 1) = mean(squeeze(attError(:, 1, :)), 2);
attErrorAvg(:, 2) = mean(squeeze(attError(:, 2, :)), 2);
attErrorAvg(:, 3) = mean(squeeze(attError(:, 3, :)), 2);

figure;
for k = 1:3
    subplot(3, 1, k), ...
    plot(tVec(2:end), rad2deg(attErrorAvg(2:end, k))), ...
    legend(['\theta_' num2str(k)]), ...
    xlabel('seconds'), ...
    ylabel('degrees'), ...
    xlim([tVec(1) tVec(end)]);
end