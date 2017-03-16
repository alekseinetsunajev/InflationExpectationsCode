xl = {'91','92','93','94','95','96','97', '98','99','00','01','02','03','04', '05','06','07','08','09','10','11', '12','13','14','15'};
figure1 = figure('Name','Data', 'NumberTitle','off', 'PaperType','<custom>','PaperSize',[8 5]);
h(1) = plot(y(:, 1),'k:',  'LineWidth', 1.5, 'DisplayName','\pi_t^{e, s}');
ax(1) = gca;
set(ax(1),  'XTickLabel', xl, 'XTick',[9:12:308 ], 'fontsize', 12);
hold on
plot(3.06*ones(309 , 1),'k:',  'LineWidth', 1);
xlim([0 308])
h(2) = plot(y(:, 2),'k',  'LineWidth', 1.5, 'DisplayName','\pi_t^{e, l}');
legend([h(1); h(2)]);
