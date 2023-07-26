function [poly_rd] = R_d(esp)
   
    ee = [1 2 4 6 8 10 12];
    rd = [1100 1200 1400 1750 1850 2000 2020];

%     figure
%     grid on;
%     grid minor;
%     plot(ee, rd, '.k', 'MarkerSize', 8);
%     xlabel('esp');
%     ylabel('rd')
%     hold on

for i = 1:6
    poly = polyfit(ee, rd, i);
    poly_rd= polyval(poly, ee);
    [m(i), k(i)] = max(abs(poly_rd - rd));
end

poly = polyfit(ee, rd, 2);

poly_rd = (polyval(poly, esp))';
% plot(ee, poly_rd,  'green');
% legend('data', 'interpolation');


end