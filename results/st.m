%% Stability calculations

MU_0 = 1.2566371e-6;
EPS_0 = 8.85418781762e-12;
C_CONST = 299792458.0;
dz = 0.01e-6;
alpha = 0.99
dt = 2^(1-1/alpha) * (sqrt(EPS_0*MU_0) * dz)^(1/alpha)
classic_dt = sqrt(EPS_0*MU_0)*dz



%% Plot xi pseudo-polynomial 
% -------------------------------
figure(1)
u = [0:0.01:2];
coeff_2a_xi = [1.8, 2, 2.5];
y = ones(size(coeff_2a_xi, 2)+1,size(u, 2));
y(1, :) = u +1;

beta = 0.8;

hold on;
grid on;
plot(u, y(1, :), 'LineWidth', 1.2);
my_legend=cell(size(coeff_2a_xi, 2)+1,1);
my_legend{1} = 'u+1';
% legend('u+1') 
for i = 1:size(coeff_2a_xi, 2)
    y(i+1, :) = coeff_2a_xi(i) .* (u).^beta;
    plot(u, y(i+1, :));
    my_legend{i+1} = strcat(num2str(coeff_2a_xi(i)), '*u^{', num2str(beta), '}');
end
yl = ylim; % current y-axis limits
plot([1 1],[yl(1) yl(2)], '--black');
xlabel('u')
legend(my_legend)

% Plot intersection points
for i = 1:size(coeff_2a_xi, 2)
    eps  = 0.001;
    idx = find(y(1, :) - y(i+1, :) < eps, 1); %// Index of coordinate in array
    px = u(idx);
    py = y(1, idx);
    plot(px,py,'.','MarkerSize',18, 'Color', 'red', 'HandleVisibility','off');
end
% -------------------------------

%% Plot smallest unstable dt with analytical stability boundary
figure(2)
alpha_arr_num = [9.950000e-001, 9.900000e-001, 9.800000e-001, 9.700000e-001, 9.500000e-001, 9.000000e-001, 8.500000e-001, 8.000000e-001, 7.500000e-001, 7.000000e-001, 6.500000e-001, 6.000000e-001, 5.500000e-001, 5.100000e-001];
unstable_dt = [5.540861e-017, 4.570280e-017, 3.091105e-017, 2.073877e-017, 9.102785e-018, 9.899597e-019, 8.292799e-020, 5.095305e-021, 2.158250e-022, 5.819434e-024, 8.999972e-026, 6.947495e-028, 2.214829e-030, 9.901580e-033];

alpha_arr_an = alpha_arr_num(1):(alpha_arr_num(end)-alpha_arr_num(1))/200:alpha_arr_num(end);
dt_arr_an = 2.^(1.0-1./alpha_arr_an) .* (sqrt(EPS_0*MU_0) * dz).^(1./alpha_arr_an);

plot(alpha_arr_an, dt_arr_an)
hold on
scatter(alpha_arr_num, unstable_dt)
set(gca,'yscale','log')
legend(["granica analityczna", "wyniki numeryczne"])
xlabel("\alpha")
ylabel("\Delta" + "t [s]")
title(sprintf("Granica stabilnoœci, dz=%e m", dz))
grid on
