L = 10e-6


fig_gauss_propag;
close;
figure()
h = plot(t4, output_signal4, "--", 'Color', [0, 0.4470, 0.7410]);
hold on
dz = 0.02e-6;
[field_at_z_002, t_arr_002] = FieldAtPositions('Ex_dz_002.bin', [int32(0.1e-6/dz)+int32(L/dz)]);
dz = 0.01e-6;
[field_at_z_001, t_arr_001] = FieldAtPositions('Ex_dz_001.bin', [int32(0.1e-6/dz)+int32(L/dz)]);
dz = 0.005e-6;
[field_at_z_0005, t_arr_0005] = FieldAtPositions('Ex_dz_0005.bin', [int32(0.1e-6/dz)+int32(L/dz)]);


plot(t_arr_002, field_at_z_002, 'Color', [0.4940, 0.1840, 0.5560]);
plot(t_arr_001, field_at_z_001, 'Color', [0.9290, 0.6940, 0.1250]);
plot(t_arr_0005, field_at_z_0005, 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('t [s]');
ylabel("Ex [V/m]");
xlim([t_arr_002(1), t_arr_002(end)]);
title(sprintf('Pole Ex w punkcie %e m', 0.1e-6+L));
grid on;
legend(["Algorytm w dziedzinie \newline"+"czêstotliwoœci", "FDTD, "+"\Delta"+"z=0.02"+"\mu"+"m", ...
        "FDTD, "+"\Delta"+"z=0.01"+"\mu"+"m", "FDTD, "+"\Delta"+"z=0.005"+"\mu"+"m"])

    
uistack(h, 'top');

