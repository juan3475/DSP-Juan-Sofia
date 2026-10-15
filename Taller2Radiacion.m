% =========================================================================
% SIMULACIÓN DE DIPOLO DE MEDIA ONDA EN BANDA ISM DE 5.8 GHz
% Antenna Toolbox - MATLAB
% =========================================================================

clear; clc; close all;

%% 1. Parámetros de Diseño
c = 3e8;                          % Velocidad de la luz [m/s]
f = 5800e6;                       % Frecuencia central: 5800 MHz (5.8 GHz)
lam = c / f;                      % Longitud de onda: lambda = 0.0517 m (51.7 mm)
L = 0.47 * lam                   % Longitud total del dipolo: L = 24.8 mm
ancho = 0.001;                    % Ancho equivalente (hilo de radio r = 0.25 mm)

%% 2. Creación del Objeto Dipolo
d = dipole('Length', L, 'Width', ancho);

% --- FIGURA 1: Geometría 3D de la Antena ---
figure(1);
set(gcf, 'Color', 'w');
show(d);
title('Geometría 3D del Dipolo de Media Onda (L = 24.8 mm)', ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje X (m)');
ylabel('Eje Y (m)');
zlabel('Eje Z (m)');

% --- FIGURA 2: Mallado Numérico (MoM) ---
figure(2);
set(gcf, 'Color', 'w');
mesh(d);
title('Mallado Numérico de la Antena (Método de los Momentos)', ...
      'FontSize', 12, 'FontWeight', 'bold');
xlabel('Eje X (m)');
ylabel('Eje Y (m)');
zlabel('Eje Z (m)');

%% 3. Análisis de Frecuencia: S11 y SWR (VSWR)

f_rango = linspace(0.8*f, 1.2*f, 101); % Rango de 4640 MHz a 6960 MHz

s = sparameters(d, f_rango);

s11_mag = abs(squeeze(s.Parameters(1,1,:)));
s11_db = 20*log10(s11_mag);

swr_rango = (1 + s11_mag) ./ (1 - s11_mag);

% --- FIGURA 3: Coeficiente de Reflexión S11 ---
figure(3);
set(gcf, 'Color', 'w');

plot(f_rango/1e6, s11_db, 'b-', 'LineWidth', 2);

grid on;
grid minor;

title('Coeficiente de Reflexión S_{11} vs. Frecuencia', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

xlabel('Frecuencia (MHz)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylabel('Pérdidas por Retorno S_{11} (dB)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

xline(f/1e6, '--r', ...
      'Frecuencia Central (5800 MHz)', ...
      'LineWidth', 1.2, ...
      'LabelOrientation', 'aligned');

ylim([-35 0]);

% --- FIGURA 4: Relación de Onda Estacionaria (SWR / VSWR) ---
figure(4);
set(gcf, 'Color', 'w');

plot(f_rango/1e6, swr_rango, 'b-', 'LineWidth', 2);

hold on;
grid on;
grid minor;

yline(2.0, '--m', ...
      'Límite Aceptable (SWR = 2.0)', ...
      'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left');

xline(f/1e6, '--r', ...
      'Frecuencia Central (5800 MHz)', ...
      'LineWidth', 1.2, ...
      'LabelOrientation', 'aligned');

title('Relación de Onda Estacionaria (SWR) vs. Frecuencia', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

xlabel('Frecuencia (MHz)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylabel('SWR (Referencia 50 \Omega)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylim([1 5]);

xlim([min(f_rango)/1e6 max(f_rango)/1e6]);

swr_f0 = vswr(d, f, 50);

text(f/1e6 + 50, ...
     swr_f0 + 0.3, ...
     sprintf('SWR a 5.8 GHz = %.2f', swr_f0), ...
     'FontSize', 10, ...
     'FontWeight', 'bold', ...
     'Color', [0 0.5 0]);

%% 4. Patrones de Radiación (3D, Plano E y Plano H)

% --- FIGURA 5: Patrón Tridimensional 3D ---
figure(5);
set(gcf, 'Color', 'w');

pattern(d, f);

title('Patrón de Radiación 3D - Directividad (dBi) a 5.8 GHz', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

% --- FIGURA 6: Corte Plano E (Elevación / Plano XZ) ---
figure(6);
set(gcf, 'Color', 'w');

patternElevation(d, f);

title('Patrón de Radiación - Plano E (Elevación / Plano XZ)', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

% --- FIGURA 7: Corte Plano H (Azimut / Plano XY) ---
figure(7);
set(gcf, 'Color', 'w');

patternAzimuth(d, f, 0);

title('Patrón de Radiación - Plano H (Azimut / Plano XY)', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

%% 5. Impedancia de Entrada de la Antena

Zin_rango = impedance(d, f_rango);

% --- FIGURA 8: Impedancia de Entrada vs. Frecuencia ---
figure(8);
set(gcf, 'Color', 'w');

plot(f_rango/1e6, ...
     real(Zin_rango), ...
     'b-', ...
     'LineWidth', 2);

hold on;

plot(f_rango/1e6, ...
     imag(Zin_rango), ...
     'r-', ...
     'LineWidth', 2);

grid on;
grid minor;

xline(f/1e6, '--k', ...
      'Frecuencia Central (5800 MHz)', ...
      'LineWidth', 1.2, ...
      'LabelOrientation', 'aligned');

yline(50, '--m', ...
      'Referencia 50 \Omega', ...
      'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left');

yline(0, '--k', ...
      'X = 0 \Omega', ...
      'LineWidth', 1);

title('Impedancia de Entrada vs. Frecuencia', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

xlabel('Frecuencia (MHz)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylabel('Impedancia (\Omega)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

legend('Resistencia R', ...
       'Reactancia X', ...
       'Location', 'best');

xlim([min(f_rango)/1e6 max(f_rango)/1e6]);

%% 6. Variación de Longitud del Dipolo

% Longitudes que se van a evaluar
% Se varía desde 0.40 lambda hasta 0.56 lambda

factores_L = 0.40:0.02:0.56;

longitudes = factores_L * lam;

% Rango de frecuencia utilizado para buscar la resonancia
f_variacion = linspace(4.0e9, 7.5e9, 151);

% Vector para almacenar las frecuencias de resonancia
f_resonancia = zeros(size(longitudes));

% Matriz para almacenar la reactancia
X_variacion = zeros(length(longitudes), length(f_variacion));

%% 7. Simulación para Cada Longitud

for k = 1:length(longitudes)

    % Longitud actual
    L_actual = longitudes(k);

    % Crear dipolo con la nueva longitud
    d_var = dipole('Length', L_actual, ...
                   'Width', ancho);

    % Calcular impedancia para el rango de frecuencias
    Z_var = impedance(d_var, f_variacion);

    % Extraer reactancia
    X_var = imag(Z_var);

    % Guardar reactancia
    X_variacion(k,:) = X_var;

    % Buscar el punto donde la reactancia es más cercana a cero
    [~, indice_res] = min(abs(X_var));

    % Guardar frecuencia de resonancia
    f_resonancia(k) = f_variacion(indice_res);

end

%% 8. Reactancia para Diferentes Longitudes

% --- FIGURA 9: Reactancia vs. Frecuencia ---
figure(9);
set(gcf, 'Color', 'w');

hold on;

for k = 1:length(longitudes)

    plot(f_variacion/1e9, ...
         X_variacion(k,:), ...
         'LineWidth', 1.5);

end

yline(0, '--k', ...
      'Condición de Resonancia X = 0 \Omega', ...
      'LineWidth', 1.5);

grid on;
grid minor;

title('Variación de la Frecuencia de Resonancia con la Longitud del Dipolo', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

xlabel('Frecuencia (GHz)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylabel('Reactancia X (\Omega)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

% Crear leyendas automáticamente
leyendas = strings(1,length(longitudes));

for k = 1:length(longitudes)

    leyendas(k) = sprintf('L = %.2f mm', ...
                          longitudes(k)*1000);

end

legend(leyendas, ...
       'Location', 'best');

xlim([min(f_variacion)/1e9 ...
      max(f_variacion)/1e9]);

%% 9. Frecuencia de Resonancia vs. Longitud

% --- FIGURA 10 ---
figure(10);
set(gcf, 'Color', 'w');

plot(longitudes*1000, ...
     f_resonancia/1e9, ...
     'o-', ...
     'LineWidth', 2, ...
     'MarkerSize', 7);

grid on;
grid minor;

title('Frecuencia de Resonancia vs. Longitud del Dipolo', ...
      'FontSize', 12, ...
      'FontWeight', 'bold');

xlabel('Longitud total del dipolo (mm)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylabel('Frecuencia de Resonancia (GHz)', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

% Marcar longitud del diseño original
xline(L*1000, '--r', ...
      'Longitud Diseño Original', ...
      'LineWidth', 1.2);

% Marcar frecuencia objetivo
yline(f/1e9, '--m', ...
      'Frecuencia Objetivo 5.8 GHz', ...
      'LineWidth', 1.2);

%% 10. Tabla de Resultados de Resonancia

fprintf('\n');
fprintf('============================================================\n');
fprintf('       VARIACIÓN DE LONGITUD Y FRECUENCIA DE RESONANCIA\n');
fprintf('============================================================\n');
fprintf(' Longitud (mm)     L/lambda      Frecuencia Resonancia (GHz)\n');
fprintf('------------------------------------------------------------\n');

for k = 1:length(longitudes)

    fprintf('    %6.2f          %.2f                %.3f\n', ...
            longitudes(k)*1000, ...
            factores_L(k), ...
            f_resonancia(k)/1e9);

end

fprintf('============================================================\n');

%% 11. Procesar Figuras

drawnow;

%% 12. Resultados e Impresión en Consola

Zin = impedance(d, f);

D = directivity(d, f);

G = gain(d, f);

bw = beamwidth(d, f);

fprintf('\n');
fprintf('===================================================\n');
fprintf('        RESULTADOS DIPOLO 5.8 GHz (BANDA 5G/ISM)   \n');
fprintf('===================================================\n');

fprintf('Longitud total (L)       : %.2f mm\n', ...
        L*1000);

fprintf('Radio equivalente (r)    : %.2f mm\n', ...
        (ancho/4)*1000);

fprintf('Impedancia de entrada    : Z_in = %.2f %+.2fj Ohm\n', ...
        real(Zin), ...
        imag(Zin));

fprintf('VSWR a 5800 MHz (50 Ohm) : %.2f\n', ...
        swr_f0);

fprintf('Directividad Máxima      : %.2f dBi\n', ...
        D);

fprintf('Ganancia Máxima          : %.2f dBi\n', ...
        G);

fprintf('Ancho de Haz (Plano H)   : %.1f°\n', ...
        bw(1));

fprintf('Ancho de Haz (Plano E)   : %.1f°\n', ...
        bw(2));

fprintf('===================================================\n');