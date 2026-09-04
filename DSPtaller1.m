%%

% TALLER 1:PROCESAMIENTO DIGITAL DE SEÑALES

% Juan Esteban Rodriguez Villamil
% 30000122719
% Ana Sofia Rodriguez Hincapie
% 30000119535

% UNIVERSIDAD DE SAN BUENVANETURA

% OSCAR DIAZ

% BOGOTA D.C

% 2026

%%
% Serie de Fourier de una onda en diente de sierra o triangular.
% Defina una una señal periódica en diente de sierra o triangular (con coeficientes analíticos conocidos, distinta de la onda cuadrada). 
% Calcule los coeficientes analíticos para k = 0, . . . , 10 e implemente el cálculo numérico. Grafique el espectro bilateral de magnitud
% y de fase y presente la comparación analítico vs numérico en una tabla. Comente cómo decrecen los coeficientes con k y en qué se diferencia 
% este espectro del de la onda cuadrada.



clc;
clear;
close all;

%% ============================================================
%              SERIE DE FOURIER - ONDA DIENTE DE SIERRA
%                    (VERSIÓN CORREGIDA)
% ============================================================

%% ============================================================
%              1. PARÁMETROS DE LA SEÑAL
% ============================================================

A = 1;
T = 2;

% Frecuencia fundamental
f0 = 1/T;

% Frecuencia angular fundamental
w0 = 2*pi*f0;


%% ============================================================
%              2. DEFINICIÓN DE LA SEÑAL
% ============================================================

N = 10000;

% Tiempo para varios periodos
t = linspace(-4, 4, N);

% Señal diente de sierra periódica
x = A * (2*mod(t + T/2, T)/T - 1);


%% ============================================================
%              3. VALORES DE k
% ============================================================

k = 0:10;


%% ============================================================
%              4. COEFICIENTES ANALÍTICOS
% ============================================================

a_analitico = zeros(size(k));
b_analitico = zeros(size(k));

for n = 1:length(k)

    if k(n) == 0

        % Componente DC
        a_analitico(n) = 0;
        b_analitico(n) = 0;

    else

        % La señal es impar:
        % por lo tanto, a_k = 0
        a_analitico(n) = 0;

        % Coeficiente b_k
        b_analitico(n) = ...
            2*(-1)^(k(n)+1)/(k(n)*pi);

    end

end


%% ============================================================
%              5. CÁLCULO NUMÉRICO
% ============================================================

% Para calcular correctamente los coeficientes se utiliza
% solamente un periodo de la señal

t_periodo = linspace(-T/2, T/2, N);

x_periodo = t_periodo;

a_numerico = zeros(size(k));
b_numerico = zeros(size(k));

for n = 1:length(k)

    kk = k(n);

    % Coeficiente a_k
    a_numerico(n) = ...
        (2/T)*trapz(t_periodo, ...
        x_periodo .* cos(kk*w0*t_periodo));

    % Coeficiente b_k
    b_numerico(n) = ...
        (2/T)*trapz(t_periodo, ...
        x_periodo .* sin(kk*w0*t_periodo));

end


%% ============================================================
%              6. ELIMINAR ERRORES NUMÉRICOS MUY PEQUEÑOS
% ============================================================

a_numerico(abs(a_numerico) < 1e-10) = 0;
b_numerico(abs(b_numerico) < 1e-10) = 0;


%% ============================================================
%              7. TABLA 1: COEFICIENTES DE FOURIER
% ============================================================

% --------------------------------------------------------------
% Error cuadrático punto a punto: (analítico - numérico)^2
% --------------------------------------------------------------

error2_a = (a_analitico - a_numerico).^2;
error2_b = (b_analitico - b_numerico).^2;

Tabla = table( ...
    k.', ...
    a_analitico.', ...
    a_numerico.', ...
    error2_a.', ...
    b_analitico.', ...
    b_numerico.', ...
    error2_b.', ...
    'VariableNames', ...
    {'k','a_Analitico','a_Numerico','Error2_a',...
    'b_Analitico','b_Numerico','Error2_b'});

disp(' ');
disp('==============================================================');
disp('TABLA 1. COMPARACIÓN DE COEFICIENTES DE FOURIER');
disp('==============================================================');

disp(Tabla);

% --------------------------------------------------------------
% Error cuadrático medio (MSE) global
%
% NOTA: la función immse() requiere el Image Processing Toolbox.
% Como no está disponible, se calcula manualmente el MSE, que es
% exactamente lo mismo que hace immse por dentro: el promedio de
% los errores cuadráticos punto a punto.
% --------------------------------------------------------------

MSE_a = mean(error2_a);
MSE_b = mean(error2_b);

fprintf('\n');
fprintf('MSE de a_k (analítico vs numérico): %.4e\n', MSE_a);
fprintf('MSE de b_k (analítico vs numérico): %.4e\n', MSE_b);


%% ============================================================
%              8. GRÁFICA DE LA SEÑAL DIENTE DE SIERRA
% ============================================================

figure;

plot(t, x, 'LineWidth', 1.5);

grid on;

xlabel('Tiempo (s)');
ylabel('Amplitud');

title('Señal diente de sierra');

xlim([-4 4]);
ylim([-1.2 1.2]);


%% ============================================================
%              9. ESPECTRO BILATERAL
% ============================================================

% Para el espectro bilateral necesitamos valores
% positivos y negativos de k

k_bilateral = -10:10;


%% ============================================================
%              10. INICIALIZAR COEFICIENTES COMPLEJOS
% ============================================================

C_analitico = zeros(size(k_bilateral));
C_numerico = zeros(size(k_bilateral));


%% ============================================================
%      11. CALCULAR COEFICIENTES COMPLEJOS (CORREGIDO)
% ============================================================
%
% CORRECCIÓN CLAVE:
% Antes, el bucle recorría k_bilateral en orden (-10,...,10) y
% al llegar a un k negativo intentaba reflejar un valor positivo
% que AÚN no había sido calculado (seguía en 0 por inicialización).
% Por eso todo el lado negativo salía en cero.
%
% La solución: separar el cálculo en dos pasadas.
%   PASO A: calcular primero TODOS los k >= 0 (positivos y DC).
%   PASO B: calcular después los k < 0 usando la simetría
%           conjugada C(-k) = conj(C(k)), ya con los positivos
%           completamente calculados.

% --------------------------------------------------------------
% PASO A: k = 0 y k > 0
% --------------------------------------------------------------

for n = 1:length(k_bilateral)

    kk = k_bilateral(n);

    if kk == 0

        % Componente DC
        C_analitico(n) = 0;
        C_numerico(n) = 0;

    elseif kk > 0

        % Índice correspondiente a k positivo dentro de
        % los arreglos a_analitico, b_analitico (k = 0:10)
        indice = kk + 1;

        % Coeficiente complejo analítico
        %
        % C_k = (a_k - j*b_k)/2

        C_analitico(n) = ...
            (a_analitico(indice) - ...
            1i*b_analitico(indice))/2;

        % Coeficiente complejo numérico

        C_numerico(n) = ...
            (a_numerico(indice) - ...
            1i*b_numerico(indice))/2;

    end

end

% --------------------------------------------------------------
% PASO B: k < 0, usando simetría conjugada
%          C(-k) = conj(C(k))
% --------------------------------------------------------------

for n = 1:length(k_bilateral)

    kk = k_bilateral(n);

    if kk < 0

        % Buscamos la posición dentro de k_bilateral que
        % corresponde a +kk (ya calculada en el PASO A)
        indice_positivo = find(k_bilateral == -kk);

        C_analitico(n) = conj(C_analitico(indice_positivo));

        C_numerico(n) = conj(C_numerico(indice_positivo));

    end

end


%% ============================================================
%              12. FRECUENCIAS DEL ESPECTRO
% ============================================================

f_bilateral = k_bilateral * f0;


%% ============================================================
%              13. MAGNITUD DEL ESPECTRO
% ============================================================

magnitud_analitico = abs(C_analitico);

magnitud_numerico = abs(C_numerico);


%% ============================================================
%              14. FASE DEL ESPECTRO
% ============================================================

fase_analitico = angle(C_analitico);

fase_numerico = angle(C_numerico);


%% ============================================================
%              15. CORREGIR FASE CUANDO C_k = 0
% ============================================================

% Cuando la magnitud es cero, la fase no está definida.
% Para mostrar la gráfica de forma limpia colocamos 0.

fase_analitico(magnitud_analitico < 1e-10) = 0;

fase_numerico(magnitud_numerico < 1e-10) = 0;


%% ============================================================
%              16. ESPECTRO BILATERAL DE MAGNITUD
% ============================================================

figure;

stem(f_bilateral, magnitud_analitico, ...
    'filled');

hold on;

stem(f_bilateral, magnitud_numerico);

grid on;

xlabel('Frecuencia (Hz)');
ylabel('|C_k|');

title('Espectro bilateral de magnitud');

legend('Analítico','Numérico');

xlim([min(f_bilateral) max(f_bilateral)]);


%% ============================================================
%              17. ESPECTRO BILATERAL DE FASE
% ============================================================

figure;

stem(f_bilateral, fase_analitico, ...
    'filled');

hold on;

stem(f_bilateral, fase_numerico);

grid on;

xlabel('Frecuencia (Hz)');
ylabel('Fase (rad)');

title('Espectro bilateral de fase');

legend('Analítico','Numérico');

xlim([min(f_bilateral) max(f_bilateral)]);


%% ============================================================
%              18. TABLA 2: ESPECTRO BILATERAL
% ============================================================

% --------------------------------------------------------------
% Error cuadrático punto a punto: (analítico - numérico)^2
% --------------------------------------------------------------

error2_magnitud = (magnitud_analitico - magnitud_numerico).^2;
error2_fase = (fase_analitico - fase_numerico).^2;

Tabla_Espectro = table( ...
    k_bilateral.', ...
    f_bilateral.', ...
    magnitud_analitico.', ...
    magnitud_numerico.', ...
    error2_magnitud.', ...
    fase_analitico.', ...
    fase_numerico.', ...
    error2_fase.', ...
    'VariableNames', ...
    {'k','Frecuencia_Hz',...
    'Magnitud_Analitica','Magnitud_Numerica','Error2_Magnitud',...
    'Fase_Analitica','Fase_Numerica','Error2_Fase'});

disp(' ');
disp('==========================================================================');
disp('TABLA 2. COMPARACIÓN DEL ESPECTRO BILATERAL: ANALÍTICO VS NUMÉRICO');
disp('==========================================================================');

disp(Tabla_Espectro);

% --------------------------------------------------------------
% Error cuadrático medio (MSE) global
%
% NOTA: se usa mean() en vez de immse() porque immse() requiere
% el Image Processing Toolbox, que no está disponible.
% --------------------------------------------------------------

MSE_magnitud = mean(error2_magnitud);
MSE_fase = mean(error2_fase);

fprintf('\n');
fprintf('MSE de la magnitud (analítico vs numérico): %.4e\n', MSE_magnitud);
fprintf('MSE de la fase (analítico vs numérico): %.4e\n', MSE_fase);


%% ============================================================
%              19. COMPARACIÓN CON LA ONDA CUADRADA
% ============================================================

% La onda cuadrada tendrá:
% - Amplitud: A = 1
% - Periodo: T = 2 s
% - Frecuencia fundamental: f0 = 0.5 Hz

% Definición de la onda cuadrada
x_cuadrada = sign(sin(w0*t));


%% ============================================================
%              20. GRÁFICA EN EL DOMINIO DEL TIEMPO
% ============================================================

figure;

plot(t, x, 'LineWidth', 1.5);

hold on;

plot(t, x_cuadrada, 'LineWidth', 1.5);

grid on;

xlabel('Tiempo (s)');
ylabel('Amplitud');

title('Comparación en el dominio del tiempo');

legend('Diente de sierra','Onda cuadrada');

xlim([-4 4]);
ylim([-1.2 1.2]);


%% ============================================================
%              21. COEFICIENTES DE LA ONDA CUADRADA
% ============================================================

% Para una onda cuadrada simétrica:
%
% a_k = 0
%
% b_k = 4/(k*pi), si k es impar
%
% b_k = 0, si k es par

b_cuadrada = zeros(size(k));

for n = 1:length(k)

    kk = k(n);

    if kk == 0

        b_cuadrada(n) = 0;

    elseif mod(kk,2) == 1

        % Armónico impar
        b_cuadrada(n) = ...
            4/(kk*pi);

    else

        % Armónico par
        b_cuadrada(n) = 0;

    end

end


%% ============================================================
%      22. COEFICIENTES COMPLEJOS DE LA ONDA CUADRADA
%                    (CORREGIDO)
% ============================================================

C_cuadrada = zeros(size(k_bilateral));

% PASO A: k = 0 y k > 0

for n = 1:length(k_bilateral)

    kk = k_bilateral(n);

    if kk == 0

        C_cuadrada(n) = 0;

    elseif kk > 0

        indice = kk + 1;

        C_cuadrada(n) = ...
            -1i*b_cuadrada(indice)/2;

    end

end

% PASO B: k < 0, simetría conjugada

for n = 1:length(k_bilateral)

    kk = k_bilateral(n);

    if kk < 0

        indice_positivo = find(k_bilateral == -kk);

        C_cuadrada(n) = conj(C_cuadrada(indice_positivo));

    end

end


%% ============================================================
%              23. MAGNITUD DEL ESPECTRO
%                  DE LA ONDA CUADRADA
% ============================================================

magnitud_cuadrada = abs(C_cuadrada);


%% ============================================================
%              24. ESPECTRO COMPARATIVO
% ============================================================

figure;

stem(f_bilateral, magnitud_analitico, ...
    'filled');

hold on;

stem(f_bilateral, magnitud_cuadrada);

grid on;

xlabel('Frecuencia (Hz)');
ylabel('|C_k|');

title('Comparación de los espectros bilaterales de magnitud');

legend('Diente de sierra','Onda cuadrada');

xlim([min(f_bilateral) max(f_bilateral)]);


%% ============================================================
%              25. TABLA 3: COMPARACIÓN DE ARMÓNICOS
% ============================================================

% Indicamos la magnitud de cada armónico
% en cada una de las señales.

armonico_diente = magnitud_analitico;

armonico_cuadrada = magnitud_cuadrada;

Tabla_Comparacion = table( ...
    k_bilateral.', ...
    f_bilateral.', ...
    armonico_diente.', ...
    armonico_cuadrada.', ...
    'VariableNames', ...
    {'k','Frecuencia_Hz',...
    'Magnitud_Diente_Sierra',...
    'Magnitud_Onda_Cuadrada'});

disp(' ');
disp('==========================================================================');
disp('TABLA 3. COMPARACIÓN DE LOS ARMÓNICOS');
disp('==========================================================================');

disp(Tabla_Comparacion);


%% ============================================================
%              26. COMENTARIO SOBRE EL DECRECIMIENTO
%                  DE LOS COEFICIENTES
% ============================================================

disp(' ');
disp('==============================================================');
disp('COMENTARIO SOBRE EL DECRECIMIENTO DE LOS COEFICIENTES');
disp('==============================================================');

disp('Los coeficientes del diente de sierra decrecen');
disp('aproximadamente como 1/k.');

disp('Por lo tanto, los primeros armónicos tienen mayor magnitud,');
disp('mientras que los armónicos de orden superior tienen menor');
disp('contribución a la señal.');

disp(' ');

disp('Los coeficientes b_k del diente de sierra alternan de signo:');
disp('+ , - , + , - , + , ...');


%% ============================================================
%              27. COMPARACIÓN CON LA ONDA CUADRADA
% ============================================================

disp(' ');
disp('==============================================================');
disp('COMPARACIÓN CON LA ONDA CUADRADA');
disp('==============================================================');

disp('El diente de sierra presenta armónicos pares e impares.');

disp('La onda cuadrada simétrica presenta únicamente');
disp('armónicos impares.');

disp(' ');

disp('Por esta razón, en el espectro del diente de sierra');
disp('aparecen componentes en todos los múltiplos enteros');
disp('de la frecuencia fundamental.');

disp(' ');

disp('En cambio, en la onda cuadrada desaparecen los armónicos');
disp('pares y solamente permanecen los armónicos impares.');

disp(' ');

disp('Ambas señales presentan un decrecimiento aproximado');
disp('de los coeficientes proporcional a 1/k.');


%% ============================================================
%              28. MOSTRAR DATOS PRINCIPALES
% ============================================================

disp(' ');
disp('==============================================================');
disp('PARÁMETROS DE LA SEÑAL');
disp('==============================================================');

fprintf('Amplitud A = %.2f\n', A);
fprintf('Periodo T = %.2f s\n', T);
fprintf('Frecuencia fundamental f0 = %.2f Hz\n', f0);
fprintf('Frecuencia angular w0 = %.2f rad/s\n', w0);
