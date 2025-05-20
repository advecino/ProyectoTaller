
%% benchmark_all_fixed.m
% Script corregido para medir tiempos de funciones del proyecto

% Lista de funciones a probar
funciones = {
    'Accel', 'AccelHarmonic', 'AccelPointMass', 'angl', 'anglesdr', 'anglesg', ...
    'AzElPa', 'Cheb3D', 'DEInteg', 'doubler', 'elements', 'EqnEquinox', ...
    'Frac', 'G_AccelHarmonic', 'gast', 'Geodetic', 'GHAMatrix', 'gibbs', ...
    'gmst', 'hgibbs', 'IERS', 'JPL_Eph_DE430', 'Legendre', 'LTC', ...
    'MeanObliquity', 'MeasUpdate', 'Mjday', 'Mjday_TDB', 'NutAngles', ...
    'NutMatrix', 'PoleMatrix', 'Position', 'PrecMatrix', 'R_x', 'R_y', ...
    'R_z', 'sign_', 'timediff', 'TimeUpdate', 'unit', 'VarEqn'};

% Variables comunes
load('DE430Coeff.mat')
global GM_Earth 
GM_Earth = 398600.4418e9;
global const
const.Arcs = 3600*180/pi;
const.MJD_J2000 = 51544.5;
const.Rad = pi/180;
eop_data = zeros(2000,13);
y = rand(6,1);
r = rand(3,1);
Mjd = 58000;

const.R_Earth = 6378137;

resultados = struct();

for i = 1:length(funciones)
    funcion = funciones{i};
    fprintf('\nProbando función: %s\n', funcion);

    try
        tic;
        switch funcion
            case 'Accel'
                aux = struct('Mjd_TT', Mjd, 'n', 10, 'm', 10);
                Accel(0, y);
            case 'AccelHarmonic'
                AccelHarmonic([7000e3; 0; 0], eye(3), 10, 10);
            case 'AccelPointMass'
                AccelPointMass(r, rand(3,1), GM_Earth);
            case 'angl'
                angl(rand(3,1), rand(3,1));
            case 'anglesdr'
                params = struct('R_Earth', 6378137);
                anglesdr(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, Mjd, Mjd+0.1, Mjd+0.2, r, r, r);
            case 'anglesg'
                anglesg(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, Mjd, Mjd+0.1, Mjd+0.2, r, r, r);
            case 'AzElPa'
                AzElPa(rand(3,1));
            case 'Cheb3D'
                Cheb3D(0.1, 13, 0, 1, rand(13,1), rand(13,1), rand(13,1));
            case 'DEInteg'
                DEInteg(@(t,y) y, 0, 1, 1e-8, 1e-8, 6, y);
            case 'doubler'
                doubler(1, 1, 7000e3, 7000e3, 6800e3, 7200e3, ...
                    [1;0;0], [0;1;0], [0;0;1], [7000e3;0;0], [0;7000e3;0], [0;0;7000e3], ...
                    0, 60, 'y');
            case 'elements'
                elements([r; rand(3,1)]);
            case 'EqnEquinox'
                EqnEquinox(Mjd);
            case 'Frac'
                Frac(5.67);
            case 'G_AccelHarmonic'
                G_AccelHarmonic([7000e3; 0; 0], eye(3), 10, 10);
            case 'gast'
                gast(Mjd);
            case 'Geodetic'
                Geodetic(r);
            case 'GHAMatrix'
                GHAMatrix(Mjd);
            case 'gibbs'
                gibbs(rand(3,1), rand(3,1), rand(3,1));
            case 'gmst'
                gmst(Mjd);
            case 'hgibbs'
                hgibbs(rand(3,1), rand(3,1), rand(3,1), Mjd, Mjd+0.1, Mjd+0.2);
            case 'IERS'
                IERS(eop_data, Mjd, 'n');
            case 'JPL_Eph_DE430'
                [~,~,r_Earth,~,~,~,~,~,~,~,~] = JPL_Eph_DE430(Mjd);
            case 'Legendre'
                Legendre(10, 10, pi/4);
            case 'LTC'
                LTC(0.5, 0.5);
            case 'MeanObliquity'
                MeanObliquity(Mjd);
            case 'MeasUpdate'
                x = zeros(6,1); z = x; P = eye(6); G = eye(6); s = ones(6,1); g = x; n = 6; K = zeros(6,6);
                MeasUpdate(x, z, g, s, G, P, n);
            case 'Mjday'
                Mjday(2020, 5, 1, 12, 0, 0);
            case 'Mjday_TDB'
                Mjday_TDB(Mjd);
            case 'NutAngles'
                NutAngles(Mjd);
            case 'NutMatrix'
                NutMatrix(Mjd);
            case 'PoleMatrix'
                PoleMatrix(0.1, 0.2);
            case 'Position'
                Position(0.5, 0.5, 100);
            case 'PrecMatrix'
                PrecMatrix(Mjd, Mjd+1);
            case 'R_x'
                R_x(pi/6);
            case 'R_y'
                R_y(pi/6);
            case 'R_z'
                R_z(pi/6);
            case 'sign_'
                sign_(3, -1);
            case 'timediff'
                timediff(0.3341, 37);
            case 'TimeUpdate'
                TimeUpdate(eye(6), eye(6), eye(6)*0.01);
            case 'unit'
                unit(rand(3,1));
            case 'VarEqn'
                VarEqn(0, zeros(42,1));

            otherwise
                continue
        end
        tiempo = toc;
        fprintf('Tiempo de ejecución de %s: %.6f s\n', funcion, tiempo);
        resultados.(funcion) = tiempo;
    catch ME
        fprintf('Error ejecutando %s: %s\n', funcion, ME.message);
        resultados.(funcion) = NaN;
    end
end

% Resumen
fprintf('\nResumen de tiempos:\n');
nombres = fieldnames(resultados);
for i = 1:length(nombres)
    fprintf('%20s: %.6f s\n', nombres{i}, resultados.(nombres{i}));
end

save('resultados_benchmark.mat', 'resultados');
