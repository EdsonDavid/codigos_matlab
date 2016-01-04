%% Programa para obtener K y f correspondientes al elemento isoparametrico 
%% lagrangiano cubico (de cuatro nodos) unidimensional 

clear all; clc
format short
%% Definicion de variables
syms xi x1 r1 r2 L E A b

%% Defino el radio del elemento finito en una posición determinada 
%% y defino el area transversal correspondiente como funcion del radio
r = (r2-r1)*xi/2 + (r2+r1)/2;     
A = pi*r^2;

%% Defino las posiciones de los nodos
x4 = x1 + L;
x3 = x1 + 2/3*L;
x2 = x1 + L/3; 

%% Funciones de forma Lagrangianas
N1 = poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),xi);
N2 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),xi);
N3 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),xi);
N4 = poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),xi);

%% Interpolacion de la geometria y sus derivadas
x   = simple(N1*x1 + N2*x2 + N3*x3 +N4*x4);       dx_dxi = diff(x,   xi);
xi_ = solve(['x = ' char(x)],xi);                 dxi_dx = diff(xi_, 'x');
% recuerde que se debe garantizar que dx_dxi>0 y dxi_dx>0

%% Definicion de la matriz de forma N y matriz de deformacion del elemento B
N = [N1 N2 N3 N4];
B = simple([diff(N1,xi) diff(N2,xi) diff(N3,xi) diff(N4,xi)]*dxi_dx);

%% "matriz constitutiva"
D = E*A;

%% Calculo la matriz de rigidez del elemento
K = int(B.'*D*B*dx_dxi,xi,-1,1);
disp('E*pi/2*L *')
disp(K/(E*pi/(2*L)))

%% Calculo la matriz de fuerzas nodales equivalentes del elemento
b = (1/5)*x^2;
f = int(N.'*b*dx_dxi,xi,-1,1);
disp('(1/200)*')
disp(200*f)
