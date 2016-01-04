clear, clc, close all

%% constantes
X = 1; Y = 2; TH = 3; Dens=78;

%% Unidades en kN y m

load cercha_ej_2 

LaG = barra(:,[1 2]);  % local a global
mat = barra(:,3);      % material

%        area      inercias_y       modulo de elasticidad
%        A(m^2)     I(m^4)          E(KPa)
props = [.04^2      .04^4/12         2e8
         pi*.05^2   pi*.05^4/4       2e8
         pi*.04^2   pi*.05^4/4       2e8];

A = props(:,1);   I = props(:,2);   E = props(:,3);

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
nbar = size(LaG,1);  % numero de EFs (numero de filas de LaG)
ngdl = 3*nno;        % numero de grados de libertad (dos por nodo)

%% gdl: grados de libertad
% fila = nodo
% col1 = gdl en direccion x
% col2 = gdl en direccion y
gdl  = [ (1:3:ngdl)' (2:3:ngdl)' (3:3:ngdl)' ]; % nodos vs grados de libertad

%% cargas aplicadas (gdl carga)

theta=90-atan2(4.5,9)*180/pi;
cargas_aplica = [ ...
                    gdl(2,X)  -10*cosd(theta)
                    gdl(2,Y)  -10*sind(theta)
                    gdl(3,Y)  -12
                    gdl(4,Y)  -20
                    gdl(5,X)  -10*cosd(theta)
                    gdl(5,Y)  -10*sind(theta)
                    gdl(6,Y)  -40
                    gdl(7,X)  -20*cosd(theta)
                    gdl(7,Y)  -20*sind(theta)];
                
dofs_cargados  = cargas_aplica(:,1);

f = zeros(ngdl, 1);
f(dofs_cargados) = cargas_aplica(:,2);

%% Se dibuja la estructura junto con su numeracion
figure(1); 
hold on;
for e = 1:nbar
   line(xnod(LaG(e,:),X), xnod(LaG(e,:),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx = (xnod(LaG(e,1),X) + xnod(LaG(e,2),X))/2;
   cgy = (xnod(LaG(e,1),Y) + xnod(LaG(e,2),Y))/2;   
   h = text(cgx, cgy, num2str(e)); set(h,'Color', [1 0 0]);
end

axis equal
grid minor
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'));
title('Numeracion de la estructura');

%% Peso propio de cada barra
peso_prop=zeros(nbar,1);
for e=1:nbar
    peso_prop(e)=Dens*A(mat(e));
end

%% Longitudes y ángulos de cada barra
L = zeros(nbar,1);
ang = zeros(nbar,1);
for e=1:nbar
   x1 = xnod(LaG(e,1), X);  x2 = xnod(LaG(e,2), X);
   y1 = xnod(LaG(e,1), Y);  y2 = xnod(LaG(e,2), Y);
   
   L(e) = sqrt((x2-x1)^2 + (y2-y1)^2);
   ang(e) = atan2((y2-y1),(x2-x1))*180/pi;
end
%% fuerzas nodales equivalentes para las diferentes barras
% (en este ejemplo las fuerzas nodales equivalentes estas siendo 
% especificadas con respecto al sistema de coordenadas globales)
fe = cell(nbar,1);
qxloc = cell(nbar,1);
qyloc = cell(nbar,1);
%            fxi         fyi                         m   
%            ton         ton                       ton-m
for e=1:nbar
    fe{e}=[   0   -peso_prop(e)*L(e)/2   -peso_prop(e)*L(e)^2/12 ...%i
              0   -peso_prop(e)*L(e)/2    peso_prop(e)*L(e)^2/12]'; %j
    qxloc{e} = @(x) -peso_prop(e)*sind(ang(e));
    qyloc{e} = @(x) -peso_prop(e)*cosd(ang(e));
end

%% ensamblo la matriz de rigidez global
Kg   = zeros(ngdl);   % separo memoria
Ke  = cell(nbar,1);
T   = cell(nbar,1);
idx = zeros(nbar,6);
for e = 1:nbar  % para cada barra
   % saco los 6 gdls de la barra e
   idx(e,:) = [gdl(LaG(e,1),:) gdl(LaG(e,2),:)];
   
   x1 = xnod(LaG(e,1), X);  x2 = xnod(LaG(e,2), X);
   y1 = xnod(LaG(e,1), Y);  y2 = xnod(LaG(e,2), Y);
   
     
   
   % matriz de transformacion de coordenadas para la barra e
   c = (x2-x1)/L(e);   s = (y2-y1)/L(e);  % sin y cos de la inclinacion
   T{e} = [ c  s  0  0  0  0;        
           -s  c  0  0  0  0;        
            0  0  1  0  0  0;
            0  0  0  c  s  0;
            0  0  0 -s  c  0;
            0  0  0  0  0  1];
         
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para la barra e
   AE = A(mat(e))*E(mat(e));       L2=L(e)^2;
   EI = E(mat(e))*I(mat(e));       L3=L(e)^3;
   Kloc = [ AE/L(e)    0         0        -AE/L(e)    0          0       
              0     12*EI/L3   6*EI/L2      0     -12*EI/L3    6*EI/L2
              0      6*EI/L2   4*EI/L(e)    0      -6*EI/L2    2*EI/L(e)
           -AE/L(e)    0         0         AE/L(e)       0          0
              0    -12*EI/L3  -6*EI/L2      0      12*EI/L3   -6*EI/L2
              0      6*EI/L2   2*EI/L(e)    0      -6*EI/L2    4*EI/L(e)];

   % matriz de rigidez local en coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};            
   Kg(idx(e,:),idx(e,:)) = Kg(idx(e,:),idx(e,:)) + Ke{e}; % sumo Ke{e} a K global
   f(idx(e,:))          = f(idx(e,:))          + fe{e}; % sumo a f global
end;

%% grados de libertad del desplazamiento conocidos (c) y desconocidos (d)
apoyos = [...
   gdl(1,Y)  0  
   gdl(2,X)  0
   gdl(2,Y)  0     
];

c = apoyos(:,1);
d = setdiff(1:ngdl, c);
%% Introduzco soporte inclinado
ang_sop = 330;
Tg = eye(ngdl);
Tg(1:2,1:2) = [ cosd(ang_sop) sind(ang_sop); -sind(ang_sop) cosd(ang_sop) ];
K = Tg*Kg*Tg';
%%
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd   |   | Kcc Kcd || ac |   | fd |  Recuerde que qc = 0
%|      | = |         ||    | - |    |
%| qc=0 |   | Kdc Kdd || ad |   | fc | 

% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = apoyos(:,2); % desplazamientos conocidos en contorno

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);        % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a  = zeros(ngdl,1);   a(c) = ac;  a(d) = ad; % desplazamientos
q  = zeros(ngdl,1);   q(c) = qd;             % fuerzas nodales equivalentes

%% retorno las fuerzas y los desplazamientos en el sistema de coordenadas
%% donde los grados de libertad son paralelos a los ejes
qg = Tg'*q;
ag = Tg'*a;

%% imprimo las fuerzas internas en cada barra referidas a las coordenadas
%% globales
qe_coord_loc = cell(nbar,1);
for e = 1:nbar % para cada barra
   fprintf('\n\n Fuerzas internas para barra %d en coord. globales = \n', e);
   qe_coord_glob = Ke{e}*ag(idx(e,:)) - fe{e};
   disp(qe_coord_glob)
   
   fprintf('\n\n Fuerzas internas para barra %d en coord. locales = \n', e);
   qe_coord_loc{e} = T{e}*qe_coord_glob;
   disp(qe_coord_loc{e});   
end;


%% imprimo los resultados
format short g
disp('Desplazamientos nodales                                                ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
vect_mov = reshape(ag,3,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g m, v = %12.4g m, theta = %12.4g m\n', ...
      i, vect_mov(i,X), vect_mov(i,Y), vect_mov(i,TH));
end;

disp(' ');
disp('Fuerzas nodales de equilibrio (solo se imprimen los diferentes de cero)  ');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
qq = reshape(qg,3,nno)';
for i = 1:nno   
   if ~isequal(qq(i,:), [0 0 0])
      fprintf('Nodo %3d qx = %12.4g kN, qy = %12.4g, mom = %12.4g kN \n', ...
         i, qq(i,X), qq(i,Y),  qq(i,TH));
   end;
end;

%% Dibujar la estructura y su deformada
esc = 50;
xdef = xnod + esc*vect_mov(:,[1 2]);

figure(2); hold on; title('Deformada exagerada');     xlabel('x, m'); ylabel('y, m'); axis equal
figure(3); hold on; title('Fuerza axial (kN)');       xlabel('x, m'); ylabel('y, m'); axis equal
figure(4); hold on; title('Fuerza cortante (kN)');    xlabel('x, m'); ylabel('y, m'); axis equal
figure(5); hold on; title('Momento flector (kN-m)');  xlabel('x, m'); ylabel('y, m'); axis equal

for e = 1:nbar
   x1 = xnod(LaG(e,1), X);  x2 = xnod(LaG(e,2), X);
   y1 = xnod(LaG(e,1), Y);  y2 = xnod(LaG(e,2), Y);

   dibujar_barra_deformada_portico(A(mat(e)), E(mat(e)), I(mat(e)), ...
      x1,y1, x2,y2, qxloc{e}, qyloc{e}, qe_coord_loc{e}, T{e}*a(idx(e,:)), ...
      esc, 0.01, 0.7, 1);
end

%%
return;
