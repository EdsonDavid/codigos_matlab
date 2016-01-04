clear, clc, close all        % borro la memoria, la pantalla y las figuras

%% definicion del problema
% Calcule los desplazamientos y las reacciones en el empotramiento 
% de la viga mostrada
% 
% | b (carga distribuida de magnitud b)                 |
% |                                                     |
% |==*==*==0==*==*==0==*==*==0==*==*==0==*==*==0==*==*==| 
% 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
% |<--------------longitud L de la barra--------------->|   
%         
%          el area transversal de la barra es A
% -----------------------------------------------------------------
% Se usaron seis elementos isoparametricos lagrangianos cubicos
% -----------------------------------------------------------------

%% defino las variables
nef = 6;                      % numero de elementos finitos (EF)
nno = 3*nef + 1;              % numero de nodos
ngdl = nno;                   % el # de grados de libertad es el mismo # de nodos
E   = 200e9;    % Pa          % modulo de elasticidad de la barra
L   = 3;        % m           % longitud de la barra
P   = 5000;      % N           % carga nodal al final de la barra

xnod = linspace(0,L,nno);     % posicion de los nodos

le   = repmat(L/nef, nef, 1); % longitud de cada EF
A = @(r) pi*r.^2;         % area transversal del elemento

LaG = [1  2  3  4                  % definicion de EFs con respecto a nodos
       4  5  6  7
       7  8  9 10
       10 11 12 13
       13 14 15 16
       16 17 18 19];
   
xnod_ele = zeros(6,4);

for e = 1:nef
    xnod_ele(e,: ) = xnod(LaG(e,:));
end
%% Relacion de cargas puntuales
f = zeros(nno,1); % vector de fuerzas nodales equivalentes global
f(10) = P;       % relaciono la carga puntual en el nodo "nno"

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = zeros(nno);   % matriz de rigidez global
 
for e = 1:nef     % ciclo sobre todos los elementos finitos
   idx = LaG(e,:);
   
   r1 = radio(xnod_ele(e,1));
   r2 = radio(xnod_ele(e,end)); 
   
   Ke = [  (E*pi*(1247*r1^2 + 257*r1*r2 + 50*r2^2))/(420*le(e)), -(3*E*pi*(353*r1^2 + 57*r1*r2 + 31*r2^2))/(280*le(e)),    (3*E*pi*(47*r1^2 - 3*r1*r2 + 19*r2^2))/(140*le(e)),  -(E*pi*(163*r1^2 - 53*r1*r2 + 163*r2^2))/(840*le(e))
 -(3*E*pi*(353*r1^2 + 57*r1*r2 + 31*r2^2))/(280*le(e)),  (27*E*pi*(32*r1^2 + 13*r1*r2 + 11*r2^2))/(140*le(e)), -(27*E*pi*(29*r1^2 + 19*r1*r2 + 29*r2^2))/(280*le(e)),    (3*E*pi*(19*r1^2 - 3*r1*r2 + 47*r2^2))/(140*le(e))
    (3*E*pi*(47*r1^2 - 3*r1*r2 + 19*r2^2))/(140*le(e)), -(27*E*pi*(29*r1^2 + 19*r1*r2 + 29*r2^2))/(280*le(e)),  (27*E*pi*(11*r1^2 + 13*r1*r2 + 32*r2^2))/(140*le(e)), -(3*E*pi*(31*r1^2 + 57*r1*r2 + 353*r2^2))/(280*le(e))
  -(E*pi*(163*r1^2 - 53*r1*r2 + 163*r2^2))/(840*le(e)),    (3*E*pi*(19*r1^2 - 3*r1*r2 + 47*r2^2))/(140*le(e)), -(3*E*pi*(31*r1^2 + 57*r1*r2 + 353*r2^2))/(280*le(e)),  (E*pi*(50*r1^2 + 257*r1*r2 + 1247*r2^2))/(420*le(e))];

   x1 = xnod_ele(e,1);
   
   fe = [(le(e)*(le(e)^2 + 4*le(e)*x1 + 15*x1^2))/600
         (3*le(e)*x1*(2*le(e) + 5*x1))/200
         (3*le(e)*(le(e) + x1)*(3*le(e) + 5*x1))/200
         (le(e)*(12*le(e)^2 + 26*le(e)*x1 + 15*x1^2))/600]; % vector de fuerzas nodales equivalentes
   
   K(idx,idx) = K(idx,idx) + Ke;
   f(idx,:)   = f(idx,:)   + fe;
end;

%% grados de libertad del desplazamiento conocidos y desconocidos
c = [1 19];    d = setdiff(1:19,c);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |  % recuerde que qc=0 (siempre)
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = [0;0];               % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(ngdl,1);  q(c) = qd;             % fuerzas nodales equivalentes

%% se realizan unos calculos intermedios que necesitaremos mas adelante
nint = 10;           % numero de puntos donde se interpolará dentro del EF
xi = linspace(-1,1,nint)'; % coordenadas naturales
N = [-(9*xi.^3)/16 + (9*xi.^2)/16 + xi/16 - 1/16 ...
            (27*xi.^3)/16 - (9*xi.^2)/16 - (27*xi)/16 + 9/16 ...
            (27*xi)/16 - (9*xi.^2)/16 - (27*xi.^3)/16 + 9/16 ...
            (9*xi.^3)/16 + (9*xi.^2)/16 - xi/16 - 1/16]; % matr. de func. de forma
xx    = cell(nef,1); % interpol de posiciones (geometria) en el elemento
uu    = cell(nef,1); % interpol desplazamientos en el elemento
axial = cell(nef,1); % fuerzas axiales en el elemento
def   = cell(nef,1); % deformaciones en el elemento
esf   = cell(nef,1); % esfuerzos en el elemento
for e = 1:nef        % ciclo sobre todas los elementos finitos
   
   r1  = radio(xnod_ele(e,1));
   r2  = radio(xnod_ele(e,end));
   r_n = (r2-r1)*xi/2 + (r2+r1)/2;
   
   De  = E*A(r_n);
   Je  = le(e)/2;     % Jacobiano del elemento ( = dx_dxi)
   Be  = (1/Je)*[(9*xi)/8 - (27*xi.^2)/16 + 1/16 ...
                   (81*xi.^2)/16 - (9*xi)/8 - 27/16 ...
                    27/16 - (81*xi.^2)/16 - (9*xi)/8 ...
                   (27*xi.^2)/16 + (9*xi)/8 - 1/16]; % matriz de deformacion del elemento

   % vector de desplazamientos nodales del elemento a^{(e)}
   ae = [a(LaG(e,1));    a(LaG(e,2));    a(LaG(e,3));   a(LaG(e,4))   ]; % = a(LaG(e,:))';
   
   % vector de posiciones de nodos locales
   xe = [xnod(LaG(e,1)); xnod(LaG(e,2)); xnod(LaG(e,3)); xnod(LaG(e,4))]; %=xnod(LaG(e,:))';

   xx{e} = N*xe; % interpola sobre la geometría (coord naturales a geométricas)
   uu{e} = N*ae; % interpola sobre los desplazamientos
   
   def{e}   = Be*ae;
   esf{e}   = E*def{e};
   axial{e} = De.*def{e}; % fuerzas axiales en elemento finito e
end

%% Grafico la solucion
%% 1) grafico los desplazamientos de la barra
figure                             % cree un nuevo lienzo
x = linspace(0,L,30);              % 30 puntos unif/ distrib. entre 0 y L
subplot(2,2,1);
hold on;                           % no borre el lienzo 
for e = 1:nef % ciclo sobre todos los elementos finitos
   plot(xx{e}, uu{e}, 'b-'); % grafico solucion por MEF
end;
title('Desplazamiento');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Desplazamiento (m)')       % titulo del eje Y

%% 2) grafico la carga axial de la barra
subplot(2,2,2);                             % cree un nuevo lienzo
hold on;                           % no borre el lienzo
for e = 1:nef % ciclo sobre todos los elementos finitos
    plot(xx{e}, axial{e}, 'c-'); % grafico solucion por MEF
end;
title('Carga axial');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Carga axial (N)')          % titulo del eje Y
%% 3) grafico el esfuerzo en la barra
subplot(2,2,3);                             % cree un nuevo lienzo
hold on;                           % no borre el lienzo
for e = 1:nef % ciclo sobre todos los elementos finitos
    plot(xx{e}, esf{e}, 'r-'); % grafico solucion por MEF
end;
title('Esfuerzo \sigma');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Esfuerzo (Pa)')          % titulo del eje Y
%% 4) grafico las deformaciones la barra
subplot(2,2,4);                             % cree un nuevo lienzo
hold on;                           % no borre el lienzo
for e = 1:nef % ciclo sobre todos los elementos finitos
    plot(xx{e}, def{e}, 'm-'); % grafico solucion por MEF
end;
title('Deformacion \epsilon');
xlabel('Eje X (m)')                % titulo del eje X
ylabel('Deformacion')          % titulo del eje Y

%%
return; % bye, bye!
