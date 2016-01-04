clear, clc, %close all,        % borro la memoria, la pantalla y las figuras
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
E   = 200e9;    % Pa          % modulo de elasticidad de la barra
L   = 3;        % m           % longitud de la barra
P   = 5000;      % N           % carga nodal al final de la barra

xnod = linspace(0,L,nno);     % posicion de los nodos

le   = repmat(L/nef, nef, 1); % longitud de cada EF

LaG = [1  2  3  4                  % definicion de EFs con respecto a nodos
       4  5  6  7
       7  8  9  10
       10 11 12 13
       13 14 15 16
       16 17 18 19];

xnod_ele = zeros(6,4);

for e = 1:nef
    xnod_ele(e,:) = xnod(LaG(e,:));
end

%% Parametros de la cuadratura de Gauss-Legendre
n_int_gl = 3;                 % orden de la cuadratura de Gauss-Legendre

% El comando:
[xi_gl,w_gl] = gausslegendre_quad(n_int_gl);
% calcula las raices (xi_gl) y los pesos (w_gl) de polinomios de Legendre

% >> [x_gl,w_gl] = gausslegendre_quad(1)
% x_gl = 0;
% w_gl = 2;
% >> [x_gl,w_gl] = gausslegendre_quad(2)
% x_gl = [  -0.577350269189626;  0.577350269189626 ];
% w_gl = [   1.000000000000000;  1.000000000000000 ];
% >> [x_gl,w_gl] = gausslegendre_quad(3)
% x_gl = [  -0.774596669241483;                  0; 0.774596669241483 ];
% w_gl = [   0.555555555555556;  0.888888888888889; 0.555555555555556 ];
% >> [x_gl,w_gl] = gausslegendre_quad(4)
% x_gl = [  -0.861136311594054; -0.339981043584857; 0.339981043584856; 0.861136311594053 ];
% w_gl = [   0.347854845137453;  0.652145154862547; 0.652145154862547; 0.347854845137453 ];

%% Relacion de cargas puntuales
f = zeros(nno,1); % vector de fuerzas nodales equivalentes global
f(10) = P;       % relaciono la carga puntual en el nodo "nno"

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = zeros(nno);   % matriz de rigidez global
A = @(r) pi*r.^2;
for e = 1:nef     % ciclo sobre todos los elementos finitos
   idx = LaG(e,:);
   Je = le(e)/2;  % Jacobiano del elemento ( = dx_dxi)
   % Calculo las matrices de rigidez y el vector de fuerzas nodales 
   % equivalentes del elemento
   
   x1 = xnod_ele(e,1);
   x2 = x1 + le(e)/3;
   x3 = x1 + (2/3)*le(e);
   x4 = x1 + le(e);
   
   r1 = radio(xnod_ele(e,1));
   r2 = radio(xnod_ele(e,end));
   
   Ke  = zeros(4);
   fe  = zeros(4,1);
   
   for m = 1:n_int_gl
      % matriz de deformacion del elemento
      xi = xi_gl(m);
      r_n = (r2-r1)*xi/2 + (r2+r1)/2;
      De = E*A(r_n);         % matriz constitutiva del elemento
      Be = (1/Je)*[(9*xi)/8 - (27*xi^2)/16 + 1/16 ...
                   (81*xi^2)/16 - (9*xi)/8 - 27/16 ...
                    27/16 - (81*xi^2)/16 - (9*xi)/8 ...
                   (27*xi^2)/16 + (9*xi)/8 - 1/16];
      Ke = Ke + w_gl(m)*Be'*De*Be*Je; % matriz de rigidez del elemento e

      % vector de fuerzas nodales equivalentes
      
      N = [-(9*xi^3)/16 + (9*xi^2)/16 + xi/16 - 1/16 ...
            (27*xi^3)/16 - (9*xi^2)/16 - (27*xi)/16 + 9/16 ...
            (27*xi)/16 - (9*xi^2)/16 - (27*xi^3)/16 + 9/16 ...
            (9*xi^3)/16 + (9*xi^2)/16 - xi/16 - 1/16]; % matr. de func. de forma
      
      x = N(1)*x1 + N(2)*x2 + N(3)*x3 + N(4)*x4;
      b = (1/5)*x^2;
      fe = fe + w_gl(m)*N'*b*Je; % vector de fuerzas nodales equivalentes
   end;

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
ac = [0 ; 0];               % desplazamientos conocidos

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(nno,1);  a(c) = ac;  a(d) = ad; % desplazamientos 
q = zeros(nno,1);  q(c) = qd;             % fuerzas nodales equivalentes

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
