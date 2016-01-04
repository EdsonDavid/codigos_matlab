function dibujar_barra_deformada_cercha(A, E, I, x1,y1, x2,y2, qxloc,qyloc, qe, ae, esc2, esc3, esc4, esc5)
% Este programa dibuja APROXIMADAMENTE el elemento de portico deformado.
% OJO no da la solucion exacta, ya se para esta se requiere un xinit muy
% fino y una tolerancia en la solucion de la ecuacion diferencial muy baja
%
% El diagrama de momento flector se grafica en el lado opuesto de la fibra a
% tracción

%{
A = area
E = E
I = Ix local
(x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
qxloc = @(x) x.^2;   % carga en la direccion del eje x local (function handle)
qyloc = @(x) 0;      % carga en la direccion del eje y local (function handle)
qe = [ 0.01          % U1, V1, M1 reacciones del nodo 1 en coordenadas locales
      -0.01
       0.04
      -0.01          % U2, V2, M2 reacciones del nodo 2 en coordenadas locales
       0.02
      -0.07 ];
ae = [ 0.01          % u1, v1, t1 desplazamientos nodo 1 en coordenadas locales
      -0.01
       0.04
      -0.01          % u2, v2, t2 desplazamientos nodo 2 en coordenadas locales
       0.02
      -0.07 ];
esc2 = 10;            % escalamiento de la deformada
esc3 = 10;            % escalamiento del diagrama de axiales
esc4 = 10;            % escalamiento del diagrama de cortantes
esc5 = 10;            % escalamiento del diagrama de momentos
%}

%% resolver la ecuacion diferencial
xinit = linspace(0, sqrt((x2-x1)^2 + (y2-y1)^2), 51);
sol   = bvpinit(xinit, zeros(6,1));
sol   = bvp5c(@f, @bc, sol);

%% Calculos intermedios
s     = sol.x;
axial = sol.y(6,:);          % Fuerza axial [kN]
V     = sol.y(4,:);          % Fuerza cortante [kN]
M     = sol.y(3,:);          % Momento flector [kN/m]
u     = sol.y(5,:);          % Desplazamiento horizontal de la viga [m]
v     = sol.y(1,:);          % Desplazamiento vertical de la viga [m]
%theta = atan(sol.y(2,:));   % Angulo de giro  [rad]

% rotacion de la solucion antes de dibujar
ang = atan2(y2-y1, x2-x1);
T   = [ cos(ang)  -sin(ang)    % matriz de rotacion
        sin(ang)   cos(ang) ];

%% Dibujar de deformada
figure(2)
pos = T*[ s + esc2*u; esc2*v ];

xx = pos(1,:) + x1;
yy = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2);

%% Dibujar los diagramas de fuerza axial 
figure(3)
pos = T*[ s; esc3*axial ]; % escalamiento del diagrama

ss = pos(1,:) + x1;
aa = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2);
text(ss(1),   aa(1),   num2str(-qe(1)));
text(ss(end), aa(end), num2str( qe(3)));

%% Dibujar los diagramas de fuerza cortante
figure(4)
pos = T*[ s; esc4*V ]; % escalamiento del diagrama

ss = pos(1,:) + x1;
vv = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 vv y2], 'r-','LineWidth',2);
text(ss(1),   vv(1),   num2str( qe(2)));
text(ss(end), vv(end), num2str(-qe(4)));


%% Dibujar los diagramas de momento flector
figure(5)
pos = T*[ s; esc5*M ]; % escalamiento del diagrama

ss = pos(1,:) + x1;
mm = pos(2,:) + y1;

plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 mm y2], 'r-','LineWidth',2);
text(ss(26),   mm(26),   num2str( M(26)));
% text(ss(end), mm(end), num2str(-qe(6)));


%% -----------------------------------------------------------------------
   function dydx = f(x,y)
      % aqui se implementa la ecuacion diferencial para vigas de material
      % homogeneo y seccion transversal constante (A, E, I, qx, qy las provee la
      % funcion exterior)
      %      d^4 v(x)
      % E I ---------- = q(x)
      %        dx^4
      %
      %      d^2 u(x)
      % A E ---------- = -b(x)
      %        dx^2

      dydx = zeros(6,1);
      %         y(1)          = v
      dydx(1) = y(2);       % = theta
      dydx(2) = y(3)/(E*I); % = M/(EI)
      dydx(3) = y(4);       % = V
      dydx(4) = qyloc(x);   % = qyloc
      dydx(5) = y(6)/(A*E); % = u
      dydx(6) = -qxloc(x);   % = faxial
   end

%% ------------------------------------------------------------------------
   function res = bc(YL,YR)
      % condiciones de frontera (externas e internas)
      res = [ % tramo 1:   YL: apoyo izquierdo, YR: apoyo derecho
              YL(5) - ae(1)          % uloc(0)
              YL(1) - ae(2)          % vloc(0)
              YL(3)                  % M(0)=0
              YR(5) - ae(3)          % uloc(L)
              YR(1) - ae(4)          % vloc(L)
              YR(3)        ];        % M(L)=0
   end
end
%%End
