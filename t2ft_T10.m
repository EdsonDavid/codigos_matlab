function ft = t2ft_T10(xnod, lado, carga, espesor)
% Esta función convierte las fuerzas superficiales aplicadas a un elemento
% finito triangular de 10 nodos a sus correspondientes cargas nodales
% equivalentes ft
%
%
% lado = 1452, 2673, 3891
%
% carga = [ t1x t1y t4x t4y t5x t5y t2x t2y ]; % si carga se aplica sobre lado 1452
%         [ t2x t2y t6x t6y t7x t7y t3x t3y ]; % si carga se aplica sobre lado 2673
%         [ t3x t3y t8x t8y t9x t9y t1x t1y ]; % si carga se aplica sobre lado 3891

%% Se definen algunas constantes
X = 1; Y = 2;

%% Parametros de la cuadratura de Gauss-Legendre
n_gl = 5;                        % orden de la cuadratura de Gauss-Legendre
[x_gl, w_gl] = gausslegendre_quad(n_gl);

%% Se definen las funciones de forma unidimensionales y sus derivadas
NN = @(xi) [ ...
-((xi - 1)*(3*xi - 1)*(3*xi + 1))/16    % N1
(9*(xi - 1)*(3*xi - 1)*(xi + 1))/16     % N2
-(9*(xi - 1)*(3*xi + 1)*(xi + 1))/16    % N3
((3*xi - 1)*(3*xi + 1)*(xi + 1))/16];   % N4

dNN_dxi = @(xi) [ ...
- (3*(3*xi - 1)*(xi - 1))/16 - (3*(3*xi + 1)*(xi - 1))/16 - ((3*xi - 1)*(3*xi + 1))/16      % dN1_dxi
(27*(xi - 1)*(xi + 1))/16 + (9*(3*xi - 1)*(xi - 1))/16 + (9*(3*xi - 1)*(xi + 1))/16         % dN2_dxi
- (27*(xi - 1)*(xi + 1))/16 - (9*(3*xi + 1)*(xi - 1))/16 - (9*(3*xi + 1)*(xi + 1))/16       % dN3_dxi
(3*(3*xi - 1)*(xi + 1))/16 + (3*(3*xi + 1)*(xi + 1))/16 + ((3*xi - 1)*(3*xi + 1))/16 ];     % dN4_dxi

%% Se definen los indices de los lados
switch lado
   case 1452,  idx = [ 1 4 5 2 ];
   case 2673,  idx = [ 2 6 7 3 ];
   case 3891,  idx = [ 3 8 9 1 ];     
   otherwise, error('Unicamente se permiten los lados 1452, 2673 ó 3891');
end

%% Se calcula el vector de fuerzas distribuidas en los nodos
te = zeros(20,1);
te(reshape([2*idx-1; 2*idx],8,1)) = carga(:);

%% Se calcula la integral
suma   = zeros(20);
N      = zeros(1,10);
dN_dxi = zeros(1,10);
for p = 1:n_gl
   N(idx)      = NN(x_gl(p));
   
   matN = [ N(1) 0    N(2) 0    N(3) 0    N(4) 0    N(5) 0    N(6) 0    N(7) 0    N(8) 0    N(9) 0    N(10) 0
            0    N(1) 0    N(2) 0    N(3) 0    N(4) 0    N(5) 0    N(6) 0    N(7) 0    N(8) 0    N(9) 0    N(10) ];
 
   dN_dxi(idx) = dNN_dxi(x_gl(p));
   
   dx_dxi = dN_dxi*xnod(:,X);
   dy_dxi = dN_dxi*xnod(:,Y);
   ds_dxi = sqrt(dx_dxi^2 + dy_dxi^2);

   suma = suma + matN'*matN*ds_dxi*w_gl(p);
end;

ft = espesor*suma*te;

%%
return;
