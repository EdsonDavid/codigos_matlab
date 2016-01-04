clear, clc, close all   % borro la memoria, la pantalla y las figuras

%% ------------------------------------------------------------------------
%% NOTA: este codigo SOLO es apropiado para TENSION PLANA usando elementos
%% triangulares isoparamétricos de 10 nodos
%% ------------------------------------------------------------------------

%% definicion del problema
% Calcule los desplazamientos y las reacciones en los empotramiento, las
% deformaciones y los esfuerzos de la estructura en TENSION PLANA mostrada 
% en la figura adjunta

%% defino las variables/constantes
X    = 1;           % un par de constantes que ayudaran en la
Y    = 2;           % lectura del codigo
Ee   = 2000e6;       % modulo de elasticidad del solido (Pa) = 200 GPa
nue  = 0.30;        % coeficiente de Poisson
te   = 0.1;        % espesor del solido (m)
rhoe = 2300;        % densidad (kg/m^3)
g    = 9.81;        % aceleracion de la gravedad (m/s^2)
be = [0; -rhoe*g];  % vector de fuerzas masicas del elemento

%% cargar
% xnod - posicion de los nodos
% LaG  - definicion de elementos finitos con respecto a nodos
% Pic  - Matriz que contiene la direccion en que se dibujará el gráfico
malla1

nno  = size(xnod,1); % numero de nodos (numero de filas de xnod)
ngdl = 2*nno;        % numero de grados de libertad (dos por nodo)
gdl  = [(1:2:ngdl)' (2:2:ngdl)']; % nodos vs grados de libertad
nef = size(LaG,1);  % numero de EFs (numero de filas de LaG)

%% Se definen las restricciones 
ngdl_res = size(restricciones,1); % numero de grados de libertad restringidos
restric = zeros(ngdl_res,2);
for i = 1:ngdl_res
%                       nodo                direccion           desplazamiento    
   restric(i,:) = [ gdl(restricciones(i,1), restricciones(i,2)) restricciones(i,3) ];
end

%% Relacion de cargas puntuales
f = zeros(ngdl,1);     % vector de fuerzas nodales equivalentes global
f(gdl(28,Y)) = -5000;  % carga puntual en el nodo 13 dir Y
f(gdl(49,Y)) = -5000;  % carga puntual en el nodo 21 dir Y

%% Se dibuja la malla de elementos finitos
figure; hold on;
cgx = zeros(1,nef); cgy = zeros(1,nef); % almacena el centro de gravedad
for e = 1:nef
   line(xnod(Pic(e,[1:9 1]),X), xnod(Pic(e,[1:9 1]),Y));
   
   % Calculo la posicion del centro de gravedad del triangulo
   cgx(e) = mean(xnod(Pic(e,[1 3 5 7]),X));
   cgy(e) = mean(xnod(Pic(e,[1 3 5 7]),Y));
   h = text(cgx(e), cgy(e), num2str(e)); 
   set(h,'Color', [1 0 0], 'FontSize',12);
end
plot(xnod(:,X), xnod(:,Y), 'r*');
text(xnod(:,X), xnod(:,Y), num2str((1:nno)'), 'FontSize',12);
axis equal tight
title('Malla de elementos finitos','FontSize',26);

%% Funciones de forma serendipitas del elemento rectangular de 8 nodos:
% NOTA estas funciones de forma y sus derivadas se encontraron con el
% programa c5_funciones_forma_lagrangianos_rect_2D_8_nodos.m

Nforma = @(alpha,beta) [ ...
-((alpha + beta - 1)*(3*alpha + 3*beta - 1)*(3*alpha + 3*beta - 2))/2%N1
(alpha*(3*alpha - 1)*(3*alpha - 2))/2                                %N2
(beta*(3*beta - 1)*(3*beta - 2))/2                                   %N3
(9*alpha*(alpha + beta - 1)*(3*alpha + 3*beta - 2))/2                %N4
-(9*alpha*(3*alpha - 1)*(alpha + beta - 1))/2                        %N5
(9*alpha*beta*(3*alpha - 1))/2                                       %N6
(9*alpha*beta*(3*beta - 1))/2                                        %N7
-(9*beta*(3*beta - 1)*(alpha + beta - 1))/2                          %N8
(9*beta*(alpha + beta - 1)*(3*alpha + 3*beta - 2))/2                 %N9
-27*(alpha + beta - 1)*alpha*beta                 ];                 %N10

%% Derivadas de N con respecto a alpha
dN_dalpha = @(alpha,beta) [ ...
-(27*alpha^2 + 54*alpha*beta - 36*alpha + 27*beta^2 - 36*beta + 11)/2%dN1_dalpha 
(27*alpha^2 - 18*alpha + 2)/2                                        %dN2_dalpha 
0                                                                    %dN3_dalpha 
(9*(9*alpha^2 + 12*alpha*beta - 10*alpha + 3*beta^2 - 5*beta + 2))/2 %dN4_dalpha 
-(9*(6*alpha*beta - beta - 8*alpha + 9*alpha^2 + 1))/2               %dN5_dalpha 
(9*beta*(6*alpha - 1))/2                                             %dN6_dalpha 
(9*beta*(3*beta - 1))/2                                              %dN7_dalpha 
-(9*beta*(3*beta - 1))/2                                             %dN8_dalpha 
(9*beta*(6*alpha + 6*beta - 5))/2                                    %dN9_dalpha 
-27*beta*(2*alpha + beta - 1)                            ];  		 %dN10_dalpha

%% Derivadas de N con respecto a beta
dN_dbeta = @(alpha,beta) [ ...
-(27*alpha^2 + 54*alpha*beta - 36*alpha + 27*beta^2 - 36*beta + 11)/2%dN1_dbeta 
0                                                                    %dN2_dbeta 
(27*beta^2 - 18*beta + 2)/2                                          %dN3_dbeta 
(9*alpha*(6*alpha + 6*beta - 5))/2                                   %dN4_dbeta 
-(9*alpha*(3*alpha - 1))/2                                           %dN5_dbeta 
(9*alpha*(3*alpha - 1))/2                                            %dN6_dbeta 
(9*alpha*(6*beta - 1))/2                                             %dN7_dbeta 
-(9*(6*alpha*beta - 8*beta - alpha + 9*beta^2 + 1))/2                %dN8_dbeta 
(9*(3*alpha^2 + 12*alpha*beta - 5*alpha + 9*beta^2 - 10*beta + 2))/2 %dN9_dbeta 
-27*alpha*(alpha + 2*beta - 1)                             ];        %dN10_dbeta

%% Parametros de la cuadratura de Gauss-Legendre
% se asumira aqui el mismo orden de la cuadratura tanto en la direccion de
% alpha como en la direccion de beta
orden_gl = 4;                 % orden de la cuadratura de Gauss-Legendre

% El comando:
xw  = TriGaussPoints(orden_gl);
% Retorna una matriz de valores que contiene alpha, beta y el peso para
% cada punto de Gauss_Legendre
n_gl = size(xw,1); % El numero de puntos de Gauss es igual al numero de
                   % filas de la matriz xw


% Constantes para legibilidad
ALPHA  = 1;
BETA   = 2;
PESO_X_2 = 3;

%% ensamblo la matriz de rigidez global y el vector de fuerzas nodales
%  equivalentes global
K = sparse(ngdl,ngdl);   % matriz de rigidez global como RALA (sparse)
N = cell(nef,n_gl); % contenedor para las matrices de forma
B = cell(nef,n_gl); % contenedor para las matrices de deformacion

% matriz constitutiva del elemento para TENSION PLANA
De = [ Ee/(1-nue^2)     Ee*nue/(1-nue^2)  0
       Ee*nue/(1-nue^2) Ee/(1-nue^2)      0
       0                0                 Ee/(2*(1+nue)) ];
idx = cell(nef,20);
for e = 1:nef          % ciclo sobre todos los elementos finitos
   idx{e} = [ gdl(LaG(e,1),:)  gdl(LaG(e,2) ,:) ...
              gdl(LaG(e,3),:)  gdl(LaG(e,4) ,:) ...
              gdl(LaG(e,5),:)  gdl(LaG(e,6) ,:) ...
              gdl(LaG(e,7),:)  gdl(LaG(e,8) ,:) ...
              gdl(LaG(e,9),:)  gdl(LaG(e,10),:)];
   
   % Calculo las matrices de rigidez y el vector de fuerzas nodales
   % equivalentes del elemento
   Ke = zeros(20);
   fe = zeros(20,1);
   det_Je = zeros(n_gl,1); % en esta matriz se almacenaran los Jacobianos

   for p = 1:n_gl
     alpha_gl = xw(p,ALPHA);
     beta_gl  = xw(p,BETA);
     w_gl     = xw(p,PESO_X_2);

     % Se evaluan las funciones de forma en los puntos de integracion
     % de Gauss-Legendre
     NNforma = Nforma(alpha_gl, beta_gl);

     % Se evaluan las derivadas de las funciones de forma en los puntos
     % de integracion de Gauss-Legendre
     ddN_dalpha = dN_dalpha(alpha_gl, beta_gl);xe = xnod(LaG(e,:),X);
     ddN_dbeta  = dN_dbeta(alpha_gl, beta_gl); ye = xnod(LaG(e,:),Y);

     dx_dalpha = sum(ddN_dalpha.*xe);dy_dalpha = sum(ddN_dalpha.* ye);
     dx_dbeta  = sum(ddN_dbeta.*xe); dy_dbeta  = sum(ddN_dbeta .* ye);

     % Se ensambla la matriz Jacobiana del elemento
     Je = [ dx_dalpha   dy_dalpha
            dx_dbeta  dy_dbeta ];

     % Se calcula el determinante del Jacobiano
     det_Je(p) = det(Je);

     N{e,p} = zeros(2,2*10);
     B{e,p} = zeros(3,2*10);
     for i = 1:10
        % Se ensambla la matriz de funciones de forma N
        N{e,p}(:,[2*i-1 2*i]) = [ NNforma(i)  0         
                                    0           NNforma(i) ];

        % Se ensambla la matriz de deformacion del elemento B
        dNi_dx = (+dy_dbeta*ddN_dalpha(i) - dy_dalpha*ddN_dbeta(i))/det_Je(p);
        dNi_dy = (-dx_dbeta*ddN_dalpha(i) + dx_dalpha*ddN_dbeta(i))/det_Je(p);
        B{e,p}(:,[2*i-1 2*i]) = [ dNi_dx 0          % aqui se ensambla
                                    0      dNi_dy     % y asigna la matriz
                                    dNi_dy dNi_dx ];  % B_i
     end;

     % se arma la matriz de rigidez del elemento e
     Ke = Ke + B{e,p}'*De*B{e,p}*det_Je(p)*te*w_gl*0.5;

     % vector de fuerzas nodales equivalentes
     fe = fe + N{e,p}'*be*det_Je(p)*te*w_gl*0.5;
   end;
   
   if any(any(det_Je <= 0))
      error('Existen elementos con det_Je negativo en el elemento %d.\n', e);
   end;

   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke;
   f(idx{e},:)      = f(idx{e},:)   + fe;
end;   

%% Muestro la configuracion de la matriz K (K es rala)
figure
spy(K);
title('Los puntos representan los elementos diferentes de cero', ...
   'FontSize', 26);

%% grados de libertad del desplazamiento conocidos y desconocidos  
c = restric(:,1);   d = setdiff(1:ngdl,c)';

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |

%% extraigo las submatrices y especifico las cantidades conocidas
Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos
ac = restric(:,2);   % desplazamientos conocidos
qc = zeros(size(d)); % cargas de equilibrio en nodos libres ( = 0 siempre)

%% resuelvo el sistema de ecuaciones
ad = Kdd\((fc+qc)-Kdc*ac);   % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd;   % calculo fuerzas de equilibrio desconocidas
a = zeros(ngdl,1);  a(c) = ac;  a(d) = ad; % desplazamientos
q = zeros(ngdl,1);  q(c) = qd;  q(d) = qc; % fuerzas nodales equivalentes

%% imprimo los resultados
format short g
disp('Nodo   Despl_x (m)   Despl_y (m) = ');     [1:nno; reshape(a,2,nno)]'
disp('Nodo Fuerzas nodales equiv. X, Y (N) = '); [1:nno; reshape(f,2,nno)]'
disp('Nodo Fuerzas nodales equil. X, Y (N) = '); [1:nno; reshape(q,2,nno)]'

%% Dibujo la malla de elementos finitos y las deformaciones de esta
delta = reshape(a,2,nno)';
escala = 500;             % factor de escalamiento de la deformada
xdef = xnod + escala*delta; % posicion de la deformada
figure
hold on
for e = 1:nef
   h1 = line(xnod(Pic(e,[1:9 1]),X), xnod(Pic(e,[1:9 1]),Y)); % original
   set(h1, 'Color', [0 0 1]); % color expresado en notacion RBG entre 0 y 1
   h2 = line(xdef(Pic(e,[1:9 1]),X), xdef(Pic(e,[1:9 1]),Y)); % deformada
   set(h2, 'Color', [1 0 0]);
end
axis equal tight;
legend('Posicion original','Posicion deformada','Location', 'SouthOutside');
title(sprintf('Deformada escalada %d veces',escala), 'FontSize', 26);

%% Se calcula para cada elemento las deformaciones y los esfuerzos
def = cell(nef,n_gl);
esf = cell(nef,n_gl);
for e = 1:nef
   ae = a(idx{e});            % desplazamientos de los gdl del elemento e
   for pp = 1:n_gl
     def{e,pp} = B{e,pp}*ae;    % calculo las deformaciones
     esf{e,pp} = De*def{e,pp};  % calculo los esfuerzos
   end
end;

%% Se extrapolan los esfuerzos y las deformaciones a los nodos
num_elem_ady = zeros(nno,1);  % numero de elementos adyacentes
sx  = zeros(nno,1);
sy  = zeros(nno,1);
sz  = zeros(nno,1);
txy = zeros(nno,1);
txz = zeros(nno,1);
tyz = zeros(nno,1);

ex  = zeros(nno,1);
ey  = zeros(nno,1);
gxy = zeros(nno,1);

%Obtenida del programa CALCULO_MATRIZ_A.m
A = [...
    874799735707088903060234567639742069766696140119888322982719471/6924130959366159398569203094474447789951355622294637289801650121,                                   -346303756053555706690259881829999232791301885945/542320188875593523061575368668780375715793743651,                                   -346303756053555706690259881829999232791301885945/542320188875593523061575368668780375715793743651,     5304922101131221743636930890505677272249686220484903435308081040/2831316331698249573294624114613230060280521042109587066458030313,                                      30726668459170822842103082545230535625021735440/221757794123770646191371894450089387718270614403,                                      30726668459170822842103082545230535625021735440/221757794123770646191371894450089387718270614403
  -4421470208598866526684732636371649005270575429282970075802226705/6924130959366159398569203094474447789951355622294637289801650121,   -4421470208599005856748482269806018140346894534696859933288600305/6924130959365963078660830129619099499667897523798628172466448459,     874799735706959326153262780774653858978472741674252749878422799/6924130959365963078660830129619099499667897523798628172466448459,      633725113677585298208291536990106318291606321647079853242248720/4573664843512557003014392800529063943530072452638563722739895121,    8238426477808268964727371626898119786139219527926516214618852560/59457642965661555236436177502425484456749332304775894102407944967,  111403364123753592502775936988408568491395919039889953202875600080/59457642965661555236436177502425484456749332304775894102407944967
  -4421470208598866526684732636371649005270575429282970075802226705/6924130959366159398569203094474447789951355622294637289801650121,     874799735706959326153262780774653858978472741674252749878422799/6924130959365963078660830129619099499667897523798628172466448459,   -4421470208599005856748482269806018140346894534696859933288600305/6924130959365963078660830129619099499667897523798628172466448459,      633725113677585298208291536990106318291606321647079853242248720/4573664843512557003014392800529063943530072452638563722739895121,  111403364123753592502775936988408568491395919039889953202875600080/59457642965661555236436177502425484456749332304775894102407944967,    8238426477808268964727371626898119786139219527926516214618852560/59457642965661555236436177502425484456749332304775894102407944967
-12949822834293858164406818455391370964913717423655895642962940057/62317178634295434587122827850270030109562200600651735608214851089,  70941675960220413115546506511554107335479035942391946640002176391/62317178634293667707947471166571895497011077714187653552198036131, -28838632667212068167172998321397423553989471690253258398349601401/62317178634293667707947471166571895497011077714187653552198036131, 250195832825517427013533814016014800637863344384405036187664829776/535118786690969169352683957661900481393018476958711955560567729157,  93882680791233745718053993102770088220305891543023283417287717712/535118786690953997127925597521829360110743990742983046921671504703, -59298980112331790993499270751472753213702447768120577969353428144/535118786690953997127925597521829360110743990742983046921671504703
-28838632667211724453641720067425544190025532131864470839317778585/62317178634295434587122827850270030109562200600651735608214851089,  70941675960219972616760439057939307296049713307676217487402521991/62317178634293667707947471166571895497011077714187653552198036131, -12949822834294613117253830623270207595442692495855649501448186489/62317178634293667707947471166571895497011077714187653552198036131, -59298980112323716205469462143228721100494240961907765960296575664/535118786690969169352683957661900481393018476958711955560567729157,  93882680791233667205728224776656335883183769600233595381333848912/535118786690953997127925597521829360110743990742983046921671504703,     2749404756324220891300227000076316929285115701373407087466625776/5880426227373120847559621950789333627590593304867945570567818733
   7882408440024516243206415635428441785428389318194786713756905455/6924130959366159398569203094474447789951355622294637289801650121,  -9612877555736603165843624224349289016197094294321181088473659603/20772392878097889235982490388857298499003692571395884517399345377,   -1438869203810212660980626391256205672290575672650022801768878833/6924130959365963078660830129619099499667897523798628172466448459,  31294226930408029744123059804661772791299758410642970163547227248/178372928896989723117561319220633493797672825652903985186855909719, -19766326704110462806463099019396281424713218873990095842596280720/178372928896984665709308532507276453370247996914327682307223834901,      916468252108075392654785344418837003082895391626080671930334800/1960142075791040282519873983596444542530197768289315190189272911
   7882408440024516243206415635428441785428389318194786713756905455/6924130959366159398569203094474447789951355622294637289801650121,   -1438869203810212660980626391256205672290575672650022801768878833/6924130959365963078660830129619099499667897523798628172466448459,  -9612877555736603165843624224349289016197094294321181088473659603/20772392878097889235982490388857298499003692571395884517399345377,  31294226930408029744123059804661772791299758410642970163547227248/178372928896989723117561319220633493797672825652903985186855909719,      916468252108075392654785344418837003082895391626080671930334800/1960142075791040282519873983596444542530197768289315190189272911, -19766326704110462806463099019396281424713218873990095842596280720/178372928896984665709308532507276453370247996914327682307223834901
-28838632667211724453641720067425544190025532131864470839317778585/62317178634295434587122827850270030109562200600651735608214851089, -12949822834294613117253830623270207595442692495855649501448186489/62317178634293667707947471166571895497011077714187653552198036131,  70941675960219972616760439057939307296049713307676217487402521991/62317178634293667707947471166571895497011077714187653552198036131, -59298980112323716205469462143228721100494240961907765960296575664/535118786690969169352683957661900481393018476958711955560567729157,     2749404756324220891300227000076316929285115701373407087466625776/5880426227373120847559621950789333627590593304867945570567818733,  93882680791233667205728224776656335883183769600233595381333848912/535118786690953997127925597521829360110743990742983046921671504703
-12949822834293858164406818455391370964913717423655895642962940057/62317178634295434587122827850270030109562200600651735608214851089, -28838632667212068167172998321397423553989471690253258398349601401/62317178634293667707947471166571895497011077714187653552198036131,  70941675960220413115546506511554107335479035942391946640002176391/62317178634293667707947471166571895497011077714187653552198036131, 250195832825517427013533814016014800637863344384405036187664829776/535118786690969169352683957661900481393018476958711955560567729157, -59298980112331790993499270751472753213702447768120577969353428144/535118786690953997127925597521829360110743990742983046921671504703,  93882680791233745718053993102770088220305891543023283417287717712/535118786690953997127925597521829360110743990742983046921671504703
 26528821251593498010868447155674864368119809231785434713698316135/62317178634295434587122827850270030109562200600651735608214851089,                                    230869170702380909340590833810227437762679894023/542320188875593523061575368668780375715793743651,                                    230869170702380909340590833810227437762679894023/542320188875593523061575368668780375715793743651, -49430558866850411534469927350157281533603573162297482129647284912/535118786690969169352683957661900481393018476958711955560567729157,                                 -3871560225855698721715462964789361133130322623024/41912223089392652130169288051066894278753146122167,                                 -3871560225855698721715462964789361133130322623024/41912223089392652130169288051066894278753146122167];


for e = 1:nef
   sx(LaG(e,:),:) = sx(LaG(e,:),:)   + A * [ esf{e,1}(1)
                                             esf{e,2}(1)
                                             esf{e,3}(1)
                                             esf{e,4}(1)
                                             esf{e,5}(1)
                                             esf{e,6}(1)];

   sy(LaG(e,:),:) = sy(LaG(e,:),:)   + A * [ esf{e,1}(2)
                                             esf{e,2}(2)
                                             esf{e,3}(2)
                                             esf{e,4}(2)
                                             esf{e,5}(2)
                                             esf{e,6}(2)];
                                        
   txy(LaG(e,:),:) = txy(LaG(e,:),:) + A * [ esf{e,1}(3)
                                             esf{e,2}(3)
                                             esf{e,3}(3)
                                             esf{e,4}(3)
                                             esf{e,5}(3)
                                             esf{e,6}(3)];                       
                                          
   ex(LaG(e,:),:) = ex(LaG(e,:),:)   + A * [ def{e,1}(1)
                                             def{e,2}(1)
                                             def{e,3}(1)
                                             def{e,4}(1)
                                             def{e,5}(1)
                                             def{e,6}(1)];

   ey(LaG(e,:),:) = ey(LaG(e,:),:)   + A * [ def{e,1}(2)
                                             def{e,2}(2)
                                             def{e,3}(2)
                                             def{e,4}(2)
                                             def{e,5}(2)
                                             def{e,6}(2)];
                                        
   gxy(LaG(e,:),:) = gxy(LaG(e,:),:) + A * [ def{e,1}(3)
                                             def{e,2}(3)
                                             def{e,3}(3)
                                             def{e,4}(3)
                                             def{e,5}(3)
                                             def{e,6}(3)];                                                                 
                                          
   num_elem_ady(LaG(e,:),:) = num_elem_ady(LaG(e,:),:) + 1;
end

%% Alisado (promedio de los esfuerzos en los nodos)
sx  =  sx./num_elem_ady;  ex  =  ex./num_elem_ady;
sy  =  sy./num_elem_ady;  ey  =  ey./num_elem_ady;
txy = txy./num_elem_ady;  gxy = gxy./num_elem_ady;

%% Se calculan las deformacion ez en tension plana
ez  = -(nue/Ee)*(sx+sy);

%% Se imprimen y grafican las deformaciones en los nodos
disp('Deformaciones: (Nodo,ex,ey,ez,gxy) = '); 
disp([(1:nno)'  ex  ey  ez  gxy])
figure
subplot(2,2,1); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),ex(Pic(e,:)))
end;
ylabel('\epsilon_x','FontSize',26); axis equal tight; colorbar; 

subplot(2,2,2); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),ey(Pic(e,:)))
end;
ylabel('\epsilon_y','FontSize',26); axis equal tight; colorbar;

subplot(2,2,3); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),ez(Pic(e,:)))
end;
ylabel('\epsilon_z','FontSize',26); axis equal tight; colorbar;

subplot(2,2,4); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),gxy(Pic(e,:)))
end;
ylabel('\gamma_{xy}','FontSize',26); axis equal tight; colorbar;

%% Se imprimen y grafican los esfuerzos en los nodos
disp('Esfuerzos (Pa):  (Nodo,sx,sy,txy) = '); 
disp([(1:nno)'  sx  sy  txy])
figure
subplot(2,2,1); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),sx(Pic(e,:)))
end;
ylabel('\sigma_x (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(2,2,2); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),sy(Pic(e,:)))
end;
ylabel('\sigma_y (Pa)','FontSize',26); axis equal tight; colorbar;

subplot(2,2,3); hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),txy(Pic(e,:)))
end;
ylabel('\tau_{xy} (Pa)','FontSize',26); axis equal tight; colorbar;

%% Se calculan y grafican para cada elemento los esfuerzos principales y
%% sus direcciones
% NOTA: esto solo es valido para el caso de TENSION PLANA).
% En caso de DEFORMACION PLANA se deben calcular los valores y vectores 
% propios de la matriz de tensiones de Cauchy
%   [dirppales{e}, esfppales{e}] = eig([sx  txy 0    % matriz de esfuerzos
%                                       txy sy  0    % de Cauchy
%                                       0   0   0]);

s1   = (sx+sy)/2 + sqrt(((sx-sy)/2).^2+txy.^2); % esfuerzo normal maximo
s2   = (sx+sy)/2 - sqrt(((sx-sy)/2).^2+txy.^2); % esfuerzo normal minimo
tmax = (s1-s2)/2;                               % esfuerzo cortante maximo
ang  = 0.5*atan2(2*txy, sx-sy); % angulo de inclinacion de s1

%% Calculo de los esfuerzos de von Mises
s3 = zeros(size(s1));
sv = sqrt(((s1-s2).^2 + (s2-s3).^2 + (s1-s3).^2)/2);

%% imprimo los resultados
disp('Nodo,s1(Pa),s2(Pa),tmax(Pa),angulo(rad) = '); 
disp([(1:nno)'  s1  s2  tmax  ang])
disp('Nodo,Esfuerzos de von Mises (Pa) = ');
disp([(1:nno)'  sv]);

%% s1, s2, taumax
esc = 0.5; % escala para graficar las flechas

figure
hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),s1(Pic(e,:)))
end;

% Grafique lineas que indican las direcciones principales de sigma_1
norma = 1; % = s1 si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...   % En el nodo grafique una flecha (linea)
   norma.*cos(ang),norma.*sin(ang),... % indicando la direccion principal de sigma_1
   esc,...                       % con una escala esc
   'k', ...                      % de color negro
  'ShowArrowHead','off',...      % una flecha sin cabeza
  'LineWidth',2,...              % con un ancho de linea 2
  'Marker','.');                 % y en el punto (x,y) poner un punto '.'
quiver(xnod(:,X),xnod(:,Y),...   % la misma flecha ahora en la otra direccion,
   norma.*cos(ang+pi),norma.*sin(ang+pi),...  % es decir girando 180 grados
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('\sigma_1 (Pa)','FontSize',26); colorbar

figure
hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),s2(Pic(e,:)))
end;
% Grafique lineas que indiquen direcciones principales de sigma_2
norma = 1; % = s2 si quiere proporcional
quiver(xnod(:,X),xnod(:,Y),...             % flecha indicando la direccion
   norma.*cos(ang+pi/2),norma.*sin(ang+pi/2),... % principal de sigma_2
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
   norma.*cos(ang-pi/2),norma.*sin(ang-pi/2),...
   esc,'k',...
   'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('\sigma_2 (Pa)','FontSize',26); colorbar

figure;
hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),tmax(Pic(e,:)))
end;
% Grafique lineas que indiquen direcciones principales de tau_max,
norma = 1; % = tmax si quiere proporcional
quiver(xnod(:,X),xnod(:,Y), ...
       norma.*cos(ang+pi/4),norma.*sin(ang+pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
       norma.*cos(ang-pi/4),norma.*sin(ang-pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
       norma.*cos(ang+3*pi/4),norma.*sin(ang+3*pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
quiver(xnod(:,X),xnod(:,Y),...
       norma.*cos(ang-3*pi/4),norma.*sin(ang-3*pi/4),'k',...
       'ShowArrowHead','off','LineWidth',2,'Marker','.');
axis equal tight;
title('\tau_{max} (Pa)','FontSize',26); colorbar

figure; hold on;
for e = 1:nef
   fill(xnod(Pic(e,:),X),xnod(Pic(e,:),Y),sv(Pic(e,:)))
end;
ylabel('\sigma_v (Pa)','FontSize',26); axis equal tight; colorbar;
title('Esfuerzos de von Mises (Pa)','FontSize',26);                

%%
return; % bye, bye!
