% Calculo las funciones de forma
syms xi
N1 = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[1 0 0 0],3),xi));
N2 = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 1 0 0],3),xi));
N3 = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 1 0],3),xi));
N4 = factor(poly2sym(polyfit([-1 -1/3 1/3 1],[0 0 0 1],3),xi));

disp('Funciones de forma');
disp(N1); disp(N2); disp(N3); disp(N4);
disp('-----------0----------0-------------0---------');
disp('Derivadas de las funciones de forma');
disp(diff(N1,xi)); disp(diff(N2,xi)); disp(diff(N3,xi)); disp(diff(N4,xi));