xw = TriGaussPoints(4);

A1 = sym('A1',[size(xw,1) 6]);
A2 = sym('A2',[10 6]);

ALPHA = 2;
BETA  = 3;
ALPHA_x_BETA = 4;
ALPHA_POW_2 = 5;
BETA_POW_2 = 6;

A1(:,1) = ones(size(xw,1),1);
A1(:,(ALPHA:BETA)) = xw(:,1:2);
A1(:,ALPHA_x_BETA) = A1(:,2).*A1(:,3);
A1(:,ALPHA_POW_2)  = A1(:,2).^2;
A1(:,BETA_POW_2)   = A1(:,3).^2;
    for i = 1:10
        A2(i,1) = 1;
        switch i
            case 1
                A2(i,ALPHA) = 0;
                A2(i,BETA)  = 0;
            case 2
                A2(i,ALPHA) = 1;
                A2(i,BETA)  = 0;
            case 3
                A2(i,ALPHA) = 0;
                A2(i,BETA)  = 1;
            case 4
                A2(i,ALPHA) = 1/3;
                A2(i,BETA)  = 0;
            case 5
                A2(i,ALPHA) = 2/3;
                A2(i,BETA)  = 0;
            case 6
                A2(i,ALPHA) = 2/3;
                A2(i,BETA)  = 1/3;
            case 7
                A2(i,ALPHA) = 1/3;
                A2(i,BETA)  = 2/3;
            case 8
                A2(i,ALPHA) = 0;
                A2(i,BETA)  = 2/3;
            case 9
                A2(i,ALPHA) = 0;
                A2(i,BETA)  = 1/3;
            case 10
                A2(i,ALPHA) = 1/3;
                A2(i,BETA)  = 1/3;
        end
        A2(i,ALPHA_x_BETA) = A2(i,ALPHA)*A2(i,BETA);
        A2(i,ALPHA_POW_2) = A2(i,ALPHA)^2;
        A2(i,BETA_POW_2) = A2(i,BETA)^2;
    end

    A = A2*inv(A1);
    disp(A)