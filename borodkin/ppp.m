clc; clear;

tt = 0;
% for T=15:15
for T=15:5:35
    jj=0;
    tt=tt+1;
    for I=400:100:2100
        jj=jj+1;
        DATA(:,:,jj,tt) = table2array(readtable("files\"+T+"\"+I+".DAT"));
    end
end
full_num = numel(DATA(:,1,1,1));



shift_left = 25;
shift_right = 35;
o = '-o';
pol = 18;

for t=1:1
    for j=1:jj
        figure;
        hold on;
        [~, index_max] = max(DATA(:,2,j,t));
        arr_picks(j,t) = index_max;
        non0 = find(index_max);
        index_first = non0(1);
        
        x_left = index_max - shift_left;
        x_right = index_max + shift_right;
        zoom = x_left:x_right;
        c_left_x = (1:x_left).';
%         plot(DATA(:,1,j,t),DATA(:,2,j,t),o);
%         plot([zoom_x(end) zoom_x(end)],[half_y 0]);
%         plot([DATA(x_left,1,j,t) DATA(x_left,1,j,t)],[0.6 0.8]);
        shift_x = DATA(index_first,1,j,t) - DATA(index_max,1,j,t);
        shift_y = (sum(DATA(:,2,j,t)) - sum(DATA(zoom,2,j,t)))/(full_num-shift_left-shift_right-1);

        zoom_x = DATA(zoom,1,j,t) + shift_x;
        zoom_y = DATA(zoom,2,j,t) - shift_y;

        c_left_y = DATA(c_left_x,2,j,t) - shift_y;
        c_deff = max(c_left_y)-min(c_left_y);
        if (c_deff > 0.02)
            fprintf('%4.4f found at I=%d T=%d !\n',c_deff,j,t);
            
            c_right_x = (x_right:full_num).';
            c_right_y = DATA(c_right_x,2,j,t) - shift_y;
            base = ((sum(c_left_y) + sum(c_right_y)) / (numel(c_left_x) + numel(c_right_x)));
            fprintf('base: %d\n',base);
%             plot_hor(zoom_x(1),zoom_x(end),base);
            fprintf('replace shift_y=%4.4f with base\n',shift_y);
            shift_y = base;
        end
       
        [p,~,mu] = polyfit(zoom_x,zoom_y,pol);

        zoom_x2 = zoom_x(1):0.01:zoom_x(end);
        zoom_y2 = polyval(p,zoom_x2,[],mu);


        plot(zoom_x,zoom_y);
        plot(zoom_x2,zoom_y2,o);

        [zmax_y, h_index] = max(zoom_y2);
        half_y = abs(zmax_y - shift_y)/2;
        plot_hor(zoom_x(1),zoom_x(end),half_y);
        
        M = find_points(zoom_x2,zoom_y2,half_y);
        plot(M(:,1),M(:,2),'*g')
        half_X(j,t) = M(2,1)-M(1,1);
        hold off;
%         pause;
    end
end

figure;
grid on;
hold on;
for t=1:tt
    I = half_X(:,t);
    name = t*5+10;
    plot((1:numel(I))*100+300, I, o,'DisplayName', num2str(name));
end
hold off;
title(legend('Location','best'),'T, C');
xlabel('I, mA');
ylabel('FWHM, nm');
        

function plot_hor(first,last,y)
    plot([first last],[y y])
end

function [a] = VW(x,y,n)
    W = vander(x);
    W = W(1:n,2:n);
    A = W'*W;
    b = W'*y;
    a = A\b;
end

function [M] = find_points(x1,y1,y2)
    id = find( abs(diff(sign(y2 - y1)))>0 );
    M = zeros(length(id),2);
    for i = 1:length(id)
        j = id(i); % индекс левой пары точек перед пересечением
        k1 = (y1(j+1)-y1(j))/(x1(j+1)-x1(j));
        k2 = 0;
        b1 = y1(j) - k1*x1(j);
        b2 = y2;
        % составляем систему и находим точку пересечения:
        % k1*x - y = -b1
        % k2*x - y = -b2
        A = [k1, -1; % матрица коэф-тов
            k2, -1];
        B = [-b1;  % столбец левой части
            -b2];
        u = A\B; % решаем систему
        % сохраняем координаты:
        M(i,:) = u;
    end
end


%         zoom_xx_yy = [zoom_x2; zoom_y2];
%         plot(zoom_x,interpol(zoom_x,zoom_y,2*V),o);
% function [iy] = interpol(x,y,l)
%     iy = 0;
%     a = VW(x,y,l+1);
%     for i=1:l
%         iy = iy + a(i)*x.^(l-i);
%     end
%     
% end

%         zoom_xx_yy(:,:,j1,t1) = [zoom_x zoom_y];
