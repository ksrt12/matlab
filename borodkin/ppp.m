clc
clear

% tt = 0;
% for T=15:5:35
%     ii=0;
%     tt=tt+1;
%     for I=400:100:2100
%         ii=ii+1;
%         DATA(:,:,ii,tt) = table2array(readtable("files\"+T+"\"+I+".DAT"));
%     end
% end

% save('DATA.mat', 'DATA' )

tt = 5; ii = 18;
DATA = matfile('DATA.mat').DATA;
full_num = length(DATA);

gauss  = false;
draw   = false;
spline = false;

shift_left  = 25;
shift_right = 35;
pol = 18;
% draw   = true;
gauss  = true;
% spline = true;

for t=1:5
    name_t = t*5+10;
    for i=1:18
        name_j = i*100+300;
        [~, index_max] = max(DATA(:,2,i,t));
        arr_picks(i,t) = index_max;
        non0 = find(index_max);
        index_first = non0(1);
        
        x_left  = index_max - shift_left;
        x_right = index_max + shift_right;
        zoom = (x_left:x_right).';

%         plot(DATA(:,1,j,t),DATA(:,2,j,t),o);
%         plot([zoom_x(end) zoom_x(end)],[half_y 0]);
%         plot([DATA(x_left,1,j,t) DATA(x_left,1,j,t)],[0.6 0.8]);
        shift_x = DATA(index_first,1,i,t) - DATA(index_max,1,i,t);
        shift_y = (sum(DATA(:,2,i,t)) - sum(DATA(zoom,2,i,t)))/(full_num-shift_left-shift_right-1);

        zoom_x = DATA(zoom,1,i,t) + shift_x; % смещение по x
        zoom_y = DATA(zoom,2,i,t) - shift_y; % смещение по y

        c_left_x = (1:x_left).';
        c_left_y = DATA(c_left_x,2,i,t) - shift_y;

        c_right_x = (x_right:full_num).';
        c_right_y = DATA(c_right_x,2,i,t) - shift_y;
        base = ((sum(c_left_y) + sum(c_right_y)) / (numel(c_left_x) + numel(c_right_x)));  % база

        [y_max, h_index] = max(zoom_y);
        half_y = abs(y_max - base)/2;

        M0 = find_points(zoom_x,zoom_y,half_y);
    if (gauss) % gauss
        [zoom_x2,zoom_y2] = interpol(zoom_x,zoom_y,pol,0.01);
        M2 = find_points(zoom_x2,zoom_y2,half_y);
      if (length(M2) < 2)
        M2 = M0;
      end
        M = M2;
    elseif (spline) % spline
        zoom_x2 = zoom_x(1):0.01:zoom_x(end);
        zoom_y2 = interp1(zoom_x,zoom_y,zoom_x2,'spline');
        M = find_points(zoom_x2,zoom_y2,half_y);
    else
        M = M0;
    end

        x_ymax = zoom_x(h_index);

        [m_min_r_index, m_min_l_index, wtf] = find_rl_index(M,x_ymax);
        if (wtf)
            fprintf("T=%d, I=%d\n",name_t,name_j);
            [m_min_r_index, m_min_l_index, wtf] = find_rl_index(M0,x_ymax);
            M = M0;
        end

        half_X(i,t) = abs(M(m_min_r_index,1)-M(m_min_l_index,1)); % FWHM

    if (draw)
        figure('Name',"T="+name_t+",I="+name_j);
        title("T="+name_t+",I="+name_j);
        hold on;
        plot_hor(zoom_x,base);
        plot(zoom_x,zoom_y);
%         findpeaks(zoom_y,zoom_x,'Annotate','extents','WidthReference','halfheight')
      if (gauss || spline)
        plot(zoom_x2,zoom_y2,'-o');
%         findpeaks(zoom_y2,zoom_x2,'Annotate','extents','WidthReference','halfheight');
      end
        plot_hor(zoom_x,half_y);
        plot(M(:,1),M(:,2),'*g');
        hold off;
    end
    end
end

if (draw == 0)
    if (gauss)
        dname = "gauss";
    elseif (spline)
        dname = "spline";
    else
        dname = "";
    end
    final(half_X,tt,dname);
end
function final(half_X,tt,dname)
    point = ["." "+" "*" "o" "x"];
    color = ["r" "c" "g" "b" "m"];
    figure('Name',dname);
    grid on;
    hold on;
    x = 1:length(half_X);
    x_fori = x(1):0.1:x(end);
    ox = x*100+300;
    for t=1:tt
        I = half_X(:,t);
        name = t*5+10;
        y = interp1(x,I,x_fori,'spline');
        plot(x_fori*100+300, y, color(t),'DisplayName', '');
        plot(ox, I, point(t)+color(t),'DisplayName', num2str(name));
    end
    hold off;

    title(legend('Location','best'),'T, C');
    xlabel('I, mA');
    ylabel('FWHM, nm');
end

function [r,l,wtf] = find_rl_index(M,x_ymax)
    m = M(:,1) - x_ymax;
	r = find(m>0, 1 );
	l = find(m<0, 1, 'last' );
    if (isempty(r) || isempty(l))
    	fprintf('wtf: ');
    	wtf = true;
    else
        wtf = false;
    end
end

function plot_hor(x,y)
    plot([x(1) x(end)],[y y])
end

function [x,y] = interpol(x_old,y_old,pol,step)
    [p,~,mu] = polyfit(x_old,y_old,pol);
    x = x_old(1):step:x_old(end);
    y = polyval(p,x,[],mu);
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
