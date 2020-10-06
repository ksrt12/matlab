% BORODKIN'S LABA 2
% created by kazakovstepan

dfname = "DATA.mat";
dfpath = "files";

DATA = checkdata(dfname,dfpath);
full_num = length(DATA{1,1});

% default vars start
gauss    = true;  % использовать Гауссову интерполяцию
spline   = false; % исgользовать интерполяцию сплайнами
draw     = false; % рисовать каждый график для текущих T,I
orig     = false; % рисовать исходные графики
pol      = 18;    % степень полинома для интерполяции
fpgauss  = 5;     % степень полинома для финального графика [0..i_last]|false
offset_l = 25;    % промежуток справа от пика
offset_r = 35;    % промежуток слева  от пика
t_first  = 1;
i_first  = 1;
[i_last, t_last] = size(DATA);
% default vars end

% fpgauss  = false;
% draw     = true;
% orig     = true;
% spline   = true;
% gauss    = false; 
% t_first  = 5;
% t_last   = 5;
% i_first  = 1;
% i_last   = 18;


if (gauss && spline)
    disp("Selected gauss instead of spline");
    spline = false;
end

if (orig)
    fig_num = t_last;
elseif (draw)
    fig_num = (t_last-t_first+1)*(i_last-i_first+1);
else
    fig_num = "final";
end

if (draw && (fig_num > 20))
    fprintf("Will be printed %d figures\n",fig_num);
    conf = input("Сontinue? (y/n): ",'s')+"";
    if ((conf ~= "y") && (conf ~= "Y"))
        draw = false;
    end
end

fprintf("Print %s figures\n",""+fig_num)
for t=t_first:t_last
    name_t = zt(t);
    if (orig)
        figure('Name',"T="+name_t);
        hold on;
    end
    for i=i_first:i_last
        name_i = zi(i);
        [~, index_max] = max(DATA{i,t}(:,2));
        
        x_left  = index_max - offset_l;
        x_right = index_max + offset_r;
        zoom = (x_left:x_right).';

        shift_x = DATA{i,t}(1,1) - DATA{i,t}(index_max,1);
        shift_y = (sum(DATA{i,t}(:,2)) - sum(DATA{i,t}(zoom,2)))/(full_num-offset_l-offset_r-1);

        zoom_x = DATA{i,t}(zoom,1) + shift_x; % смещение по x
        zoom_y = DATA{i,t}(zoom,2) - shift_y; % смещение по y

        c_left_x = (1:x_left).';
        c_left_y = DATA{i,t}(c_left_x,2) - shift_y;

        c_right_x = (x_right:full_num).';
        c_right_y = DATA{i,t}(c_right_x,2) - shift_y;
        base = ((sum(c_left_y) + sum(c_right_y)) / (length(c_left_x) + length(c_right_x)));  % база

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
            fprintf("T=%d, I=%d\n",name_t,name_i);
            [m_min_r_index, m_min_l_index, wtf] = find_rl_index(M0,x_ymax);
            M = M0;
        end

        half_X(i,t) = abs(M(m_min_r_index,1)-M(m_min_l_index,1)); % FWHM

    if (orig)
        plot(DATA{i,t}(:,1),DATA{i,t}(:,2),'DisplayName',""+name_i);
%         plot([zoom_x(end) zoom_x(end)],[half_y 0]);
%         plot([DATA{i,t}(x_left,1) DATA{i,t}(x_left,1)],[0.6 0.8]);
    elseif (draw)
        figure('Name',"T="+name_t+",I="+name_i);
        title("T="+name_t+",I="+name_i);
        hold on;
        plot_hor(zoom_x,base,'base');
        plot(zoom_x,zoom_y,'DisplayName','zoomed');
%         findpeaks(zoom_y,zoom_x,'Annotate','extents','WidthReference','halfheight')
        if (gauss || spline)
          plot(zoom_x2,zoom_y2,'-.','DisplayName','apprx');
%         findpeaks(zoom_y2,zoom_x2,'Annotate','extents','WidthReference','halfheight');
        end
        plot_hor(zoom_x,half_y,'half y');
        plot(M([m_min_l_index,m_min_r_index],1),M([m_min_l_index,m_min_r_index],2),'*g','DisplayName','L R points');
        legend('Location','best');
        hold off;
    end
    end
    if (orig)
        hold off;
        legend('Location','best');
        title("T="+name_t);
    end
end

if (~draw && ~orig)
    if (gauss)
        dname = "gauss"+pol;
    elseif (spline)
        dname = "spline";
    else
        dname = "";
    end
    if (fpgauss)
        dname = dname+":gauss"+fpgauss;
    else
        dname = dname+":spline";
    end
    final(half_X,dname,fpgauss);
end

function final(half_X,dname,p)
    [i_end,t_end] = size(half_X);
    point = ["." "+" "*" "o" "x"];
    color = ["r" "c" "g" "b" "m"];
    figure('Name',dname);
    disp(dname);
    grid on;
    hold on;
    i = 1:i_end;
    i2 = 1:0.1:i_end;
    for t=1:t_end
        fwhm = half_X(:,t);
        if (p)
            [~,y] = interpol(i,fwhm,p,0.1);
        else
            y = interp1(i,fwhm,i2,'spline');
        end
        plot(zi(i2), y, color(t),'DisplayName', '');
        plot(zi(i), fwhm, point(t)+color(t),'DisplayName', num2str(zt(t)));
    end
    hold off;

    title(legend('Location','best'),'T, °C');
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

function [z_i] = zi(i)
    z_i = i*100+300;
end

function [z_t] = zt(t)
    z_t = t*5+10;
end

function plot_hor(x,y,dname)
    plot([x(1) x(end)],[y y],'DisplayName',dname);
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

function [v] = varname(~)
	v=inputname(1);
end

function [data] = checkdata(fname,fpath)
if exist((fname),'file')
    data = matfile(fname).data;
    disp("Reading "+fname);
else
    disp(fname+" is not found!");
    disp("Creating it...");
    t = 0;
    for T=15:5:35
        i=0;
        t=t+1;
        for I=400:100:2100
            i=i+1;
            data{i,t} = table2array(readtable(fpath+"\"+T+"\"+I+".DAT"));
        end
    end
    save(fname,'data');
    disp(fname+" created and saved!");
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
