close all;
clear all;
clc;

%% Load the data
pos=load('sample'); % x y ...
net=load('network');% x y fn ft nx ny
ctr=load('control');% ...

% bidouille infernale !!! A REGLER !!! ****
%dec=1; % a supprimer
%for i=1:size(net,1)
%    net(i,7)=net(i,7)+dec;
%    net(i,8)=net(i,8)+dec;       
%end

%% Display parameters
nbound = 4; % number of bodies witch are not bodies at the begin (corners and fix points)

plot_disk_centers       = 0;
plot_disks              = 0;
plot_disk_velocities    = 0;
plot_normal_vectors     = 0;
plot_walls              = 1;
plot_contact_points     = 0;
plot_center2center      = 0;
plot_all_forces         = 0;
plot_strong_forces      = 1;

graph_fn_ft             = 1;
graph_contact_dir       = 0;

%% THE PLOT
figure;
hold on;

% Disk centers
if (plot_disk_centers)
    plot (pos(:,1),pos(:,2),'or');
end

% Disk radii
if (plot_disks)
    for i = 1:size(pos,1)
        [xc,yc,zc] = cylinder(pos(i,12),12);
        plot(xc(1,:)+pos(i,1),yc(1,:)+pos(i,2),'-k');
    end
end

% Disk velocities
if (plot_disk_velocities)
    fact = 1.0;
    quiver(pos(:,1),pos(:,2),pos(:,4),pos(:,5),fact);
end

% The walls
% 3---2
% |   |
% 4---1
if (plot_walls)
    plot ([pos(1,1) pos(2,1)],[pos(1,2) pos(2,2)],'-r');
    plot ([pos(2,1) pos(3,1)],[pos(2,2) pos(3,2)],'-r');
    plot ([pos(3,1) pos(4,1)],[pos(3,2) pos(4,2)],'-r');
    plot ([pos(4,1) pos(1,1)],[pos(4,2) pos(1,2)],'-r');
end


% The contact points
if (plot_contact_points)
    plot (net(:,1),net(:,2),'xk');
end

% Normal vectors
if (plot_normal_vectors)
    L=0.004; % length of line segment (for display)
    for i=1:size(net,1)
        plot ([net(i,1)-L*net(i,5) net(i,1)+L*net(i,5)],[net(i,2)-L*net(i,6) net(i,2)+L*net(i,6)],'-k')
    end
end

% Center to center (to check numbering)
if (plot_center2center)
    for i=1:size(net,1)
        plot ([pos(net(i,7),1) pos(net(i,8),1)],[pos(net(i,7),2) pos(net(i,8),2)],'-k')       
    end
end

% Normal forces
if (plot_all_forces)
    fmean=0.0;
    for i=1:size(net,1)
        fmean = fmean + net(i,3);
    end
    fmean = fmean / double(size(net,1))
    
    if (0) 
        L=0.002;
        for i=1:size(net,1)    
            plot ([net(i,1)-L*net(i,3)*net(i,5)/fmean  net(i,1)+L*net(i,3)*net(i,5)/fmean],...
                  [net(i,2)-L*net(i,3)*net(i,6)/fmean  net(i,2)+L*net(i,3)*net(i,6)/fmean],...
                  '-b')
        end
    else
        R=0.002;
        for i=1:size(net,1)
            tx = -net(i,6)*R*min(1.0,net(i,3)/fmean);
            ty =  net(i,5)*R*min(1.0,net(i,3)/fmean);
            if (net(i,7)>nbound)
                line ([pos(net(i,7),1)+tx  pos(net(i,8),1)+tx  pos(net(i,8),1)-tx  pos(net(i,7),1)-tx  pos(net(i,7),1)+tx],...
                      [pos(net(i,7),2)+ty  pos(net(i,8),2)+ty  pos(net(i,8),2)-ty  pos(net(i,7),2)-ty  pos(net(i,7),2)+ty])
            end
        end
    end
    
end

% strong/weak network
if (plot_strong_forces)
    factor = 0.8;
    fmean=0.0;
    for i=1:size(net,1)
        fmean = fmean + net(i,3);
    end
    fmean = fmean / double(size(net,1))
    
    L=0.006;
    for i=1:size(net,1)
        if (net(i,3) >= factor*fmean & net(i,7) > nbound & net(i,8) > nbound)
            plot ([pos(net(i,7),1) pos(net(i,8),1)],[pos(net(i,7),2) pos(net(i,8),2)],'LineWidth',3)
        else
            if (net(i,7) > nbound & net(i,8) > nbound)
            plot ([pos(net(i,7),1) pos(net(i,8),1)],[pos(net(i,7),2) pos(net(i,8),2)],'-k')
            end
        end
    end
end



%grid on;
axis equal;
set(gca,'Fontsize',20,'Fontweight','demi','Fontname', 'Helvetica');


%% en test
if (0)
    figure;
    fact = 1.0;
    quiver(pos(:,1),pos(:,2),pos(:,4),pos(:,5),fact);
    hold on;
    i=2
    plot ([pos(i,1) pos(i,1)+fact*ctr(i,4)],[pos(i,2) pos(i,2)+fact*ctr(i,5)],'- or');
    grid on;
    axis equal;
end

%% fn vs ft
if (graph_fn_ft)
    figure;
    plot (net(:,3),net(:,4),'.k');
    axis equal;
end

%% contact direction
if (graph_contact_dir)
    figure;
    plot (net(:,5),net(:,6),'.k');
    axis equal;
end

%% PDFs
if (1)
[n,xout] = hist(net(:,3),60);
%[n_ft,xout_ft] = hist(net(:,4),30);
n=n./size(net,1);
%n_ft=n_ft./size(net,1);
figure; 
semilogy (xout,n,'- .k');hold on;
%semilogy (xout_ft,n_ft,'- .r');

figure;
loglog(xout,xout.*n,'-r');
end




%% TEST

ft_up=0.0;
fn_up=0.0;
for i=1:size(net,1)
        if (net(i,7) == 2)
            fn_up = fn_up + net(i,3);
            ft_up = ft_up + net(i,4);
        end
end
fn_up
ft_up

