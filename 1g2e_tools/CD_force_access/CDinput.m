% Program for the determination of the contact network of a serie of image
% from the position of the center of mass of the particle and the size of
% their radius.

% The CONF file should come with the position of the corner but without
% the fixed points. And without the first line with the number of
% particles
% 
disp('========== BEGIN');
%% Cleaning
disp('Cleanning memory');
clear all;
close all;


%% Parameters and input of the program

% Photo information
 
nphoto = 101; % the photo number

% Quantity add to the radii
addradius = 1.23; % [pixels]

% Folder in which there is the output file (with)
%folder = strcat('/home/3S-LAB/rmaurin/Stage/20120404experiment/DICfile/sample3/PostProcess/');
folder = strcat('');

radCorner = 142;  % [pixels]

dt      = 5.0;    % [sec/photo]
gamma   = 4.8e-3; % [rad/sec]
density = 2900.0; % [kg/m3]

% Conversion factors
pxl2m = 1.45e-4;


%% Initialization 
disp('Initialization');
evaluatemaxcontactnumber = 5000; 
nbcontact = 0;
Overlap = 0;
Vect = zeros(2,1);  
Cornerposition = zeros(4,1);
a = 0;
p = 0;
mmtolcontact = 2; % (5)tolerance on the contact with the side
g = 0;
k = 0;
f = 0;
fixedpoint = 0;
contactwall = 0;

%% Load of the file and caracterisation of the sample and machine
disp('Load file');  

if nphoto<10
   file = strcat('CONF000');
elseif nphoto<100
   file = strcat('CONF00');
elseif nphoto<1000
   file = strcat('CONF0');
elseif nphoto>=1000
   file = strcat('CONF');
end

u=num2str(nphoto);  % convert the number of the photo to a string 

path = strcat(folder,file,u); % path = name of the file corresponding to a data we want to analyze
B = importdata(path);

nbpart = length(B);

% De-reverse the y axis and redefine the origin
for i = 1:nbpart
    B(i,2) = -B(i,2)+4032;
end

mass = zeros(nbpart,1);
mom  = zeros(nbpart,1);


% Put the real value of the radius, and compute mass and inertial moment
for i = 1:nbpart 
    if (B(i,4)==1)    % the corner are considered as particle with radius 142 pxl which correspond approximatelly to the reality
        f=f+1;
        B(i,4) = radCorner+addradius;
        Cornerposition(f,1) = B(i,1);
        Cornerposition(f,2) = B(i,2);
        if f>4
            display('Error: more than 4 Corner!!');
            break
        end
        area = pi*(B(i,4)*pxl2m)^2;
        mass(i) = area*0.06*density*100;  % *100 to get a mass really bigger than the particles for the wall
        mom(i) = 0.5*(B(i,4)*pxl2m)^2*mass(i);
    elseif B(i,4) == 48
        B(i,4)=48.21+addradius;
        area = pi*(B(i,4)*pxl2m)^2;
        mass(i) = area*0.06*density;        
        mom(i) = 0.5*(B(i,4)*pxl2m)^2*mass(i);
    elseif  B(i,4)== 55
        B(i,4) = 55.10+addradius;
        area = pi*(B(i,4)*pxl2m)^2;
        mass(i) = area*0.06*density;
        mom(i) = 0.5*mass(i)*(B(i,4)*pxl2m)^2;
    elseif B(i,4) == 76
        B(i,4) = 75.75 + addradius;
        area = pi*(B(i,4)*pxl2m)^2;
        mass(i) = area*0.06*density;
        mom(i) = 0.5*mass(i)*(B(i,4)*pxl2m)^2;
    end
end


% Determine the parameters necessary for the condition of contact
% with the machine 
% 3--2
% |  |
% 4--1
for h = 1:4

    if h==4    % couple of points: 1-2,2-3,3-4,4-1, to have the last one well
        a=1;
    else
        a=h+1;
    end

    px(h) = Cornerposition(h,1);   %position of the corner i
    py(h) = Cornerposition(h,2);
    lx(h) = Cornerposition(a,1)-px(h);    % vector linking the corner i-i+1
    ly(h) = Cornerposition(a,2)-py(h);
    L(h) = sqrt(lx(h)*lx(h)+ly(h)*ly(h));   
    lx(h) = lx(h)/L(h);ly(h)=ly(h)/L(h);       
    nx(h) = -ly(h);ny(h)=lx(h);           %unit vector perpendicular to the vector linking the corner i-i+1 
end



%% Contact network determination
display('Tracking of contacts')

for i=5:nbpart
    % position and radius of the particle i 
    xi = B(i,1);     % center and radius of grain i
    yi = B(i,2);
    ri = B(i,4);  
    posi = [xi; yi]; % position of the center of the grain i

    % Contact with the plane part of the walls, and then the circular one    
    for h=1:4
        % plane
        ddx = xi-px(h); % vector between the corner and the particle considered
        ddy = yi-py(h);
        dn = abs(ddx*nx(h)+ddy*ny(h))-ri;    % Projection on the perpendicular direction to the wall - radius of the particle
        if (dn <= mmtolcontact)&&(i>4) %
            nbcontact = nbcontact+1;      
            Contact(nbcontact,1) = xi-(ri+dn/2)*nx(h); % store the position of the contact
            Contact(nbcontact,2) = yi-(ri+dn/2)*ny(h);
            Contact(nbcontact,3) = nx(h);            % store its orientation, from the wall to the rod
            Contact(nbcontact,4) = ny(h);
            Contact(nbcontact,5) = h;                % store the number of the "wall" associated with the contact
            Contact(nbcontact,6) = i;                % wall
            contactwall = contactwall+1;
        end
        
        % circle
        % position and radius of the particle i 
        xi = B(h,1);     % center and radius of grain i
        yi = B(h,2);
        ri = B(h,4);  
        posi = [xi; yi]; % position of the center of the grain i
        xj = B(i,1); % center and radius of grain j
        yj = B(i,2);
        rj = B(i,4);

        % Vector associated with the considered couple of particle (i,j)
        Vect(1) = xj-xi; 
        Vect(2) = yj-yi; 
        normVect = norm(Vect);
        unitVect = Vect/normVect;

        % sumradius = sum of the radius + a contribution added
        % to the radius in order to take into account the
        % imprecision of the radius measurement 

        sumradius = ri+rj;   
        posj = [xj; yj]; % position of the center of the grain j 
        dn = normVect - sumradius;


        % if it is not a fixed point and the length of the vector is smaller or equal to the sum of the 2 radius + a distance in pxl added to both radius, it means that there is contact

        if (ri>2) && (rj>2) && (dn <= 0)  
            nbcontact = nbcontact + 1;  % count a contact
            contactwall = contactwall+1;
            Overlap = Overlap + sumradius - normVect;  % quantify the overlap in length (not area)

            % Position of the contact and the point
            % corresponding to it
            Contact(nbcontact,1) = posi(1) + (ri+dn/2)*unitVect(1);   % Calculate the contact network and put it in a matrix
            Contact(nbcontact,2) = posi(2) + (ri+dn/2)*unitVect(2);
            Contact(nbcontact,3) = unitVect(1);            % Store also the number of the particles involved in the contact
            Contact(nbcontact,4) = unitVect(2);            % from i to j
            Contact(nbcontact,5) = h;
            Contact(nbcontact,6) = i;
        end
        
        xi = B(i,1);     % center and radius of grain i
        yi = B(i,2);
        ri = B(i,4);  
        posi = [xi; yi];

    end

    % Contact between particles

    for j=i+1:nbpart % sum over all the couple of particle (i,j)
        % Position and radius of the couple of particle (i,j)
        xj = B(j,1); % center and radius of grain j
        yj = B(j,2);
        rj = B(j,4);

        % Vector associated with the considered couple of particle (i,j)
        Vect(1) = xj-xi; 
        Vect(2) = yj-yi; 
        normVect = norm(Vect);
        unitVect = Vect/normVect;

        % sumradius = sum of the radius + a contribution added
        % to the radius in order to take into account the
        % imprecision of the radius measurement 

        sumradius = ri+rj;   
        posj = [xj; yj]; % position of the center of the grain j 
        dn = normVect - sumradius;


        % if it is not a fixed point and the length of the vector is smaller or equal to the sum of the 2 radius + a distance in pxl added to both radius, it means that there is contact

        if (ri>2) && (rj>2) && (dn <= 0)  
            nbcontact = nbcontact + 1;  % count a contact
            Overlap = Overlap + sumradius - normVect;  % quantify the overlap in length (not area)

            % Position of the contact and the point
            % corresponding to it
            Contact(nbcontact,1) = posi(1) + (ri+dn/2)*unitVect(1);   % Calculate the contact network and put it in a matrix
            Contact(nbcontact,2) = posi(2) + (ri+dn/2)*unitVect(2);
            Contact(nbcontact,3) = unitVect(1);            % Store also the number of the particles involved in the contact
            Contact(nbcontact,4) = unitVect(2);            % from i to j
            Contact(nbcontact,5) = i;
            Contact(nbcontact,6) = j;
        end
    end % end of the loop over the particle j
end % end of the over the particle i

%% Plot the contacts
figure(1);hold on;
quiver(Contact(:,1),Contact(:,2),Contact(:,3),Contact(:,4),0.3,'-k');
%quiver(Contact(:,1),Contact(:,2),-Contact(:,3),-Contact(:,4),0.3,'-k');
% the walls
[x,y,z] = cylinder(radCorner,20);
for h = 1:4
     if h==4    % couple of points: 1-2,2-3,3-4,4-1, to have the last one well
        a=1;
    else
        a=h+1;
    end
    plot(x(1,:)+px(h),y(1,:)+py(h),'-r');
    plot([px(h) px(a)],[py(h) py(a)],'-r');
end
% the rods
for i = 5:nbpart
    [x,y,z] = cylinder(B(i,4),12);
    plot(x(1,:)+B(i,1),y(1,:)+B(i,2),'-g');
end
axis equal;

% Calculation of the coordination number and mean overlap
% associated with the current alpha radius factor
nbcoordination = (contactwall + nbcontact*2)/(nbpart-(fixedpoint+4));   % remove the fixed point and corner (4) from the nbpart
fprintf('Coordination number   = %f\n',nbcoordination);
fprintf('nb contacts           = %d\n',nbcontact);
fprintf('nb contacts with wall = %d\n',contactwall);

Overlap(nphoto) = Overlap/(nbpart-(fixedpoint+4));        

%% Computation of particle velocities

Particlespeed = zeros(nbpart,2);
%importdata of the previous picture
if nphoto-1<10
   file = strcat('CONF000');
elseif nphoto-1<100
   file = strcat('CONF00');
elseif nphoto-1<1000
   file = strcat('CONF0');
elseif nphoto-1>=1000
   file = strcat('CONF');
end

u = num2str(nphoto-1); % convert the number of the photo to a string 
path2 = strcat(folder,file,u); % path = name of the file corresponding to a data we want to analyze
Av = importdata(path2);

for i = 1:length(Av)
    Av(i,2) = - Av(i,2)+4032;
end

for i = 1:nbpart
    Particlespeed(i,1) = (B(i,1) - Av(i,1))/dt;   % speed of the particle
    Particlespeed(i,2) = (B(i,2) - Av(i,2))/dt;
    Particlespeed(i,3) = (B(i,3) - Av(i,3))/dt;   % speed of rotation of the particle
    
    
    X(i) = B(i,1);
    Y(i) = B(i,2);
    w(i) = Particlespeed(i,1);
    v(i) = Particlespeed(i,2);   
end

%% Plot of velocities
figure(2);hold on;
quiver(X,Y,w,v,'-k');
% the walls
[x,y,z] = cylinder(radCorner,20);
for h = 1:4
     if h==4    % couple of points: 1-2,2-3,3-4,4-1, to have the last one well
        a=1;
    else
        a=h+1;
    end
    plot(x(1,:)+px(h),y(1,:)+py(h),'-r');
    plot([px(h) px(a)],[py(h) py(a)],'-r');
end
axis equal;


%% Save file
fprintf('Save\n');

fid = fopen('1g2e_data.txt', 'w');

fprintf(fid, 'xgrav 0.0\n');
fprintf(fid, 'ygrav 0.0\n');
fprintf(fid, 'epsf 0.00001\n');
fprintf(fid, 'niterconv 500\n');
fprintf(fid, 'nitermx 100000\n');
fprintf(fid, 'nitermn 300\n');
fprintf(fid, 'en 0.0\n');
fprintf(fid, 'mu 0.3\n'); % mu different entre parts et murs ???
fprintf(fid, 'et 0.0\n');
fprintf(fid, 'dt %12.8f\n',dt);
fprintf(fid, 'END_PARAMETERS\n\n');
    
fprintf(fid, '%i\n', size(B,1));
for i=1:size(B,1)
    % grain mass mom x y rot vx vy vrot
    fprintf(fid, '%12.8e %12.8e %12.8f   %12.8f %12.8f %12.8f   %12.8f %12.8f %12.8f\n',mass(i),mom(i),pxl2m*B(i,4),...
        pxl2m*B(i,1),pxl2m*B(i,2),B(i,3),...
        pxl2m*Particlespeed(i,1),pxl2m*Particlespeed(i,2),Particlespeed(i,3));
end
fprintf(fid, '\n');

fprintf(fid, '%i\n', size(Contact,1));
for i=1:size(Contact,1)
    % contact i j x y nx ny
    fprintf(fid, '%i %i   %12.8f %12.8f   %12.8f %12.8f\n',Contact(i,5),Contact(i,6),pxl2m*Contact(i,1),pxl2m*Contact(i,2),Contact(i,3),Contact(i,4));
end
fprintf(fid, '\n');
    
% ctr.vxImposed >> ctr.vyImposed >> ctr.vrotImposed;
% ctr.xvalue >> ctr.yvalue >> ctr.rotvalue;
fprintf(fid, '4\n');
    
    
    
% Value measured from the outputfile of 1gamma2epsilon for the speed
% and 
    
    
%     %import the file
%     I = dlmread('/home/3S-LAB/rmaurin/Stage/20120404experiment/Macroscopicstressandstrain/2012-04-04/04-cis50.dat','',5,0);
%     t = 5*nphoto;  % time corresponding to the picture considered
%     a = round(t/6); % The measurement in the 1g2e file is made every 6 sec. So to know the line, take the integer part of the number obtained to get the file as near as possible to the photo.
%     
% 
%     measured_gammadot = (I(a,6) - I(a+1,6))/dt;   % the angle is counted in the wrong sens, so inverse the relation...
    
    
    measured_gammadot = 0.026694/dt;
    measuredFtopy = -10^3*(0.722643 + 0.71576)/2;% force imposed on the top plate, - to get it reverse
    
    
    
    % Determination of the speed of deformation (angle)
    
    l14(1)= B(4,1) -  B(1,1);
    l14(2)= B(4,2) -  B(1,2);
    n14 = l14/norm(l14);
    
    l12(1)= (B(2,1) -  B(1,1));
    l12(2)= (B(2,2) -  B(1,2));
    norml12 = norm(l12);
    
    l2 = l12*n14';

    gamma1B = asin(l2/norml12) ;

    
    l14(1)= Av(4,1) -  Av(1,1);
    l14(2)= Av(4,2) -  Av(1,2);
    n14 = l14/norm(l14);
    
    l12(1)= (Av(2,1) -  Av(1,1));
    l12(2)= (Av(2,2) -  Av(1,2));
    norml12 = norm(l12);
    l2 = l12*n14';

    gamma1Av = asin(l2/norml12);
    
    
    gammadot1DIC = (180/pi)*(gamma1B - gamma1Av)/dt;
    
 
    % verification 
    
    
    l32(1)= B(2,1) -  B(3,1);
    l32(2)= B(2,2) -  B(3,2);
    n32 = l32/norm(l32);
    
    l34(1)= B(4,1) -  B(3,1);
    l34(2)= B(4,2) -  B(3,2);
    norml34 = norm(l34);
    
    l2 = l34*n32';

    gamma2B = asin(l2/norml34);

    l32(1)= Av(2,1) -  Av(3,1);
    l32(2)= Av(2,2) -  Av(3,2);
    n32 = l32/norm(l32);
    
    l34(1)= Av(4,1) -  Av(3,1);
    l34(2)= Av(4,2) -  Av(3,2);
    norml134 = norm(l34);
    
    l2 = l34*n32';


    gamma2Av = asin(l2/norml12);
    
    
    gammadot2DIC = (180/pi)*(gamma2B - gamma2Av)/dt;
        
    
%     fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(Ap(1,1)-B(1,1))*pxl2m,(Ap(1,2)-B(1,2))*pxl2m,measured_gammadot); % surf should be 'measured'
%     fprintf(fid, '1 0 1    %12.8f %12.8f %12.8f\n',(Ap(2,1)-B(2,1))*pxl2m,measuredFtopy,0); 
%     fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(Ap(3,1)-B(3,1))*pxl2m,(Ap(3,2)-B(3,2))*pxl2m,measured_gammadot);
%     fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(Ap(4,1)-B(4,1))*pxl2m,(Ap(4,2)-B(4,2))*pxl2m,0); 

    fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(B(1,1)-Av(1,1))*pxl2m,(B(1,2)-Av(1,2))*pxl2m,gammadot1DIC); % surf should be 'measured'
    %fprintf(fid, '1 0 1    %12.8f %12.8f %12.8f\n',(B(2,1)-Av(2,1))*pxl2m,-2*measuredFtopy,0);
    fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(B(2,1)-Av(2,1))*pxl2m,(B(2,2)-Av(2,2))*pxl2m,0);
    fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(B(3,1)-Av(3,1))*pxl2m,(B(3,2)-Av(3,2))*pxl2m,gammadot2DIC);
    fprintf(fid, '1 1 1    %12.8f %12.8f %12.8f\n',(B(4,1)-Av(4,1))*pxl2m,(B(4,2)-Av(4,2))*pxl2m,0); 


    fprintf(fid, '\n');
    fclose(fid);


disp('========== END');
