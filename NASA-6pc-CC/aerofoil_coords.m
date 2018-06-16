%
%   See NASA-TM-2005-213545, Trailig edge blowing on a two-dimensional
%   six-percent thick elliptical circulation control airfoil up to
%   transonic conditions. MG Alexander, SG Anders, SK Johnson, JP Florance
%   & DF Keller
%
%
%
%   M. Forster, 2013
%

clear all
format long

xmax=0.5;           % half chord length
np=61;              % number points (odd)
cc=0.9;             % truncation for CC wing, cc=1 for pure ellipse.
cc=double(1-2*(1-cc));      % correction for using half chord

h=4;                % h/c = 0.0012, 0.0020, 0.0026, 0.0007
r=3;                % ratio = 1.78, 2.38, 2.98

ratio=[1.78, 2.38, 2.98];              % ratio = coanda ellipse major/minor

disp('chord length =');
disp(2*xmax);
ymax=0.06*xmax;      % 6% thickness
cmax=0.0075*2*xmax;     %max camber = 0.75% chord
z=zeros(1,np);        % zero matrix - for z coordinate, also for x axis on plot

x=linspace(-4.5,4.5,np-2);
x=[-xmax,double(xmax*tanh(x(1:(np+1)/2))),double(cc*xmax*tanh(x((np+1)/2+1:np-2))),cc*xmax];    % exponential distribution of points about LE and TE
y=ymax*sqrt(1-(x/xmax).^2);             % eqn of ellipse

%elliptical camber
% camber=cmax*sqrt(1-(x/xmax).^2);        % eqn of camber
%
%circular camber
y1=(cmax^2-xmax^2)/(2*cmax);
% x^2+(y+y1)^2=r^2 -> (cmax+y1)^2=r^2 && xmax^2+y1^2=r^2
% -> cmax^2 + 2*cmax*y1 = xmax^2
rsq=xmax^2+y1^2;
camber=sqrt(rsq-x.^2)+y1;
theta=atand(-xmax/y1)

u=[x;camber+y;z];
l=[x;camber-y;z];

plot(x,u(2,:))
hold on
plot(x,l(2,:))
plot(x,camber,'--g')
plot(x,z,'k')

% 151124 new
%beta=linspace(0,pi,np);
%xc=-cos(beta)/2;
xc=x;

yt=ymax*sqrt(1-(xc/xmax).^2);
y1=(cmax^2-xmax^2)/(2*cmax);
rsq=xmax^2+y1^2;
camber=sqrt(rsq-xc.^2)+y1;

theta = atan(xc./(y1+camber)); % angle wrt y axis
unew = [xc-yt.*sin(theta);camber + yt.*cos(theta);z];
lnew = [xc+yt.*sin(theta);camber - yt.*cos(theta);z];

plot(unew(1,:),unew(2,:),'r')
plot(lnew(1,:),lnew(2,:),'r')

%
axis equal


fileID = fopen('151124upper_points.dat','w');
fprintf(fileID,'%3.0f %1s\n',np,'2');
fprintf(fileID,'%10f %10f %10f\n',unew);
fclose(fileID);

fileID = fopen('151124lower_points.dat','w');
fprintf(fileID,'%3.0f %1s\n',np,'2');
fprintf(fileID,'%10f %10f %10f\n',lnew);
fclose(fileID);

if isequal(cc,1)
    disp('unblown aerofoil')
else
    disp('cc aerofoil')

%%%%  COANDA SURFACE  %%%%%


rs=0.91*xmax/30;
rte=ratio*rs;


sh=0.03*2*xmax/30;            %[.001246,.001994,.002599,.000748]
thick=y(np)-rs-sh;

xe=linspace(x(np)-rte(r),x(np)+rte(r),np);
ye=rs*sqrt(1-((xe-x(np))/rte(r)).^2);             % eqn of ellipse

% ellipse profile

% plot(xe,camber(np)+ye)
% plot(xe,camber(np)-ye)

ue=[xe;camber(np)+ye;z];
le=[xe;camber(np)-ye;z];

% plenum wall profile

xp=linspace(xmax*(0.72*2-1),x(np),np);

plot(xp,camber(np)+y(np)-thick*ones(1,np),'g')
plot(xp,camber(np)-y(np)+thick*ones(1,np),'g')

up=[xp;camber(np)+y(np)-thick*ones(1,np);z];
lp=[xp;camber(np)-y(np)+thick*ones(1,np);z];

% fileID = fopen('upper_skin.dat','w');
% fprintf(fileID,'%3.0f %1s\n',np,'2');
% fprintf(fileID,'%10f %10f %10f\n',up);
% fclose(fileID);
%
% fileID = fopen('lower_skin.dat','w');
% fprintf(fileID,'%3.0f %1s\n',np,'2');
% fprintf(fileID,'%10f %10f %10f\n',lp);
% fclose(fileID);


% inner plenum wall

F=[0.85,0.8,0.75];

% plot(xp,camber(np)+y(np)-thick*ones(1,np)-0.009*2*xmax,'m')
% plot(xp,camber(np)-y(np)+thick*ones(1,np)+0.009*2*xmax,'m')

u1=[xp(1:F(r)*np),xe(0.15*np:np+1);camber(np)+y(np)-thick*ones(1,F(r)*np)-0.009*2*xmax,camber(np)+ye(0.15*np:np+1);z(1:F(r)*np),z(0.15*np:np+1)];
u2=[xp(1:F(r)*np),xe(0.15*np:np+1);camber(np)-y(np)+thick*ones(1,F(r)*np)+0.009*2*xmax,camber(np)-ye(0.15*np:np+1);z(1:F(r)*np),z(0.15*np:np+1)];

plot(u1(1,:),u1(2,:),'m')
plot(u2(1,:),u2(2,:),'m')

% fileID = fopen('upper_coanda.dat','w');
% fprintf(fileID,'%3.0f %1s\n',(0.85+F(r))*np,'2');
% fprintf(fileID,'%10f %10f %10f\n',u1);
% fclose(fileID);
%
% fileID = fopen('lower_coanda.dat','w');
% fprintf(fileID,'%3.0f %1s\n',(0.85+F(r))*np,'2');
% fprintf(fileID,'%10f %10f %10f\n',u2);
% fclose(fileID);

end
