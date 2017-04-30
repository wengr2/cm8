%% Continuum Mechanics, Exercise 8
% Paul Kulyk
% Raphael Wenger
%
% paul.kulyk@students.unibe.ch
% raphael.wenger@students.unibe.ch
%
% Due Mai 2, 2017

%% Include the predefined tensor math functions
addpath('../../Matlab/');
if exist('imported','var') ~= 1
    def_symbols;
    imported = 1;
end
%% 8.1 Stress vectors
% define the geometry and build the tetrahedron
Xi = [[1/10 0 0]; [0 1/10 0]; [0 0 1/10]; [0 0 0]; ]';
% density
rho = 1000; % kg/m^3
g = 10; % m/s^2 (close enough...)
syms T Tmax t

% Rotation around e1,e2,e3 respectively
R1 = @(theta) [ [ 1 0 0 ]; [ 0 cos(theta) -sin(theta) ]; [ 0 sin(theta) cos(theta) ];];
R2 = @(theta) [ [ cos(theta) 0 -sin(theta) ]; [ 0 1 0 ]; [ sin(theta) 0 cos(theta) ];];
R3 = @(theta) [ [ cos(theta) -sin(theta) 0 ]; [ sin(theta) cos(theta) 0 ]; [ 0 0 1 ];];

%Motion of the tetrahedron as in exercice 7
Tmax = 1/2;
T = t; %Back to symbolics as the master said
Rt = R3(2*pi*T/Tmax);
bt = [ 0 0 3/20*T/Tmax]';
y =@(R,x,b) R*x + b;

%Contact force defined with handles
F_contact  =  @(theta) [-pi^2/15*(cos(theta)-sin(theta)); -pi^2/15*(cos(theta)+sin(theta));5/3];
F_contact0 =  @(fact,theta) fact*[cos(theta)+sin(theta);cos(theta)-sin(theta);0];

%Apply the given values
F_con  = F_contact(4*pi*T);
F_con0 = F_contact0(125+pi^2*(4+60*T)/3000,4*pi*T);

% Function handle to clean up all these triple integrals
TripInt =@(fun,v1,l1,u1,v2,l2,u2,v3,l3,u3) int( int( int( fun, v1, l1, u1), v2, l2, u2), v3, l3, u3);
 
%% Forces vectors on the verticles  
COG     =@(xx) 6*V/M*TripInt(rho*xx,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
faceArea= @(curNode1,curNode2,curNode3) cm.norm(cm.cross_product((curNode1-curNode2),(curNode3-curNode2)));  


yt = y(Rt,Xi(1:3,1:3)*b,bt); % Simple transform of b into x into y...
V = 1/6000; % same result as above but keep it symbolic
M = 6*V*TripInt(rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);          
yc = COG(yt);

%Final positions
yi(:,1) = y(Rt,Xi(:,1),bt);
yi(:,2) = y(Rt,Xi(:,2),bt);
yi(:,3) = y(Rt,Xi(:,3),bt);
yi(:,4) = y(Rt,Xi(:,4),bt);
yi

%Get the normals to the faces, centers and area
Ai(1)=faceArea(yi(:,3),yi(:,2),yi(:,4));
Ai(2)=faceArea(yi(:,4),yi(:,3),yi(:,1));
Ai(3)=faceArea(yi(:,1),yi(:,4),yi(:,2));
Ai(4)=faceArea(yi(:,2),yi(:,1),yi(:,3));
Ai

%Use the providen function for the normals to the surface
[ynormi ycenti] = cm.get_tetra_normal(yi(:,1),yi(:,2),yi(:,3),yi(:,4))

%Computation of the Area-weighted normas
for i = 1:4
    Aini(:,i) = Ai(i)*ynormi(:,i);
end

%Only non-singular matrix is @i=4
i=4;
wt = -4*cm.invert(cm.scalar_product(ycenti(:,i),Aini(:,i))*I+cm.dyadic_product11(ycenti(:,i),Aini(:,i)))*(F_con0-cm.cross_product(yc,F_con));

                          
% Calculate the linear momentum
lt  = 6*V*TripInt(dydt*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);

%% 7.1.(1) verify linear momentum is velocity of the center of mass times the total mass
M           = 6*V*TripInt(rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
% Center of gravity for the velocity of y
COG =@(xx) 6*V/M*TripInt(rho*xx,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
dydt_cog = COG(dydt);
cm.show1(simplify(lt - M*dydt_cog));

% Compute angular momentum with respect to the origin
l_hat_0 = simplify(6*V*TripInt(cm.cross_product(yt,dydt)*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1));

% Calculate the angular velocity
dRdt = diff(Rt,T);
% Put the epsilon tensor together
epsilon = levi;
omega = simplify(-1/2*cm.transform_32(epsilon,cm.composition22(dRdt,cm.transpose2(Rt))));

% Know the tensor of inertia in case of a rotation
X = Xi(1:3,1:3)*b;
Xcog = COG(X);
%Xcog    = 6*V/M*TripInt(rho*X,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);

Jc = 6*V*TripInt(((cm.norm(X - Xcog)^2*I)-cm.dyadic_product11((X-Xcog),(X-Xcog)))*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
Jc_t   = simplify(Rt*Jc*cm.transpose2(Rt));
yt_cog = COG(yt);

%% 7.1.(2)  Verify the angular momentum about the origin
cm.show1(simplify(l_hat_0 - (Jc_t*omega + M*cm.cross_product(yt_cog,dydt_cog))));

%% 7.2 Inertial Forces and Moments

%%% Calculate the inertial force
f_ine = diff(lt,T);
%f_ine_alt  = 6*V*TripInt(d2yd2t*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
%cm.show1(simplify(f_ine_alt-f_ine));

%%% Calculate the inertial moment
fhat_ine_0 = diff(l_hat_0,T);
%fhat_ine_0_alt = 6*V*TripInt(cm.cross_product(yt,d2yd2t)*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
%cm.show1(simplify(fhat_ine_0_alt-fhat_ine_0));

%% 7.2.(3)
d2yd2t_cog = diff(diff(yt_cog,T),T);
cm.show1(simplify(M*d2yd2t_cog - f_ine));

%% 7.2.(4)
cm.show1(simplify(fhat_ine_0 - (diff(Jc_t,T)*omega + Jc_t*diff(omega,T) + M*cm.cross_product(yt_cog ,d2yd2t_cog))));

%% 7.3 Volume forces and moments

% Calculate the gravitational force on the tetrahedron
f_vol = -6*V*TripInt(rho*g*e3,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);

%% 7.3.(5) Gravitational volume force at center of mass
cm.show1(simplify(f_vol + M*g*e3));

%% 7.3.(6) Moment due to gravity
fhat_vol_0 = -6*V*TripInt(cm.cross_product(yt,rho*g*e3),b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
cm.show1(simplify(fhat_vol_0 + cm.cross_product(yt_cog,M*g*e3)));

%% 7.4 Plots

%Tetrahedron coordinates
Yi = zeros(3,4);
for i = 1:4
    Yi(1:3,i) = eval(subs(y(Rt,Xi(:,i),bt),T,Tmax));
end
Yt_cog = eval(subs(yt_cog,T,Tmax));
F_ine = eval(subs(f_ine,T,Tmax)); % seems wrong
F_vol = eval(subs(f_vol,T,Tmax)); 
Fhat_ine_0 = eval(subs(fhat_ine_0,T,Tmax)); 
Fhat_vol_0 = eval(subs(fhat_vol_0,T,Tmax)); 

scale = 1/10;



%% 7.4.(7) Inertial force on tetrahedron
% Show me the inertia
figure();
% Show me the tetrahedrae
cm.plot_tetra(Yi(:,1),Yi(:,2),Yi(:,3),Yi(:,4));
hold on;
cm.plot_tetra(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
% Show me some Centers of Gravity
cm.plot_tetra_point(Yt_cog(1),Yt_cog(2),Yt_cog(3),'blue');
cm.plot_tetra_point(Xcog(1),Xcog(2),Xcog(3),'blue');
cm.plot_vector(Yt_cog,Yt_cog+F_ine*scale,2,'m'  );

%% 7.4.(8) Inertial force about the origin
% Show me the moment (about the origin)
figure();
% Show me the tetrahedrae
cm.plot_tetra(Yi(:,1),Yi(:,2),Yi(:,3),Yi(:,4));
hold on;
cm.plot_tetra(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
% Show me some Centers of Gravity
cm.plot_tetra_point(Yt_cog(1),Yt_cog(2),Yt_cog(3),'blue');
cm.plot_tetra_point(Xcog(1),Xcog(2),Xcog(3),'blue');
cm.plot_vector(zeros(3,1),Fhat_ine_0*scale,2,'m'  );

%% 7.4.(9) Gravity force on tetrahedron
% Show me the gravity
figure();
% Show me the tetrahedrae
cm.plot_tetra(Yi(:,1),Yi(:,2),Yi(:,3),Yi(:,4));
hold on;
cm.plot_tetra(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
% Show me some Centers of Gravity
cm.plot_tetra_point(Yt_cog(1),Yt_cog(2),Yt_cog(3),'blue');
cm.plot_tetra_point(Xcog(1),Xcog(2),Xcog(3),'blue');
cm.plot_vector(Yt_cog,Yt_cog+F_vol*scale,2,'m'  );

%% 7.4.(10) Gravity moment on tetrahedron about origin
% Show me the moment due to gravity about the origin
figure();
scale = 1/2; % Scale this up just to visualize better
% Show me the tetrahedrae
cm.plot_tetra(Yi(:,1),Yi(:,2),Yi(:,3),Yi(:,4));
hold on;
cm.plot_tetra(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
% Show me some Centers of Gravity
cm.plot_tetra_point(Yt_cog(1),Yt_cog(2),Yt_cog(3),'blue');
cm.plot_tetra_point(Xcog(1),Xcog(2),Xcog(3),'blue');
cm.plot_vector(zeros(3,1),Fhat_vol_0*scale,2,'m'  );
scale = 1/10;
%% 7.4.(v2) Plotall

% Show me the inertia
figure();
% Show me the tetrahedrae
cm.plot_tetra(Yi(:,1),Yi(:,2),Yi(:,3),Yi(:,4));
hold on;
cm.plot_tetra(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
% Show me some Centers of Gravity
cm.plot_tetra_point(Yt_cog(1),Yt_cog(2),Yt_cog(3),'blue');
cm.plot_tetra_point(Xcog(1),Xcog(2),Xcog(3),'blue');
cm.plot_vector(Yt_cog,Yt_cog+F_ine*scale,2,'m'  );

% Show me the moment (about the origin)
hold on;
cm.plot_vector(zeros(3,1),Fhat_ine_0*scale,2,'g'  );

% Gravity force on tetrahedron
hold on;
cm.plot_vector(Yt_cog,Yt_cog+F_vol*scale,2,'r'  );

% Show me the moment due to gravity about the origin
hold on;
cm.plot_vector(zeros(3,1),Fhat_vol_0*scale,2,'c'  );

%% 7.5 Homogenous Deformation
% If volume is not preserved may need to recalculate V

F =@(T,Tmax) [[ (1+T/Tmax) 3^0.5*T/Tmax*(1+T/Tmax) 0 ]; [0 (1+T/Tmax) 0]; [0 0 (1+T/Tmax)*(1+2*T/Tmax)]; ];
Ft = F(T,Tmax);
y2 =@(R,F,x,b) R*F*x+b;
yt = y2(Rt,Ft,Xi(1:3,1:3)*b,bt);
dydt = diff(yt,T);
Yi = zeros(3,4);
for i = 1:4
    Yi(1:3,i) = eval(subs(y2(Rt,Ft,Xi(:,i),bt),T,Tmax));
end 
lt  = 6*V*TripInt(dydt*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
%Yt_cog = eval(subs(subs(subs(subs(COG(yt),T,Tmax),b1,0.1),b2,0.1),b3,0.1));
Yt_cog = eval(subs(COG(yt),T,Tmax));
F_ine = eval(subs(diff(lt,T),T,Tmax));
l_hat_0 = simplify(6*V*TripInt(cm.cross_product(yt,dydt)*rho,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1));
Fhat_ine_0 = eval(subs(diff(l_hat_0,T),T,Tmax));
F_vol = -6*V*TripInt(rho*g*e3,b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);% Probably not correct
fhat_vol_0 = -6*V*TripInt(cm.cross_product(yt,rho*g*e3),b1,0,1-b2-b3,b2,0,1-b3,b3,0,1);
Fhat_vol_0 = eval(subs(fhat_vol_0,T,Tmax));


%%
figure(2);
% Show me the tetrahedrae
hold on;
cm.plot_tetra(Yi(:,1),Yi(:,2),Yi(:,3),Yi(:,4));
hold on;
cm.plot_tetra(Xi(:,1),Xi(:,2),Xi(:,3),Xi(:,4));
% Show me some Centers of Gravity
hold on;
cm.plot_tetra_point(Yt_cog(1),Yt_cog(2),Yt_cog(3),'blue');
% Show me the inertia
hold on;
cm.plot_vector(Yt_cog,Yt_cog+F_ine*scale,2,'m'  );
% Show me the moment (about the origin)
hold on;
cm.plot_vector(zeros(3,1),Fhat_ine_0*scale,2,'m'  );
% Show me the gravity
hold on;
cm.plot_vector(Yt_cog,Yt_cog+F_vol*scale,2,'m'  );
% Show me the moment due to gravity about the origin
hold on;
cm.plot_vector(zeros(3,1),Fhat_vol_0*scale,2,'m'  );
% Show me the moment due to gravity about the origin
