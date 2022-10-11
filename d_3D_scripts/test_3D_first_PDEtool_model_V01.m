%test
% 3D


clear all;
%close all;
clc;

% 1. create pde system
model = createpde;

% 2. Import CAD geometry
g1 = importGeometry(model,'Halfbreast_tumor_seeds_V05.stl');
rotate(g1,90,[0 -40 0],[0 40 0])
figure(10);
%subplot(1,2,1)
pdegplot(model,'FaceLabels','on','FaceAlpha',0.2);

%subplot(1,2,2)
%pdegplot(model,'CellLabels','on','FaceAlpha',0.2);

% 3. create mesh 
mesh2 = generateMesh(model,'GeometricOrder','linear','Hmin',0.1); %Hmin defines a lighter model, less detailed while Hmax increases the level of analysis
figure(11)
pdemesh(model,'FaceAlpha',0.3);

%%
% 4. define BC and IC
%Boundary conditions 
% are either Neumann ('g') = flux across the boundaries + 'q' = convective
% flux % if not defined = 0
%or Dirichlet ('u') = surfaces held at fixed concentration/temperatures

% -- NOT POSSIBLE TO SET (??) DUE TO HIGHER TEMP. BODY (??) -- % Possible,
% maximum amount of heat able to conduct out, due to parameters
flux_out = 1700; % 85% of metabolic breast heat generation heat per area W/m2     % Heat flux out of skin found [source]
flux_body = 600; % 30% from metabolic breast heat generation heat generation body
flux_seed = -10000; % heat per area       % randomly chosen ?
flux_tumor = 1700;                       % randomly chosen ?
%T_body = 37 ; %body temperature 37 *C

% tumor
applyBoundaryCondition(model,'neumann','face',1,'g',flux_tumor); % 'face' 1 for 3-D % 20 degrees
% breast: 1. outside, 2. attached to body
applyBoundaryCondition(model,'neumann','face',2,'g',flux_body); % 'face' 2 for 3-D % 20 degrees
applyBoundaryCondition(model,'neumann','face',3,'g',flux_out); % 'face' 2 for 3-D % 20 degrees
% seeds: 1. wall, 2. top, 3. bottom
no_seeds = 6;
no_faces_seeds = no_seeds*3; % each seed is cilinder with 3 faces: wall, top and bottom
applyBoundaryCondition(model,'neumann','face',[4:no_faces_seeds],'g',flux_seed); % 'face' 3 for 3-D
% applyBoundaryCondition(model,'neumann','face',4,'g',flux_seed); % 'face' 4 for 3-D
% applyBoundaryCondition(model,'neumann','face',5,'g',flux_seed); % 'face' 5 for 3-D 
% % seed 4: 1. wall, 2. top, 3. bottom
% applyBoundaryCondition(model,'neumann','face',6,'g',flux_seed); % 'face' 5 for 3-D 
% applyBoundaryCondition(model,'neumann','face',7,'g',flux_seed); % 'face' 2 for 3-D
% applyBoundaryCondition(model,'neumann','face',8,'g',flux_seed); % 'face' 3 for 3-D
% % seed 5: 1. wall, 2. top, 3. bottom
% applyBoundaryCondition(model,'neumann','face',9,'g',flux_seed); % 'face' 4 for 3-D
% applyBoundaryCondition(model,'neumann','face',10,'g',flux_seed); % 'face' 5 for 3-D 
% applyBoundaryCondition(model,'neumann','face',11,'g',flux_seed); % 'face' 5 for 3-D 
% % seed 6: 1. wall, 2. top, 3. bottom
% applyBoundaryCondition(model,'neumann','face',12,'g',flux_seed); % 'face' 1 for 3-D
% applyBoundaryCondition(model,'neumann','face',13,'g',flux_seed); % 'face' 2 for 3-D
% applyBoundaryCondition(model,'neumann','face',14,'g',flux_seed); % 'face' 1 for 3-D
% % seed 7: 1. wall, 2. top, 3. bottom
% applyBoundaryCondition(model,'neumann','face',15,'g',flux_seed); % 'face' 2 for 3-D
% applyBoundaryCondition(model,'neumann','face',16,'g',flux_seed); % 'face' 1 for 3-D
% applyBoundaryCondition(model,'neumann','face',17,'g',flux_seed); % 'face' 2 for 3-D
% % seed 8: 1. wall, 2. top, 3. bottom
% applyBoundaryCondition(model,'neumann','face',18,'g',flux_seed); % 'face' 2 for 3-D
% applyBoundaryCondition(model,'neumann','face',19,'g',flux_seed); % 'face' 1 for 3-D
% applyBoundaryCondition(model,'neumann','face',20,'g',flux_seed); % 'face' 2 for 3-D

% Initial condition
% find out if time dependent or not: if not --> elliptic eq.
                                    % if    --> parabolic eq.
Tini_breast = 32.8; % degrees: initial body temperature of breast
Tini_tumour = 37; % degrees: initial body temperature of breast
T_TA = 50; % degrees: temperature of seeds after activation magnetic field
setInitialConditions(model,Tini_breast);
setInitialConditions(model,T_TA,'Cell',(3:8));
%setInitialConditions(singeldomain,'cell',2,T0);

% 5. put in properties
T_end = 3*60*60; 
tlist = 0:3*60:T_end; % times at which to compute a solution specified in the array tlist
                    % between 0 and 60 minutes in steps of 60 secs
% specify the coefficients for Poissons equation in the model
% m * d2u/dt2 + d * du/dt - nabla(c*nabla(u) + a * u )= f

%possible to define per face or per cell the physical properties
%breast 
f_b = 2000; % Q = heat source/sink / metabolic heat generation
%f_b = 0; % Q = heat source/sink / metabolic heat generation
a_b = 0; % convective flux % rho*Cp*v = a = 0
c_b = 0.48 ; % c = k = diffusion/conduction coefficient % approximately body
d_b = (1050*3770); % 0 if time independent % d = rho*Cp
specifyCoefficients(model,'Cell',2,'m',0,'d',d_b,'c',c_b,'a',a_b,'f',f_b);

% tumor
f_t = 10000; % Q = heat source/sink / metabolic heat generation
%f_t = 0; % Q = heat source/sink / metabolic heat generation
a_t = 0; % convective flux % rho*Cp*v = a = 0
c_t = 0.48 ; % c = k = diffusion/conduction coefficient % approximately body
d_t = (1050*3852); % 0 if time independent % d = rho*Cp
specifyCoefficients(model,'Cell',1,'m',0,'d',d_t,'c',c_t,'a',a_t,'f',f_t);

% seeds
f_s = 785000; % Q = heat source/sink / metabolic heat generation
a_s = 0; % convective flux % rho*Cp*v = a = 0
c_s = 0.48; % c = k = diffusion/conduction coefficient % approximately body
d_s = (900*1500000); % 0 if time independent % d = rho*Cp
specifyCoefficients(model,'Cell',(3:8),'m',0,'d',d_s,'c',c_s,'a',a_s,'f',f_s);


% 6. solve PDE system
results = solvepde(model,tlist);
T = results.NodalSolution;

%% Temperature plots
clc;

% Check intervals of time
% for i = 0:4
%     for j = 1:4 
%        column(i+1,j) = (1+(j-1)+(4*i))*3
%        time(i+1,j) = (1+(j-1)+(4*i))*3*3
%     end
% end


for i = 0:4
    for j = 1:4 
        figure(i+1+5);
        subplot(2,2,j);
        pdeplot3D(model,"ColorMapData",T(:,((1+(j-1)+(4*i))*3)),'FaceAlpha',0.4);
        xlim([5 15]);
        ylim([5 15]);
        zlim([-15 -5]);
        caxis([30 65]);
        title(['Temperature at time =', num2str((1+(j-1)+(4*i))*3*3),' min'])
        axis on
    end
end

%% Practice with 2D visualisation of results
% clc;
% [cgradx,cgrady,cgradz] = evaluateCGradient(results);
% 
% for i = 1:10:60
%     figure
%     pdeplot3D(model,'FlowData',[cgradx(:,i) cgrady(:,i) cgradz(:,i)]) %shows the direction en magnitude of flux, larger flux, larger longer arrow
%     xlabel('x')
%     ylabel('y')
%     xlim([5 15]);
%     ylim([5 15]);
%     zlim([-15 -5]);
% end

%% Trying to get a 2D slice

[X_resultsXY,Y_resultsXY] = meshgrid(mesh2.Nodes(1,:),mesh2.Nodes(2,:));
[X_resultsXZ,Z_resultsXZ] = meshgrid(mesh2.Nodes(1,:),mesh2.Nodes(3,:));
[Y_resultsYZ,Z_resultsYZ] = meshgrid(mesh2.Nodes(2,:),mesh2.Nodes(3,:));
x = 1;
y = 2;
z = 3;

for j = 1:3
    figure(j)
    for i = 1:10:61
        %surf(X_resultsXZ,Z_resultsXZ,Temperature_per_element)
        plot(mesh2.Nodes(j,:),results.NodalSolution(:,i),'.')
        xlim([-5 55])
        xlabel('X direction')
        ylabel('Temperature (*C)')
        hold on
    end
end


% figure;
% plot3(mesh2.Nodes(1,:),mesh2.Nodes(3,:),results.NodalSolution(:,40),'-')
% xlabel("x axis");
% ylabel("z axis");
% zlabel("Temperature");
% 
% figure;
% plot3(mesh2.Nodes(1,:),mesh2.Nodes(2,:),results.NodalSolution(:,40),'-')
% xlabel("x axis");
% ylabel("y axis");
% zlabel("Temperature");
% 
% figure;
% plot3(mesh2.Nodes(2,:),mesh2.Nodes(3,:),results.NodalSolution(:,40),'-')
% xlabel("y axis");
% ylabel("z axis");
% zlabel("Temperature");
% xlim([5 15]);
% %ylim([5 15]);
% zlim([5 15]);
% caxis([20 50]);
% figure(12);
% %subplot(2,2,j);
% pdeplot3D(model,"ColorMapData",T(:,20),'FaceAlpha',0.2);
% xlim([5 15]);
% ylim([5 15]);
% zlim([5 15]);
% caxis([20 50]);
% axis on
% 
% figure(13);
% %subplot(2,2,j);
% pdeplot3D(model,"ColorMapData",T(:,1),'FaceAlpha',0.2);
% xlim([5 15]);
% ylim([5 15]);
% zlim([5 15]);
% caxis([20 50]);
% axis on
% subplot(2,2,4)
% pdeplot3D(model,"ColorMapData",T(:,18),'FaceAlpha',0.3);

%% Old version
%Boundary conditions are either Neumann ('g') = flux across the boundaries 
%or Dirichlet ('u') = surfaces held at fixed temperatures
% applyBoundaryCondition(model,'face',1,'g',-0.1); % 'face' 1 for 3-D
% applyBoundaryCondition(model,'face',2,'g',-1); % 'face' 2 for 3-D
% applyBoundaryCondition(model,'face',3,'g',-0.1); % 'face' 3 for 3-D
% applyBoundaryCondition(model,'face',4,'g',-0.1); % 'face' 4 for 3-D
% applyBoundaryCondition(model,'face',5,'g',-1); % 'face' 5 for 3-D 
% applyBoundaryCondition(model,'face',6,'g',1); % 'face' 5 for 3-D 
% 
% applyBoundaryCondition(model,'face',7,'g',-0.1); % 'face' 2 for 3-D
% applyBoundaryCondition(model,'face',8,'g',-1); % 'face' 3 for 3-D
% applyBoundaryCondition(model,'face',9,'g',-0.1); % 'face' 4 for 3-D
% applyBoundaryCondition(model,'face',10,'g',1); % 'face' 5 for 3-D 
% applyBoundaryCondition(model,'face',11,'g',1); % 'face' 5 for 3-D 
% 
% f = 0;
% a = 0;
% c = 1;
% 
% % specify the coefficients for Poissons equation in the model
% % m * d2u/dt2 + d * du/dt - nabla(c*nabla(u) + a * u = f
% specifyCoefficients(model,'m',0,'d',1,'c',c,'a',a,'f',f);
% setInitialConditions(model,350);
% 
% 
% tlist = 0:1:20;
% results = solvepde(model,tlist);
% 
% u = results.NodalSolution
% 
% 
% subplot(2,2,2);
% pdeplot3D(model,"ColorMapData",u(:,18));
% 
% % pdegplot(model,'FaceLabels','on') % 'FaceLabels' for 3-D 
% % %pdegplot(model,'EdgeLabels','on') % 'EdgeLabels' for 2-D
% 
% % changes made to practice with Git
% %y = linspace(2,6,100)
% 
% % lets see what this does if it is the only thi g changed
