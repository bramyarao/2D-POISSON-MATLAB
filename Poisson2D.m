%=========================================================
% NUMERICAL ANALYSIS: USING MESHFREE METHOD

% USING THE REPRODUCING KERNEL COLLOCATION METHOD TO SOLVE 
% THE 2D POISSONS PROBLEM
%=========================================================
clc
clc
close all 

tic;


%-------------------------
%INPUT PARAMETERS
%-------------------------
marker_size = 60; % For scatter plot
showPlot = true; %Plotting is done if true
printStatements = false; %Printing is done if true

%Domain
xdim1=0;
xdim2=1;
ydim1=0;
ydim2=1;

num_pts = [10 10 40 40];%[10 10 40 40]  [5 5 5 5]
% No. of Source points in the each direction
NS_x = num_pts(1); %No. of Source points in the x-direction
NS_y = num_pts(2); %No. of Source points in the y-direction

% No. of Collocation points in the each direction
CP_x  = num_pts(3); %No. of Collocation points in the x-direction
CP_y  = num_pts(4); %No. of Collocation points in the y-direction

%-------------------------
% SOURCE POINTS
%-------------------------
[NS] = forming_NS_NC.source(xdim1, xdim2, ydim1, ydim2, NS_x, NS_y);

%-------------------------
% COLLOCATION POINTS
%-------------------------
[NC,NI_c,NEB] = forming_NS_NC.collocation(xdim1, xdim2, ydim1, ydim2, CP_x, CP_y);

%-------------------------------------------------------------------------
% ----------------------------POISSONS----------------------------------
%-------------------------------------------------------------------------
basis = 2;   % Code only works for quadratic basis
% No. of source points and total number of collocation points
no_NS = size(NS,1);
no_NC = size(NC,1);
no_NEB = size(NEB,1);
h = 1/(sqrt(no_NS)-1);
ss = (basis+1)*h; % Support size for the RK SF

sq_alphag = no_NS;  % Weight for the essential boundary
sq_alphah = 1; % Weight for the natural boundary
%-------------------------------------------------------------------------

% Solving the differntial equation
% the A matriz will be of size no_NC x no_NS since u(x,y) is a scalar

%-------------------------
% Forning A matrix
%-------------------------
[A1] = part_of_NI(NC,NS,ss);
[A2] = part_of_NEB(NEB,NS,ss,sq_alphag);

A = [A1;A2]; % A matrix is formed

%-------------------------
% Forming b vector 
%-------------------------
 
b = zeros(no_NC+no_NEB,1);
no_int = size(NC,1); % No. of interior points
int_1 = 1; % Position counter for the b matrix

% INTERIOR points force term
for int_2 = 1:no_int
    xtemp = NC(int_2,1);
    ytemp = NC(int_2,2);        
    b(int_1) = b(int_1)+ (xtemp^2 + ytemp^2)*exp(xtemp*ytemp);         
    int_1 = int_1+1;
end
% int_1 is getting incremented
% EB points
for int_2 = 1:size(NEB,1)
    xtemp = NEB(int_2,1);
    ytemp = NEB(int_2,2);        
    b(int_1) = b(int_1)+ sq_alphag*exp(xtemp*ytemp);         
    int_1 = int_1+1;
end
clear int_1 int_2

%-------------------------
% Solving the system
%------------------------
a = A\b;


%-----------------------------------------------------------------------
% Comparing Solutions Along the Diagonal through (0,0) & (1,1)
%-----------------------------------------------------------------------
x_con = xdim1:0.05:xdim2;
y_con = ydim1:0.05:ydim2;

u_con = zeros(length(x_con),1);
u_exact_con = zeros(length(x_con),1);

for int1 = 1:length(x_con) 
    x = x_con(int1);
    y = y_con(int1);
    [P] = required_nodes(x,y,NS,ss);
    [SI] = SF2D.SF_2D(x,y,NS,P,ss);
    
    %for getting uh value of u_approx at each of the points (x,y)

    u_con(int1) = u_con(int1) + SI*a;
        
    % Finding u_exact at the point (x,y)
    u_exact_con(int1) = u_exact_con(int1) + exp(x*y);
end

% Diagonal length
DL = zeros(length(x_con),1);
for int6 = 1:length(x_con)
    DL(int6,1)  = sqrt((x_con(int6))^2 + (y_con(int6))^2);
end


%-----------------------------------------------------------------------
% Plotting the numerical solution
%-----------------------------------------------------------------------
uScatter = zeros(size(NC,1),1);
uScatter_exact = zeros(size(NC,1),1);

for int1 = 1:size(NC,1)
    x = NC(int1,1);
    y = NC(int1,2);
    [P] = required_nodes(x,y,NS,ss);
    [SI] = SF2D.SF_2D(x,y,NS,P,ss);
    
    %for getting uh value of u_approx at each of the points (x,y)
    uScatter(int1) = uScatter(int1) + SI*a;
        
    % Finding u_exact at the point (x,y)
    uScatter_exact(int1) = uScatter_exact(int1) + exp(x*y);
end



%-------------------------
%PRINTING
%-------------------------
if (printStatements == true)
    fprintf('Source points %d \n',size(NS,1));
    fprintf('Collocation points %d \n',size(NC,1));
    fprintf('Interior collocation points %d \n',size(NI_c,1));
    fprintf('EB collocation points %d \n',size(NEB,1));
end

%-------------------------
%PLOTTING
%-------------------------
if (showPlot == true)
    
    %1
    plot1 = figure(1); hold on %Source points
    scatter(NS(:,1),NS(:,2),marker_size,'r','filled')
    xlabel('x','FontSize',20,'FontName','Times New Roman')
    ylabel('y','FontSize',20,'FontName','Times New Roman')
    title('Source Points','FontSize',16)
    axis equal
    hold off
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')
    
    %2
    plot2 = figure(2); hold on %Collocation points
    legd = zeros(1,2);
    hold on
    legd(1,1) = scatter(NI_c(:,1),NI_c(:,2),marker_size,'b','filled');     
    legd(1,2) = scatter(NEB(:,1),NEB(:,2),marker_size,'black','filled');  

    xlabel('x','FontSize',20,'FontName','Times New Roman')
    ylabel('y','FontSize',20,'FontName','Times New Roman')
    title('Collocation Points','FontSize',16)
    axis equal
    hold off
    legend(legd,{'Interior','Essential boundary'})
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')
    
    %3
    plot3 = figure(3);
    plot(DL,u_con,'ok',DL,u_exact_con,'-k','LineWidth',2,'MarkerSize',10)
    % title('Comparision of u along diagonal','FontSize', 16)
    xlabel('Diagonal length','FontSize', 20,'FontName','Times New Roman')
    ylabel('u(x,y)','FontSize', 20,'FontName','Times New Roman')
    legend('RKCM','Analytical')
    set(gca,'FontSize',16)
    set(gca,'FontName','Times New Roman')
    
    %4
    plot4 = figure(4);
    plot(DL,abs(u_con-u_exact_con),'-r','LineWidth',3)
    title('Error of u along diagonal','FontSize', 16,'FontName','Times New Roman')
    xlabel('Diagonal length','FontSize', 20,'FontName','Times New Roman')
    ylabel('u_ex-u^h','FontSize', 20,'FontName','Times New Roman')
    set(gca,'FontSize',16)
    set(gca,'FontName','Times New Roman')
    
    %5
    plot5 = figure(5);
    scatter(NC(:,1), NC(:,2), marker_size, uScatter, 'filled');
    xlabel('x','FontSize',20,'FontName','Times New Roman')
    ylabel('y','FontSize',20,'FontName','Times New Roman')
    title('u(x,y) Numerical (RKCM)','FontSize',16)
    axis equal
    colorbar
    colormap hsv
    %caxis([0 3])
    hold off
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')
    
    %6
    plot6 = figure(6);
    scatter(NC(:,1), NC(:,2), marker_size, uScatter_exact, 'filled');
    xlabel('x','FontSize',20,'FontName','Times New Roman')
    ylabel('y','FontSize',20,'FontName','Times New Roman')
    title('u(x,y) Analytical (Exact)','FontSize',16)
    axis equal
    colorbar
    colormap hsv
    %caxis([0 3])
    hold off
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')
    
end


%-------------------------
%TIME CALC
%-------------------------
toc;


