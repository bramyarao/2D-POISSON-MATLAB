%=========================================================
% USING THE REPRODUCING KERNEL COLLOCATION METHOD TO SOLVE 
% THE 2D POISSONS PROBLEM
%=========================================================
clc
clc
close all 
tic
timerval = tic;

%-------------------------
%INPUT PARAMETERS
%-------------------------
marker_size = 60; % For scatter plot
showPlot = true; %Plotting is done if true
printStatements = true; %Printing is done if true

%Domain
xdim1=0;
xdim2=1;
ydim1=0;
ydim2=1;

% No. of Source points in the each direction
NS_x = 10; %No. of Source points in the x-direction
NS_y = 10; %No. of Source points in the y-direction

% No. of Collocation points in the each direction
CP_x  = 20; %No. of Collocation points in the x-direction
CP_y  = 20; %No. of Collocation points in the y-direction

%-------------------------
% SOURCE POINTS
%-------------------------
[NS] = forming_NS_NC.source(xdim1, xdim2, ydim1, ydim2, NS_x, NS_y);

%-------------------------
% COLLOCATION POINTS
%-------------------------
[NC,NI_c,NEB] = forming_NS_NC.collocation(xdim1, xdim2, ydim1, ydim2, CP_x, CP_y);


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
    for int1 = 1:size(NS,1)
        scatter(NS(int1,1),NS(int1,2),marker_size,'r','filled')    
    end

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
    
end

