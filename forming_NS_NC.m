% Functions for creating the source and collocation points
classdef forming_NS_NC
   methods(Static)
       
       function [NS] = source(xdim1, xdim2, ydim1, ydim2, NS_x, NS_y)
           % Constructing the nodes in the x and y directions
            dx_s = (xdim2-xdim1)/(NS_x-1);
            dy_s = (ydim2-ydim1)/(NS_y-1);

            x_s = xdim1:dx_s:xdim2;
            y_s = ydim1:dy_s:ydim2;            

            % Forming all the x,y coordinates of the Source points in the domain
            NP_s = length(x_s)*length(y_s); % Total number of points
            NS = zeros(NP_s,2);        
            int_1=1;
            for count_1 = 1:length(x_s)
                for count_2 = 1:length(y_s)
                    NS(int_1,1) = NS(int_1,1)+ x_s(count_1);
                    NS(int_1,2) = NS(int_1,2)+ y_s(count_2);
                    int_1 = int_1+1;
                end    
            end
       end % func source
       
       function [NC,NI_c,NEB] = collocation(xdim1, xdim2, ydim1, ydim2, CP_x, CP_y)
       
            % Nodal distance of collocation points
            dx_c = (xdim2-xdim1)/(CP_x-1);
            dy_c = (ydim2-ydim1)/(CP_y-1);

            x_c = xdim1:dx_c:xdim2;
            y_c = ydim1:dy_c:ydim2;

            % Forming all the x,y coordinates of the Collocation points in the domain
            NP_c = length(x_c)*length(y_c); % Total no. of Collocation points no_NI
            NC_total = zeros(NP_c,2);        
            int_1=1;
            for count_1 = 1:length(x_c)
                for count_2 = 1:length(y_c)
                    NC_total(int_1,1) = NC_total(int_1,1)+ x_c(count_1);
                    NC_total(int_1,2) = NC_total(int_1,2)+ y_c(count_2);
                    int_1=int_1+1;
                end    
            end

            % Splittng the collocation points into interior and boundary points
            % Constructing the INTERIOR COLLOCATION POINTS NI_c
            NP_in = (length(x_c)-2)*(length(y_c)-2);
            NI_c = zeros(NP_in,2);
            % Constructing Nodes on the Essential Boundary
            NP_EB = (length(x_c)+length(y_c)+length(x_c)+length(y_c)-4);
            NEB = zeros(NP_EB,2);

            int_2 = 1; int_3 = 1;
            for int_1 = 1:NP_c % NP_c = Total no. of Collocation points no_NI
                x_temp = NC_total(int_1,1);
                y_temp = NC_total(int_1,2);
                if (x_temp == xdim1 || x_temp == xdim2 || y_temp == ydim1 || y_temp == ydim2)
                    NEB(int_3,1) = NEB(int_3,1) + x_temp;
                    NEB(int_3,2) = NEB(int_3,2) + y_temp;
                    int_3=int_3+1;
                else
                    NI_c(int_2,1)= NI_c(int_2,1) + x_temp;
                    NI_c(int_2,2)= NI_c(int_2,2) + y_temp;
                    int_2=int_2+1;
                end
            end
            
            % FORMING THE COLLOCATION POINTS MATRIX 
            NC = [NI_c;NEB];
       end % func collocation
       
   end % method
end %class


