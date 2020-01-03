classdef DSFyy   
    methods(Static)
        
        function [SIyy] = DSF_yy(x,y,NS,P,ss)
            %this function gives the second derivative of the Shape function w.r.t 'x' 
            % at any point (x,y) centered at node (xI,yI)
            % For linear basis we need 3x3 M matrix and for quadratic basis we need 6x6
            % M matrix in 2D

            P = transpose(P);
            P_len = size(P,2);

            M   =zeros(6,6);
            M_y =zeros(6,6);
            M_yy=zeros(6,6);

            %Evaluation of the Moment matrix, here we take the node positions
            %as xxI and yyI for the summation process to be carried out easily

            for integer_2=1:P_len
                xxI = P(1,integer_2);
                yyI = P(2,integer_2);

                h = [1 ;x-xxI; y-yyI;(x-xxI)^2; (y-yyI)^2; (x-xxI)*(y-yyI)];
                h_y = [0;0;1;0;(2*(y-yyI));(x-xxI)];
                h_yy = [0;0;0;0;2;0];

                %For FINDING PHI
                clear zz
                zz=(sqrt((x-xxI)^2+(y-yyI)^2))/ss;  

                [phi] = DSFyy.phi_eval(zz);
                [dphi] = DSFyy.dphi_eval(x,y,xxI,yyI,zz,ss); 
                [ddphi] = DSFyy.ddphi_eval(x,y,xxI,yyI,zz,ss);      

                M = M + (h*transpose(h)*phi);
                M_y = M_y + (h_y*transpose(h)*phi) + (h*transpose(h_y)*phi) + (h*transpose(h)*dphi);
                M_yy = M_yy + (h_yy*transpose(h)*phi) + (h_y*transpose(h_y)*phi) + (h_y*transpose(h)*dphi)  + (h_y*transpose(h_y)*phi) + (h*transpose(h_yy)*phi) + (h*transpose(h_y)*dphi)+ (h_y*transpose(h)*dphi) + (h*transpose(h_y)*dphi) + (h*transpose(h)*ddphi); 

            end

            %--------------------------------------------------------------------

            %Evaluating inv(M),x i.e derivative of inv of M
            InvM_y = -1*inv(M)*M_y*inv(M);
            InvM_yy = -inv(M)*((M_yy*inv(M))+(2*M_y*InvM_y));

            %--------------------------------------------------------------------

            %MAIN matrices for computing Derivative of Shape function
            SIyy = zeros(1,size(NS,1));
            % Derivatives of H's
            for int_1 = 1:size(NS,1)
                xI = NS(int_1,1);
                yI = NS(int_1,2);

                Ho = [1 ;0;0;0;0;0];
                H = [1 ;x-xI; y-yI;(x-xI)^2; (y-yI)^2; (x-xI)*(y-yI)];  
                H_y = [0;0;1;0;(2*(y-yI));(x-xI)];
                H_yy = [0;0;0;0;2;0];    

                % For finding PHI
                clear z
                z=(sqrt((x-xI)^2+(y-yI)^2))/ss;

                [PHI] = DSFyy.phi_eval(z);
                [DPHI] = DSFyy.dphi_eval(x,y,xI,yI,z,ss); 
                [DDPHI] = DSFyy.ddphi_eval(x,y,xI,yI,z,ss); 

                SIyy(int_1) = SIyy(int_1) + transpose(Ho)*((InvM_yy*H*PHI)+(inv(M)*H_yy*PHI)+(inv(M)*H*DDPHI)...
                                     + (2*InvM_y*H_y*PHI)+(2*InvM_y*H*DPHI)+(2*inv(M)*H_y*DPHI));    
            end
        end
        
        %================= SUB FUNCTIONS ==============================
        function [phi] = phi_eval(zz)
            %USING Quintic B-SPLINE
            if zz>=0 && zz<(1/3)
                phi = (11/20)-(9/2)*zz^2 + (81/4)*zz^4 -(81/4)*zz^5;
            else if zz>=1/3 && zz<2/3
                    phi = (17/40) + (15/8)*zz - (63/4)*zz^2 + (135/4)*zz^3 - (243/8)*zz^4 + (81/8)*zz^5;
                else if zz>=2/3 && zz<1
                        phi = (81/40) - (81/8)*zz + (81/4)*zz^2 - (81/4)*zz^3 + (81/8)*zz^4 - (81/40)*zz^5;
                    else if zz>=1
                            phi = 0;
                        end
                    end
                end
            end
        end

        %Phi derivative wrt 'y'
        function [dphi] = dphi_eval(x,y,xxI,yyI,zz,ss)
            if zz>=0 && zz<(1/3)
                dphi = (81*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2))/(2*ss^4) - (((405*y)/4 - (405*yyI)/4)*((x - xxI)^2 + (y - yyI)^2)^(3/2))/ss^5 - (9*y - 9*yyI)/ss^2;
            else if zz>=1/3 && zz<2/3
                    dphi = (15*(2*y - 2*yyI))/(16*ss*((x - xxI)^2 + (y - yyI)^2)^(1/2)) - (63*(2*y - 2*yyI))/(4*ss^2) + (405*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2)^(1/2))/(8*ss^3) + (405*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2)^(3/2))/(16*ss^5) - (243*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2))/(4*ss^4);
                else if zz>=2/3 && zz<1
                        dphi = (81*(2*y - 2*yyI))/(4*ss^2) - (81*(2*y - 2*yyI))/(16*ss*((x - xxI)^2 + (y - yyI)^2)^(1/2)) - (243*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2)^(1/2))/(8*ss^3) - (81*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2)^(3/2))/(16*ss^5) + (81*(2*y - 2*yyI)*((x - xxI)^2 + (y - yyI)^2))/(4*ss^4);
                    else if zz>=1
                            dphi = 0;
                        end
                    end
                end
            end
        end

        % Phi 2nd derivative wrt 'y'
        function  [ddphi] = ddphi_eval(x,y,xxI,yyI,zz,ss)
            if zz>=0 && zz<(1/3)
                ddphi = (81*(2*y - 2*yyI)^2)/(2*ss^4) - (405*((x - xxI)^2 + (y - yyI)^2)^(3/2))/(4*ss^5) + (81*((x - xxI)^2 + (y - yyI)^2))/ss^4 - 9/ss^2 - (1215*(2*y - 2*yyI)^2*((x - xxI)^2 + (y - yyI)^2)^(1/2))/(16*ss^5);
            else if zz>=1/3 && zz<2/3
                    ddphi = 15/(8*ss*((x - xxI)^2 + (y - yyI)^2)^(1/2)) + (405*((x - xxI)^2 + (y - yyI)^2)^(1/2))/(4*ss^3) + (405*((x - xxI)^2 + (y - yyI)^2)^(3/2))/(8*ss^5) - (243*(y - yyI)^2)/ss^4 - (243*(x^2 - 2*x*xxI + xxI^2 + y^2 - 2*y*yyI + yyI^2))/(2*ss^4) - 63/(2*ss^2) - (15*(y - yyI)^2)/(8*ss*((x - xxI)^2 + (y - yyI)^2)^(3/2)) + (405*(y - yyI)^2)/(4*ss^3*((x - xxI)^2 + (y - yyI)^2)^(1/2)) + (1215*((x - xxI)^2 + (y - yyI)^2)^(1/2)*(y - yyI)^2)/(8*ss^5);
                else if zz>=2/3 && zz<1
                        ddphi = (81*(y - yyI)^2)/ss^4 - (243*((x - xxI)^2 + (y - yyI)^2)^(1/2))/(4*ss^3) - (81*((x - xxI)^2 + (y - yyI)^2)^(3/2))/(8*ss^5) - 81/(8*ss*((x - xxI)^2 + (y - yyI)^2)^(1/2)) + (81*(x^2 - 2*x*xxI + xxI^2 + y^2 - 2*y*yyI + yyI^2))/(2*ss^4) + 81/(2*ss^2) + (81*(y - yyI)^2)/(8*ss*((x - xxI)^2 + (y - yyI)^2)^(3/2)) - (243*(y - yyI)^2)/(4*ss^3*((x - xxI)^2 + (y - yyI)^2)^(1/2)) - (243*((x - xxI)^2 + (y - yyI)^2)^(1/2)*(y - yyI)^2)/(8*ss^5);
                    else if zz>=1
                            ddphi = 0;
                        end
                    end
                end
            end
        end

    end %mthods
end %class

