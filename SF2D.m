classdef SF2D
    methods(Static)
        
        function [si] = SF_2D(x,y,NS,P,ss)
            %this function gives Shape function at any point (x,y) for node (xI,yI)
            % For linear basis we need 3x3 M matrix and for quadratic basis we need 6x6
            % M matrix in 2D

            P = transpose(P); % 2 Rows
            P_len = size(P,2);

            M=zeros(6,6); % basis = 2

            %Evaluation of the Moment matrix, here we take the node positions
            %as xxI and yyI for the summation process to be carried out easily

            for integer_1=1:P_len
                xxI = P(1,integer_1);
                yyI = P(2,integer_1);

                HB = [1;x-xxI; y-yyI; (x-xxI)^2; (y-yyI)^2; (x-xxI)*(y-yyI)]; 

                zz=(sqrt((x-xxI)^2+(y-yyI)^2))/ss;
                %phy = PHI is the weight function

                [phi] = SF2D.phi_eval(zz);
                M = M + (HB*transpose(HB)*phi);
            end

            %After we get the Moment matrix we construct the SF
            si = zeros(1,size(NS,1));
            for int_1 = 1:size(NS,1)
                xI = NS(int_1,1);
                yI = NS(int_1,2);

                Ho = [1 ;0;0;0;0;0];
                H = [1 ;x-xI; y-yI;(x-xI)^2; (y-yI)^2; (x-xI)*(y-yI)];

                z=(sqrt((x-xI)^2+(y-yI)^2))/ss;
                [PHI] = SF2D.phi_eval(z);

                si(int_1) = si(int_1) + transpose(Ho)*inv(M)*H*PHI;
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


        
    end %methods
end %class

