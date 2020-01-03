function [A] = part_of_NI(NC,NS,ss)

A = zeros(size(NC,1),size(NS,1));

int_row = 1; % Row counter for A
int_col = 1; % Column counter for A
for int_1 = 1:size(NC,1) % Looping over the no. of collocation pts. in the domain
    % arranging rows of A
    % We need the RK shape function at the collocation point x,y centered
    % at the Source point

    x = NC(int_1,1);
    y = NC(int_1,2);    
    [P] = required_nodes(x,y,NS,ss);
    [SIxx] = DSFxx.DSF_xx(x,y,NS,P,ss);
    [SIyy] = DSFyy.DSF_yy(x,y,NS,P,ss);
        
    for int_2 = 1:size(NS,1)  %Looping over the no. of source pts.        

        a_inter =SIxx(int_2)+SIyy(int_2);

        % Arranging in in the matrix
        A(int_row,int_col) = A(int_row,int_col) + a_inter;
        int_col = int_col + 1;
    end
    int_col=1;
    int_row = int_row+1;
end
