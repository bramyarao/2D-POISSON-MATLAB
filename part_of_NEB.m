function [A] = part_of_NEB(NEB,NS,ss,sq_alphag)

A = zeros(size(NEB,1),size(NS,1));

int_row = 1; % Row counter for A
int_col = 1; % Column counter for A

for int_1 = 1:size(NEB,1) % Looping over the no. of collocation pts. in the domain    

    x = NEB(int_1,1);
    y = NEB(int_1,2);    
    [P] = required_nodes(x,y,NS,ss);
    [SI] = SF2D.SF_2D(x,y,NS,P,ss);
        
    for int_2 = 1:size(NS,1)         
        a_inter = sq_alphag *SI(int_2);  % weight on the EB

        % Arranging in in the matrix
        A(int_row,int_col) = A(int_row,int_col) + a_inter;

        int_col = int_col + 1;
    end
    int_col=1;
    int_row = int_row+1;
end
