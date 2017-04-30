classdef cm
% ***********************************************************
% Marc Stadelmann - marc.stadelmann@istb.unibe.ch
% April 2017
%
% Continuum Mechanics class - University of Bern.
%
% History
% -------
% 1.3.17 - V3 - composition44 did not work when combining symbolic and
%               numeric input.
% 7.3.17   V4 - Added gradient02 and gradient22 functions.
% 12.3.17  V5 - Added tetrahedron + gradientTransform functions.
% 19.3.17  V6 - Added tetra face normal functions
% 21.3.17  V7 - cm.* was missing in get normal function.
% 26.3.17     - Added round decimals function
% 4.4.17   V8 - Bugfix + visualizeTensorAtP() 
%               and plot_tetra_point
% 11.4.17  V9   Bugfix / improved version of plot_vector()
% 27.4.17  V10  Bugfix: get_tetra_normal() returned the normals in a not
%               very smart order..
%
% ***********************************************************    

    methods (Static)
    % All functions are defined as static..
    % **********************************************************
%===============================================================      
    function s = cross_product(a,b)
    % **********************************************************
    % computes the cross product of two first order tensors.
    % index notation: si = Eijk aj bk
    % Where Eijk = Levi-Civita Permutation Tensor
    % lecture notes p. 42
    % **********************************************************
    if strcmp(class(a),'sym') == 1 || strcmp(class(b),'sym') == 1
        s = vpa(zeros(3,1)); % Variable-precision arithmetic:
    else
        s = zeros(3,1);
    end
    % define Levi-Civita Permutation Tensor
    levi = zeros(3,3,3);
    levi(1,2,3) = 1;
    levi(1,3,2) = -1;
    levi(2,3,1) = 1;
    levi(2,1,3) = -1;
    levi(3,1,2) = 1;
    levi(3,2,1) = -1;
    for i=1:3
        for j=1:3
            for k = 1:3
               s(i) = s(i) + levi(i,j,k)*a(j)*b(k);          
            end
        end 
    end
    end
%===============================================================       
    function S = dyadic_product11(a,b)
    % **********************************************************
    % Computes the dyadic product of two first order tensors.
    % Sij = ai*bj
    % lecture notes p. 45
    % **********************************************************
    if strcmp(class(a),'sym') == 1 || strcmp(class(b),'sym') == 1
        S = vpa(zeros(3,3));
    else
        S = zeros(3,3);
    end

    for i=1:3
        for j=1:3
            S(i,j) = a(i)*b(j); 
        end 
    end
    end
%===============================================================
    function K4 = dyadic_product22(A,B)
    % ***********************************************************
    % Computes the dyadic product of two second order tensors.
    % K4ijkl = Aij*Bkl
    % lecture notes p. 48
    % ***********************************************************
    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3
                    K4(i,j,k,l) = A(i,j)*B(k,l); 
                end
            end
        end 
    end
    end
%===============================================================
    function s = frobenius22(A,B)
    % ***********************************************************
    % Frobenius inner product of two second order tensors.
    % (Also second order scalar product or double contraction) 
    % s = A:B
    % s = Aij * Bij
    % lecture notes p. 46
    % Holzapfel p. 14
    % ***********************************************************
        if strcmp(class(A),'sym') == 1 || strcmp(class(B),'sym') == 1
            s = vpa(0); % Variable-precision arithmetic:
        else
            s = 0;
        end

        for i=1:3
            for j=1:3
                s = s + A(i,j)*B(i,j);
            end
        end
    end      
%===============================================================
    function s = frobenius44(A,B)
    % ***********************************************************
    % Frobenius inner product of two fourth order tensors.
    % s = A4 :: B4
    % s = Aijkl * Bijkl
    % lecture notes p. 49
    % ***********************************************************
        if strcmp(class(A),'sym') == 1 || strcmp(class(B),'sym') == 1
            s = vpa(0); % Variable-precision arithmetic:
        else
            s = 0;
        end

        for i=1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        s = s + A(i,j,k,l) * B(i,j,k,l);
                    end
                end
            end
        end
    end
%===============================================================
    function C = composition22(A,B)
    % ***********************************************************
    % Composition of two second order tensors.
    % Cij = Aik Bkj
    % -> matlab function: A * B
    % lecture notes p. 47
    % ***********************************************************
    if strcmp(class(A),'sym') == 1 || strcmp(class(B),'sym') == 1
        C = vpa(zeros(3,3));
    else
        C = zeros(3,3);
    end

    for i=1:3
        for j=1:3
            for k = 1:3
            	C(i,j) = C(i,j) + A(i,k)*B(k,j);    
            end
        end 
    end
    end   
%===============================================================
    function C4 = composition44(A4,B4)
    % ***********************************************************
    % Composition of two fourth order tensors.
    % Cijkl = Aijmn Bmnkl
    % ***********************************************************
    if strcmp(class(A4),'sym') == 1 || strcmp(class(B4),'sym') == 1
        C4 = vpa(zeros(3,3,3,3));
    else
        C4 = zeros(3,3,3,3);
    end

    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3
                    for m=1:3
                        for n=1:3
                            C4(i,j,k,l) = C4(i,j,k,l) + A4(i,j,m,n)*B4(m,n,k,l); 
                        end
                    end   
                end
            end
        end 
    end
    end    
%===============================================================
    function C = invert(A)
    % ***********************************************************
    % Inverts a second order tensor (3*3 matrix only!).
    % A is only invertible if det(A) != 0 (non-singular matrix)
    % ***********************************************************
    if length(A) == 3
        determinant = cm.det2(A);
        if determinant == 0
            disp('Singular Matrix - no Inverse!')
            C = 'NaN';
        else
            v11 =  A(2,2)*A(3,3) - A(2,3) * A(3,2);
            v12 = -A(1,2)*A(3,3) + A(1,3) * A(3,2);
            v13 =  A(1,2)*A(2,3) - A(1,3) * A(2,2);
            v21 = -A(2,1)*A(3,3) + A(2,3) * A(3,1);
            v22 =  A(1,1)*A(3,3) - A(1,3) * A(3,1);
            v23 = -A(1,1)*A(2,3) + A(1,3) * A(2,1);
            v31 =  A(2,1)*A(3,2) - A(2,2) * A(3,1);
            v32 = -A(1,1)*A(3,2) + A(1,2) * A(3,1);
            v33 =  A(1,1)*A(2,2) - A(1,2) * A(2,1);    
            C = (1/determinant)  * [v11 v12 v13;v21 v22 v23;v31 v32 v33];
        end
    else
        disp('Wrong dimension!')
        C = 'NaN';
    end
    end    
%===============================================================    
    function K4 = symmetric_product22(A,B)
    % ***********************************************************
    % Computes the symmetric poruct of two second order tensors.
    % K4ijkl = 0.5 * (Aik*Bjl + Ail*Bjk).
    % lecture notes p. 48
    % ***********************************************************
    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3
                    K4(i,j,k,l) = (1/2)*(A(i,k)*B(j,l)+A(i,l)*B(j,k)); 
                end
            end
        end 
    end
    end
%===============================================================
    function K4 = tensor_product22(A,B) 
    % ***********************************************************
    % Computes the tensor product of two second order tensors.
    % K4ijkl = Aik*Bjl
    % lecture notes p. 48
    % ***********************************************************
    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3
                    K4(i,j,k,l) = A(i,k)*B(j,l); 
                end
            end
        end 
    end
    end
%===============================================================
    function K4 = tensor_product22_2D(A,B)
    % ***********************************************************
    % Computes the tensor product of two second order tensors in 2D.
    % K4ijkl = Aik*Bjl
    % ***********************************************************
    for i=1:2
        for j=1:2
            for k = 1:2
                for l=1:2
                    K4(i,j,k,l) = A(i,k)*B(j,l); 
                end
            end
        end 
    end
    end
%===============================================================
    function K4 = transposed_product22(A,B)
    % ***********************************************************
    % Computes the transposed product of two second order tensors.
    % K4ijkl = Ail*Bjk
    % lecture notes p. 48
    % ***********************************************************    
    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3             
                    K4(i,j,k,l) = A(i,l)*B(j,k);                 
                end
            end
        end 
    end
    end  
%===============================================================
    function rad = deg2rad(a)
    % ***********************************************************
    % Converts degrees to radians.
    % ***********************************************************
    rad = a/360*2*pi;
    end
%===============================================================
    function deg = rad2deg(a)
    % ***********************************************************
    % Converts radians to degrees.
    % ***********************************************************
    deg = a/(2*pi)*360;
    end
%===============================================================
    function show1(a1)
    % ***********************************************************
    % Shows a first order tensor in a pretty way.
    % ***********************************************************
    if strcmp(class(a1),'sym') == 1
            outstr ={
            '/'  char(a1(1)) '\'
            '|'  char(a1(2)) '|'
            '\'  char(a1(3)) '/'
             };
    else
            outstr ={
            '/'  num2str(a1(1)) '\'
            '|'  num2str(a1(2)) '|' 
            '\'  num2str(a1(3)) '/'
             };
    end
    
    outstr1 = strcat(outstr,{' '});  %// add whitespace
    %// Convert to char array
    outstr_char = char(outstr1{:});
    %// Get size parameters
    [m,n] = size(outstr1);
    p = size(outstr_char,2);
    %// Reshape + Permute Magic to create a 
    %// char array "replica" of input cell array
    out = reshape(permute(reshape(outstr_char.',p,m,[]),[1 3 2]),n*p,m).';
    %// Display the char array
    disp(out)
        
    end
%===============================================================
    function show2(A33)
    % ***********************************************************
    % Shows a second order tensor in a pretty way.
    % ***********************************************************
    if strcmp(class(A33),'sym') == 1
            outstr ={
            '/'  char(A33(1,1))  char(A33(1,2))  char(A33(1,3))  '\'
            '|'  char(A33(2,1))  char(A33(2,2))  char(A33(2,3))  '|'
            '\'  char(A33(3,1))  char(A33(3,2))  char(A33(3,3))  '/'
             };
    else
            outstr ={
            '/'  num2str(A33(1,1))  num2str(A33(1,2))  num2str(A33(1,3))  '\'
            '|'  num2str(A33(2,1))  num2str(A33(2,2))  num2str(A33(2,3))  '|' 
            '\'  num2str(A33(3,1))  num2str(A33(3,2))  num2str(A33(3,3))  '/'
             };
    end
    
    outstr1 = strcat(outstr,{' '});  %// add whitespace
    %// Convert to char array
    outstr_char = char(outstr1{:});
    %// Get size parameters
    [m,n] = size(outstr1);
    p = size(outstr_char,2);
    %// Reshape + Permute Magic to create a 
    %// char array "replica" of input cell array
    out = reshape(permute(reshape(outstr_char.',p,m,[]),[1 3 2]),n*p,m).';
    %// Display the char array
    disp(out)
        
    end
%===============================================================
    function show3(A333)
    % ***********************************************************
    % Shows a thrid order tensor in a pretty way.
    % ***********************************************************
    if strcmp(class(A333),'sym') == 1
            outstr ={
            '/' '/'  char(A333(1,1,1))  char(A333(1,1,2))  char(A333(1,1,3))  '\' '\'
            '|' '|'  char(A333(1,2,1))  char(A333(1,2,2))  char(A333(1,2,3))  '|' '|'    
            '|' '\'  char(A333(1,3,1))  char(A333(1,3,2))  char(A333(1,3,3))  '/' '|' 
            '|' ' '  ''  ''  ''  ' ' '|'
            '|' '/'  char(A333(2,1,1))  char(A333(2,1,2))  char(A333(2,1,3))  '\' '|'
            '|' '|'  char(A333(2,2,1))  char(A333(2,2,2))  char(A333(2,2,3))  '|' '|'    
            '|' '\'  char(A333(2,3,1))  char(A333(2,3,2))  char(A333(2,3,3))  '/' '|' 
            '|' ' '  ''  ''  ''  ' ' '|'
            '|' '/'  char(A333(3,1,1))  char(A333(3,1,2))  char(A333(3,1,3))  '\' '|'
            '|' '|'  char(A333(3,2,1))  char(A333(3,2,2))  char(A333(3,2,3))  '|' '|'    
            '\' '\'  char(A333(3,3,1))  char(A333(3,3,2))  char(A333(3,3,3))  '/' '/'  
             };
    else
            outstr ={
            '/' '/'  num2str(A333(1,1,1))  num2str(A333(1,1,2))  num2str(A333(1,1,3))  '\' '\'
            '|' '|'  num2str(A333(1,2,1))  num2str(A333(1,2,2))  num2str(A333(1,2,3))  '|' '|'    
            '|' '\'  num2str(A333(1,3,1))  num2str(A333(1,3,2))  num2str(A333(1,3,3))  '/' '|' 
            '|' ' '  ''  ''  ''  ' ' '|'
            '|' '/'  num2str(A333(2,1,1))  num2str(A333(2,1,2))  num2str(A333(2,1,3))  '\' '|'
            '|' '|'  num2str(A333(2,2,1))  num2str(A333(2,2,2))  num2str(A333(2,2,3))  '|' '|'    
            '|' '\'  num2str(A333(2,3,1))  num2str(A333(2,3,2))  num2str(A333(2,3,3))  '/' '|' 
            '|' ' '  ''  ''  ''  ' ' '|'
            '|' '/'  num2str(A333(3,1,1))  num2str(A333(3,1,2))  num2str(A333(3,1,3))  '\' '|'
            '|' '|'  num2str(A333(3,2,1))  num2str(A333(3,2,2))  num2str(A333(3,2,3))  '|' '|'    
            '\' '\'  num2str(A333(3,3,1))  num2str(A333(3,3,2))  num2str(A333(3,3,3))  '/' '/'  
             };
    end
    
    outstr1 = strcat(outstr,{' '});  %// add whitespace
    %// Convert to char array
    outstr_char = char(outstr1{:});
    %// Get size parameters
    [m,n] = size(outstr1);
    p = size(outstr_char,2);
    %// Reshape + Permute Magic to create a 
    %// char array "replica" of input cell array
    out = reshape(permute(reshape(outstr_char.',p,m,[]),[1 3 2]),n*p,m).';
    %// Display the char array
    disp(out)
        
    end
%===============================================================
    function show4(A3333)
    % ***********************************************************
    % Shows a fourth order tensor in a pretty way.
    % ***********************************************************
    if strcmp(class(A3333),'sym') == 1
        outstr ={
        '/' '/'  char(A3333(1,1,1,1))  char(A3333(1,1,1,2))  char(A3333(1,1,1,3))  '\' '/'  char(A3333(1,2,1,1))  char(A3333(1,2,1,2))  char(A3333(1,2,1,3))  '\' '/'  char(A3333(1,3,1,1))  char(A3333(1,3,1,2))  char(A3333(1,3,1,3))  '\' '\'
        '|' '|'  char(A3333(1,1,2,1))  char(A3333(1,1,2,2))  char(A3333(1,1,2,3))  '|' '|'  char(A3333(1,2,2,1))  char(A3333(1,2,2,2))  char(A3333(1,2,2,3))  '|' '|'  char(A3333(1,3,2,1))  char(A3333(1,3,2,2))  char(A3333(1,3,2,3))  '|' '|'    
        '|' '\'  char(A3333(1,1,3,1))  char(A3333(1,1,3,2))  char(A3333(1,1,3,3))  '/' '\'  char(A3333(1,2,3,1))  char(A3333(1,2,3,2))  char(A3333(1,2,3,3))  '/' '\'  char(A3333(1,3,3,1))  char(A3333(1,3,3,2))  char(A3333(1,3,3,3))  '/' '|' 
        '|' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '|'
        '|' '/'  char(A3333(2,1,1,1))  char(A3333(2,1,1,2))  char(A3333(2,1,1,3))  '\' '/'  char(A3333(2,2,1,1))  char(A3333(2,2,1,2))  char(A3333(2,2,1,3))  '\' '/'  char(A3333(2,3,1,1))  char(A3333(2,3,1,2))  char(A3333(2,3,1,3))  '\' '|'
        '|' '|'  char(A3333(2,1,2,1))  char(A3333(2,1,2,2))  char(A3333(2,1,2,3))  '|' '|'  char(A3333(2,2,2,1))  char(A3333(2,2,2,2))  char(A3333(2,2,2,3))  '|' '|'  char(A3333(2,3,2,1))  char(A3333(2,3,2,2))  char(A3333(2,3,2,3))  '|' '|'    
        '|' '\'  char(A3333(2,1,3,1))  char(A3333(2,1,3,2))  char(A3333(2,1,3,3))  '/' '\'  char(A3333(2,2,3,1))  char(A3333(2,2,3,2))  char(A3333(2,2,3,3))  '/' '\'  char(A3333(2,3,3,1))  char(A3333(2,3,3,2))  char(A3333(2,3,3,3))  '/' '|' 
        '|' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '|'
        '|' '/'  char(A3333(3,1,1,1))  char(A3333(3,1,1,2))  char(A3333(3,1,1,3))  '\' '/'  char(A3333(3,2,1,1))  char(A3333(3,2,1,2))  char(A3333(3,2,1,3))  '\' '/'  char(A3333(3,3,1,1))  char(A3333(3,3,1,2))  char(A3333(3,3,1,3))  '\' '|'
        '|' '|'  char(A3333(3,1,2,1))  char(A3333(3,1,2,2))  char(A3333(3,1,2,3))  '|' '|'  char(A3333(3,2,2,1))  char(A3333(3,2,2,2))  char(A3333(3,2,2,3))  '|' '|'  char(A3333(3,3,2,1))  char(A3333(3,3,2,2))  char(A3333(3,3,2,3))  '|' '|'    
        '\' '\'  char(A3333(3,1,3,1))  char(A3333(3,1,3,2))  char(A3333(3,1,3,3))  '/' '\'  char(A3333(3,2,3,1))  char(A3333(3,2,3,2))  char(A3333(3,2,3,3))  '/' '\'  char(A3333(3,3,3,1))  char(A3333(3,3,3,2))  char(A3333(3,3,3,3))  '/' '/'  
         };
    else
        outstr ={
        '/' '/'  num2str(A3333(1,1,1,1))  num2str(A3333(1,1,1,2))  num2str(A3333(1,1,1,3))  '\' '/'  num2str(A3333(1,2,1,1))  num2str(A3333(1,2,1,2))  num2str(A3333(1,2,1,3))  '\' '/'  num2str(A3333(1,3,1,1))  num2str(A3333(1,3,1,2))  num2str(A3333(1,3,1,3))  '\' '\'
        '|' '|'  num2str(A3333(1,1,2,1))  num2str(A3333(1,1,2,2))  num2str(A3333(1,1,2,3))  '|' '|'  num2str(A3333(1,2,2,1))  num2str(A3333(1,2,2,2))  num2str(A3333(1,2,2,3))  '|' '|'  num2str(A3333(1,3,2,1))  num2str(A3333(1,3,2,2))  num2str(A3333(1,3,2,3))  '|' '|'    
        '|' '\'  num2str(A3333(1,1,3,1))  num2str(A3333(1,1,3,2))  num2str(A3333(1,1,3,3))  '/' '\'  num2str(A3333(1,2,3,1))  num2str(A3333(1,2,3,2))  num2str(A3333(1,2,3,3))  '/' '\'  num2str(A3333(1,3,3,1))  num2str(A3333(1,3,3,2))  num2str(A3333(1,3,3,3))  '/' '|' 
        '|' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '|'
        '|' '/'  num2str(A3333(2,1,1,1))  num2str(A3333(2,1,1,2))  num2str(A3333(2,1,1,3))  '\' '/'  num2str(A3333(2,2,1,1))  num2str(A3333(2,2,1,2))  num2str(A3333(2,2,1,3))  '\' '/'  num2str(A3333(2,3,1,1))  num2str(A3333(2,3,1,2))  num2str(A3333(2,3,1,3))  '\' '|'
        '|' '|'  num2str(A3333(2,1,2,1))  num2str(A3333(2,1,2,2))  num2str(A3333(2,1,2,3))  '|' '|'  num2str(A3333(2,2,2,1))  num2str(A3333(2,2,2,2))  num2str(A3333(2,2,2,3))  '|' '|'  num2str(A3333(2,3,2,1))  num2str(A3333(2,3,2,2))  num2str(A3333(2,3,2,3))  '|' '|'    
        '|' '\'  num2str(A3333(2,1,3,1))  num2str(A3333(2,1,3,2))  num2str(A3333(2,1,3,3))  '/' '\'  num2str(A3333(2,2,3,1))  num2str(A3333(2,2,3,2))  num2str(A3333(2,2,3,3))  '/' '\'  num2str(A3333(2,3,3,1))  num2str(A3333(2,3,3,2))  num2str(A3333(2,3,3,3))  '/' '|' 
        '|' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '' '|'
        '|' '/'  num2str(A3333(3,1,1,1))  num2str(A3333(3,1,1,2))  num2str(A3333(3,1,1,3))  '\' '/'  num2str(A3333(3,2,1,1))  num2str(A3333(3,2,1,2))  num2str(A3333(3,2,1,3))  '\' '/'  num2str(A3333(3,3,1,1))  num2str(A3333(3,3,1,2))  num2str(A3333(3,3,1,3))  '\' '|'
        '|' '|'  num2str(A3333(3,1,2,1))  num2str(A3333(3,1,2,2))  num2str(A3333(3,1,2,3))  '|' '|'  num2str(A3333(3,2,2,1))  num2str(A3333(3,2,2,2))  num2str(A3333(3,2,2,3))  '|' '|'  num2str(A3333(3,3,2,1))  num2str(A3333(3,3,2,2))  num2str(A3333(3,3,2,3))  '|' '|'    
        '\' '\'  num2str(A3333(3,1,3,1))  num2str(A3333(3,1,3,2))  num2str(A3333(3,1,3,3))  '/' '\'  num2str(A3333(3,2,3,1))  num2str(A3333(3,2,3,2))  num2str(A3333(3,2,3,3))  '/' '\'  num2str(A3333(3,3,3,1))  num2str(A3333(3,3,3,2))  num2str(A3333(3,3,3,3))  '/' '/'  
         };
    end

    outstr1 = strcat(outstr,{' '});  %// add whitespace
    %// Convert to char array
    outstr_char = char(outstr1{:});
    %// Get size parameters
    [m,n] = size(outstr1);
    p = size(outstr_char,2);
    %// Reshape + Permute Magic to create a 
    %// char array "replica" of input cell array
    out = reshape(permute(reshape(outstr_char.',p,m,[]),[1 3 2]),n*p,m).';
    %// Display the char array
    disp(out)
        
    end
%===============================================================
    function det = det2(A)
    % ***********************************************************
    % Computes the determinant of a second order tensor.
    % det = levi(i,j,k)*A(1,i)*A(2,j)*A(3,k)    
    % (Matlab command: det(A)
    % ***********************************************************
    det = 0;

    levi = zeros(3,3,3);
    levi(1,2,3) = 1;
    levi(1,3,2) = -1;
    levi(2,3,1) = 1;
    levi(2,1,3) = -1;
    levi(3,1,2) = 1;
    levi(3,2,1) = -1;

    for i=1:3
        for j = 1:3
            for k =1:3            
                det = det + levi(i,j,k)*A(1,i)*A(2,j)*A(3,k);
            end
        end
    end
    end
%===============================================================
    function grad = gradient02(s,A)
    % ***********************************************************
    % Gradient of a scalar versus a second order tensor.
    % The result is a second order tensor.
    % e.g. nabla_A(trace(A)) = I is cm.gradient02(cm.trace(A),A).
    % ***********************************************************
    grad=vpa(zeros(3,3));

    for i=1:3
        for j=1:3
            grad(i,j) = diff(s,A(i,j));
        end 
    end
    end
%===============================================================
    function grad_mat = gradient22(A,B)
    % ***********************************************************
    % Gradient of a 2nd order tensor vs a second order tensor.
    % The result is a fourth order tensor.
    % e.g. nabla_A(A^2) is cm.gradient22(A^2,A).
    % ***********************************************************
    grad_mat=vpa(zeros(3,3,3,3));
    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3
                    grad_mat(i,j,k,l) = diff(A(i,j),B(k,l));
                end
            end
        end 
    end
    end
%===============================================================
    function s = mixed_product(a,b,c)
    % ***********************************************************
    % Computes the mixed product of three first order tensors.
    % s = Eijk aj bk ci
    % Eijk = Levi-Civita Permutation Tensor
    % ***********************************************************
        if strcmp(class(a),'sym') == 1 || strcmp(class(b),'sym') == 1 || strcmp(class(c),'sym') == 1
            % symbolic or numeric input?
            crossp = vpa(zeros(3,1)); % Variable-precision arithmetic:
            s = vpa(0);
        else
            crossp = zeros(3,1);
            s = 0;
        end
        % define Levi-Civita Permutation Tensor
        levi = zeros(3,3,3);
        levi(1,2,3) = 1;
        levi(1,3,2) = -1;
        levi(2,3,1) = 1;
        levi(2,1,3) = -1;
        levi(3,1,2) = 1;
        levi(3,2,1) = -1;
        for i=1:3
            for j=1:3
                for k = 1:3
                   crossp(i) = crossp(i) + levi(i,j,k)*a(j)*b(k);          
                end
            end 
            s = s + crossp(i) * c(i);
        end
    end
%===============================================================
    function R = rot(adeg,d)
    % ***********************************************************
    % Creates rotation matrix.
    % R = rot(angle [deg], axis [x,y,z]).
    % (e.g. R = cm.rot(30,'z')).
    % ***********************************************************
    a = cm.deg2rad(adeg);
    switch d
        case 'x'
            % positive rotation around X axis
            R=[1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];        
        case 'y'
            % positive rotation around Y axis
            R=[cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];        
        case 'z'
            % positive rotation around Z axis
            R=[cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    end
    end
%===============================================================
    function B = transpose2(A)
    % ***********************************************************
    % Transposes a second order tensor
    % Bij = Aji
    % ***********************************************************
        for i=1:3
            for j=1:3
                B(i,j) = A(j,i); 
            end 
        end
    end
%===============================================================
    function B4 = transpose4(A4)
    % ***********************************************************
    % Transposes a fourth order tensor
    % Bijkl = Aklij
    % ***********************************************************
        for i=1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                        B4(i,j,k,l) = A4(k,l,i,j); 
                    end
                end
            end 
        end
    end    
%===============================================================
    function bool = check_symmetry(A)
    % ***********************************************************
    % Returns 1 if 3*3 input matrix is symmetric, otherwise 0.
    % ***********************************************************
    if A(2,1) == A(1,2) && A(3,1) == A(1,3) && A(2,3) == A(3,2)
        bool = 1;
    else
        bool = 0;
    end
    end
%===============================================================
    function symm = check_minor_symmetry(A4)
    % ***********************************************************
    % Returns 1 if 4th order tensor has minor symmetry.
    % Minor symmetry = Sijkl = Sjikl, Sijkl = Sijlk
    % ***********************************************************
    symm = 1;
    for i = 1:3
        for j = 1:3
            temp = A4(:,:,i,j);
            
            if temp(2,1) == temp(1,2) && temp(3,1) == temp(1,3) && temp(2,3) == temp(3,2)
                symm = 1;
            else
                symm = 0;
                break;
            end
        end
    end
    if symm == 1
        for i = 1:3
            for j = 1:3
                
                temp = squeeze(A4(i,j,:,:));

                if temp(2,1) == temp(1,2) && temp(3,1) == temp(1,3) && temp(2,3) == temp(3,2)
                    symm = 1;
                else
                    symm = 0;
                    break;
                end
            end
        end
    end
    end    
%===============================================================
    function b = isomorph_33_91(A)
    % ***********************************************************
    % Converts 3*3 matrix into 9*1 vector.
    % If A is symmetric, the outputvector will be 6*1.
    % lecture notes p. 33 / 34
    % ***********************************************************
        if cm.check_symmetry(A) == 0 % check if input matrix is symmetric
            if strcmp(class(A),'sym') == 1
                b = vpa(zeros(9,1));
            else
                b = zeros(9,1);
            end
            b = reshape(A.',[9,1]);
        else % symmetric matrix
            sprintf('\n Tensor is symmetric!')
            if strcmp(class(A),'sym') == 1
                b = vpa(zeros(6,1));
            else
                b = zeros(6,1);
            end
            b(1) = A(1,1);
            b(2) = A(2,2);
            b(3) = A(3,3);
            b(4) = sqrt(2) * A(2,3);
            b(5) = sqrt(2) * A(3,1);
            b(6) = sqrt(2) * A(1,2);
        end
    end
%===============================================================
    function B = isomorph_3333_99(A4)
    % ***********************************************************
    % Isomorphism between 4th order tensor (3*3*3*3 matrix)
    % and 9*9 matrix.
    % If the input tensor posses minor symmetries, the resulting
    % matrix is further reduced to a 6*6 matrix.
    % lecture notes p. 35 / 36
    % ***********************************************************
        if cm.check_minor_symmetry(A4) == 0 
            if strcmp(class(A4),'sym') == 1
                B = vpa(zeros(9,9));
            else
                B = zeros(9,9);
            end

            for i=1:3
                for j=1:3
                    ii = (i-1)*3 + j;
                    B(:,ii) = reshape(A4(:,:,i,j).',[9,1]);
                end
            end
        else % symmetric matrix
            
            sprintf('\n Tensor is symmetric!')

            if strcmp(class(A4),'sym') == 1
                B = vpa(zeros(6,6));
            else
                B = zeros(6,6);
            end

            B(1,1) = A4(1,1,1,1);
            B(2,1) = A4(2,2,1,1);
            B(3,1) = A4(3,3,1,1);
            B(4,1) = sqrt(2) * A4(2,3,1,1);
            B(5,1) = sqrt(2) * A4(3,1,1,1);
            B(6,1) = sqrt(2) * A4(1,2,1,1);

            B(1,2) = A4(1,1,2,2);
            B(2,2) = A4(2,2,2,2);
            B(3,2) = A4(3,3,2,2);
            B(4,2) = sqrt(2) * A4(2,3,2,2);
            B(5,2) = sqrt(2) * A4(3,1,2,2);
            B(6,2) = sqrt(2) * A4(1,2,2,2);

            B(1,3) = A4(1,1,3,3);
            B(2,3) = A4(2,2,3,3);
            B(3,3) = A4(3,3,3,3);
            B(4,3) = sqrt(2) * A4(2,3,3,3);
            B(5,3) = sqrt(2) * A4(3,1,3,3);
            B(6,3) = sqrt(2) * A4(1,2,3,3);

            B(1,4) = sqrt(2) * A4(1,1,2,3);
            B(2,4) = sqrt(2) * A4(2,2,2,3);
            B(3,4) = sqrt(2) * A4(3,3,2,3);
            B(4,4) = 2 * A4(2,3,2,3);
            B(5,4) = 2 * A4(3,1,2,3);
            B(6,4) = 2 * A4(1,2,2,3);

            B(1,5) = sqrt(2) * A4(1,1,3,1);
            B(2,5) = sqrt(2) * A4(2,2,3,1);
            B(3,5) = sqrt(2) * A4(3,3,3,1);
            B(4,5) = 2 * A4(2,3,3,1);
            B(5,5) = 2 * A4(3,1,3,1);
            B(6,5) = 2 * A4(1,2,3,1);

            B(1,6) = sqrt(2) * A4(1,1,1,2);
            B(2,6) = sqrt(2) * A4(2,2,1,2);
            B(3,6) = sqrt(2) * A4(3,3,1,2);
            B(4,6) = 2 * A4(2,3,1,2);
            B(5,6) = 2 * A4(3,1,1,2);
            B(6,6) = 2 * A4(1,2,1,2);

        end
    end
%===============================================================
    function visualizeTensor(D)
    % ***********************************************************
    % Visualizes a tensor using a sphere / ellipsoid.
    % ***********************************************************
    currentD = findall(0,'Type','figure','Tag','D');
    if ~isempty(currentD)
        p = get(currentD, 'Position');
        figure(currentD);
        [az,el] = view;
        delete(currentD);
        set(0, 'DefaultFigurePosition', p);
        figure;
    else
        figure;
    end

    sz=size(D);
    if length(sz)==2
        nx=1;ny=1;
    elseif length(sz)==3
        nx=sz(3);ny=1;
    elseif length(sz)==4
        nx=sz(3);ny=sz(4);
    end

    n=size(D,3);
    for i=1:nx
        for j=1:ny
            [v,l]=eig(D(:,:,i,j));
            [X,Y,Z]=ellipsoid(0,0,0,l(1,1),l(2,2),l(3,3),20);
            sz=size(X);
            for x=1:sz(1)
                for y=1:sz(2)
                    A=[X(x,y) Y(x,y) Z(x,y)]';
                    A=v*A;
                    X(x,y)=A(1);Y(x,y)=A(2);Z(x,y)=A(3);
                end
            end
            X=X+(i-1)*2;
            Y=Y+(j-1)*2;
            h(i)=surf(X,Y,Z);
            if i==1 & j==1
                hold;
            end
        end
    end

    if exist('az') > 0
        view(az,el); %
    end

    grid on;
    xlabel('X') % x-axis label
    ylabel('Y') % y-axis label
    zlabel('Z') % z-axis label
    axis equal
    colormap([0.2 0.2 0.7])
    set(gcf,'Tag','D');
    end
%===============================================================
    function visualizeF(F)
    % ***********************************************************
    % Visualizes a tensor using cube.
    % This is especially usefull to visualize a deformation 
    % gradient F.
    % ***********************************************************
    currentF = findall(0,'Type','figure','Tag','F');

    if ~isempty(currentF)
        p = get(currentF, 'Position');
        figure(currentF);
        [az,el] = view;
        delete(currentF);
        set(0, 'DefaultFigurePosition', p);
        figure
    else
        figure
    end

    plot3(0,0,0,'kx','MarkerSize',10); hold on;

    if exist('az') > 0
        view(az,el); %
    end

    pC=[0; 0; 0];
    pX=[1; 0; 0];
    pY=[0; 1; 0];
    pZ=[0; 0; 1];
    pYZ=[0; 1; 1];
    pXZ=[1; 0; 1];
    pXY=[1; 1; 0];
    pXYZ=[1; 1; 1];

    quiver3(pC(1),pC(2),pC(3),2*pX(1),2*pX(2),2*pX(3),'r')
    quiver3(pC(1),pC(2),pC(3),2*pY(1),2*pY(2),2*pY(3),'b')
    quiver3(pC(1),pC(2),pC(3),2*pZ(1),2*pZ(2),2*pZ(3),'g')
    axis equal;

    cm.plotCube(pC,pX,pY,pZ,pYZ,pXZ,pXY,pXYZ,'k')
    min = -2;
    max = 2;

    pC=F*pC;
    pX=F*pX;
    pY=F*pY;
    pZ=F*pZ;
    pYZ=F*pYZ;
    pXZ=F*pXZ;
    pXY=F*pXY;
    pXYZ=F*pXYZ;

    quiver3(pC(1),pC(2),pC(3),2*pX(1),2*pX(2),2*pX(3),'r')
    quiver3(pC(1),pC(2),pC(3),2*pY(1),2*pY(2),2*pY(3),'b')
    quiver3(pC(1),pC(2),pC(3),2*pZ(1),2*pZ(2),2*pZ(3),'g')
    axis equal;

    cm.plotCube(pC,pX,pY,pZ,pYZ,pXZ,pXY,pXYZ,'r')
    min = -2;
    max = 2;
    grid on;
    xlabel('X') % x-axis label
    ylabel('Y') % y-axis label
    zlabel('Z') % z-axis label
    axis([min, max, min, max, min, max]) 
    set(gcf,'Tag','F')
    end
%===============================================================
    function visualizeF_2D(F)
        % ***********************************************************
        % Visualizes a tensor using square.
        % This is especially usefull to visualize a 2D deformation 
        % gradient F.
        % ***********************************************************
    currentF = findall(0,'Type','figure','Tag','F2d');

    if ~isempty(currentF)
        p = get(currentF, 'Position');
        figure(currentF);
        [az,el] = view;
        delete(currentF);
        set(0, 'DefaultFigurePosition', p);
        figure
    else
        figure
    end

    if exist('az') > 0
        view(az,el); 
    end

    pC=[0; 0];
    pX=[1; 0];
    pY=[0; 1];
    pXY=[1; 1];

    quiver(pC(1),pC(2),2*pX(1),2*pX(2),'r'); hold on;
    quiver(pC(1),pC(2),2*pY(1),2*pY(2),'b'); hold on;

    axis equal;
    line([pC(1),pX(1),pXY(1),pY(1),pC(1)],[pC(2),pX(2),pXY(2),pY(2),pC(2)],'color','k')
    
    pC=F*pC;
    pX=F*pX;
    pY=F*pY;
    pXY=F*pXY;

    line([pC(1),pX(1),pXY(1),pY(1),pC(1)],[pC(2),pX(2),pXY(2),pY(2),pC(2)],'color','r')

    min = -3;
    max = 3;
    grid on;
    xlabel('X') % x-axis label
    ylabel('Y') % y-axis label

    axis([min, max, min, max, min, max]) 
    set(gcf,'Tag','F2d')
    end
%===============================================================
    function tr = trace(A)
    % ***********************************************************
    % computes the trace of a second order tensor.
    % trace(A) = Aij * dij = Aii
    % (dij = kroenecker delta, Holzapfel p.5)
    % ***********************************************************
    tr = 0;
    for i=1:3
        tr = tr + A(i,i);
    end
    end
%===============================================================
    function u = transform_21(Q,v)
    % ***********************************************************
    % Transforms a first order tensor by the application of a 
    % second order tensor.
    % ui = Qij vj
    % lecture notes p. 44
    % ***********************************************************
    if strcmp(class(Q),'sym') == 1 || strcmp(class(v),'sym') == 1
        u = vpa(zeros(3,1));
    else
        u = zeros(3,1);
    end

    for i=1:3
        for j=1:3
            u(i)=u(i)+Q(i,j)*v(j);
        end 
    end
    end
%===============================================================
    function s = transform_32(A3,B)
    % ***********************************************************
    % Computes the multiplication of a second order tensor
    % by a third order tensor.
    % The third order tensor has to be in a 3*3*3 format!
    % si = Aijk Bjk
    % ***********************************************************
    
    % symbolic or numeric? -> used for zero initialisation
    if strcmp(class(A3),'sym') == 1 || strcmp(class(B),'sym') == 1
        s = vpa(zeros(3,1));
    else
        s = zeros(3,1);
    end
    for i=1:3
        for j=1:3
            for k = 1:3
               s(i) = s(i) + A3(i,j,k)*B(j,k);          
            end
        end 
    end
    end
%===============================================================
    function Y = transform_42(A4,A)
    % ***********************************************************
    % Transformation of a second order tensor by the application
    % of a fourth order tensor.
    % Yij = A4ijkl * Xkl
    % lecture notes p. 47
    % ***********************************************************
    if strcmp(class(A),'sym') == 1
        Y = vpa(zeros(3,3));
    else
        Y = zeros(3,3);
    end
    for i=1:3
        for j=1:3
            for k = 1:3
                for l=1:3
                    Y(i,j) = Y(i,j) + A4(i,j,k,l)*A(k,l);
                end
            end
        end 
    end
    end
%===============================================================
    function Y = transform_42_2D(A4,A)
    % ***********************************************************
    % Transformation of a second order tensor by the application
    % of a fourth order tensor - same but in 2D!
    % Yij = A4ijkl * Xkl
    % ***********************************************************
    if strcmp(class(A4),'sym') == 1 && strcmp(class(A),'sym') == 1
        Y = vpa(zeros(2,2));
    else
        Y = zeros(2,2);
    end

    for i=1:2
        for j=1:2
            for k = 1:2
                for l=1:2
                    Y(i,j) = Y(i,j) + A4(i,j,k,l)*A(k,l);
                end
            end
        end 
    end
    end
%===============================================================
    function nrm = norm(a)
    % ***********************************************************
    % Computes the norm of a tensor.
    % norm = sqrt(ai*ai)
    % lecture notes p. 41, 49
    % ***********************************************************
    nrm = 0;
    if size(a,2) == 1
        % first order    
        for i=1:3
        nrm = nrm + a(i)*a(i);
        end
    else
        % second order
        for i=1:3
            for j = 1:3
                nrm = nrm + a(i,j)*a(i,j);
            end
        end
    end
    nrm = sqrt(nrm);
    end
%===============================================================
    function s = first_inv(A)
    % ***********************************************************
    % computes the first invariant of a second order tensor.
    % first invariant = trace!
    % ***********************************************************
    s = cm.trace(A);
    end
%===============================================================
    function s = sec_inv(A)
    % ***********************************************************
    % computes the second invariant of a second order tensor.
    % sec(A) = 0.5 (Aii Ajj - Aji Aij)
    % sec(A) = 0.5 ( tr(A)^2 - tr(A^2) )
    % ***********************************************************
    s = 0;
    for i=1:3
        for j = 1:3
            s = s + (1/2)*(A(i,i)*A(j,j)-A(j,i)*A(i,j));
        end
    end
    end
%===============================================================
    function s = third_inv(A)
    % ***********************************************************
    % computes the third invariant of a second order tensor.
    % first invariant = trace!
    % ***********************************************************
    s = cm.det2(A);
    end
%===============================================================
    function plot_vector(p0,p1,lineWidth,color)
    % ***********************************************************
    %   plot_vector(p0,p1, line thickness, color) plots a line 
    %   vector with arrow pointing from point p0 to point p1.
    %
    %   Example:
    %       p0 = [1 2 3].';   % Coordinate of the first point p0
    %       p1 = [4 5 6].';   % Coordinate of the second point p1
    %       cm.plot_vector(p0,p1,2,'m')
    % ***********************************************************
    p0 = double(p0);
    p1=double(p1);
      if max(size(p0))==3
          if max(size(p1))==3
              x0 = p0(1);
              y0 = p0(2);
              z0 = p0(3);
              x1 = p1(1);
              y1 = p1(2);
              z1 = p1(3);
              plot3([x0;x1],[y0;y1],[z0;z1],'LineWidth',lineWidth,'Color',color);   % Draw a line between p0 and p1

              p = p1-p0;
              alpha = 0.1;  % Size of arrow head relative to the length of the vector
              beta = 0.1;  % Width of the base of the arrow head relative to the length

              hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
              hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
              hw = [z1-alpha*p(3);z1;z1-alpha*p(3)];

              hold on
              plot3(hu(:),hv(:),hw(:),'LineWidth',lineWidth,'Color',color)  % Plot arrow head
              hold off
          else
              error('p0 and p1 must have the same dimension')
          end
      elseif max(size(p0))==2
          if max(size(p1))==2
              x0 = p0(1);
              y0 = p0(2);
              x1 = p1(1);
              y1 = p1(2);
              plot([x0;x1],[y0;y1],'LineWidth',lineWidth,'Color',color);   % Draw a line between p0 and p1

              p = p1-p0;
              alpha = 0.1;  % Size of arrow head relative to the length of the vector
              beta = 0.1;  % Width of the base of the arrow head relative to the length

              hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
              hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];

              hold on
              plot(hu(:),hv(:),'LineWidth',lineWidth,'Color',color)  % Plot arrow head
              hold off
          else
              error('p0 and p1 must have the same dimension')
          end
      else
          error('this function only accepts 2D or 3D vector')
      end
    end
%===============================================================
    function plotCube(pC,pX,pY,pZ,pYZ,pXZ,pXY,pXYZ,col)
    % ***********************************************************
    % Plots a cube using 8 input points.
    % ***********************************************************
    x=[pC(1),pX(1),pXY(1),pY(1),pC(1),pZ(1),pXZ(1),pX(1),pXZ(1),pXYZ(1),pXY(1),pXYZ(1),pYZ(1),pY(1),pYZ(1),pZ(1)];
    y=[pC(2),pX(2),pXY(2),pY(2),pC(2),pZ(2),pXZ(2),pX(2),pXZ(2),pXYZ(2),pXY(2),pXYZ(2),pYZ(2),pY(2),pYZ(2),pZ(2)];
    z=[pC(3),pX(3),pXY(3),pY(3),pC(3),pZ(3),pXZ(3),pX(3),pXZ(3),pXYZ(3),pXY(3),pXYZ(3),pYZ(3),pY(3),pYZ(3),pZ(3)];
    plot3(x,y,z,col)
    end
%===============================================================
    function plotCircle3D(center,normal,radius)
    % ***********************************************************
    % Plots a circle in 3D using center, normal and radius as input.
    % ***********************************************************
    theta=0:0.01:2*pi;
    v=null(normal);
    points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
    plot3(points(1,:),points(2,:),points(3,:),'r-');
    end
%===============================================================
    function s = scalar_product(a,b)
    % ***********************************************************
    % computes the scalar product of two first order tensors.
    % s = ai * bi
    % lecture notes p. 40
    % ***********************************************************
        % dimension of first order tensor
        dim = size(a,1);
        % symbolic or numeric input?
        if strcmp(class(a),'sym') == 1 || strcmp(class(b),'sym') == 1
            s = vpa(0); % Variable-precision arithmetic:
        else
            s = 0;
        end
        for i=1:dim
            s = s + a(i)*b(i);
        end
    end
%===============================================================
    function S = deviator(A)
    % ***********************************************************
    % computes the deviator
    % dev(A) = A - (1/3) * tensor_trace(A) * I
    % ***********************************************************
    I = eye(3);
    S = A - (1/3) * tensor_trace(A) * I;
    end
%===============================================================
    function [] = mohrplot(S)
    % ***********************************************************
    % Mohrs cricle
    % ***********************************************************

    Center = (S(1,1)+S(2,2))/2;
    R = sqrt(((S(1,1)-S(2,2))/2)^2+S(1,2)^2);
    sigmaMax = Center + R
    sigmaMin = Center - R

    % Plot Mohr's circle
    c_handle = circle([Center,0],R,100,'-k');
    hold on
    plot([S(1,1), S(2,2)],[-S(1,2), S(1,2)],'k-','LineWidth',2)
    plot(Center,0,'ko','MarkerFaceColor','k')   % Center
    plot(S(1,1),-S(1,2),'ko','MarkerFaceColor','k') %A
    plot(S(2,2),S(1,2),'ko','MarkerFaceColor','k') %B
    plot(Center,R,'ko','MarkerFaceColor','k')   % Tau_max
    plot(Center,-R,'ko','MarkerFaceColor','k')  % Tau_max (negative)
    plot(Center+R, 0, 'ko','MarkerFaceColor','k')   % sigma_1
    plot(Center-R, 0, 'ko','MarkerFaceColor','k')	% sigma_2

    % Plot options
    axis equal
    set(c_handle,'Color','k','LineWidth',2)
    ylim([-1.25*R, 1.25*R])
    xlim([Center-1.25*R, Center+1.25*R])
    plot([Center-1.25*R, Center+1.25*R],[0, 0],'k-') % line through 0 tau

    % Annotate graph
    %text(S(1,1)+R/12,-S(1,2),sprintf('A (%d, %d)',S(1,1),-S(1,2)),...
    %     'HorizontalAlignment','left','FontSize',18)
    %text(S(2,2)-R/12,S(1,2),sprintf('B (%d, %d)',S(2,2),S(1,2)),...
    %     'HorizontalAlignment','right','FontSize',18)
    text(Center,1.1*R,'\tau_{max}',...
         'HorizontalAlignment','center','FontSize',18)
    text(Center,-1.1*R,'\tau_{min}',...
         'HorizontalAlignment','center','FontSize',18)
    text(Center+1.1*R,R/12,'\sigma_{max}',...
         'HorizontalAlignment','center','FontSize',18)
    text(Center-1.1*R,R/12,'\sigma_{min}',...
         'HorizontalAlignment','center','FontSize',18)

    % More plot options
    set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
    set(get(gca,'XLabel'),'String','\sigma (MPa)','FontSize',24,'FontWeight','bold','FontName','Times')
    set(get(gca,'YLabel'),'String','\tau (MPa)','FontSize',24,'FontWeight','bold','FontName','Times')

    hold off

    end
%===============================================================
    function plot_tetra(node1, node2, node3, node4)
    % ***********************************************************
    % Draws a tetrahedron using the four nodes
    % (node1, node2, node3, node4).
    % ***********************************************************

    node1 = double(node1);
    node2 = double(node2);
    node3 = double(node3);
    node4 = double(node4);
    
    X = [node1(1) node2(1) node3(1) node4(1)]';
    Y = [node1(2) node2(2) node3(2) node4(2)]';
    Z = [node1(3) node2(3) node3(3) node4(3)]';    
    T = [1 2 3; 1 2 4; 2 3 4; 1 3 4]; 

    xxmin=min([node1(1),node2(1),node3(1),node4(1)]);
    yymin=min([node2(2),node2(2),node3(2),node4(2)]);
    zzmin=min([node1(3),node2(3),node3(3),node4(3)]);

    xxmax=max([node1(1),node2(1),node3(1),node4(1)]);
    yymax=max([node2(2),node2(2),node3(2),node4(2)]);
    zzmax=max([node1(3),node2(3),node3(3),node4(3)]);

    xrange = 0.1*(xxmax - xxmin);
    yrange = 0.1*(yymax - yymin);
    zrange = 0.1*(zzmax - zzmin);

    th1 = trisurf(T,X,Y,Z); hold on;
    alpha(th1,0.5)

    scatter3(node1(1),node1(2),node1(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node2(1),node2(2),node2(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node3(1),node3(2),node3(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node4(1),node4(2),node4(3),'filled','MarkerFaceColor','r'); hold on;

    text(node1(1) + 0.5 * xrange,node1(2) + 0.5 * yrange,node1(3) + 0.5 * zrange,'1','FontSize', 15,'Color','red')
    text(node2(1) + 0.5 * xrange,node2(2) + 0.5 * yrange,node2(3) + 0.5 * zrange,'2','FontSize', 15,'Color','red')
    text(node3(1) + 0.5 * xrange,node3(2) + 0.5 * yrange,node3(3) + 0.5 * zrange,'3','FontSize', 15,'Color','red')
    text(node4(1) - 0.5 * xrange,node4(2) + 0.5 * yrange,node4(3) - 0.5 * zrange,'4','FontSize', 15,'Color','red')


    set(th1, 'FaceColor', 'b')
    xlabel('e1') % x-axis label
    ylabel('e2') % y-axis label
    zlabel('e3') % z-axis label
    xlim([xxmin-xrange xxmax+xrange])
    ylim([yymin-yrange yymax+yrange])
    zlim([zzmin-zrange zzmax+zrange])
    daspect([1, 1, 1])
    axis equal

    az = 17;
    el = 64;
    view(az, el);

    end
%===============================================================
    function face = plot_tetra_face(n1, n2, n3, b1, b2)
    % ***********************************************************
    % Computes and draws a point at BARYCENTRIC COORDINATES b1,b2 
    % on the tetrahedron face formed by node1, node2 and node3.
    % ***********************************************************
    face = b1 * n1 + b2 * n2 + (1 - b1 - b2) * n3;

    h = scatter3(face(1),face(2),face(3),'filled','MarkerFaceColor','m'); hold on;
    set(h,'SizeData',100);
    end
%===============================================================
    function edge = plot_tetra_edge(n1, n2, b)
    % ***********************************************************
    % Computes and draws point at b on the tetrahedron edge
    % formed by node1 and node2.
    % b = 1 -> at node1
    % b = 0 -> at node2
    % ***********************************************************
    edge = b*n1+(1-b)*n2;

    h = scatter3(edge(1),edge(2),edge(3),'filled','MarkerFaceColor','g'); hold on;
    set(h,'SizeData',100);
    
    end
%===============================================================
    function plot_tetra_point(x, y, z, col)
    % ***********************************************************
    % Draws a point at [x y z] in the color "col"
    % ***********************************************************

    h = scatter3(x, y, z,'filled','MarkerFaceColor',col); hold on;
    set(h,'SizeData',100);
    
    end    
%===============================================================
    function plot_tetra_dual(node1A, node2A, node3A, node4A, node1B, node2B, node3B, node4B)
    % ***********************************************************
    % Draws two tetrahedrons using the eight nodes
    % (node1A, node2A, node3A, node4A, node1B node2B ....).
    %
    % This is usefull to visualize the transformation between two
    % tetrahedronds.
    % ***********************************************************
    figure

    node1A = double(node1A);
    node2A = double(node2A);
    node3A = double(node3A);
    node4A = double(node4A);

    node1B = double(node1B);
    node2B = double(node2B);
    node3B = double(node3B);
    node4B = double(node4B);
    
    XA = [node1A(1) node2A(1) node3A(1) node4A(1)]';
    YA = [node1A(2) node2A(2) node3A(2) node4A(2)]';
    ZA = [node1A(3) node2A(3) node3A(3) node4A(3)]';    
    TA = [1 2 3; 1 2 4; 2 3 4; 1 3 4]; 

    XB = [node1B(1) node2B(1) node3B(1) node4B(1)]';
    YB = [node1B(2) node2B(2) node3B(2) node4B(2)]';
    ZB = [node1B(3) node2B(3) node3B(3) node4B(3)]';    
    TB = [1 2 3; 1 2 4; 2 3 4; 1 3 4]; 

    xxminA=min([node1A(1),node2A(1),node3A(1),node4A(1)]);
    yyminA=min([node2A(2),node2A(2),node3A(2),node4A(2)]);
    zzminA=min([node1A(3),node2A(3),node3A(3),node4A(3)]);

    xxmaxA=max([node1A(1),node2A(1),node3A(1),node4A(1)]);
    yymaxA=max([node2A(2),node2A(2),node3A(2),node4A(2)]);
    zzmaxA=max([node1A(3),node2A(3),node3A(3),node4A(3)]);

    xrangeA = 0.1*(xxmaxA - xxminA);
    yrangeA = 0.1*(yymaxA - yyminA);
    zrangeA = 0.1*(zzmaxA - zzminA);

    xxminB=min([node1B(1),node2B(1),node3B(1),node4B(1)]);
    yyminB=min([node2B(2),node2B(2),node3B(2),node4B(2)]);
    zzminB=min([node1B(3),node2B(3),node3B(3),node4B(3)]);

    xxmaxB=max([node1B(1),node2B(1),node3B(1),node4B(1)]);
    yymaxB=max([node2B(2),node2B(2),node3B(2),node4B(2)]);
    zzmaxB=max([node1B(3),node2B(3),node3B(3),node4B(3)]);

    xrangeB = 0.1*(xxmaxB - xxminB);
    yrangeB = 0.1*(yymaxB - yyminB);
    zrangeB = 0.1*(zzmaxB - zzminB);

    th2 = trisurf(TB,XB,YB,ZB); hold on;
    set(th2, 'FaceColor', 'g')
    alpha(th2,0.7)
 
    th1 = trisurf(TA,XA,YA,ZA); hold on;
    set(th1, 'FaceColor', 'b')
    alpha(th1,0.7)
 
    xlabel('e1') % x-axis label
    ylabel('e2') % y-axis label
    zlabel('e3') % z-axis label
    axis equal

    scatter3(node1A(1),node1A(2),node1A(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node2A(1),node2A(2),node2A(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node3A(1),node3A(2),node3A(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node4A(1),node4A(2),node4A(3),'filled','MarkerFaceColor','r'); hold on;

    text(node1A(1) + 0.5 * xrangeA,node1A(2) + 0.5 * yrangeA,node1A(3) + 0.5 * zrangeA,'1','FontSize', 15,'Color','red')
    text(node2A(1) + 0.5 * xrangeA,node2A(2) + 0.5 * yrangeA,node2A(3) + 0.5 * zrangeA,'2','FontSize', 15,'Color','red')
    text(node3A(1) + 0.5 * xrangeA,node3A(2) + 0.5 * yrangeA,node3A(3) + 0.5 * zrangeA,'3','FontSize', 15,'Color','red')
    text(node4A(1) - 0.5 * xrangeA,node4A(2) + 0.5 * yrangeA,node4A(3) - 0.5 * zrangeA,'4','FontSize', 15,'Color','red')

    text(node1B(1) + 0.5 * xrangeB,node1B(2) + 0.5 * yrangeB,node1B(3) + 0.5 * zrangeB,'1','FontSize', 15,'Color','red')
    text(node2B(1) + 0.5 * xrangeB,node2B(2) + 0.5 * yrangeB,node2B(3) + 0.5 * zrangeB,'2','FontSize', 15,'Color','red')
    text(node3B(1) + 0.5 * xrangeB,node3B(2) + 0.5 * yrangeB,node3B(3) + 0.5 * zrangeB,'3','FontSize', 15,'Color','red')
    text(node4B(1) - 0.5 * xrangeB,node4B(2) + 0.5 * yrangeB,node4B(3) - 0.5 * zrangeB,'4','FontSize', 15,'Color','red')

    scatter3(node1B(1),node1B(2),node1B(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node2B(1),node2B(2),node2B(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node3B(1),node3B(2),node3B(3),'filled','MarkerFaceColor','r'); hold on;
    scatter3(node4B(1),node4B(2),node4B(3),'filled','MarkerFaceColor','r'); hold on;

    quiver3(node1A(1),node1A(2),node1A(3),node1B(1)-node1A(1),node1B(2)-node1A(2),node1B(3)-node1A(3),'m')
    quiver3(node2A(1),node2A(2),node2A(3),node2B(1)-node2A(1),node2B(2)-node2A(2),node2B(3)-node2A(3),'m')
    quiver3(node3A(1),node3A(2),node3A(3),node3B(1)-node3A(1),node3B(2)-node3A(2),node3B(3)-node3A(3),'m')
    quiver3(node4A(1),node4A(2),node4A(3),node4B(1)-node4A(1),node4B(2)-node4A(2),node4B(3)-node4A(3),'m')

    az = 17;
    el = 64;
    view(az, el);
    axis equal
    end

%===============================================================
    function gradM = gradientTransform(a,b)
    % ***********************************************************
    % Computes the gradient of the motion:
    % F = nambla_x y(x,t)
    % ***********************************************************
    gradM=vpa(zeros(3,3));

    for i=1:3
        for j=1:3
            gradM(i,j) = diff(a(i),b(j));
        end 
    end
    end
%===============================================================
    function [normal, centeroid] = get_tetra_normal(curNode1,curNode2,curNode3, curNode4)
    % ***********************************************************
    % Returns two 3*4 matrices. The columns of the first matrix 
    % correspond to the normals of the tetrahedron faces defined
    % by the four input nodes. The columns of the second matrix
    % correspond to the location of the face centeroids.
    % ***********************************************************
    faceNormal = @(curNode1,curNode2,curNode3) cm.cross_product((curNode1-curNode2),(curNode3-curNode2)) / cm.norm(cm.cross_product((curNode1-curNode2),(curNode3-curNode2)));  
    faceCenteroid = @(nodeA,nodeB,nodeC)  [(nodeA(1) + nodeB(1) + nodeC(1)) / 3;(nodeA(2) + nodeB(2) + nodeC(2)) / 3;(nodeA(3) + nodeB(3) + nodeC(3)) / 3];
    
    if strcmp(class(curNode1),'sym') == 1
        disp('symbolic!')
        %*********************************************************************   
        normal1 = -1 * faceNormal(curNode1,curNode2,curNode3);
        centeroid1 = faceCenteroid(curNode1,curNode2,curNode3);
        %*********************************************************************
        normal2 = faceNormal(curNode1,curNode2,curNode4);
        centeroid2 = faceCenteroid(curNode1,curNode2,curNode4);
        %*********************************************************************    
        normal3 = -1 * faceNormal(curNode3,curNode2,curNode4);
        centeroid3 = faceCenteroid(curNode3,curNode2,curNode4); 
        %*********************************************************************    
        normal4 = -1 * faceNormal(curNode1,curNode3,curNode4);
        centeroid4 = faceCenteroid(curNode1,curNode3,curNode4);
    else
        %*********************************************************************   
        normal1 = -1 * faceNormal(curNode1,curNode2,curNode3);
        centeroid1 = faceCenteroid(curNode1,curNode2,curNode3);
        remaining = curNode4;
        if norm(remaining - centeroid1 + normal1) > norm(remaining - centeroid1 - normal1)
            normal1 = normal1 * (-1);
        end
        %*********************************************************************
        normal2 = faceNormal(curNode1,curNode2,curNode4);
        centeroid2 = faceCenteroid(curNode1,curNode2,curNode4);
        remaining = curNode3;
        if norm(remaining - centeroid2 + normal2) > norm(remaining - centeroid2 - normal2)
            normal2 = normal2 * (-1);
        end
        %*********************************************************************    
        normal3 = -1 * faceNormal(curNode3,curNode2,curNode4);
        centeroid3 = faceCenteroid(curNode3,curNode2,curNode4); 
        remaining = curNode1;
        if norm(remaining - centeroid3 + normal3) > norm(remaining - centeroid3 - normal3)
            normal3 = normal3 * (-1);
        end
        %*********************************************************************    
        normal4 = -1 * faceNormal(curNode1,curNode3,curNode4);
        centeroid4 = faceCenteroid(curNode1,curNode3,curNode4);
        remaining = curNode2;
        if norm(remaining - centeroid4 + normal4) > norm(remaining - centeroid4 - normal4)
            normal4 = normal4 * (-1);
        end
        %********************************************************************* 
    end   
    %normal = [normal1 normal2 normal3 normal4];
    %centeroid = [centeroid1 centeroid2 centeroid3 centeroid4];   
    normal = [normal3 normal4 normal2 normal1];
    centeroid = [centeroid3 centeroid4 centeroid2 centeroid1];   
    %axis auto 
    end
%===============================================================
 function [curNormals, curCenteroids] = plot_tetra_normal(curNode1,curNode2,curNode3,curNode4)
    % ***********************************************************
    % Returns two 3*4 matrices. The columns of the first matrix 
    % correspond to the normals of the tetrahedron faces defined
    % by the four input nodes. The columns of the second matrix
    % correspond to the locations of the face centeroids.
    % Additionally the four nodes are potted into the current fig.
    % ***********************************************************
    vectorLength = 0.1;
    
    [curNormals, curCenteroids] = cm.get_tetra_normal(curNode1,curNode2,curNode3, curNode4);
    h = scatter3(curCenteroids(1,1),curCenteroids(2,1),curCenteroids(3,1),'filled','MarkerFaceColor','k'); hold on;
    set(h,'SizeData',100);
    quiver3(curCenteroids(1,1),curCenteroids(2,1),curCenteroids(3,1),vectorLength*curNormals(1,1),vectorLength*curNormals(2,1),vectorLength*curNormals(3,1),'m','linewidth',2);hold on;
    %*********************************************************************    
    h = scatter3(curCenteroids(1,2),curCenteroids(2,2),curCenteroids(3,2),'filled','MarkerFaceColor','k'); hold on;
    set(h,'SizeData',100);
    quiver3(curCenteroids(1,2),curCenteroids(2,2),curCenteroids(3,2),vectorLength*curNormals(1,2),vectorLength*curNormals(2,2),vectorLength*curNormals(3,2),'m','linewidth',2);hold on;
    %*********************************************************************    
    h = scatter3(curCenteroids(1,3),curCenteroids(2,3),curCenteroids(3,3),'filled','MarkerFaceColor','k'); hold on;
    set(h,'SizeData',100);
    quiver3(curCenteroids(1,3),curCenteroids(2,3),curCenteroids(3,3),vectorLength*curNormals(1,3),vectorLength*curNormals(2,3),vectorLength*curNormals(3,3),'m','linewidth',2);hold on;
    %*********************************************************************        
    h = scatter3(curCenteroids(1,4),curCenteroids(2,4),curCenteroids(3,4),'filled','MarkerFaceColor','k'); hold on;
    set(h,'SizeData',100);
    quiver3(curCenteroids(1,4),curCenteroids(2,4),curCenteroids(3,4),vectorLength*curNormals(1,4),vectorLength*curNormals(2,4),vectorLength*curNormals(3,4),'m','linewidth',2);hold on;
    %********************************************************************* 
    axis auto
    end
%===============================================================
    function dataOut = roundDecimals(data,n)
    % ***********************************************************
    % Rounds data to n digits after the coma.
    % Same as round(number,n)
    % ***********************************************************
    rndFact = 10^n;
    dataOut = vpa(round(data*rndFact) / rndFact);

    end
%===============================================================
    function visualizeTensorAtP(D,p)
    % ***********************************************************
    % Plots a tensor D in form of an ellipsoid at a given location p.
    % D is a 3*3 matrix, p is a 3*1 vector.
    % ***********************************************************
    
    D = double(D);
    p = double(p);
    
    sz=size(D);
    if length(sz)==2
        nx=1;ny=1;
    elseif length(sz)==3
        nx=sz(3);ny=1;
    elseif length(sz)==4
        nx=sz(3);ny=sz(4);
    end

    n=size(D,3);
    for i=1:nx
        for j=1:ny
            [v,l]=eig(D(:,:,i,j));
            [X,Y,Z]=ellipsoid(0,0,0,l(1,1),l(2,2),l(3,3),20);
            sz=size(X);
            for x=1:sz(1)
                for y=1:sz(2)
                    A=[X(x,y) Y(x,y) Z(x,y)]';
                    A=v*A;
                    X(x,y)=A(1);Y(x,y)=A(2);Z(x,y)=A(3);
                end
            end
            X=X+(i-1)*2 + p(1);
            Y=Y+(j-1)*2 + p(2);
            Z = Z+p(3);
            h(i)=surf(X,Y,Z);
            if i==1 & j==1
                hold;
            end
        end
    end

    if exist('az') > 0
        view(az,el); %
    end

    xlabel('X') % x-axis label
    ylabel('Y') % y-axis label
    zlabel('Z') % z-axis label
    axis equal
    colormap([0.2 0.2 0.7])
    end
%===============================================================
%===============================================================
%===============================================================




    % *******************************************************
    % END OF FUNCTION DEFINITION
    % *******************************************************
    end
end




