function T = RecMesh(nx, ny, L, R, L0, R0)
    % This function is designed for the Direct Boundary Problem in order to
    % move the whole mesh
    % 
    
    % Assign default arguments:
    if nargin < 6; R0 = 0; end
    if nargin < 5; L0 = 0;end
    if nargin == 3
        R = L;
    elseif nargin < 3 
        L = 1; R = 1;
    end
    if nargin < 2; ny = nx; end
    
    % Compute the number of nodes and allocate space for Nodes.
    
    Nv=(nx+1)*(ny+1);
    T.N = Nv;
    T.Nx = nx;
    T.Ny = ny;
    T.Nodes=zeros(Nv,2);
    k=0;    dx=L/nx;    dy=R/ny;
    
    for j=0:ny
       y = R0 + j*dy; 
       for i=0:nx
          x = L0 + i*dx;
          k = k+1;
          T.Nodes(k,:) = [x,y];
       end
    end
    
    
    % Compute the number of triangles and allocate space for Elements:
    
    T.Ne = 3*nx*ny + nx + ny;
    T.Nt = 2*nx*ny;
    T.Edge = zeros(T.Ne, 2);
    T.Tri = zeros(T.Nt, 3);
    T.TrEg = zeros(T.Nt, 3);
    T.EdgeFlag = zeros(T.Ne, 1);
    ind = 1;
    for j = 1:ny+1
         if j==1;      T.EdgeFlag = 1;   end
         if j==(ny+1); T.EdgeFlag = 3;   end
         for i = 1:nx
             T.Edge(ind, :) = [(j-1)*(nx+1)+i (j-1)*(nx+1)+i+1];
             ind = ind+1;
         end
    end
    ind1 = ind-1;
    for i = 1:nx+1
        if i == 1;      T.EdgeFlag = 4;   end
        if i == (nx+1); T.EdgeFlag = 2;   end
        for j = 1:ny
            T.Edge(ind, :) = [(nx+1)*(j-1)+i (nx+1)*j+i];
            ind = ind+1;
        end
    end
    ind2 = ind-1;
    for j = 1:ny
        for i = 1:nx
            T.Edge(ind, :) = [(nx+1)*(j-1)+i+1 (nx+1)*j+i];
            ind = ind+1;
        end
    end
    ind = 1;
    for j = 1:ny
        for i = 1:nx
            T.Tri(ind:ind+1,:) = [(j-1)*(nx+1)+i (j-1)*(nx+1)+i+1 j*(nx+1)+i; j*(nx+1)+i (j-1)*(nx+1)+i+1 j*(nx+1)+i+1];
            T.TrEg(ind:ind+1,:) = [(j-1)*nx+i (j-1)*nx+i+ind2 (i-1)*ny+j+ind1;
                (j-1)*nx+i+ind2 i*ny+j+ind1 j*nx+i];
            ind = ind+2;
        end
    end
end