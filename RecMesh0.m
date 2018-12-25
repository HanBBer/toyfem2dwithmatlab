function T = RecMesh0(nx, ny, L, R, L0, R0)
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
    
    % define NodePtrs:
    
    T.NodePtrs=zeros(Nv,1);
    
    % Compute the number of free nodes and define T.FNodePtrs:

    Nf = Nv - 2*nx - 2*ny;
    T.FNodePtrs = zeros(Nf,1);
    for j = 1:ny-1
       T.FNodePtrs((j-1)*(nx-1)+1:j*(nx-1))=(j*(nx+1)+2:(j+1)*(nx+1)-1)';
    end
    
    % Compute the number of constrained nodes and define CNodePtrs:

    Nc = 2*nx + 2*ny;
    T.CNodePtrs=[ (1:nx+1)'; (nx+2:nx+1:(nx+1)*(ny-1)+1)'; (2*nx+2:nx+1:(nx+1)*ny)'; (nx*(ny+1)+1:(nx+1)*(ny+1))' ];
    T.NodePtrs(T.FNodePtrs)=(1:Nf)';
    T.NodePtrs(T.CNodePtrs)=(1:Nc)';
    
    % Compute the number of triangles and allocate space for Elements:
    
    T.Tri = delaunay(T.Nodes(:,1), T.Nodes(:,2));
    T.Nt = 2*nx*ny;
end