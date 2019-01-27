%dbstop if error;

% This is a script for testing

N = 1;
Nx = N; Ny = N;
T = RecMesh([Nx, Ny], [1,1], [0,0]);
figure(1); ShowMesh(T, [1,1,1]);

T2 = Refine(T);
figure(2); ShowMesh(T2, [1,1,1]);

% Verify the validity of Refine.m
%T3 = Refine(T2);
%figure(3); ShowMesh(T3, [1,1,1]);
