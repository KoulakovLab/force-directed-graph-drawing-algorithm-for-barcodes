%   (c) Alex Koulakov (akula@cshl.edu) 2019 
%

function [ C, N1, N2 ] = build_adjacency_matrix( ClusterC1 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

CCC = ClusterC1;

[N1,N2] = size(CCC);
C = zeros(N1+N2, N1+N2);

C(1:N1,(N1+1):(N1+N2))=CCC;
C((N1+1):(N1+N2),1:N1)=CCC';


end

