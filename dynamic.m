function [ out] = dynamic( input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
 input = zeros(3,100);
 input(1,[5:20:100]) = 1;
 input(2,[30:70]) = 1;
 %input(1,20:25) = 1;
end

X = zeros(3,size(input,2));
inhib_weight   = -0.4;
input_weight   =  0.5;
forward_weight =  0.3;
mod_weight     =  0.3;
back_weight    = -0.1;
A = diag(ones(1,3)) * inhib_weight;
A(2,1) = forward_weight;
A(3,2) = forward_weight;
A(1,2) = back_weight;
A(2,3) = back_weight;

B = zeros(3);
%B(2,1) = mod_weight;
B(2,1) = mod_weight;

C = zeros(3);
C(1,1) = input_weight;
%C(3,2) = input_weight;

for step = 2:size(input,2)
 X(:,step) = X(:,step-1) + A*X(:,step-1) + input(2,step)*B*X(:,step-1) + C*input(:,step);
end
 figure(301);
 subplot(2,1,2);
 plot(X');
 legend('Region A','Region B','Region C')
 subplot(2,1,1);
 plot(input');
 legend('Input 1','Input 2','Input 3')
 %pause
end

