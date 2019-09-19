v2noatt = load('PPI_V2xNoAtt');
v2att   = load('PPI_V2xAtt.mat');
v5noatt = load('PPI_V5xNoAtt.mat');
v5att   = load('PPI_V5xAtt.mat');
figure(101);
plot(v2noatt.PPI.ppi,v5noatt.PPI.ppi,'k.');
hold on;
plot(v2att.PPI.ppi,v5att.PPI.ppi,'r.');

%x is design matrix
x  = v2noatt.PPI.ppi(:);
x  = [x, ones(size(x))];

y  = v5noatt.PPI.ppi(:);

%least squares fit
B  = x\y;
y1 = B(1)*x(:,1)+B(2);
plot(x(:,1),y1,'k-');



x = v2att.PPI.ppi(:);
x = [x, ones(size(x))];
y = v5att.PPI.ppi(:);
B = x\y;
y1 = B(1)*x(:,1)+B(2);
plot(x(:,1),y1,'r-');
legend('No Attention','Attention')
xlabel('V2 activity')
ylabel('V5 response')
title('Psychophysiologic Interaction')