linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));

n = 10;

MarkerEdgeColors = jet(n);
Markers = ['o','x','+','*','s','d','v','^','<','>','p','h','.',...
           '+','*','o','x','^','<','h','.','>','p','s','d','v',...
           'o','x','+','*','s','d','v','^','<','>','p','h','.'];

x = 0:0.01:2;

X = repmat(x, n, 1);
Y = exp(-0.5*X.*repmat((1:10)', 1, length(x)));

figure
hold on
for i=1:n
  plot(X(i,:), Y(i,:),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
end
hold off

figure
hold on
for i=1:n
  plot(X(i,:), Y(i,:),[linestyles{i}], 'Color',MarkerEdgeColors(i,:));
end
hold off
