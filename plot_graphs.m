function plot_graphs(G)

numGraphs = length(G);
for i = 1:numGraphs
  figure
  plot(graph(G(i).W))
end

end
