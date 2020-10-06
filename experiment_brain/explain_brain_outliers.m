load('data/brain_data_66')

for i = 1:size(CC, 3)
  brain_graph = CC(:,:,i)*100;
  CC(:,:,i) = double(brain_graph >= min(max(brain_graph)));
  
  fprintf('%d: %d edges\n', i, sum(sum(triu(CC(:,:,i)))));
end

brain_pairs_similarity = [];
for i = 1:size(CC, 3)
  for j = (i+1):size(CC, 3)
    pair_similarity = sum(sum(triu(CC(:,:,i), 1).*triu(CC(:,:,j), 1)));
    brain_pairs_similarity = [brain_pairs_similarity pair_similarity];
  end
end

[max_vals, max_idxs] = maxk(brain_pairs_similarity, 3);
figure
stem(brain_pairs_similarity, 'LineWidth', 2)
hold on
stem(max_idxs, max_vals, 'LineWidth', 2)
hold off
grid on
xlabel('Brain graph pairs')
ylabel('Pair similarity')
