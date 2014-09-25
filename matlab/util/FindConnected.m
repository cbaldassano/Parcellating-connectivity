function visited = FindConnected(A, seed)
visited = [];
to_visit = seed;
while(~isempty(to_visit))
    if (mod(length(visited),100)==1)
        disp(['V: ' num2str(length(visited)) ', TV: ' num2str(length(to_visit)) ', SUM: ' num2str(length(visited) + length(to_visit))]);
    end
    visited = [visited to_visit(1)];
    suc = find(A(to_visit(1),:));
    to_visit = [to_visit suc(~ismember(suc, to_visit) & ~ismember(suc, visited))];
    to_visit(1) = [];
end
end