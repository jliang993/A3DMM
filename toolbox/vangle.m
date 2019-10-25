function t = vangle(v1, v2)

if size(v1, 2)==1
    t = ((v1')*v2) / (norm(v1)*norm(v2));
else
    t = (v1*(v2')) / (norm(v1)*norm(v2));
end

% t = abs(t);