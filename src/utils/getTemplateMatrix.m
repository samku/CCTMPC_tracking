function F = getTemplateMatrix(nx,f)
% Generate a matrix template in nx-dimensions with a given f complexity.

% Define a 2D polar representation of the polytope template
% (according to the original paper)
F_perm = [];
angle = 2*pi/f;
for i = 1:f
    F_perm = [F_perm; cos((i-1)*angle) sin((i-1)*angle)];
end

% Project 2D template into higher dimensions (is necessary)
if nx > 2
    x_permutations = permn(1:nx,2);
    F = [eye(nx);-eye(nx)];
    for i=1:size(x_permutations,1)
        if x_permutations(i,1)~=x_permutations(i,2)
            F = [F; ones(floor(f/2),nx)];
            F = [F; -ones(f-floor(f/2),nx)];
            F(end-f+1:end,x_permutations(i,:)) = F_perm;
        end
    end
elseif nx==1
    F = [1;-1]; % only two possible directions.
else
    F = F_perm; % return the 2D template
end

end


function perms = permn(items, r)
% permutations with repetitions

% Create a cell array with 'r' copies of 'items'
itemsCell = repmat({items}, 1, r);

% Generate all combinations (permutations with repetition)
[itemsGrid{1:r}] = ndgrid(itemsCell{:});

% Convert the result to a matrix where each row is a permutation
perms = reshape(cat(r+1, itemsGrid{:}), [], r);
end

