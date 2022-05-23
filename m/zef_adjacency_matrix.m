function [A, I, J] = zef_adjacency_matrix(nodes, tetra)

    % Constructs a sparse square adjacency matrix or stensil A for a given set
    % of nodes and tetrahedra constructed from them. Also returns the nonzero
    % indices I and J of A.

    N = size(nodes,1);
    A = spalloc(N,N,0);

    for i = 1 : 4

        for j = i + 1 : 4

            A = A + sparse(            ...
                tetra(:,i),            ...
                tetra(:,j),            ...
                ones(size(tetra,1),1), ...
                N,                     ...
                N                      ...
            );

        end
    end

    % Stensils are symmetric, as they describe an undirected graph.

    A = A + A';

    % Take care of the diagonal.

    for i = 1 : 4

        A = A + sparse(            ...
            tetra(:,i),            ...
            tetra(:,i),            ...
            ones(size(tetra,1),1), ...
            N,                     ...
            N                      ...
        );

    end

    % Find indices and values of nonzero elements and force them into ones.

    [I,J,K] = find(A);

    K = ones(size(K));

    A = sparse(I,J,K);

end
