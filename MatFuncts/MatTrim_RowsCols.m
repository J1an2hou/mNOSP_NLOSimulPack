function A_out = MatTrim_RowsCols(A, rowsToRemove, colsToRemove)
% MatTrim_RowsCols Removes specified rows and columns from any matrix.
%
%   A_out = removeRowsCols_general(A, rowsToRemove, colsToRemove)
%
% Inputs:
%   A              - Input matrix (m x n)
%   rowsToRemove   - Row indices to remove (a row vector)
%   colsToRemove   - Column indices to remove (a row vector)
%   % For example, to remove the row indices 1, 2, 4
%   % then, make rowsToRemove = [1,2,4];
%   % Similarly for colsToRemove

% Output:
%   A_out          - Resulting matrix after removals

    sz = size(A);
    m = sz(1);
    n = sz(2);
    if nargin == 2
        colsToRemove = rowsToRemove;
    end
    % Validate row/column indices
    if any(rowsToRemove < 1 | rowsToRemove > m)
        error('Some row indices are out of bounds.');
    end
    if any(colsToRemove < 1 | colsToRemove > n)
        error('Some column indices are out of bounds.');
    end

    % Determine indices to keep
    rowsToKeep = setdiff(1:m, rowsToRemove);
    colsToKeep = setdiff(1:n, colsToRemove);

    % Handle 2D or 3D input
    if ismatrix(A)
        A_out = A(rowsToKeep, colsToKeep);
    elseif ndims(A) == 3
        A_out = A(rowsToKeep, colsToKeep, :);
    else
        error('Only 2D or 3D matrices are supported.');
    end
end
