function Mat = Inverse_vectorization(vec,len)
%INVERSE_VECTORIZATION convert a vector into a matrix
%   vec: vector, len: column length
numColumns = length(vec)/len;
Mat = zeros(len,numColumns);
for i_column = 1:numColumns
    Mat(:,i_column) = vec((i_column-1)*len+1:i_column*len);
end

end

