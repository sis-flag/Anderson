function mark = my_watershed(w)

N = size(w,1);
mark = zeros(size(w));
num_mark = 10;

% 给所有像素排序，从低到高浸没
[~, ind] = sort(w(:));
for k = 1:length(ind)
    indi = mod(ind(k)-1, N) +1;
    indj = ceil(ind(k) / N);
    
    % 边界上的不算
    if indi == 1 || indi == N || indj == 1 || indj == N
        continue
    end
    
    % 找出当前点附近的标记
    mloc = mark(indi-1:indi+1, indj-1:indj+1);

    % 删去重复的标记，看看标记有几种
    mloc = unique(mloc(:));
    mloc(mloc == 0) = [];
    mloc(isnan(mloc)) = [];
    
    % 如果这点周围全没被浸没过，它就是一种新的水
    if isempty(mloc)
        mark(indi, indj) = num_mark + 1;
        num_mark = num_mark + 1;
        % 如果这点周围只被一种水浸没，它就是也是这种
    elseif length(mloc) == 1
        mark(indi, indj) = mloc;
        % 如果这点周围被很多种水浸没，它就是分水岭
    else
        mark(indi, indj) = nan;
    end
    
%     if mod(k, floor(N*N/10)) == 0
%         figure();
%         hold on
%         s = pcolor(x1,x2, mark);
%         s.LineStyle = 'none';
%         colorbar;
%     end
end