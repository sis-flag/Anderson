function mark = my_watershed(w)

N = size(w,1);
mark = zeros(size(w));
num_mark = 10;

% �������������򣬴ӵ͵��߽�û
[~, ind] = sort(w(:));
for k = 1:length(ind)
    indi = mod(ind(k)-1, N) +1;
    indj = ceil(ind(k) / N);
    
    % �߽��ϵĲ���
    if indi == 1 || indi == N || indj == 1 || indj == N
        continue
    end
    
    % �ҳ���ǰ�㸽���ı��
    mloc = mark(indi-1:indi+1, indj-1:indj+1);

    % ɾȥ�ظ��ı�ǣ���������м���
    mloc = unique(mloc(:));
    mloc(mloc == 0) = [];
    mloc(isnan(mloc)) = [];
    
    % ��������Χȫû����û����������һ���µ�ˮ
    if isempty(mloc)
        mark(indi, indj) = num_mark + 1;
        num_mark = num_mark + 1;
        % ��������Χֻ��һ��ˮ��û��������Ҳ������
    elseif length(mloc) == 1
        mark(indi, indj) = mloc;
        % ��������Χ���ܶ���ˮ��û�������Ƿ�ˮ��
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