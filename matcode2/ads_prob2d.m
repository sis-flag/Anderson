function prob = ads_prob2d(V, B, ph, pg)

% boundary codition type
% 'D' for Dirichlet, 'R' for Robin, 'P' for periodic
bd = 'R';
[Mx, My] = size(V);

    function a = a11(~, ~)
        a = 1;
    end

    function a = a12(~, ~)
        a = 0;
    end

    function a = a22(~, ~)
        a = 1;
    end

    function b = b1(x, y)
        b = B(1) * exp(x);
    end

    function b = b2(x, y)
        b = B(2) * sin(y);
    end

    function c = c(x, y)
        indx = floor(x*Mx) + 1;
        indy = floor(y*My) + 1;
        c = V(indx, indy);
    end

    function f = f(~,~)
        f = 1;
    end

    function h = h(~,~)
        h = ph;
    end

    function g = g(~,~)
        g = pg;
    end

prob = struct(...
    'a11', @a11, 'a12', @a12, 'a21', @a12, 'a22', @a22,...
    'b1', @b1, 'b2', @b2, 'c', @c, 'f', @f, ...
    'h', @h, 'g', @g, 'bd', bd);
end
