function prob = ads_prob1d(V, B, ph, pg)

% boundary codition type
% 'D' for Dirichlet, 'R' for Robin, 'P' for periodic
bd = 'R';
M = length(V);

    function a = a(~)
        a = 1;
    end

    function b = b(x)
        b = B * (x-1)*x;
    end

    function c = c(x)
        ind = floor(x*M) + 1;
        c = V(ind);
    end

    function f = f(~)
        f = 1;
    end

    function h = h(~)
        h = ph;
    end

    function g = g(~)
        g = pg;
    end

prob = struct('a', @a, 'b', @b, 'c', @c, 'f', @f, ...
              'h', @h, 'g', @g, 'bd', bd);
end
