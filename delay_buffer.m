function varargout = delay_buffer(cmd, varargin)
%DELAY_BUFFER  Fixed-horizon buffer returning the most recent command within the age bound.

switch lower(cmd)

case 'init'
    p  = varargin{1};
    x0 = varargin{2};
    Ts = p.spec.num.Ts_plant;
    N  = ceil(p.cyber.tau_buffer/Ts) + 2;

    assert(p.cyber.tau_buffer > p.cyber.tau_dos_hi, ...
        ['delay_buffer: horizon %.3f s does not exceed the worst-case age of ' ...
         'data %.3f s. This is the exact condition the legacy code violated.'], ...
         p.cyber.tau_buffer, p.cyber.tau_dos_hi);

    buf.X = repmat(x0(:), 1, N);
    buf.T = -inf(1, N);
    buf.T(1) = 0;
    buf.head = 1;
    buf.count = 1;
    buf.N = N;
    buf.Ts = Ts;
    buf.horizon = Ts*(N-1);
    buf.maxAgeSeen = 0;
    buf.nClamped = 0;
    buf.nx = numel(x0);
    varargout{1} = buf;

case 'push'
    buf = varargin{1};
    t   = varargin{2};
    x   = varargin{3};
    buf.head = mod(buf.head, buf.N) + 1;
    buf.X(:, buf.head) = x(:);
    buf.T(buf.head) = t;
    buf.count = min(buf.count + 1, buf.N);
    varargout{1} = buf;

case 'fetch'
    buf   = varargin{1};
    t_now = varargin{2};
    age   = varargin{3};
    if nargin >= 5, allow_clamp = varargin{4}; else, allow_clamp = false; end

    buf.maxAgeSeen = max(buf.maxAgeSeen, age);

    k = round(age/buf.Ts);
    clamped = false;
    if k > buf.count - 1
        if allow_clamp
            k = buf.count - 1;
            clamped = true;
            buf.nClamped = buf.nClamped + 1;
        else
            error('delay_buffer:horizon', ...
                ['Requested age %.4f s exceeds the stored history %.4f s ' ...
                 '(capacity %d samples at Ts = %.4f s). Increase ' ...
                 'p.cyber.tau_buffer. Silently clamping here is the defect ' ...
                 'that voided the legacy stochastic study.'], ...
                 age, (buf.count-1)*buf.Ts, buf.N, buf.Ts);
        end
    end
    k = max(k, 0);

    idx = mod(buf.head - k - 1, buf.N) + 1;
    varargout{1} = buf.X(:, idx);
    varargout{2} = t_now - buf.T(idx);
    varargout{3} = struct('index',idx,'k',k,'clamped',clamped, ...
                          'stamp',buf.T(idx),'buf',buf);

case 'fetch_at'
    buf    = varargin{1};
    tstamp = varargin{2};
    valid  = isfinite(buf.T);
    [~, j] = min(abs(buf.T - tstamp) + 1e9*(~valid));
    varargout{1} = buf.X(:, j);
    varargout{2} = buf.T(buf.head) - buf.T(j);
    varargout{3} = struct('index',j,'k',NaN,'clamped',false, ...
                          'stamp',buf.T(j),'buf',buf);

case 'stats'
    buf = varargin{1};
    varargout{1} = struct('capacity',buf.N,'count',buf.count, ...
        'horizon',buf.horizon,'maxAgeSeen',buf.maxAgeSeen, ...
        'nClamped',buf.nClamped,'Ts',buf.Ts);

otherwise
    error('delay_buffer:cmd','Unknown command "%s"', cmd);
end
end
