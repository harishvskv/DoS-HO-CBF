function [shed, st] = ufls_relay(df, t, st, spec)
%UFLS_RELAY  Staged under-frequency load-shedding relay, external to the vector field.

if nargin < 4 || isempty(spec), spec = cpmg_spec(); end
nS = numel(spec.ufls.f_Hz);

if isempty(st)
    st = struct();
    st.tripped    = false(1,nS);
    st.armed      = false(1,nS);
    st.arm_time   = nan(1,nS);
    st.trip_time  = nan(1,nS);
    st.shed_total = 0;
    st.latching   = true;
    st.min_df     = 0;
end

st.min_df = min(st.min_df, df);

df_reset = spec.hz2pu(spec.ufls.reset_Hz);

for k = 1:nS
    pickup = spec.ufls.df_pu(k);

    if ~st.tripped(k)
        if df <= pickup
            if ~st.armed(k)
                st.armed(k)    = true;
                st.arm_time(k) = t;
            elseif (t - st.arm_time(k)) >= spec.ufls.delay_s(k)
                st.tripped(k)   = true;
                st.trip_time(k) = t;
            end
        else
            st.armed(k)    = false;
            st.arm_time(k) = NaN;
        end
    else
        if ~st.latching && df >= df_reset
            st.tripped(k)   = false;
            st.armed(k)     = false;
            st.arm_time(k)  = NaN;
        end
    end
end

shed = sum(spec.ufls.shed_frac(st.tripped));
shed = min(max(shed,0),1);
st.shed_total = shed;

end
