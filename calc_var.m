function R = calc_var( y )

[p,T] = size(y);
range_ff =[ 0.25,0.5];
        [psd_Y,ff]= pwelch(y',round(T/8),[],1000,1);
        psd_Y = psd_Y';
        ind = double(ff>range_ff(1));
        ind(ff>range_ff(2))= 0;
        ind(ind == 0) = NaN;
        ind = repmat(ind',p,1);
        R = diag(exp(nanmean(log(psd_Y.*ind/2),2)));

    end