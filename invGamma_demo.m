%% load contour w from file, and generate contour f as a cubic spline
% w and f do not have to be splines, but they need to be continuse and differntiable

datafile = 'fig4_data';  % fig5_data


load(datafile);

fR2C = @(x) complex(x(:,1), x(:,2));
fC2R = @(x) [real(x) imag(x)];

n = numel(w);
fGenPP = @(c) cscvn(fC2R(c([1:end 1]))');
fEvalCurve = @(pp) fR2C( ppval(pp, linspace(pp.breaks(1), pp.breaks(end), n))' );

f = fEvalCurve(fGenPP(cf1));

dw = w([2:end 1])-w;
df = f([2:end 1])-f;


%%
B = fC2R( polygonOffset(w(1:ceil(n/100):end), 0.2, true) );
[x, t] = cdt(B, [], nv0);
x = fR2C(x);
nv = numel(x);


%%
fGammaL = @(z) 2./( (w-z) - dw.*conj((w-z)./dw) );
fGammaW = @(z) 2./conj( (w-z) - dw.*conj((w-z)./dw) )./(z-w).^2;
fGammaMV = @(z) 1./( abs((z-w)).*(w-z) );
fNormalize = @(z) z/sum(z.*dw);
fNormalizedGamma = @(z) fNormalize( fGammaMV(z) );
fMapGamma = @(z) sum( fNormalizedGamma(z).*(f.*dw+df.*(z-w)) );
fMapGammaA = @(z) sum( conj(fNormalizedGamma(z)).*(f.*conj(dw)+df.*conj(z-w)) );


fz = arrayfun(fMapGamma, x);


figure;
fSetMeshColor = @(h, cdata) set(h, 'EdgeAlpha', 0.1, 'CData', cdata, 'FaceColor', 'interp', 'linewidth', 1);
%% draw forward map
iknots = 1:ceil(n/30):n;
iknots = iknots(1:end-1);
subplot(331);
h0 = drawmesh(t, x);  title('Source');
hold on; plot(w(iknots), 'o', 'MarkerFaceColor', [0 0 0.7], 'MarkerEdgeColor', 'none');
plot(w, '-k');
colormap(jet(8192));

subplot(332); 
h = drawmesh(t, fz);  title('Target');
hold on; plot(f(iknots), 'o', 'MarkerFaceColor', [0 0 0.7], 'MarkerEdgeColor', 'none');
plot(f, '-k');
colormap(jet(8192));

fSetMeshColor([h; h0], abs(x-mean(x)));


%% init using reverse map
fGammaL2 = @(z) 2./( (f-z) - df.*conj((f-z)./df) );
fGammaW2 = @(z) 2./conj( (f-z) - df.*conj((f-z)./df) )./(z-f).^2;
fGammaMV2 = @(z) 1./( abs((z-f)).*(f-z) );
fNormalize2 = @(z) z/sum(z.*df);
fNormalizedGamma2 = @(z) fNormalize2( fGammaMV2(z) );
fMapGamma2 = @(z) sum( fNormalizedGamma2(z).*(w.*df+dw.*(z-f)) );

%% inverse map
z = arrayfun(fMapGamma2, fz);
subplot(333); h = drawmesh(t, z); title(0);
fSetMeshColor(h, abs(z-mean(x)));
colormap(jet(8192));

%%
converged = false(nv, 1);
zIters = {};
stats = [ norm( z-x ) mean( abs(z-x) ) minmax( abs(z-x)' ) ];

nIter = 6;
for it=1:nIter
    for i=1:nv
        if converged(i), continue; end

        w_f = fz(i) - fMapGamma(z(i));
        r = w-z(i);

        %% mean value
        gamma = fGammaMV(z(i));
        dgammadz = 3/2*gamma./r;
        EE = sum( gamma.*dw );
        
        a = (sum(dgammadz.*(f.*dw-r.*df))+sum(gamma.*df))/EE ...  
            -sum(gamma.*(f.*dw-r.*df))*sum(dgammadz.*dw) / EE^2;

        b = (sum(conj(dgammadz).*(f.*conj(dw)-conj(r).*df))+sum(conj(gamma).*df))/conj(EE) ...
            -sum(conj(gamma).*(f.*conj(dw)-conj(r).*df))*conj(sum(dgammadz.*dw)) /conj(EE)^2;

        dz = (w_f*conj(a) - conj(w_f)*b)/(abs(a)^2-abs(b)^2);

        stepsize = 1;
        z(i) = z(i) + dz*stepsize;

%         if abs(dz)<1e-5, converged(i) = true; end
    end
    
%     set(h, 'Vertices', fC2R(z), 'CData', abs(z-mean(x))); drawnow;

    zIters{end+1} = z;
    stats(end+1,:) = [ norm( z-x ) mean( abs(z-x) ) minmax( abs(z-x)' ) ];
    if all(converged), break; end
end


for i=1:nIter
    subplot(3,3,i+3); h=drawmesh(t, zIters{i}); title(i);
    fSetMeshColor(h, abs(zIters{end}-mean(x)));
    colormap(jet(8192));
end
