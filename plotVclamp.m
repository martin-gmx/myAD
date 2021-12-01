function plotVclamp

for i = 1:3
    Vclamp = getVclamp(i-1);
    Vclamp(:,1) = cumsum(Vclamp(:,1));
    Vclamp = [0 Vclamp(1,2); Vclamp];
    figure(i); clf;
    plot(Vclamp(:,1)/1e3,Vclamp(:,2));
    xlabel('t (s)'); ylabel('V (mV)');
    box off; xlim(Vclamp([1 end],1)/1e3);
end
