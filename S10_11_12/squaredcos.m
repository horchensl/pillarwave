function x = squaredcos(ax,freqs)

% ax - axis
% freqs - parameter frequencies in Hz

if isa(ax,'ax')
    do_ifft=false;
    if ~isfrequencydomain(ax)     % transform axis to frequency domain if needed
        do_ifft=true;
        ax=fft(ax);
    end
    f = ax.values;
else
    f=abs(ax);
end

m=zeros(size(f));
range = find(abs(f)>freqs(1) & abs(f)<freqs(2));
m(range) = cos((f(range)-freqs(2)).*pi./2./(freqs(1)-freqs(2))).^2;
range = find(abs(f)>=freqs(2) & abs(f)<=freqs(3));
m(range) = ones(size(range));
range = find(abs(f)>freqs(3) & abs(f)<freqs(4));
m(range) = cos((f(range)-freqs(3)).*pi./2./(freqs(3)-freqs(4))).^2;

if isa(ax,'ax')
    x=data(m,ax);
    if do_ifft
        x=real(ifft(x));
    end
else
    x=m;
end


