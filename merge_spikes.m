function merged = merge_spikes(spikes)
spikes(end) = 0;
T = length(spikes);
merged = zeros(size(spikes));
t = 1;
while t <T+1
    index = t;
    while spikes(index) ~= 0
    index = index+1;
    end
    merged(t) = sum(spikes(t:index));
    t = index+1;
end
end