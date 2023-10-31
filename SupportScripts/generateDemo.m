function demoStruct = generateDemo(D, delta_t)
    demoStruct.t = delta_t * [0:1:length(D)-1];
    demoStruct.pos = D;
    demoStruct.vel = [diff(D,1,2)/delta_t zeros(size(demoStruct.pos,1), 1)];
    demoStruct.acc = [diff(demoStruct.vel,1,2)/delta_t zeros(size(demoStruct.pos,1), 1)];
end