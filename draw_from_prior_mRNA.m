function th = draw_from_prior_mRNA(n,num_params)
proposal = [abs(randn(n,num_params-1)),(rand(n,1)+1)/2];
%want to transform this to avoid excessively long jumps
%could do so by rejection sampling instead
max_jump_length = 20;
a = atan2(1,max_jump_length);
gradient = 1 - 2*a/pi;
r = sqrt(proposal(:,1).^2 + proposal(:,4).^2); %radius in appropriate plane
alpha = atan2(proposal(:,4),proposal(:,1));
%squish_transform = @(x) gradient*x + a; %define this linear transform in
%radial coordinates
transformed = [r.*cos(gradient*alpha + a), r.*sin(gradient*alpha + a)];

th = proposal;
th(:,[1,4]) = transformed;
end
