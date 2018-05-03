function b = is_outside_prior(theta,prior_width,params)
b = max(max((theta>prior_width/2+params.ref),(theta<params.ref-prior_width/2)));
