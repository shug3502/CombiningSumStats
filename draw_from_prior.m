function samples = draw_from_prior(n,prior_width,params)
samples = prior_width*rand(n,params.num_params)-prior_width/2 +repmat(params.ref,n,1);
