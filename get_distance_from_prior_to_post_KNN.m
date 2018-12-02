function d_h = get_distance_from_prior_to_post_KNN(current_weights,sstar,ss,theta_star,t,prior_sample,params)
    M = ceil(params.N/params.alpha);
    dist = weighted_distance(sstar,repmat(ss,M,1,1),10.^current_weights);
    [~,sort_ind] = sort(dist);
    selected = sort_ind(1:params.N);
    %once decided distances, then can select samples
    theta_selected = theta_star(selected,:);  
    d_h = abs(hellinger_knn_estimator(prior_sample,theta_selected,5)); %based on theoretical argument, switch the order here so distance from posterior to prior
