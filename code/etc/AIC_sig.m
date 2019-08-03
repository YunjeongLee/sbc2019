function a = AIC_sig(cost, num_data)

% error_MLE=((num_data-num_par)/num_data)*sigma_est_square;

a = num_data*log(cost);

end