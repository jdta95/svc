svclm = function(
    Y,
    X,
    s,
    Y_knots,
    X_knots,
    knots,
    w_knots_start = matrix(0, nrow = nrow(X_knots), ncol = ncol(X_knots)),
    sigmasq_start = rep(1, ncol(X_knots)),
    a_sigmasq = rep(0.001, ncol(X_knots)),
    b_sigmasq = rep(0.001, ncol(X_knots)),
    tausq_start = 1,
    a_tausq = 0.001,
    b_tausq = 0.001,
    phi_start = phi_lower + (phi_upper - phi_lower) / 2,
    phi_lower,
    phi_upper,
    phi_proposal_sd_start = rep(1, ncol(X_knots)),
    phi_target_acceptance = 0.234,
    mcmc = 1000
) {
  # Call the C++ function
  return(
    svclm_cpp(
      Y,
      X,
      s,
      Y_knots,
      X_knots,
      knots,
      w_knots_start,
      sigmasq_start,
      a_sigmasq,
      b_sigmasq,
      tausq_start,
      a_tausq,
      b_tausq,
      phi_start,
      phi_lower,
      phi_upper,
      phi_proposal_sd_start,
      phi_target_acceptance,
      mcmc
    )
  )
}