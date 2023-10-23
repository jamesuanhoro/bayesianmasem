methods::setValidity("bmasempriors", function(object) {
  if (any(
    c(
      object@lkj_shape, object@sl_par, object@rc_par,
      object@sr_par, object@br_par, object@rm_par
    ) <= 0
  )) {
    paste0(
      "@lkj_shape, @sl_par, @rc_par, @sr_par, @br_par and ",
      "@rm_par must all be greater than 0"
    )
  } else if (any(c(object@lkj_shape, object@rc_par) < 1)) {
    "@lkj_shape and @rc_par must all be at least 1"
  } else {
    TRUE
  }
})
