model.frame.spbp <-
  function(formula, data = get(as.character(formula$call$data)), ...){
  survival:::model.frame.coxph(formula, data = data, ...)
}

