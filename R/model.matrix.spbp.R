model.matrix.spbp <-  function(object, data = get(as.character(object$call$data)), ...){
   survival:::model.matrix.coxph(object, data = data, ...)
}

model.frame.spbp <-  function(formula, data = get(as.character(formula$call$data)), ...){
  survival:::model.frame.coxph(formula, data = data, ...)
}

