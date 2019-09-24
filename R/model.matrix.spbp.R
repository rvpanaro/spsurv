model.matrix.spbp <-
  function(object, data = get(as.character(object$call$data)), ...){
   survival:::model.matrix.coxph(object, data = data, ...)
}
