model.matrix.spbp <-  function(object = spbp, data = get(as.character(object$call$data)), ...){
   model.matrix.default(object = object, data = data, ...)[,-1]
}
