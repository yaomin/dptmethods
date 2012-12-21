replace.names <-
function(names, target, replace) {  
  names[match(target, names)] <- replace
  names
}
