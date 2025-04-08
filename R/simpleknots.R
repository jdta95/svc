# generate knots
## 1 in every k^2 point on a grid is a knot
## ensures knot data does not become mismatched
simpleknots = function(
  Y,
  X,
  coords,
  k = 2
) {
  coord_names = colnames(coords)
  
  if(is.null(coord_names)){
    colnames(coords) = c("lat", "lon")
  }
  
  lat = sort(unique(coords[, 1]))
  lon = sort(unique(coords[, 2]))
  
  lat_knots = lat[seq(k, length(lat), by = k)]
  lon_knots = lon[seq(k, length(lon), by = k)]
  
  knots = as.matrix(expand.grid(lat_knots, lon_knots))
  
  full_df = data.frame(coords, Y, X)
  
  knots_df = data.frame(knots)
  
  colnames(knots_df) = colnames(coords)
  
  knots_df = merge(knots_df, full_df, by = colnames(coords), all.x = TRUE, all.y = FALSE)
  
  # if Y, X, or coords had names, assign them to the cols of knots_df
  colnames(knots_df)[3] = "Y"
  if (!is.null(colnames(X))) {
    colnames(knots_df)[4:ncol(knots_df)] = colnames(X)
  }
  colnames(knots_df)[1:2] = coord_names
  
  knots = as.matrix(knots_df[, 1:2])
  Y_knots = knots_df[, 3]
  X_knots = as.matrix(knots_df[, 4:ncol(knots_df)])
  
  return(list(
    Y_knots = Y_knots,
    X_knots = X_knots,
    knots = knots,
    knots_data_frame = knots_df
  ))
}