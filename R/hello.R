bifun <- function(mod_in, fish_in, fish_col = NULL, axes = c('1', '2'),
                  arrow = 0.4, ext = 1.5, size = 4, txt = 6, xlims = NULL,
                  ylims = NULL, var_sub = NULL){

  # for MCA
  if('MCA' %in% class(mod_in)){

    # data to plot
    axes <- paste0('Dim.', axes)
    obs <- data.frame(mod_in$ind$coord)
    obs <- obs[, names(obs) %in% axes]
    obs$Life_Stage <- fish_in
    vecs <- data.frame(mod_in$var$coord)
    vecs <- vecs[, names(vecs) %in% axes]

  }

  # for pca
  if('prcomp' %in% class(mod_in)){

    # data to plot
    axes <- paste0('PC', axes)
    obs <- data.frame(mod_in$x[, axes])
    obs$Life_Stage <- fish_in
    vecs <- data.frame(mod_in$rotation[, axes])

  }


  # a little more tweaks to vecs data frame
  # create vecs label for vector labels
  names(vecs) <- c('one', 'two')
  vecs_lab <- ext * vecs
  vecs_lab$labs <- row.names(vecs_lab)
  vecs$lab <- row.names(vecs)

  # remove vectors for easier viz
  if(!is.null(var_sub)){

    var_sub <- paste(var_sub, collapse = '|')
    vecs <- vecs[grepl(var_sub, vecs$lab), ]
    vecs_lab <- vecs_lab[grepl(var_sub, vecs_lab$lab), ]

  }

  ## plots

  # individual points
  p <- ggplot(obs, aes_string(x = axes[1], y = axes[2])) +
    geom_point(aes(colour = Life_Stage), size = size) +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims) +
    scale_colour_discrete('Life Stage') +
    theme_bw()

  # add vectors
  if(!is.null(arrow))
    p <- p + geom_segment(
      data = vecs,
      aes(x = 0, y = 0, xend = one, yend = two),
      arrow = arrow(length = unit(arrow, "cm"))
    )

  # add labels
  if(!is.null(txt))
    p <- p + geom_text(data = vecs_lab, aes(x = one, y = two, label = labs),
                       size = txt)

  return(p)

}
