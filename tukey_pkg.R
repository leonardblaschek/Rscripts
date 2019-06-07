tukey <- function(df, y, x, pos, alpha) {
  aov1 <- aov(data = df, y ~ x)
  # groups <- HSD.test(aov1, x, alpha = alpha)
  # groups$groups[[x]] <- rownames(groups$groups)
  # groups$groups[[y]] <- pos
  # return(groups[["groups"]])
  return(aov1)
}

tukey_groups <- function(df, y_var, x_var, pos_var, alpha_var, grp_var){
  # group_by_() .dots argument
  group_dots <- interp(~ grp_var, grp_var = as.name(grp_var))
  # do_() .dots argument
  do_dots = interp( ~ tukey(df = .,
                            y = y_var,
                            x = x_var, 
                            pos = pos_var,
                            alpha = alpha_var),
                    y_var = as.name(y_var),
                    x_var = as.name(x_var),
                    pos_var = as.numeric(pos_var),
                    alpha = as.numeric(alpha_var))
  # Operations
  out <- df %>%
    group_by_(.dots = group_dots) %>%
    do_(.dots = do_dots)
  return(out)
}