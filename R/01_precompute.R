# ==============================================================================
# 1_PRECOMPUTE_DATA.R
# Pré-calcul des modèles EVT-CART et du Krigeage pour la Web App
# ==============================================================================

library(data.table)
library(sf)
library(rpart)
library(partykit)
library(gstat)
library(sp)
library(dplyr)
library(POT)
library(ggplot2)
library(patchwork)

# --- 1. FONCTIONS UTILITAIRES ET EVT-CART ---
clean_county_name <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- toupper(x)
  x <- gsub(" COUNTY| PARISH| CITY| BOROUGH| MUNICIPALITY", "", x)
  x <- gsub("[^A-Z0-9 ]", "", x)
  x <- trimws(gsub("\\s+", " ", x))
  x
}

gpd_init <- function(y, offset, parms = NULL, wt) {
  list(
    y = y, parms = parms, numresp = 2, numy = 1,
    summary = function(yval, dev, wt, ylevel, digits)
      paste0("gamma=", format(yval[1], digits), ", sigma=", format(yval[2], digits)),
    text = function(yval, dev, wt, ylevel, digits, n, col, ...)
      paste0("g=", format(yval[1], digits))
  )
}

gpd_eval <- function(y, wt, parms) {
  if (length(y) < 10) return(list(label = c(0, 1), deviance = 1e9))
  fit <- tryCatch(fitgpd(y, threshold = 0, est = "mle"), error = function(e) NULL)
  if (is.null(fit)) return(list(label = c(0, 1), deviance = 1e9))
  list(label = c(fit$fitted.values["shape"], fit$fitted.values["scale"]), deviance = fit$deviance)
}

gpd_split <- function(y, wt, x, parms, continuous) {
  n <- length(y)
  goodness  <- rep(-Inf, n - 1)
  direction <- rep(1, n - 1)
  
  fit_p <- tryCatch(fitgpd(y, threshold = 0, est = "mle"), error = function(e) NULL)
  if (is.null(fit_p)) return(list(goodness = goodness, direction = direction))
  ll_p <- -fit_p$deviance / 2
  
  if (continuous) {
    for (i in 1:(n - 1)) {
      if (x[i] == x[i+1]) next
      if (i < 15 || (n - i) < 15) next
      
      fit_l <- tryCatch(fitgpd(y[1:i], threshold = 0, est = "mle"), error = function(e) NULL)
      fit_r <- tryCatch(fitgpd(y[(i+1):n], threshold = 0, est = "mle"), error = function(e) NULL)
      
      if (is.null(fit_l) || is.null(fit_r)) next
      goodness[i]  <- (-fit_l$deviance / 2 + -fit_r$deviance / 2) - ll_p
    }
  }
  list(goodness = goodness, direction = direction)
}

gpd_model <- list(init = gpd_init, eval = gpd_eval, split = gpd_split)

# --- 2. CHARGEMENT DES DONNÉES ---
message("── Chargement du dataset depuis CSV...")
df_model <- data.table::fread("data/processed/df_final.csv")

if (!"density" %in% names(df_model) && all(c("population", "area_km2") %in% names(df_model)))
  df_model[, density := population / area_km2]

df_model[, res_county := clean_county_name(res_county)]
df_model[, res_state  := toupper(trimws(res_state))]

message("── Chargement shapefile...")
county_sf <- sf::st_read("data/raw/cb_2018_us_county_500k", quiet = TRUE)
data("fips_codes", package = "tidycensus")
fips_lookup <- unique(as.data.table(fips_codes)[, .(state_code, state)])

suppressWarnings({ centroids <- sf::st_centroid(sf::st_transform(county_sf, 4326)) })
coords     <- sf::st_coordinates(centroids)
geo_coords <- data.table(
  res_state  = fips_lookup$state[match(county_sf$STATEFP, fips_lookup$state_code)],
  res_county = clean_county_name(county_sf$NAME),
  lon = coords[, 1], lat = coords[, 2]
)

# ── SUPPRESSION DES ANCIENNES COORDONNÉES ──
if ("lon" %in% names(df_model)) df_model[, lon := NULL]
if ("lat" %in% names(df_model)) df_model[, lat := NULL]

df_model <- merge(df_model, geo_coords, by = c("res_state", "res_county"), all.x = TRUE)
df_model <- na.omit(df_model[!is.na(lon) & !is.na(lat)], cols = c("med_age", "med_income", "density", "pct_male", "pct_65p"))

all_counties <- geo_coords[res_state %in% fips_lookup$state & !res_state %in% c("AK", "HI")]
in_dataset <- paste(df_model$res_state, df_model$res_county)
missing_counties <- all_counties[!paste(res_state, res_county) %in% in_dataset]

# --- 3. BOUCLE DE PRÉ-CALCUL SUR LES SEUILS ---
u_values <- c(180, 190, 200, 210, 220, 230) # Tu peux ajouter ou modifier les seuils ici
shiny_results <- list()

for (u in u_values) {
  message(paste("\n========== CALCUL POUR u =", u, "=========="))
  
  df_extremes <- copy(df_model[Y > u])
  df_extremes[, Y_excess := Y - u]
  
  # EVT-CART
  fit_evt <- rpart(Y_excess ~ med_age + med_income + density + pct_male + pct_65p + lon + lat,
                   data = df_extremes, method = gpd_model,
                   control = rpart.control(minsplit = 20, minbucket = 10, maxdepth = 5, cp = 0.001))
  
  df_extremes$gamma <- fit_evt$frame$yval2[fit_evt$where, 1]
  df_extremes$sigma <- fit_evt$frame$yval2[fit_evt$where, 2]
  df_extremes$leaf_idx <- as.integer(rownames(fit_evt$frame)[fit_evt$where])
  leaf_params <- unique(df_extremes[, .(leaf_idx, gamma, sigma)])[order(leaf_idx)]
  
  # DIAGNOSTICS (QQ & RL)
  leaf_ids <- sort(unique(df_extremes$leaf_idx))
  p_u <- nrow(df_extremes) / nrow(df_model)
  qq_list <- list()
  rl_list <- list()
  
  for (lid in leaf_ids) {
    obs <- sort(df_extremes$Y_excess[df_extremes$leaf_idx == lid])
    n_obs <- length(obs)
    xi_h <- df_extremes$gamma[df_extremes$leaf_idx == lid][1]
    sig_h <- df_extremes$sigma[df_extremes$leaf_idx == lid][1]
    
    p_emp <- (1:n_obs - 0.5) / n_obs
    q_theo <- if (abs(xi_h) < 1e-6) -sig_h * log(1 - p_emp) else (sig_h / xi_h) * ((1 - p_emp)^(-xi_h) - 1)
    
    qq_list[[as.character(lid)]] <- ggplot(data.frame(emp = obs, theo = q_theo), aes(x = theo, y = emp)) +
      geom_abline(slope = 1, intercept = 0, color = "#e74c3c", linetype = "dashed") +
      geom_point(color = "#2c3e50", alpha = 0.7) +
      labs(title = paste("Feuille", lid), subtitle = paste("n =", n_obs, "| xi =", round(xi_h, 3)), x = "Théoriques", y = "Empiriques") +
      theme_minimal()
    
    T_vals <- c(2, 5, 10, 20, 50, 100, 200)
    valid_T <- T_vals[T_vals * p_u > 1]
    if (length(valid_T) > 0) {
      p_vals <- 1 - 1 / (valid_T * p_u)
      rl_vals <- if (abs(xi_h) < 1e-6) u - sig_h * log(1 - p_vals) else u + (sig_h / xi_h) * ((1 - p_vals)^(-xi_h) - 1)
      rl_list[[as.character(lid)]] <- ggplot(data.frame(T = valid_T, RL = rl_vals), aes(x = T, y = RL)) +
        geom_line(color = "#e74c3c") + geom_point(color = "#e74c3c") + scale_x_log10() +
        labs(title = paste("Feuille", lid), x = "T (comtés)", y = "Y (décès/100k)") + theme_minimal()
    }
  }
  qq_final <- if (length(qq_list) > 0) wrap_plots(qq_list, ncol = min(length(qq_list), 3)) else NULL
  rl_final <- if (length(rl_list) > 0) wrap_plots(rl_list, ncol = min(length(rl_list), 3)) else NULL
  
  # PRÉDICTION ARBRE
  tree_party <- as.party(fit_evt)
  nodes_train <- predict(tree_party, newdata = df_extremes, type = "node")
  mapping_dt <- unique(data.table(node_party = as.integer(nodes_train), frame_row = as.integer(fit_evt$where)))
  node_to_frame <- setNames(mapping_dt$frame_row, mapping_dt$node_party)
  
  node_ids <- predict(tree_party, newdata = df_model, type = "node")
  frame_rows <- node_to_frame[as.character(node_ids)]
  
  df_model_u <- copy(df_model)
  df_model_u$gamma <- fit_evt$frame$yval2[frame_rows, 1]
  df_model_u$sigma <- fit_evt$frame$yval2[frame_rows, 2]
  df_model_u$leaf_idx <- as.integer(rownames(fit_evt$frame)[frame_rows])
  df_model_u$source <- "arbre"
  
  # KRIGEAGE
  pts_krig <- df_model_u[!is.na(gamma) & !is.na(sigma) & !is.na(lon) & !is.na(lat)]
  pts_sp <- as.data.frame(pts_krig[, .(lon, lat, gamma, sigma)])
  coordinates(pts_sp) <- ~ lon + lat; proj4string(pts_sp) <- CRS("+proj=longlat +datum=WGS84")
  pts_sp_albers <- spTransform(pts_sp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))
  pts_sp_albers <- pts_sp_albers[!duplicated(coordinates(pts_sp_albers)), ]
  
  missing_sp <- as.data.frame(missing_counties[!is.na(lon) & !is.na(lat), .(lon, lat)])
  coordinates(missing_sp) <- ~ lon + lat; proj4string(missing_sp) <- CRS("+proj=longlat +datum=WGS84")
  missing_sp_albers <- spTransform(missing_sp, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))
  
  vgm_gamma <- variogram(gamma ~ 1, pts_sp_albers)
  fit_vgm_gamma <- tryCatch(fit.variogram(vgm_gamma, vgm(var(pts_sp_albers$gamma), "Sph", 800000, 0)), error = function(e) vgm(var(pts_sp_albers$gamma), "Sph", 800000, 0))
  if (any(fit_vgm_gamma$range < 0)) fit_vgm_gamma <- vgm(var(pts_sp_albers$gamma), "Sph", 800000, 0)
  krig_gamma <- krige(gamma ~ 1, pts_sp_albers, missing_sp_albers, model = fit_vgm_gamma, nmax = 50, debug.level = 0)
  
  vgm_sigma <- variogram(sigma ~ 1, pts_sp_albers)
  fit_vgm_sigma <- tryCatch(fit.variogram(vgm_sigma, vgm(var(pts_sp_albers$sigma), "Sph", 800000, 0)), error = function(e) vgm(var(pts_sp_albers$sigma), "Sph", 800000, 0))
  if (any(fit_vgm_sigma$range < 0)) fit_vgm_sigma <- vgm(var(pts_sp_albers$sigma), "Sph", 800000, 0)
  krig_sigma <- krige(sigma ~ 1, pts_sp_albers, missing_sp_albers, model = fit_vgm_sigma, nmax = 50, debug.level = 0)
  
  missing_counties_valid <- copy(missing_counties)[!is.na(lon) & !is.na(lat)]
  missing_counties_valid$gamma <- krig_gamma$var1.pred
  missing_counties_valid$sigma <- krig_sigma$var1.pred
  missing_counties_valid$leaf_idx <- NA_integer_
  missing_counties_valid$source <- "krigeage"
  
  # FUSION CARTOGRAPHIE
  cols_communs <- c("res_state", "res_county", "lon", "lat", "gamma", "sigma", "leaf_idx", "source")
  df_complet <- rbind(df_model_u[, ..cols_communs], missing_counties_valid[, ..cols_communs])
  
  county_sf_cont <- subset(county_sf, as.numeric(STATEFP) < 60 & !STATEFP %in% c("02", "15"))
  county_sf_cont$res_state <- fips_lookup$state[match(county_sf_cont$STATEFP, fips_lookup$state_code)]
  county_sf_cont$res_county <- clean_county_name(county_sf_cont$NAME)
  
  county_sf_final <- left_join(county_sf_cont, as.data.frame(df_complet), by = c("res_state", "res_county"))
  county_sf_final <- sf::st_transform(county_sf_final, 4326)
  
  # SAUVEGARDE
  shiny_results[[as.character(u)]] <- list(
    sf_data = county_sf_final,
    pts_extremes = sf::st_as_sf(as.data.frame(df_extremes), coords = c("lon", "lat"), crs = 4326),
    tree = tree_party,
    qq_plots = qq_final,
    rl_plots = rl_final
  )
}

# --- 4. EXPORT ---
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(shiny_results, "data/processed/precomputed_results.rds")
message("── Terminé ! Résultats sauvegardés dans data/processed/precomputed_results.rds.")