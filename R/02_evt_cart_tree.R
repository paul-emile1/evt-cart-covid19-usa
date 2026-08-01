# ==============================================================================
# EVT-CART + PRÉDICTION HYBRIDE (AVEC MRL PLOT ET HISTOGRAMME)
# ==============================================================================

library(data.table)
library(leaflet)
library(sf)
library(rpart)
library(partykit)
library(htmltools)
library(gstat)
library(sp)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(POT)
library(patchwork)

# ==============================================================================
# 1. UTILITAIRE : NETTOYAGE UNIFIÉ DES NOMS DE COMTÉS
# ==============================================================================

clean_county_name <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- toupper(x)
  x <- gsub(" COUNTY| PARISH| CITY| BOROUGH| MUNICIPALITY", "", x)
  x <- gsub("[^A-Z0-9 ]", "", x)
  x <- trimws(gsub("\\s+", " ", x))
  x
}

# ==============================================================================
# 2. FONCTIONS EVT-CART (AVEC LE PACKAGE POT OPTIMISÉ)
# ==============================================================================

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
  list(label = c(fit$fitted.values["shape"], fit$fitted.values["scale"]), 
       deviance = fit$deviance)
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
      
      ll_l <- -fit_l$deviance / 2
      ll_r <- -fit_r$deviance / 2
      goodness[i]  <- (ll_l + ll_r) - ll_p
    }
  }
  list(goodness = goodness, direction = direction)
}

gpd_model <- list(init = gpd_init, eval = gpd_eval, split = gpd_split)

# ==============================================================================
# 3. CHARGEMENT DEPUIS LE CSV EXPORTÉ
# ==============================================================================

message("── Chargement du dataset depuis CSV...")
# /!\ ATTENTION : Modifie ce chemin si le nouveau dataset a un nom différent
df_model <- data.table::fread("data/processed/df_final.csv")

if (!"density" %in% names(df_model) && all(c("population", "area_km2") %in% names(df_model)))
  df_model[, density := population / area_km2]

df_model[, res_county := clean_county_name(res_county)]
df_model[, res_state  := toupper(trimws(res_state))]


# ==============================================================================
# 4. CHARGEMENT DU SHAPEFILE ET RÉCUPÉRATION DES COORDONNÉES
# ==============================================================================

message("── Chargement shapefile...")
county_sf <- sf::st_read(
  "data/raw/cb_2018_us_county_500k",
  quiet = TRUE
)

data("fips_codes", package = "tidycensus")
fips_lookup <- unique(as.data.table(fips_codes)[, .(state_code, state)])

suppressWarnings({
  centroids <- sf::st_centroid(sf::st_transform(county_sf, 4326))
})
coords     <- sf::st_coordinates(centroids)
geo_coords <- data.table(
  res_state  = fips_lookup$state[match(county_sf$STATEFP, fips_lookup$state_code)],
  res_county = clean_county_name(county_sf$NAME),
  lon = coords[, 1], lat = coords[, 2]
)

# ── SUPPRESSION DES ANCIENNES COORDONNÉES SI ELLES EXISTENT ──
if ("lon" %in% names(df_model)) df_model[, lon := NULL]
if ("lat" %in% names(df_model)) df_model[, lat := NULL]

# ── Fusion des coordonnées géographiques avec df_model ──
df_model <- merge(df_model, geo_coords, by = c("res_state", "res_county"), all.x = TRUE)
df_model <- na.omit(
  df_model[!is.na(lon) & !is.na(lat)],
  cols = c("med_age", "med_income", "density", "pct_male", "pct_65p")
)

message(paste("── Comtés chargés :", nrow(df_model)))

# ==============================================================================
# 5. MRL PLOT (CHOIX DU SEUIL u)
# ==============================================================================
message("── Génération du MRL Plot pour le choix de u...")

mrl_plot <- function(y, n_thresholds = 200) {
  u_seq <- quantile(y, seq(0.05, 0.97, length.out = n_thresholds))
  u_seq <- unique(sort(u_seq))
  
  mrl   <- sapply(u_seq, function(u) {
    exc <- y[y > u] - u
    if (length(exc) < 5) NA else mean(exc)
  })
  
  se <- sapply(u_seq, function(u) {
    exc <- y[y > u] - u
    if (length(exc) < 5) NA else sd(exc) / sqrt(length(exc))
  })
  
  df_m <- data.frame(u = u_seq, mrl = mrl,
                     lower = mrl - 1.96*se,
                     upper = mrl + 1.96*se) |> filter(!is.na(mrl))
  
  q_lab <- quantile(y, c(0.80, 0.85, 0.90, 0.95))
  df_ann <- data.frame(
    u   = q_lab,
    mrl = sapply(q_lab, function(u) {
      exc <- y[y > u] - u
      if (length(exc) < 5) NA else mean(exc)
    }),
    lab = paste0("q", c(80,85,90,95), "%\nn=", sapply(q_lab, function(u) sum(y > u)))
  ) |> filter(!is.na(mrl))
  
  ggplot(df_m, aes(x = u)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2196F3", alpha = 0.15) +
    geom_line(aes(y = mrl), color = "#1565C0", linewidth = 0.9) +
    geom_vline(xintercept = q_lab, linetype = "dashed", color = "#E53935", alpha = 0.4) +
    geom_point(data = df_ann, aes(x = u, y = mrl), color = "#E53935", size = 3, shape = 18) +
    geom_text(data = df_ann, aes(x = u, y = mrl, label = lab), vjust = -0.8, size = 3, color = "#E53935") +
    labs(title    = "Mean Residual Life (MRL) Plot",
         subtitle = "Recherchez le seuil à partir duquel la courbe devient approximativement linéaire",
         x = "Seuil u candidates",
         y = "Espérance de l'excès moyen") +
    theme_minimal()
}

# 1. On affiche le plot
print(mrl_plot(df_model$Y))

# Définition du seuil
u <- 210

message(paste("── Seuil u sélectionné :", round(u, 2)))

# ==============================================================================
# 5.5 HISTOGRAMME DE DISTRIBUTION (MISE EN ÉVIDENCE DES EXTRÊMES)
# ==============================================================================
message("── Génération de l'histogramme de distribution...")

p_hist <- ggplot(df_model, aes(x = Y, fill = Y > u)) +
  geom_histogram(bins = 60, color = "white", alpha = 0.85) +
  scale_fill_manual(
    name = "Classification (Approche PoT)",
    values = c("FALSE" = "#3498db", "TRUE" = "#e74c3c"),
    labels = c(paste("Comportement normal (≤", round(u,1), ")"), 
               paste("Extrêmes analysés (>", round(u,1), ")"))
  ) +
  geom_vline(xintercept = u, linetype = "dashed", color = "black", linewidth = 1) +
  annotate(
    "text", x = u, y = Inf, label = paste("Seuil u =", round(u,1)),
    vjust = 2, hjust = -0.1, fontface = "bold"
  ) +
  labs(
    title = "Distribution de la mortalité par comté (Y)",
    subtitle = "Séparation entre le corps de la distribution et la queue (valeurs extrêmes)",
    x = "Taux de mortalité (Décès)",
    y = "Nombre de comtés"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_hist)

# ==============================================================================
# 6. ARBRE EVT-CART (entraîné sur les extrêmes)
# ==============================================================================

df_extremes <- df_model[Y > u]
df_extremes[, Y_excess := Y - u]
message(paste("── Comtés extrêmes :", nrow(df_extremes)))

fit_evt <- rpart(
  Y_excess ~ med_age + med_income + density + pct_male + pct_65p + lon + lat,
  data    = df_extremes,
  method  = gpd_model,
  control = rpart.control(minsplit = 20, minbucket = 10, maxdepth = 5, cp = 0.001) 
)

message(paste("── Feuilles de l'arbre :", sum(fit_evt$frame$var == "<leaf>")))

if (nrow(fit_evt$frame) > 1) {
  par(xpd = TRUE)
  plot(fit_evt, margin = 0.1, main = paste("Arbre EVT-CART (u =", round(u, 1), ")"))
  text(fit_evt, cex = 0.8)
}

# Assignation directe via $where
df_extremes$gamma    <- fit_evt$frame$yval2[fit_evt$where, 1]
df_extremes$sigma    <- fit_evt$frame$yval2[fit_evt$where, 2]
df_extremes$leaf_idx <- as.integer(rownames(fit_evt$frame)[fit_evt$where])

leaf_params <- unique(df_extremes[, .(leaf_idx, gamma, sigma)])[order(leaf_idx)]
message("Paramètres GPD par feuille :")
print(leaf_params)

# ==============================================================================
# 6.5 DIAGNOSTICS PAR FEUILLE (QQ-PLOTS ET RETURN LEVELS)
# ==============================================================================
message("── Génération des QQ-Plots et Return Levels...")

leaf_ids <- sort(unique(df_extremes$leaf_idx))
p_u <- nrow(df_extremes) / nrow(df_model)

qq_list <- list()
rl_list <- list()

for (lid in leaf_ids) {
  obs <- sort(df_extremes$Y_excess[df_extremes$leaf_idx == lid])
  n_obs <- length(obs)
  
  xi_h <- df_extremes$gamma[df_extremes$leaf_idx == lid][1]
  sig_h <- df_extremes$sigma[df_extremes$leaf_idx == lid][1]
  
  # 1. QQ-PLOT
  p_emp <- (1:n_obs - 0.5) / n_obs
  if (abs(xi_h) < 1e-6) {
    q_theo <- -sig_h * log(1 - p_emp)
  } else {
    q_theo <- (sig_h / xi_h) * ((1 - p_emp)^(-xi_h) - 1)
  }
  
  df_qq <- data.frame(emp = obs, theo = q_theo)
  p_qq <- ggplot(df_qq, aes(x = theo, y = emp)) +
    geom_abline(slope = 1, intercept = 0, color = "#e74c3c", linetype = "dashed") +
    geom_point(color = "#2c3e50", alpha = 0.7) +
    labs(title = paste("QQ-Plot - Feuille", lid), 
         subtitle = paste("n =", n_obs, "| xi =", round(xi_h, 3)),
         x = "Quantiles théoriques", y = "Quantiles empiriques") +
    theme_minimal()
  qq_list[[as.character(lid)]] <- p_qq
  
  # 2. RETURN LEVELS
  T_vals <- c(2, 5, 10, 20, 50, 100, 200)
  valid_T <- T_vals[T_vals * p_u > 1]
  
  if(length(valid_T) > 0) {
    p_vals <- 1 - 1 / (valid_T * p_u)
    if (abs(xi_h) < 1e-6) {
      rl_vals <- u - sig_h * log(1 - p_vals)
    } else {
      rl_vals <- u + (sig_h / xi_h) * ((1 - p_vals)^(-xi_h) - 1)
    }
    
    df_rl <- data.frame(T = valid_T, RL = rl_vals)
    p_rl <- ggplot(df_rl, aes(x = T, y = RL)) +
      geom_line(color = "#e74c3c") + geom_point(color = "#e74c3c") +
      scale_x_log10() +
      labs(title = paste("Return Level - Feuille", lid), 
           x = "Période de retour T (comtés)", y = "Niveau Y") +
      theme_minimal()
    rl_list[[as.character(lid)]] <- p_rl
  }
}

if (length(qq_list) > 0) print(wrap_plots(qq_list, ncol = min(length(qq_list), 3)))
if (length(rl_list) > 0) print(wrap_plots(rl_list, ncol = min(length(rl_list), 3)))

# ==============================================================================
# 7. PRÉDICTION ARBRE SUR TOUS LES COMTÉS DU DATASET
# ==============================================================================

message("── Prédiction arbre sur tous les comtés du dataset...")

tree_party <- as.party(fit_evt)
nodes_train   <- predict(tree_party, newdata = df_extremes, type = "node")
frame_rows_train <- fit_evt$where 

mapping_dt <- unique(data.table(
  node_party  = as.integer(nodes_train),
  frame_row   = as.integer(frame_rows_train)
))
node_to_frame <- setNames(mapping_dt$frame_row, mapping_dt$node_party)

node_ids <- predict(tree_party, newdata = df_model, type = "node")
frame_rows <- node_to_frame[as.character(node_ids)]

df_model$gamma    <- fit_evt$frame$yval2[frame_rows, 1]
df_model$sigma    <- fit_evt$frame$yval2[frame_rows, 2]
df_model$leaf_idx <- as.integer(rownames(fit_evt$frame)[frame_rows])
df_model$source   <- "arbre"

message(paste("── Comtés classifiés par l'arbre :", sum(!is.na(df_model$gamma))))

# ==============================================================================
# 8. KRIGEAGE POUR LES COMTÉS MANQUANTS
# ==============================================================================

message("── Identification des comtés manquants...")

all_counties <- geo_coords[
  res_state %in% fips_lookup$state &
    !res_state %in% c("AK", "HI") 
]

in_dataset <- paste(df_model$res_state, df_model$res_county)
all_counties[, key_ := paste(res_state, res_county)]
missing_counties <- all_counties[!key_ %in% in_dataset]
missing_counties[, key_ := NULL]

message(paste("── Comtés manquants (à interpoler) :", nrow(missing_counties)))

pts_krig <- df_model[!is.na(gamma) & !is.na(sigma) & !is.na(lon) & !is.na(lat)]
pts_sp <- as.data.frame(pts_krig[, .(lon, lat, gamma, sigma)])
coordinates(pts_sp) <- ~ lon + lat
proj4string(pts_sp) <- CRS("+proj=longlat +datum=WGS84")

albers_crs      <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"
pts_sp_albers   <- spTransform(pts_sp, CRS(albers_crs))

missing_sp <- as.data.frame(missing_counties[!is.na(lon) & !is.na(lat), .(lon, lat)])
coordinates(missing_sp) <- ~ lon + lat
proj4string(missing_sp) <- CRS("+proj=longlat +datum=WGS84")
missing_sp_albers       <- spTransform(missing_sp, CRS(albers_crs))

dup_pts <- duplicated(coordinates(pts_sp_albers))
if (any(dup_pts)) {
  message(paste("── Doublons supprimés :", sum(dup_pts)))
  pts_sp_albers <- pts_sp_albers[!dup_pts, ]
}

# Krigeage Gamma
vgm_gamma <- variogram(gamma ~ 1, pts_sp_albers)
fit_vgm_gamma <- tryCatch(
  fit.variogram(vgm_gamma, model = vgm(psill = var(pts_sp_albers$gamma), model = "Sph", range = 800000, nugget = 0)),
  error = function(e) vgm(var(pts_sp_albers$gamma), "Sph", 800000, 0)
)
if (any(fit_vgm_gamma$range < 0)) fit_vgm_gamma <- vgm(var(pts_sp_albers$gamma), "Sph", 800000, 0)
krig_gamma <- krige(gamma ~ 1, pts_sp_albers, missing_sp_albers, model = fit_vgm_gamma, nmax = 50)

# Krigeage Sigma
vgm_sigma <- variogram(sigma ~ 1, pts_sp_albers)
fit_vgm_sigma <- tryCatch(
  fit.variogram(vgm_sigma, model = vgm(psill = var(pts_sp_albers$sigma), model = "Sph", range = 800000, nugget = 0)),
  error = function(e) vgm(var(pts_sp_albers$sigma), "Sph", 800000, 0)
)
if (any(fit_vgm_sigma$range < 0)) fit_vgm_sigma <- vgm(var(pts_sp_albers$sigma), "Sph", 800000, 0)
krig_sigma <- krige(sigma ~ 1, pts_sp_albers, missing_sp_albers, model = fit_vgm_sigma, nmax = 50)


missing_counties_valid <- missing_counties[!is.na(lon) & !is.na(lat)]
missing_counties_valid$gamma    <- krig_gamma$var1.pred
missing_counties_valid$sigma    <- krig_sigma$var1.pred
missing_counties_valid$leaf_idx <- NA_integer_
missing_counties_valid$source   <- "krigeage"

message(paste("── Comtés interpolés par krigeage :", nrow(missing_counties_valid)))

# ==============================================================================
# 9. TABLE COMPLÈTE : ARBRE + KRIGEAGE
# ==============================================================================

cols_communs <- c("res_state", "res_county", "lon", "lat", "gamma", "sigma",
                  "leaf_idx", "source")

df_arbre   <- df_model[, ..cols_communs]
df_krig    <- missing_counties_valid[, ..cols_communs]
df_complet <- rbind(df_arbre, df_krig)

message(paste("── Total comtés sur la carte :", nrow(df_complet)))

# ==============================================================================
# 10. JOINTURE AVEC LE SHAPEFILE POUR LA CARTE
# ==============================================================================

county_sf_cont  <- subset(county_sf,
                          as.numeric(STATEFP) < 60 &
                            !STATEFP %in% c("02", "15"))
county_sf_cont$res_state  <- fips_lookup$state[match(county_sf_cont$STATEFP, fips_lookup$state_code)]
county_sf_cont$res_county <- clean_county_name(county_sf_cont$NAME)

county_sf_final <- left_join(
  county_sf_cont,
  as.data.frame(df_complet),
  by = c("res_state", "res_county")
)

county_sf_final <- sf::st_transform(county_sf_final, 4326)

n_total      <- nrow(county_sf_final)
n_classified <- sum(!is.na(county_sf_final$gamma))

# ==============================================================================
# 11. CARTOGRAPHIE LEAFLET
# ==============================================================================

gamma_vals <- county_sf_final$gamma[!is.na(county_sf_final$gamma)]
pal_gamma  <- colorNumeric(
  palette  = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
  domain   = gamma_vals,
  na.color = "#cccccc"
)

tooltips <- sprintf(
  paste0(
    "<strong>%s, %s</strong><br/>",
    "Source : <em>%s</em><br/>",
    "Gamma (forme) : %s<br/>",
    "Sigma (échelle) : %s<br/>",
    "Feuille EVT : %s"
  ),
  county_sf_final$NAME,
  county_sf_final$res_state,
  ifelse(!is.na(county_sf_final$source), county_sf_final$source, "N/A"),
  ifelse(!is.na(county_sf_final$gamma),  round(county_sf_final$gamma, 4), "N/A"),
  ifelse(!is.na(county_sf_final$sigma),  round(county_sf_final$sigma, 4), "N/A"),
  ifelse(!is.na(county_sf_final$leaf_idx), county_sf_final$leaf_idx, "interpolé")
) |> lapply(htmltools::HTML)

pts_extremes_sf <- sf::st_as_sf(
  as.data.frame(df_extremes),
  coords = c("lon", "lat"), crs = 4326
)
tooltips_ext <- sprintf(
  "<strong>%s, %s</strong><br/>Y = %s décès/100k<br/>Feuille : %s | Gamma : %s | Sigma : %s",
  df_extremes$res_county, df_extremes$res_state,
  round(df_extremes$Y, 0),
  df_extremes$leaf_idx,
  round(df_extremes$gamma, 4),
  round(df_extremes$sigma, 4)
) |> lapply(htmltools::HTML)

ma_carte <- leaflet(options = leafletOptions(zoomControl = TRUE)) |>
  addProviderTiles("CartoDB.Positron") |>
  
  addPolygons(
    data             = county_sf_final,
    fillColor        = ~pal_gamma(gamma),
    fillOpacity      = 0.8,
    color            = "#444444",
    weight           = 0.3,
    highlightOptions = highlightOptions(weight = 2, color = "black", bringToFront = TRUE),
    label            = tooltips,
    labelOptions     = labelOptions(direction = "auto"),
    group            = "Xi GPD (arbre + krigeage)"
  ) |>
  
  addPolygons(
    data        = subset(county_sf_final, source == "krigeage"),
    fill        = FALSE,
    color       = "#888888",
    weight      = 0.8,
    dashArray   = "3,4",
    opacity     = 0.6,
    group       = "Comtés interpolés (krigeage)"
  ) |>
  
  addCircleMarkers(
    data         = pts_extremes_sf,
    radius       = 4,
    color        = "#1a1a2e",
    fillColor    = "#e94560",
    fillOpacity  = 0.9,
    weight       = 1,
    label        = tooltips_ext,
    labelOptions = labelOptions(direction = "auto"),
    group        = "Comtés extrêmes (entraînement)"
  ) |>
  
  addLegend(
    position = "bottomleft",
    pal      = pal_gamma,
    values   = gamma_vals,
    title    = "Xi GPD",
    opacity  = 0.85
  ) |>
  
  addLayersControl(
    overlayGroups = c(
      "Xi GPD (arbre + krigeage)",
      "Comtés interpolés (krigeage)",
      "Comtés extrêmes (entraînement)"
    ),
    options = layersControlOptions(collapsed = FALSE)
  )

ma_carte

# ==============================================================================
# 12. VARIOGRAMMES (diagnostic)
# ==============================================================================

par(mfrow = c(1, 2))
plot(vgm_gamma, fit_vgm_gamma, main = "Variogramme - Gamma")
plot(vgm_sigma, fit_vgm_sigma, main = "Variogramme - Sigma")
par(mfrow = c(1, 1))

# ==============================================================================
# 13. RÉSUMÉ FINAL
# ==============================================================================

message("\n========== RÉSUMÉ ==========")
message(paste("Seuil u                       :", u))
message(paste("Comtés dans le dataset        :", nrow(df_model)))
message(paste("Comtés extrêmes (train arbre) :", nrow(df_extremes)))
message(paste("Feuilles de l'arbre           :", sum(fit_evt$frame$var == "<leaf>")))
message(paste("Comtés classifiés par arbre   :", sum(df_complet$source == "arbre")))
message(paste("Comtés interpolés par krigeage:", sum(df_complet$source == "krigeage")))
message(paste("Couverture carte              :",
              round(100 * n_classified / n_total, 1), "%"))
message("Paramètres GPD par feuille :")
print(leaf_params)
message("=============================\n")