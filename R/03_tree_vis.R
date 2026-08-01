# ==============================================================================
# SCRIPT INDÉPENDANT : GÉNÉRATION DE L'ARBRE COLORÉ (u = 210)
# ==============================================================================

library(data.table)
library(sf)
library(rpart)
library(rpart.plot)
library(POT)

# --- 1. FONCTIONS GPD ---
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
    # C'est ici la correction pour rpart.plot :
    text = function(yval, dev, wt, ylevel, digits, n, use.n)
      paste0("xi=", round(yval[1], 3), "\nn=", n) 
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


# --- 2. CHARGEMENT ET PRÉPARATION (Inclus la géométrie) ---
# Correction du nom du fichier
df_model <- data.table::fread("data/processed/df_final.csv")

if (!"density" %in% names(df_model) && all(c("population", "area_km2") %in% names(df_model))) {
  df_model[, density := population / area_km2]
}

df_model[, res_county := clean_county_name(res_county)]
df_model[, res_state  := toupper(trimws(res_state))]

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

# ── SUPPRESSION DES ANCIENNES COORDONNÉES AVANT FUSION ──
if ("lon" %in% names(df_model)) df_model[, lon := NULL]
if ("lat" %in% names(df_model)) df_model[, lat := NULL]

df_model <- merge(df_model, geo_coords, by = c("res_state", "res_county"), all.x = TRUE)
df_model <- na.omit(df_model[!is.na(lon) & !is.na(lat)], cols = c("med_age", "med_income", "density", "pct_male", "pct_65p"))


# --- 3. ENTRAÎNEMENT POUR U = 210 ---
u <- 210  # Nouveau seuil
df_extremes <- copy(df_model[Y > u])
df_extremes[, Y_excess := Y - u]

fit_evt <- rpart(Y_excess ~ med_age + med_income + density + pct_male + pct_65p + lon + lat,
                 data = df_extremes, method = gpd_model,
                 control = rpart.control(minsplit = 20, minbucket = 10, maxdepth = 5, cp = 0.001))

# --- 4. LE HACK DÉFINITIF POUR LE DESIGN CARRÉ, COLORÉ & RÈGLES DE SPLIT ---
fit_evt_plot <- fit_evt 
fit_evt_plot$method <- "anova" 
class(fit_evt_plot) <- "rpart"

# On force la valeur principale pour que le dégradé de couleur s'applique sur xi (feuilles)
fit_evt_plot$frame$yval <- fit_evt_plot$frame$yval2[, 1]

# Avec type = 0, cette fonction ne s'appliquera qu'aux FEUILLES terminales
node_fun <- function(x, labs, digits, varlen) {
  # Extraction de Xi (colonne 1) et Sigma (colonne 2)
  xi_arrondi <- round(x$frame$yval2[, 1], 3)
  sig_arrondi <- round(x$frame$yval2[, 2], 1)
  n_vals <- x$frame$n
  
  # Format pour les feuilles avec les vraies lettres grecques
  paste0("ξ = ", xi_arrondi, "\nσ = ", sig_arrondi, "\nn = ", n_vals)
}

# --- 5. AFFICHAGE AVEC LE DESIGN EXACT (100% Carré & Bleu) ---
par(mar = c(1, 1, 3, 1))

rpart.plot(fit_evt_plot, 
           type = 0,                   # La règle de split entière DANS le noeud interne
           extra = 0,                  # Désactive tout texte parasite
           nn = FALSE,                 # Retire les petits numéros au-dessus
           fallen.leaves = FALSE,  
           branch.lty = 1,             # Lignes pleines
           box.palette = "Blues",      # Colore les FEUILLES avec un dégradé bleu (selon le risque)
           split.box.col = "#C6DBEF",  # COLORIE LES NOEUDS INTERNES avec un beau bleu clair
           round = 0,                  # Paramètre global pour les angles droits
           leaf.round = 0,             # FORCE les feuilles à être strictement CARRÉES
           split.round = 0,            # FORCE les noeuds internes à être strictement CARRÉS
           shadow.col = "gray",        # Ombre légère pour le volume
           node.fun = node_fun,        # Applique notre texte personnalisé aux feuilles
           cex = 0.78,              
           main = paste("Arbre EVT-CART Spatio-Démographique (u =", u, ")") # Titre dynamique
)