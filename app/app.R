# ==============================================================================
# APP.R - Interface RShiny Moderne avec Police Locale (.ttf)
# ==============================================================================

library(shiny)
library(bslib)
library(shinyWidgets)
library(leaflet)
library(ggplot2)
library(partykit)
library(sf)
library(RColorBrewer)
library(htmltools)

# --- Chargement instantané des données pré-calculées ---
results <- readRDS("../data/processed/precomputed_results.rds")
u_choices <- as.numeric(names(results))

# --- Configuration du Thème (Design Dashboard Bleu Roi + Ta Police Locale) ---
my_theme <- bs_theme(
  version = 5,
  bootswatch = "lumen", # Base épurée et lumineuse
  primary = "#1A4789"   # Bleu Roi profond
) |> 
  bs_add_rules("
    /* 1. Déclaration de ta police locale située dans le dossier www/ */
    @font-face {
      font-family: 'MozillaCustomText';
      src: url('MozillaText-VariableFont_wght.ttf') format('truetype'); /* <--- Change le nom si besoin */
    }

    /* 2. Application globale de ta police à toute l'application */
    body, h1, h2, h3, h4, h5, h6, .nav-link, .card-header, .control-label, .irs {
      font-family: 'MozillaCustomText', sans-serif !important;
    }

    /* Le grand bandeau du haut */
    .custom-banner {
      background: linear-gradient(135deg, #002366 0%, #1A4789 100%);
      color: white;
      padding: 30px 40px;
      display: flex;
      align-items: center;
      box-shadow: 0 4px 12px rgba(0,0,0,0.15);
      z-index: 1050;
      position: relative;
    }
    .custom-banner h2 { 
      margin: 0; 
      font-weight: 600; 
      font-size: 1.6rem; 
      letter-spacing: 0.5px; 
    }
    
    /* Design des onglets en dessous du titre */
    .nav-underline .nav-link {
      color: #6c757d;
      font-size: 16px;
      font-weight: 500;
      padding: 15px 25px;
      transition: all 0.3s ease;
    }
    .nav-underline .nav-link:hover {
      color: #1A4789;
      background-color: #f8f9fa;
      border-radius: 5px 5px 0 0;
    }
    .nav-underline .nav-link.active { 
      color: #1A4789 !important; 
      font-weight: 600; 
      border-bottom: 3px solid #1A4789 !important; 
    }
    
    /* En-tête des cartes */
    .card-header { 
      background-color: white; 
      font-weight: 600; 
      color: #1A4789; 
      border-bottom: 2px solid #f0f0f0; 
      font-size: 1.1rem;
      padding: 15px 20px;
    }
  ")

# ==============================================================================
# UI
# ==============================================================================
ui <- page_fillable(
  padding = 0, 
  theme = my_theme,
  
  # --- LE GRAND BANDEAU HEADER ---
  div(class = "custom-banner",
      tags$img(src = "logo.jpg", height = "110px", style = "margin-right: 25px;"), 
      h2("Modélisation Spatiale des Extrêmes de Mortalité COVID-19 aux États-Unis")
  ),
  
  # --- LE CORPS DE L'APPLICATION ---
  layout_sidebar(
    padding = 20,
    fillable = TRUE,
    
    # Barre latérale des filtres
    sidebar = sidebar(
      width = 320,
      bg = "#f4f6f9",
      h5("Paramètres du modèle", style = "color: #1A4789; font-weight: 600; margin-bottom: 20px;"),
      
      sliderTextInput(
        inputId = "u_val",
        label = "Seuil d'extrêmes (u) :",
        choices = sort(u_choices),
        selected = sort(u_choices)[1], # Sélectionne dynamiquement la première valeur disponible (ex: 180)
        grid = TRUE
      ),
      br(),
      helpText(
        "Explorez la segmentation spatiale du risque épidémique.",
        br(), br(),
        "L'application charge dynamiquement les modèles pré-calculés pour garantir une navigation instantanée."
      )
    ),
    
    # La zone principale avec les onglets de diagnostics
    navset_card_underline(
      nav_panel("Carte Interactive", 
                leafletOutput("map", height = "100%")
      ),
      nav_panel("Arbre CART", 
                plotOutput("tree_plot", height = "100%")
      ),
      nav_panel("QQ-Plots", 
                plotOutput("qq_plot", height = "100%")
      ),
      nav_panel("Return Levels", 
                plotOutput("rl_plot", height = "100%")
      )
    )
  )
)

# ==============================================================================
# SERVER (Optimisé avec leafletProxy pour zéro latence)
# ==============================================================================
server <- function(input, output, session) {
  
  # --- 1. Récupération réactive des données ---
  current_res <- reactive({
    req(input$u_val)
    results[[as.character(input$u_val)]]
  })
  
  # --- 2. Rendu INITIAL de la carte (S'exécute UNE SEULE FOIS) ---
  output$map <- renderLeaflet({
    # On initialise une carte vide avec juste le fond de carte Positron
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |> 
      addProviderTiles("CartoDB.Positron") |>
      setView(lng = -96, lat = 38, zoom = 4) |>
      addLayersControl(
        overlayGroups = c("Xi GPD", "Zones Interpolées"), 
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  # --- 3. Mise à jour DYNAMIQUE de la carte (Sans la recréer) ---
  observe({
    # Cette section s'exécute à chaque fois que 'u_val' change
    data_sf <- current_res()$sf_data
    
    gamma_vals <- data_sf$gamma[!is.na(data_sf$gamma)]
    pal_gamma <- colorNumeric(palette = rev(brewer.pal(11, "RdBu")), domain = gamma_vals, na.color = "#E0E0E0")
    
    tooltips <- sprintf(
      "<div style='font-size: 13px;'>
         <strong style='color:#1A4789; font-size: 15px;'>%s, %s</strong><br/>
         <hr style='margin: 5px 0;'>
         Source : <em>%s</em><br/>
         <strong>ξ (Forme)</strong> : %s<br/>
         <strong>σ (Échelle)</strong> : %s<br/>
         Feuille : %s
       </div>",
      data_sf$NAME, data_sf$res_state,
      ifelse(!is.na(data_sf$source), data_sf$source, "N/A"),
      ifelse(!is.na(data_sf$gamma), round(data_sf$gamma, 4), "N/A"),
      ifelse(!is.na(data_sf$sigma), round(data_sf$sigma, 1), "N/A"),
      ifelse(!is.na(data_sf$leaf_idx), data_sf$leaf_idx, "Interpolé")
    ) |> lapply(HTML)
    
    # On utilise leafletProxy pour modifier la carte existante
    leafletProxy("map", data = data_sf) |>
      clearShapes() |>   # Efface uniquement les anciens polygones
      clearControls() |> # Efface l'ancienne légende
      addPolygons(
        fillColor = ~pal_gamma(gamma), fillOpacity = 0.85, 
        color = "white", weight = 0.5,
        highlightOptions = highlightOptions(weight = 2, color = "#1A4789", bringToFront = TRUE),
        label = tooltips, group = "Xi GPD"
      ) |>
      addPolygons(
        data = subset(data_sf, source == "krigeage"),
        fill = FALSE, color = "#444444", weight = 1, dashArray = "3,4", opacity = 0.5,
        group = "Zones Interpolées"
      ) |>
      addLegend(position = "bottomleft", pal = pal_gamma, values = gamma_vals, 
                title = "Paramètre ξ (Forme)", opacity = 1)
  })
  
  # --- 4. Rendu de l'arbre ---
  output$tree_plot <- renderPlot({
    tree_obj <- current_res()$tree
    pts_ext  <- current_res()$pts_extremes
    
    # 1. On recrée un tableau classique et on récupère lon/lat depuis la géométrie sf
    pts_df <- as.data.frame(pts_ext)
    coords <- sf::st_coordinates(pts_ext)
    pts_df$lon <- coords[, 1]
    pts_df$lat <- coords[, 2]
    
    # 2. Création du dictionnaire de correspondance (Faux numéros -> Vrais numéros rpart)
    party_nodes <- predict(tree_obj, newdata = pts_df, type = "node")
    map_df <- unique(data.frame(
      party_id = party_nodes, 
      rpart_id = pts_df$leaf_idx, 
      gamma = pts_df$gamma, 
      sigma = pts_df$sigma
    ))
    
    # 3. Fonction pour dessiner DIRECTEMENT chaque feuille
    my_terminal_panel <- function(node) {
      info <- map_df[map_df$party_id == node$id, ]
      
      if (nrow(info) > 0) {
        rpart_id <- info$rpart_id[1]
        xi_val   <- round(info$gamma[1], 3)
        sig_val  <- round(info$sigma[1], 1)
      } else {
        rpart_id <- "?"
        xi_val   <- "?"
        sig_val  <- "?"
      }
      
      txt <- paste0("Feuille ", rpart_id, "\n\n",
                    "ξ = ", xi_val, "\n",
                    "σ = ", sig_val)
      
      grid::grid.rect(gp = grid::gpar(fill = "white", col = "#1A4789", lwd = 1.5), 
                      width = grid::unit(0.9, "npc"), height = grid::unit(0.85, "npc"))
      grid::grid.text(txt, gp = grid::gpar(fontsize = 12, fontfamily = "sans", col = "#1A4789", fontface = "bold"))
    }
    
    # 4. On trace l'arbre avec notre boîte sur mesure
    plot(tree_obj, 
         main = paste("Arbre EVT-CART pour u =", input$u_val),
         type = "extended",
         gp = grid::gpar(fontsize = 12, fontfamily = "sans"), 
         ep_args = list(justmin = 15),
         inner_panel = node_inner(tree_obj, pval = FALSE, id = FALSE), # Masque les faux numéros internes
         terminal_panel = my_terminal_panel
    )
  })
  
  # --- 5. Rendu des QQ-Plots ---
  output$qq_plot <- renderPlot({ current_res()$qq_plots })
  
  # --- 6. Rendu des Return Levels ---
  output$rl_plot <- renderPlot({ current_res()$rl_plots })
}

shinyApp(ui, server)