#NOTA: El codigo presenta problema en la reescritura de los archivos, 
#por lo que es necesario borrar la region de GRASS cada vez que se ejecute 
#el codigo completo o ejecutarlo por partes.
ejecutar_leaflet <- F
reproducir_desde_cero <- T
if(reproducir_desde_cero) {
  system('rm -r grass-data-test')
}
#Paquetes
library(rgrass7)
library(sp)
library(sf)
library(raster)

gisdbase <- 'grass-data-test' #Base de datos de GRASS GIS
wd <- getwd() #Directorio de trabajo
wd
loc <- initGRASS(gisBase = "/usr/lib/grass78/",
                 home = wd,
                 gisDbase = paste(wd, gisdbase, sep = '/'),
                 location = 'brujuelas',
                 mapset = "PERMANENT",
                 override = TRUE)

#Muestra la definición de la región
gmeta()
#Definir ruta del DEM
dem <- 'dem-cuencas-brujuelas-cumayasa.tif'

#Definir la proyección de la región basada en DEM
execGRASS(
  cmd = 'g.proj',
  flags = c('t','c'),
  georef=dem)

#Importar mapa raster
#r.in.gdal importa la fuente a GRASS
execGRASS(
  cmd = 'r.in.gdal',
  flags=c('overwrite','quiet'),
  parameters=list(
    input=dem,
    output='dem'
  )
)
#Actualizar la extensión de la región al DEM, sólo por precaución
execGRASS(
  cmd = 'g.region',
  parameters=list(
    raster = 'dem',
    align = 'dem'
  )
)
par(mfrow=c(1,1)) #para que se muestre un solo plot a la vez

#Tallado del DEM
execGRASS(
  'v.in.ogr',
  flags = 'overwrite',
  parameters = list(
    input = 'red_mtn50k_cleaned_largos.gpkg',
    output = 'red_mtn50k_cleaned_largos'
  )
)

execGRASS(
  'r.univar',
  flags = 't',
  parameters = list(
    map = 'dem'  
  )
)

execGRASS(
  'v.to.rast',
  flags = 'overwrite',
  parameters = list(
    input = 'red_mtn50k_cleaned_largos',
    output = 'red_mtn50k_cleaned_largos',
    use = 'val',
    value = 1 
  )
)

execGRASS( 
  'r.null',
  parameters = list( 
    map = 'red_mtn50k_cleaned_largos',
    null = 0   
  )
)

#r.null map=red_mtn50k_cleaned_largos null=0

execGRASS(
  'r.mapcalc',
  flags = 'overwrite',
  parameters = list( 
    expression = 'stddem = (dem - -7.10952377319336) / (517.373901367188 - -7.10952377319336)' 
  )
)

execGRASS(
  'r.mapcalc',
  flags = 'overwrite',
  parameters = list( 
    expression = 'stddemburn = stddem - red_mtn50k_cleaned_largos'
  )
)

use_sp()
stddemburn <- readRAST('stddemburn')
plot(stddemburn)

execGRASS(
  'r.mapcalc',
  flags = 'overwrite',
  parameters = list( 
    expression = 'dem_tallado = (stddemburn * (517.373901367188 - -7.10952377319336)) -7.10952377319336'
  )
)
dem_tallado_para_ver <- readRAST('dem_tallado')
plot(dem_tallado_para_ver)

#Declarando depresiones.
execGRASS(
  'v.in.ogr',
  flags = 'overwrite',
  parameters = list(
    input = 'depresiones_para_editar.gpkg',
    output = 'depresiones_para_editar'
  )
)

execGRASS(
  'v.to.rast',
  flags = 'overwrite',
  parameters = list(
    input = 'depresiones_para_editar',
    output = 'depresiones_todas',
    use = 'val',
    value = 1
  )
)
#Hidrologia computacional
# Calcular parámetros hidrográficos de interés usando `r.watershed`
execGRASS(
  "r.watershed",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = "dem_tallado",
    accumulation = "accum-de-rwshed",
    depression = 'depresiones_todas',
    stream = "stream-de-rwshed",
    drainage = "drainage-dir-de-rwshed",
    basin = 'basins',
    half_basin = 'half-basins',
    threshold = 80 #umbral de acumulacion
  )
)

execGRASS(
  'r.water.outlet',
  flags = 'overwrite',
  parameters = list(
    input = 'drainage-dir-de-rwshed',
    output = 'cuenca_brujuelas',
    coordinates = c(436668.43711353285,	2046734.011679014)
  )
)

execGRASS(
  "r.to.vect",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'cuenca_brujuelas',
    output = 'cuenca_brujuelas',
    type = 'area'
  )
)
plot(readVECT('cuenca_brujuelas'))

execGRASS(
  'v.out.ogr',
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'cuenca_brujuelas',
    output = 'cuenca_brujuelas',
    format = 'GeoJSON'
  )
)

#Importar un mapa vectorial también
demext <- 'cuenca_brujuelas.geojson'   #Es necesario para las variables de terreno
execGRASS(
  cmd = 'v.in.ogr',
  flags=c('overwrite','quiet'),
  parameters=list(
    input=demext,
    output='dem_extent'
  )
)
#plot(st_read(demext))
#Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

#Cargar en R el DEM (mapa ráster)
use_sp()
dem_sp <- readRAST('dem')
op <- par()
plot(dem_sp)

#Cargar a R el mapa vectorial de una cuenca que se encuentra alojado fuera de GRASS, hacer el plot y representar la cuenca del rio brujuelas superpuesta
rutabrujuelas <- demext
brujuelas <- st_read(rutabrujuelas)
plot(dem_sp)
plot(brujuelas, add=T, col='transparent', border='black', lwd=5);par(op[c('mfrow','mar')])

#Analizar el DEM dentro de la cuenca del rio brujuelas
dem_r0 <- raster(dem_sp)
dem_r1 <- crop(dem_r0, brujuelas)
dem_brujuelas <- mask(dem_r1, brujuelas)
plot(dem_brujuelas)

summary(dem_brujuelas)
hist(dem_brujuelas)

#Obtener variables de terreno básicas con el paquete raster dentro de R
pend_brujuelas <- terrain(x = dem_brujuelas, opt = 'slope', unit = 'degrees')
plot(pend_brujuelas)

# Traer capas a R
#Usar Spatial*
library(sp)
use_sp()
#Paquete manejo de los raster
library(raster)
#DEM
dem_para_leaflet <- raster(readRAST('dem'))
#Basins
basins <- raster(readRAST('basins'))
#Stream network
stream <- raster(readRAST('stream-de-rwshed'))
stream3857 <- projectRaster(stream, crs = CRS("+init=epsg:3857"), method = 'ngb')
#Generar un vectorial de extensión de capa en EPSG:4326
e <- extent(stream)
e <- as(e, 'SpatialPolygons')
proj4string(e) <- CRS("+init=epsg:32619")
e <- spTransform(e, CRSobj = CRS("+init=epsg:4326"))

# Visualizar capas con `leaflet`
library(leaflet)
library(leafem)
if(ejecutar_leaflet) {
  leaflet() %>%
    addProviderTiles(providers$OpenStreetMap, group = 'terrain') %>%
    addRasterImage(dem_para_leaflet, group='DEM', opacity = 0.5) %>%
    addRasterImage(
      ratify(basins),
      group='basins', opacity = 0.7,
      colors = sample(rep(RColorBrewer::brewer.pal(12, 'Set3'),1000))) %>% 
    addRasterImage(stream3857, project = F, group='str', opacity = 0.7, method = 'ngb', colors = 'blue') %>% 
    addLayersControl(
      overlayGroups = c('terrain','DEM','basins','str'),
      options = layersControlOptions(collapsed=FALSE)) %>% 
    addHomeButton(extent(e), 'Ver todo')
}

# Obtener las coordenadas de la desembocadura de la cuenca de interés
library(mapview)
mapa_red_coord_desemb <- mapview(
  stream3857, method='ngb', col.regions = 'blue',
  legend = FALSE, label = FALSE,
  maxpixels =  29743620 #ANTES: maxpixels =  910425
)
mapa_red_coord_desemb

#Convertir las coordenadas lat/lon a EPSG:32619
my_trans <- function(coords = NULL) {
  require(sp)
  pt <- SpatialPoints(matrix(coords, ncol = 2), CRS("+init=epsg:4326"))
  foo <- spTransform(pt, CRSobj = CRS("+init=epsg:32619"))
  bar <- as.vector(coordinates(foo))
  return(bar)
}
brujuelas_out <- my_trans(coords = c(-69.599987, 18.510156)) # Originalmente, c(436668.43711353285,	2046734.011679014)
brujuelas_out  #436668.4 2046734.0

## Extraer la cuenca de interés
execGRASS(
  "r.water.outlet",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'drainage-dir-de-rwshed',
    output = 'cuenca_brujuelas',
    coordinates = brujuelas_out
  )
)

## Convertir la cuenca a vectorial en GRASS
execGRASS(
  "r.to.vect",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'cuenca_brujuelas',
    output = 'cuenca_brujuelas',
    type = 'area'
  )
)

execGRASS(
  "v.buffer",
  flags = 'overwrite',
  parameters = list(
    input = 'cuenca_brujuelas',
    output = 'cuenca_brujuelas_buffer',
    distance = 100,
    type = 'area'
  )
)

execGRASS(
  "r.mask",
  flags = c('verbose','overwrite','quiet'),
  parameters = list(
    vector = 'cuenca_brujuelas_buffer' #NOTA: USANDO BUFFER!!!!!!!!!!!!!
  )
)

## Mostrar lista nuevamente
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Traer a R la cuenca
brujuelas_bas <- readVECT('cuenca_brujuelas')
brujuelas_bas
plot(brujuelas_bas)
brujuelas_bas4326 <- spTransform(brujuelas_bas, CRSobj = CRS("+init=epsg:4326"))
if(ejecutar_leaflet) {
  leaflet() %>% 
    addProviderTiles(providers$OpenStreetMap) %>%
    addPolygons(data = brujuelas_bas4326) %>% 
    leafem::addHomeButton(extent(brujuelas_bas4326), 'Ver cuenca')
}

#Usar la cuenca del rio brujuelas como máscara
# Extraer la red de drenaje de la cuenca de interés
execGRASS(
  'r.stream.extract',
  flags = 'overwrite',
  parameters = list(
    elevation= 'dem_tallado',
    accumulation= 'accum-de-rwshed',
    depression= 'depresiones_todas',
    threshold= 160,
    stream_raster= 'rstrm_brujuelas',
    stream_vector= 'rstrm_brujuelas',
    direction= 'rstrm_direccion'
  )
)

#Crear mapas de órdenes de red
execGRASS(
  "r.stream.order",
  flags = c('overwrite','quiet'),
  parameters = list(
    stream_rast = 'rstrm_brujuelas',
    direction = 'rstrm_direccion',
    elevation = 'dem_tallado',
    accumulation = 'accum-de-rwshed',
    stream_vect = 'order_all',
    strahler = 'order-strahler',
    horton = 'order-horton',
    shreve = 'order-shreve',
    hack = 'order-hack-gravelius',
    topo = 'order-topology'
  )
)
# Mostrar lista nuevamente
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
# Visualizar la red con leaflet
#Simbología única
order <- readVECT('order_all')
order4326 <- spTransform(order, CRSobj = CRS("+init=epsg:4326"))
if(ejecutar_leaflet) {
  leaflet() %>% 
    addProviderTiles(providers$OpenStreetMap, group = 'terrain') %>%
    addPolylines(
      data = order4326, weight = 3, opacity = 0.7, group = 'order',
      label = ~as.character(strahler),
      highlightOptions = highlightOptions(color = "white",
                                          weight = 5, bringToFront = F, opacity = 1),
      labelOptions = labelOptions(noHide = T,
                                  style = list(
                                    "font-size" = "8px",
                                    "background" = "rgba(255, 255, 255, 0.5)",
                                    "background-clip" = "padding-box",
                                    "padding" = "1px"))) %>% 
    leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
    addLayersControl(
      overlayGroups = c('terrain','order'),
      options = layersControlOptions(collapsed=FALSE))
}

order4326
plot(order4326)

#Simbología aplicando grosor según orden de red
if(ejecutar_leaflet) {
  leaflet() %>% 
    addProviderTiles(providers$OpenStreetMap, group = 'terrain') %>%
    addPolylines(
      data = order4326, weight = order4326$strahler*1.1, opacity = 0.7, group = 'order',
      label = ~as.character(strahler),
      highlightOptions = highlightOptions(color = "white",
                                          weight = 5, bringToFront = F, opacity = 1),
      labelOptions = labelOptions(noHide = F)) %>% 
    leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
    addLayersControl(
      overlayGroups = c('terrain','order'),
      options = layersControlOptions(collapsed=FALSE))
}

#Delimitar cuencas según orden de red de Strahler
#Obtener órdenes de red mínimo y máximo
#Estadísticas para obtener los valores mínimo y máximo del orden de red de Strahler
rinfo.ordstra <- execGRASS(
  'r.info',
  flags = 'r',
  parameters = list(
    map = 'order-strahler'
  )
)
#Órdenes de red mínimo y máximo
minmaxord <- as.numeric(
  stringr::str_extract_all(
    attributes(rinfo.ordstra)$resOut,
    "[0-9]+"
  )
)
minmaxord

### Delimitar cuencas, convertirlas de ráster a vectorial
sapply(
  min(minmaxord):max(minmaxord),
  function(x){
    execGRASS(
      "r.stream.basins",
      flags = c('overwrite','c','quiet'),
      parameters = list(
        direction = 'rstrm_direccion',
        stream_rast = 'order-strahler',
        cats = as.character(x),
        basins = paste0('cuenca_brujuelas',x)
      )
    )
    execGRASS(
      "r.to.vect",
      flags=c('overwrite','quiet'),
      parameters = list(
        input = paste0('cuenca_brujuelas',x),
        output = paste0('cuenca_brujuelas',x),
        type = 'area'
      )
    )
  }
)

#Representar las cuencas con leaflet
if(ejecutar_leaflet) {
  sapply(
    min(minmaxord):max(minmaxord),
    function(x){
      assign(
        paste0('orden', x),
        spTransform(readVECT(paste0('cuenca_brujuelas',x)), CRSobj = CRS("+init=epsg:4326")),
        envir = .GlobalEnv)
    }
  )
  paleta <- RColorBrewer::brewer.pal(12, 'Set3')
  leaflet() %>% 
    addProviderTiles(providers$OpenStreetMap, group = 'terrain') %>%
    addPolygons(data = orden4, stroke = T, weight = 2,
                color = ~paleta, fillOpacity = 0.4, group = 'O4') %>% 
    addPolygons(data = orden3, stroke = T, weight = 2,
                color = ~paleta, fillOpacity = 0.4, group = 'O3') %>%
    addPolygons(data = orden2, stroke = T, weight = 2,
                color = ~paleta, fillOpacity = 0.4, group = 'O2') %>%
    addPolygons(data = orden1, stroke = T, weight = 2,
                color = ~paleta, fillOpacity = 0.4, group = 'O1') %>%
    addPolylines(
      data = order4326, weight = order4326$strahler*1.1,
      opacity = 0.7, group = 'str_order') %>%
    leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
    addLayersControl(
      overlayGroups = c('terrain','O1','O2','O3','O4','str_order'),
      options = layersControlOptions(collapsed=FALSE))
}


#Estadísticas de red resumidas por orden de red.
execGRASS(
  "r.stream.stats",
  flags = c('overwrite','quiet','o'),
  parameters = list(
    stream_rast = 'order-strahler',
    direction = 'rstrm_direccion',
    elevation = 'dem_tallado',
    output = 'brujuelas_stats.txt'
  )
)
file.show('brujuelas_stats.txt')
d <- read.csv("brujuelas_stats.txt", skip=1, header=TRUE)
plot(num_of_streams~order, data=d, log="y")
mod <- lm(log10(num_of_streams)~order, data=d)
abline(mod)
text(2.5, 20,
     paste0('logN = ',
            round(mod$coefficients[[1]], 2),
            ' ',
            ifelse(
              mod$coefficients[[1]]<0,
              round(mod$coefficients[[1]], 2),
              paste0('+', round(mod$coefficients[[1]], 2))),
            ' ',
            'u'))
rb <- 1/10^mod$coefficients[[2]]
rb

#Estadísticas de red ampliadas
execGRASS(
  "r.stream.stats",
  flags = c('overwrite','quiet'),
  parameters = list(
    stream_rast = 'order-strahler',
    direction = 'rstrm_direccion',
    elevation = 'dem_tallado',
    output = 'brujuelas_stats_expanded.txt'
  )
)
file.show('brujuelas_stats_expanded.txt')


#mapview(order, col.regions = 'blue', legend = FALSE) ----
mapview(order, col.regions = 'blue', legend = FALSE)

# Obtener cursos más largos (cargar función propia)
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/lfp_network.R') #Cargada como función "LfpNetwork"
LfpNetwork(
  xycoords = c(436668.43711353285,	2046734.011679014),
  suffix = 'brujuelas',
  stream_vect = 'order_all',
  direction = 'rstrm_direccion'
)

# Imprimir lista de mapas ráster y vectoriales
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
# Representar con leaflet
library(leaflet)
lfp <- readVECT('LfpNetwork_lfp_all_final_brujuelas')
lfp4326 <- spTransform(lfp, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>%
  addProviderTiles(providers$OpenTopoMap, group = 'terrain') %>%
  addPolylines(
    data = lfp4326, weight = 3, opacity = 0.7, group = 'order',
    label = ~as.character(cat),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = T,
                                style = list(
                                  "font-size" = "8px",
                                  "background" = "rgba(255, 255, 255, 0.5)",
                                  "background-clip" = "padding-box",
                                  "padding" = "1px"))) %>% 
  leafem::addHomeButton(extent(lfp4326), 'Ver todo')

# Exportar a KML
execGRASS(
  'v.out.ogr',
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'LfpNetwork_lfp_all_final_brujuelas',
    output = 'lfp_brujuelas_kml.kml',
    format = 'KML',
    dsco = 'NameField=cat'
  )
)

# Obtención de perfiles longitudinales e índices de concavidad
library(sp)
use_sp()
source('https://raw.githubusercontent.com/geofis/rgrass/master/lfp_profiles_concavity.R') #Cargado como función "LfpProfilesConcavity"
brujuelas_conv_prof <- LfpProfilesConcavity(
  xycoords = my_trans(c(-69.599987, 18.510156)),
  network = 'LfpNetwork_lfp_all_final_brujuelas',
  prefix = 'brujuelas',
  dem = 'dem',
  direction = 'rstrm_direccion',
  crs = '+init=epsg:32619',
  smns = 1,
  nrow = 6)

## Mostrar resultados
#brujuelas_conv_prof$profiles
jpeg('brujuelas_conv_prof_profiles.jpg', width = 4000, height = 2000, res = 100)
brujuelas_conv_prof$profiles
dev.off()

layers <- brujuelas_conv_prof$dimensionlessprofiles$layers
filtered_layers <- lapply(layers, function(layer) {
  if (!inherits(layer$geom, c("GeomText", "GeomLabel"))) {
    return(layer)
  } else {
    return(NULL)
  }
})

brujuelas_conv_prof$dimensionlessprofiles$layers <- Filter(Negate(is.null), filtered_layers)
jpeg('brujuelas_conv_prof_dimensionlessprofiles.jpg', width = 4000, height = 2000, res = 100)
brujuelas_conv_prof$dimensionlessprofiles 
dev.off()

perfiles_a_quitar <- c("brujuelas-13", "brujuelas-16", "brujuelas-19", "brujuelas-29", "brujuelas-32", 
            "brujuelas-41", "brujuelas-44", "brujuelas-52", "brujuelas-54", "brujuelas-59", 
            "brujuelas-64", "brujuelas-65", "brujuelas-89", "brujuelas-98", "brujuelas-99", 
            "brujuelas-121", "brujuelas-16", "brujuelas-78", "brujuelas-95", "brujuelas-87",
            "brujuelas-67", "brujuelas-90")

filtered_data <- brujuelas_conv_prof$dimensionlessprofiles$data[
  !brujuelas_conv_prof$dimensionlessprofiles$data$stream %in%
    perfiles_a_quitar, ]

# Convertir el resultado a un dataframe, si no lo es ya
filtered_data <- as.data.frame(filtered_data)

# Ahora puedes trabajar con filtered_data o volver a asignarlo si necesitas
brujuelas_conv_prof$dimensionlessprofiles$data <- filtered_data

brujuelas_conv_prof$dimensionlessprofiles$data$stream <- factor(
  brujuelas_conv_prof$dimensionlessprofiles$data$stream,
  levels = levels(brujuelas_conv_prof$dimensionlessprofiles$data$stream)
)

# Verifica el contenido de filtered_data para asegurarte de que los perfiles fueron eliminados correctamente
print(head(filtered_data))
str(filtered_data)

jpeg('brujuelas_conv_prof$dimensionlessprofiles.jpg', width = 4000, height = 2000, res = 100)
brujuelas_conv_prof$dimensionlessprofiles 
dev.off()

brujuelas_conv_prof$concavityindex

## Tabla dx/dy, tanto en metros como adimensional. Útiles para construir perfiles por cuenta propia
brujuelas_conv_prof$lengthzdata %>% tibble::as.tibble()
brujuelas_conv_prof$lengthzdatadmnls %>% tibble::as.tibble()

# Clasificar los perfiles con agrupamiento jerarquico
brujuelas_conv_prof$concavityindex_limpio <- brujuelas_conv_prof$concavityindex[
  !brujuelas_conv_prof$concavityindex$stream %in% perfiles_a_quitar, ]

rownames(brujuelas_conv_prof$concavityindex_limpio) <-
  brujuelas_conv_prof$concavityindex_limpio$stream
brujuelas_conv_prof$concavityindex_limpio <- 
  brujuelas_conv_prof$concavityindex_limpio[, -1, drop = F]
brujuelas_conv_dist <- dist(brujuelas_conv_prof$concavityindex_limpio)
brujuelas_conv_dist_ward <- hclust(brujuelas_conv_dist, method = "ward.D2")
brujuelas_conv_dist_ward %>% plot()
brujuelas_conv_dist_ward_k4 <- cutree(tree = brujuelas_conv_dist_ward, k = 4)
rect.hclust(brujuelas_conv_dist_ward, k = 4)

# Perfiles colorizados por grupo
brujuelas_conv_prof$dimensionlessprofiles$data <- brujuelas_conv_prof$dimensionlessprofiles$data %>% 
  dplyr::inner_join(data.frame(stream = names(brujuelas_conv_dist_ward_k4), grupo = brujuelas_conv_dist_ward_k4))
brujuelas_conv_prof$dimensionlessprofiles$data$grupo <- factor(brujuelas_conv_prof$dimensionlessprofiles$data$grupo)
brujuelas_conv_prof$dimensionlessprofiles <- brujuelas_conv_prof$dimensionlessprofiles +
  geom_line(aes(color = factor(grupo)), size = 2) +
  facet_wrap(~ stream, ncol = 15) +
  scale_color_manual(values = c("1" = "red", "2" = "blue",
                                "3" = "green", "4" = "purple"),
                     name = "Grupo") +
  labs(color = "Grupo", fill = "Grupo") +
  theme(legend.position = "bottom",
        legend.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title.position = "top", nrow = 1))

jpeg('brujuelas_conv_prof$dimensionlessprofiles_colorizado.jpg', width = 3500, height = 2500, res = 150)

brujuelas_conv_prof$dimensionlessprofiles

dev.off()

brujuelas_conv_prof$dimensionlessprofiles2 <- brujuelas_conv_prof$dimensionlessprofiles
brujuelas_conv_prof$dimensionlessprofiles2$data$stream <- factor(
  brujuelas_conv_prof$dimensionlessprofiles2$data$stream,
  levels = unique(
    brujuelas_conv_prof$dimensionlessprofiles2$data$stream) %>%
    gtools::mixedsort(decreasing = T)
)

jpeg('brujuelas_conv_prof$dimensionlessprofiles2_colorizado.jpg', width = 3500, height = 2500, res = 150)

brujuelas_conv_prof$dimensionlessprofiles2

dev.off()



#Parámetros de cuenca con r.basin ----
# Convertir a números enteros la extensión y la resolución del DEM
dem_usado_en_region <- readRAST('dem')
library(raster)
dem_usado_en_region_r <- raster(dem_usado_en_region)
rawextent <-  extent(dem_usado_en_region_r)
rawextent
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/integerextent.R')
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/xyvector.R')
newextent <- intext(e = rawextent, r = 90, type = 'inner')
newextent
gdalUtils::gdalwarp(
  srcfile = 'dem-cuencas-brujuelas-cumayasa.tif',
  dstfile = 'demint.tif',
  te = xyvector(newextent),
  tr = c(90,90),
  r = 'bilinear',
  overwrite = T
)

## Importar a sesión de GRASS
rutademint <- 'demint.tif'
execGRASS(
  "g.proj",
  flags = c('t','c'),
  georef=rutademint)
gmeta()
execGRASS(
  "r.in.gdal",
  flags='overwrite',
  parameters=list(
    input=rutademint,
    output="demint"
  )
)
execGRASS(
  "g.region",
  parameters=list(
    raster = "demint",
    align = "demint"
  )
)
gmeta()
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
## Generar red de drenaje para obtener coordenada posteriormente
execGRASS(
  "r.stream.extract",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = 'demint',
    threshold = 80,
    stream_raster = 'stream-de-rstr',
    stream_vector = 'stream_de_rstr'
  )
)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
## Obtener coordenada
library(sp)
use_sp()
library(mapview)
netw <- spTransform(
  readVECT('stream_de_rstr'),
  CRSobj = CRS("+init=epsg:4326"))
mapview(netw, col.regions = 'blue', legend = FALSE)

## Transformar coordenada a EPSG:32619 como número entero
source('my-trans.R')
outlet <- as.integer(my_trans(c(-69.62048, 18.53288))) #no se puede usar las mismas coordenadas que el LFP

## Ejecutar `r.basin`
pref <- 'rbasin_brujuelas'
execGRASS(
  "r.basin",
  flags = 'overwrite',
  parameters = list(
    map = 'demint',
    prefix = pref,
    coordinates = outlet,
    threshold = 80,
    dir = 'salidas-rbasin/brujuelas'
  )
)
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
#> Si `r.basin` arrojara error (sólo en el caso de error, no en caso de advertencia), ejecutar este bloque para borrar las salidas anteriores y reejecutar el `r.basin`:
#execGRASS(
# "g.remove",
#flags = 'f',
#parameters = list(
# type = c('raster','vector'),
#pattern = paste0(pref, '*')
#)
#)

## Cargar los vectoriales transformados a EPSG:4326 para visualizar en leaflet
rbnetw <- spTransform(
  readVECT('rbasin_brujuelas_demint_network'),
  CRSobj = CRS("+init=epsg:4326"))
rbnetw
rbmain <- spTransform(
  readVECT('rbasin_brujuelas_demint_mainchannel'),
  CRSobj = CRS("+init=epsg:4326"))
rbmain
rbbasin <- spTransform(
  readVECT('rbasin_brujuelas_demint_basin'),
  CRSobj = CRS("+init=epsg:4326"))
rbbasin

library(leaflet)
leaflet() %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(data = rbnetw, weight = 3, opacity = 0.7) %>% 
  addPolylines(data = rbmain, weight = 3, opacity = 0.7, color = 'red') %>% 
  addPolygons(data = rbbasin) %>% 
  leafem::addHomeButton(extent(rbbasin), 'Ver cuenca')

## Explorar los parámetros de cuenca
library(readr)
rbbrujuelaspar1 <- read_csv("salidas-rbasin/brujuelas/rbasin_brujuelas_demint_parametersT.csv")
rbbrujuelaspar1 %>% tibble::as_tibble()
rbbrujuelaspar2 <- read_csv(
  "salidas-rbasin/brujuelas/rbasin_brujuelas_demint_parameters.csv",
  skip=2, col_names = c('Parameter', 'Value'))
rbbrujuelaspar2 %>% print(n=Inf)


#Curva e integral hipsométrica ----
# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
#* Nótese que los paquetes requeridos en esta sessión (`rgrass7`, `raster`, `leaflet`, `leafem`), fueron en el bloque anterior al ejecutarse el código contenido en el archivo `orden-de-red.Rmd`. Igualmente, dicho bloque de código creó todos los objetos necesarios para realizar este tutorial.
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)
## Representar cuencas
library(sp)
use_sp()
library(mapview)
bas2 <- readVECT('cuenca_brujuelas2')
bas3 <- readVECT('cuenca_brujuelas3')

## Curva e integral hipsométrica
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/integral_hypsometric_curve.R')
# source('integral_hypsometric_curve.R') #Cargada como función "HypsoIntCurve"
HypsoBasinsOrder2 <- HypsoIntCurve(
  basins = 'cuenca_brujuelas2',
  dem = 'dem',
  labelfield = 'cat',
  nrow = 1,
  labelsize = 4
)

HypsoBasinsOrder2$HypsoInt
write.table(HypsoBasinsOrder2$HypsoInt, file = "Brujuelas_HypsoBasinsOrder2$HypsoInt.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

jpeg('HypsoBasinsOrder2$HypsoCurve.jpg', width = 3500, height = 2500, res = 150)
HypsoBasinsOrder2$HypsoCurve
dev.off()

mapview(bas2, zcol='cat', col.regions = 'blue', legend = FALSE) %>%
  addStaticLabels(label = bas2$cat)

HypsoBasinsOrder3 <- HypsoIntCurve(
  basins = 'cuenca_brujuelas3',
  dem = 'dem',
  labelfield = 'cat',
  nrow = 1,
  labelsize = 4
)

HypsoBasinsOrder3$HypsoInt
write.table(HypsoBasinsOrder3$HypsoInt, file = "Brujuelas_HypsoBasinsOrder3$HypsoInt.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

HypsoBasinsOrder3$HypsoCurve
mapview(bas3, zcol='cat', col.regions = 'blue', legend = FALSE) %>%
  addStaticLabels(label = bas3$cat)

#knitr::spin("brujuelas.R", knit = T) #para tejer en HTM y Rmd