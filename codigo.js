/* ******************************************************************
CALCULO DE LA TEMPERATURA SUPERFICIAL DEL SUELO CON LANDSAT 8 TOA
EN EL ÁREA URBANA DEL MUNICIPIO DE SANTA CRUZ DE LA SIERRA - BOLIVIA
PARA ESTE EJEMPLO SE TOMÓ EN CUENTA LA FECHA DE MARZO DE 2024
********************************************************************* 
Para la elaboración de este script, se tomó como referencia el siguiente
estudio: https://waleedgeo.com/papers/waleed2022_pululc.pdf  (Mirza Waleed)
Para el enmascaramiento de nubes, de este tutorial: https://www.youtube.com/watch?v=UmqaqKIQ2gk  (Kenneth Ekpetere)
******************************************************************* */ 

// Creamos un transecto con punto de patida en el CEAM y termine en Las Brisas
var transect = ee.Geometry.LineString(
    [[-63.2273, -17.8010],// Área Munipal Curiche La Madre
     [-63.2158, -17.7972],// Rotonda 4to anillo Av. Piraí
     [-63.2046, -17.7949],// Mercado Abasto
     [-63.1885, -17.7872],// 1er anillo Av. Cañoto
     [-63.1822, -17.7833],// Plaza 24 de Septiembre
     [-63.1729, -17.7942],// Parque Urbano
     [-63.1519, -17.7893],// Feria Barrio Lindo
     [-63.0965, -17.7766],// Laguna guapilo (norte)
     [-63.0946, -17.7834],// Laguna guapilo (sur)
     [-63.0637, -17.7784]]);// Jardín Botánico
  
  // Llamamos a la colección Landsat 8
  var dataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .filterDate('2024-01-01', '2024-06-01')
    .filterBounds(table)
    .filterMetadata('CLOUD_COVER', 'less_than', 3)
  print('Colección Landsat 8 TOA',dataset)
  
  // Función para enmascarar la nube
  var getQABits = function(image, start, end, newName) {
      // Calcula los bits que necesitamos extraer.
      var pattern = 0;
      for (var i = start; i <= end; i++) {
         pattern += Math.pow(2, i);
      }
      // Devuelve una imagen de banda única de los bits de QA extraídos
      return image.select([0], [newName])
                    .bitwiseAnd(pattern)
                    .rightShift(start);
  };
  
  // Función para enmascarar los píxeles nublados.
  var cloud_shadows = function(image) {
    // Seleccione la banda QA.
    var QA = image.select(['QA_PIXEL']);
    // Obtiene el bit internal_cloud_algorithm_flag.
    return getQABits(QA, 10,11, 'Cloud_shadows').eq(1);
    // Devuelve una imagen enmascarando las zonas nubladas.
  };
    
  // Función para enmascarar los píxeles nublados.
  var clouds = function(image) {
    // Seleccione la banda QA.
    var QA = image.select(['QA_PIXEL']);
    // Obtiene el bit internal_cloud_algorithm_flag.
    return getQABits(QA, 3,4, 'Cloud').eq(0);
    // Devuelve una imagen enmascarando las zonas nubladas.
  };
  
  var maskClouds = function(image) {
    var cs = cloud_shadows(image);
    var c = clouds(image);
    image = image.updateMask(cs);
    return image.updateMask(c);
  };
  
  // Aplicamos la máscara de nubes
  var L8_mask = dataset.map(maskClouds).median()
  
  // Recortamos la imagen con el área de estudio
  var L8_urban = L8_mask.clip(table)
  
  // Creamos la simbología de visualización
  var vis_palet = {
    bands: ['B4', 'B3', 'B2'],
    min: 0.0,
    max: 0.2,
  };
  
  // Centramos el mapa y visualizamos
  Map.centerObject(table, 11);
  Map.addLayer(L8_urban, vis_palet, 'Landsat 8 (432)', false);
  
  // *************************************************************** //
  
  // Calculamos el NDVI
  var ndvi = L8_urban.normalizedDifference(['B5', 'B4']).rename('NDVI');
  
  // Visualizamos la imagen ndvi en el mapa
  Map.addLayer(ndvi, {min: -0.0994, max: 0.7591, 
    palette: ['blue', 'white', 'green']}, 'NDVI', false);
  
  // Obtenemos el valor máximo de NDVI
  var min = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: table,
    scale: 30,
    maxPixels: 1e10
  }).values().get(0))
  
  // Obtenemos el valor mínimo de NDVI
  var max = ee.Number(ndvi.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: table,
    scale: 30,
    maxPixels: 1e10
  }).values().get(0))
  
  // imprimimos los valores min y max
  print('Valor máximo de NDVI', max)
  print('Valor mínimo de NDVI', min)
  
  // ************************************************************** //
  
  // Calculamos la proporción de la vegetación
  var pv = (ndvi.subtract(min).divide(max.subtract(min))).pow(ee.Number(2))
  
  // Visualizamos la proporción de la vegetación en el mapa
  Map.addLayer(pv, {}, 'PV', false)
  
  // Definimos las constantes como a y b
  var a = ee.Number(0.004)
  var b = ee.Number(0.986)
  
  // Aplicamos la formula para EM
  var EM = pv.multiply(a).add(b).rename('EM')
  
  // Seleccionamos la "B10"
  var thermal = L8_urban.select('B10')
  
  // Aplicamos la formula para LST
  var LST = thermal.expression(
    '(T/(1+(0.00115*(T/0.48359547432))*log(e)))-273.15',{
      'T':thermal.select('B10'),
      'e':EM.select('EM')
    }).rename('LST')
  
  // Creamos la simbología para visualizar la temperatura
  var visLST = {min: 25, max:41, palette: [
  '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
  '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
  '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
  'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
  'ff0000', 'de0101', 'c21301', 'a71001', '911003'
   ]}
  
  // Visualizamos LST en el mapa
  Map.addLayer(LST, visLST,'LST');
  
  // Extraemos la temperatura mínima
  var tem_min = ee.Number(LST.reduceRegion({
    reducer: ee.Reducer.min(),
    geometry: table,
    scale: 30,
    maxPixels: 1e10
  }).values().get(0))
  
  // Extraemos la temperatura máxima
  var tem_max = ee.Number(LST.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: table,
    scale: 30,
    maxPixels: 1e10
  }).values().get(0))
  
  // imprimimos los valores min y max
  print('Temperatura máxima (°C):', tem_max)
  print('Temperatura mínima (°C):', tem_min)
  
  // ************************************************************** //
  
  // Llamamos a la colección Landsat 8 nuevamente
  var col8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")
    .filterDate('2024-03-26', '2024-03-28')
  
  // Seleccionamos la "B10" de la colección
  var temperature = col8.filterBounds(transect)
      .select(['B10'], ['temp'])
      .map(function(image) {
        return image.subtract(273.15)
            .set('system:time_start', image.get('system:time_start'));
      });
  
  // Obtenemos de las temperatura
  var temperatura = temperature.reduce(ee.Reducer.mean())
      .select([0], ['Temperatura']);
  
  // Obtenemos el NDVI
  var Ndvi = ndvi.reduce(ee.Reducer.mean())
      .select([0], ['NDVI']);
  
  //Delimitamos el transecto a partir del punto 1
  var startingPoint = ee.FeatureCollection(ee.Geometry.Point(-63.2273, -17.8010));
  var distance = startingPoint.distance(500000);
  var image = distance.addBands(Ndvi).addBands(temperatura);
  
  // Extraemos los valores de la temperatura para el transecto
  var array = image.reduceRegion(ee.Reducer.toList(), transect, 100)
                   .toArray(image.bandNames());
  
  // Orden de los puntos para el transecto.
  var distances = array.slice(0, 0, 1);
  array = array.sort(distances);
  
  // Creamos la matriz (Array) para el eje Y
  var ndviAndTemp = array.slice(0, 1);
  var distance = array.slice(0, 0, 1).project([1]);
  print('***** Gráfico NDVI vs LST (°C)*****')
  
  // Configuramos el gráfico.
  var chart = ui.Chart.array.values(ndviAndTemp, 1, distance)
      .setChartType('LineChart')
      .setSeriesNames(['NDVI', 'Temperatura'])
      .setOptions({
        title: 'Comportamiento del NDVI vs la Temperatura en el recorrido del transecto',
        titleTextStyle: {
          color: '#01625D',
          bold: true,
          fontSize: 18,
          position: 'right'
        },
        vAxes: {
          0: {
            title: 'Temperatura (°C)',
            titleTextStyle: {
              bold: true,
              color: '#01625D'
            }
          },
          1: {
            title: 'NDVI',
            titleTextStyle: {
              bold: true,
              color: '#01625D'
            },
            baselineColor: '#D9DAD9'
          }
        },
        hAxis: {
          title: 'Distancia en metros (Curiche La Madre - Jardín Botánico)',
          titleTextStyle: {
            bold: true,
            color: '#01625D'
          }
        },
        backgroundColor: '#F4F8F8', 
        interpolateNulls: true,
        pointSize: 2,
        lineWidth: 2,
        series: {
          0: {targetAxisIndex: 1},
          1: {targetAxisIndex: 0},
          2: {targetAxisIndex: 0}
        },
        chartArea: {
          width: '88%',
          height: '70%'
        },
        legend: {
          position: 'up', 
          textStyle: {fontSize: 16}
        }
      });
  
  // Imprimimos en consola el gráfico
  print(chart);
  
  // visualizamos el área de estudio
  var palet = {color: '#4E045F ', fillColor:'00000000', width: 3}
  Map.addLayer(table.style(palet), {}, 'Area Urbana SC') 
  
  // Mostramos el transecto en el mapa
  Map.addLayer(transect, {color: '#7F08B2'}, 'Transecto');
  
  //  **************************************************************************** //
  
  // Función para crear la barra de colores
  function createColorBar(palette, min, max) {
    var colorBar = ui.Panel({
      style: {
        position: 'bottom-right',
        padding: '10px 10px'
      }
    });
  
    var title = ui.Label('Temperatura (°C)', {fontWeight: 'bold'});
    var colorBarLabels = ui.Panel({
      widgets: [
        ui.Label(min),
        ui.Label((max - min) / 2 + min, {textAlign: 'center', stretch: 'horizontal'}),
        ui.Label(max)
      ],
      layout: ui.Panel.Layout.flow('horizontal')
    });
  
    var colorBarGradient = ui.Thumbnail({
      image: ee.Image.pixelLonLat().select(0),
      params: {
        bbox: [0, 1, 1, 0.1],
        dimensions: '200x30',
        format: 'png',
        min: 0,
        max: 1,
        palette: palette
      },
      style: {stretch: 'vertical', margin: '0px 8px', maxHeight: '100px'}
    });
    var footer = ui.Label('Created by: Jorge A. Zampieri', {
      fontSize: '12px',
      textAlign: 'left',
      color: '#0C6675',
      margin: '10px 0'
    });
    
    colorBar.add(title);
    colorBar.add(colorBarGradient);
    colorBar.add(colorBarLabels);
    colorBar.add(footer);
    return colorBar;
  }
  
  // Crear la barra de colores y agregarla al mapa
  var colorBar = createColorBar(visLST.palette, visLST.min, visLST.max);
  Map.add(colorBar);
  
  //Creamos el título del panel
  var title = ui.Label({
    value: 'Temperatura Superficial del Suelo (LST) - Santa Cruz de la Sierra',
    style:{
    fontWeight: 'bold',
    fontSize: '16px',
    color: '#052EA4'
    }});
  
  // Agregamos el ´titulo
  title.style().set('position', 'top-center');
  Map.add(title);
  