from qgis.core import QgsVectorLayer, QgsFeature, QgsGeometry, QgsPointXY, QgsProject
from qgis.core import QgsLineSymbol, QgsSingleSymbolRenderer
from qgis.core import QgsCategorizedSymbolRenderer, QgsSymbol, QgsRendererCategory, QgsWkbTypes
from qgis.PyQt.QtGui import QColor
from pathlib import Path
import csv
import sys

working_folder = Path(r'C:\temp\aziza\Direct_Indirect_Interface_V02\Direct_Indirect_Interface\Interface_Github\Output')  # <-- Update this path
filename='interface_K5_KF5_limiar105_maxtheta60_QT10_test-AR2019.csv'
csv_path= folder/filename

# 1. Clear all layers from the project
project = QgsProject.instance()
project.clear()  # Clears the project, including all layers[1][3][7]


########################################################################
# urb & flam
# 2. Find the most recent geopackage starting with "urb"
file_pattern = "urb*.gpkg"
gpkg_files = [os.path.join(working_folder, f) for f in os.listdir(working_folder) if f.startswith("urb") and f.endswith(".gpkg")]

if not gpkg_files:
    print("No geopackage files found starting with 'urb'")
else:
    # Get the most recent file by modification time
    latest_gpkg = max(gpkg_files, key=lambda x: os.path.getmtime(x))
    print("Loading:", latest_gpkg)

# 3. Load the geopackage
layer = QgsVectorLayer(latest_gpkg, os.path.basename(latest_gpkg), "ogr")
if not layer.isValid():
    print("Failed to load layer:", latest_gpkg)
else:
    project.addMapLayer(layer)
    # 4. Zoom to the layer extent
    canvas = iface.mapCanvas()
    canvas.setExtent(layer.extent())
    canvas.refresh()
    # 5. Set polygon fill color to green for all features
    symbol = QgsFillSymbol.createSimple({'color': 'gray', 'outline_style': 'no'})
    layer.setRenderer(QgsSingleSymbolRenderer(symbol))
    layer.triggerRepaint()

########################################################################
# flam
# find the most recent geopackage starting with "flam"
file_pattern = "urb*.gpkg"
gpkg_files = [os.path.join(working_folder, f) for f in os.listdir(working_folder) if f.startswith("flam") and f.endswith(".gpkg")]

if not gpkg_files:
    print("No geopackage files found starting with 'urb'")
else:
    # Get the most recent file by modification time
    latest_gpkg = max(gpkg_files, key=lambda x: os.path.getmtime(x))
    print("Loading:", latest_gpkg)

# 3. Load the geopackage
layer = QgsVectorLayer(latest_gpkg, os.path.basename(latest_gpkg), "ogr")
if not layer.isValid():
    print("Failed to load layer:", latest_gpkg)
else:
    project.addMapLayer(layer)
    # 4. Zoom to the layer extent
    canvas = iface.mapCanvas()
    canvas.setExtent(layer.extent())
    canvas.refresh()
    # 5. Set polygon fill color to green for all features
    symbol = QgsFillSymbol.createSimple({'color': 'green', 'outline_style': 'no'})
    layer.setRenderer(QgsSingleSymbolRenderer(symbol))
    layer.triggerRepaint()


##############################################################
# interface
# read csv
file_pattern = "interface*.csv"
csv_files = [os.path.join(working_folder, f) for f in os.listdir(working_folder) if f.startswith("interface") and f.endswith(".csv")]

if not csv_files:
    print("No csv files found starting with 'interface'")
else:
    # Get the most recent file by modification time
    latest_csv = max(csv_files, key=lambda x: os.path.getmtime(x))
    print("Loading:", latest_csv)
    
print(latest_csv)

sequences = []
current_sequence = []
current_type = None
types=[]

with open(latest_csv, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        x = float(row['x'])
        y = float(row['y'])
        linkL = int(row['linkL'])
        linkR = int(row['linkR'])
        vert_type = row['vert_type']

        # If current_type is None, start a new sequence with this type
        if current_type is None:
            current_type = vert_type

        # If the current row matches the current sequence type and both links are 1
        if vert_type == current_type and linkL == 1 and linkR == 1:
            current_sequence.append((x, y))
        else:
            # If the sequence is broken, add it if it has at least 2 points
            if len(current_sequence) >= 2:
                sequences.append(current_sequence)
                types.append(str(current_type))
            current_sequence = []
            current_type = None
            # If the new row is valid, start a new sequence
            if linkL == 1 and linkR == 1:
                current_sequence.append((x, y))
                current_type = vert_type

# Add the last sequence if it's valid
if len(current_sequence) >= 2:
    sequences.append(current_sequence)
    types.append(str(current_type))

# 
xyd = QgsVectorLayer("LineString?crs=EPSG:3763", "linestrings", "memory")
provider = xyd.dataProvider()

# Add the "length" field to the layer
xyd.startEditing()
provider.addAttributes([
    #QgsField("length", QVariant.Int),
    QgsField("type", QVariant.String)
    ])  # Use QVariant.Double if lengths are not integers
xyd.updateFields()

#Add features with their corresponding lengths
for i, seq in enumerate(sequences):
    feat = QgsFeature()
    points = [QgsPointXY(x, y) for x, y in seq]
    geom = QgsGeometry.fromPolylineXY(points)
    feat.setGeometry(geom)
    feat.setAttributes([types[i]])  # Set the "type" attribute
    provider.addFeature(feat)

xyd.commitChanges()  # Save changes

QgsProject.instance().addMapLayer(xyd)

categories = []
for i in range(1, 6):  # types 1 to 5
    # Interpolate color: red (1) -> orange (3) -> yellow (5)
    if i == 1:
        color = QColor(255, 0, 0)      # Red
    elif i == 3:
        color = QColor(255, 165, 0)    # Orange
    elif i == 5:
        color = QColor(255, 255, 0)    # Yellow
    else:
        if i < 3:
            r = 255
            g = int(165 * (i-1)/2)
            b = 0
        else:
            r = 255
            g = 165 + int(90 * (i-3)/2)
            b = 0
        color = QColor(r, g, b)
    # Create a line symbol for this category
    symbol = QgsSymbol.defaultSymbol(QgsWkbTypes.LineGeometry)
    symbol.setColor(color)
    # Optionally, set line width or style if needed
    symbol.setWidth(2)
    category = QgsRendererCategory(i, symbol, str(i))
    categories.append(category)

# Create and set the categorized renderer
renderer = QgsCategorizedSymbolRenderer("type", categories)
xyd.setRenderer(renderer)
xyd.triggerRepaint()


xyd.setName('Interface')
