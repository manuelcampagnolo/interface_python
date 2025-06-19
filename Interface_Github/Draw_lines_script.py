from qgis.PyQt.QtCore import QVariant
from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
)

# Set the maximum allowed distance between consecutive points.
# Adjust this value based on your data's spatial scale (e.g., 100 for 100 meters).
max_distance = 50.0

# Load your points layer by its name ('xydDT')
points_layer = QgsProject.instance().mapLayersByName('Interface_Points')[0]

# Create a new memory layer for line features using the same CRS as the points layer
crs = points_layer.crs().authid()
line_layer = QgsVectorLayer(f'LineString?crs={crs}', 'Ordered_Lines', 'memory')
provider = line_layer.dataProvider()

# Optionally, copy attribute fields from the points layer
fields = points_layer.fields()
provider.addAttributes(fields)
line_layer.updateFields()

# Group features by the 'vert_type' attribute
groups = {}
for feat in points_layer.getFeatures():
    vt = feat['vert_type']
    groups.setdefault(vt, []).append(feat)

# Process each group separately
for vt, feats in groups.items():
    if len(feats) < 2:
        continue  # Need at least two points to draw a line

    # Sort features arbitrarily (e.g., by feature id)
    feats.sort(key=lambda f: f.id())
    
    # Build chains using a greedy nearest neighbor approach with a distance threshold
    chains = []
    feats_remaining = feats.copy()
    while feats_remaining:
        # Choose a starting point (here, the one with the smallest x value)
        start_feat = min(feats_remaining, key=lambda f: f.geometry().asPoint().x())
        chain = [start_feat]
        feats_remaining.remove(start_feat)
        current_feat = start_feat
        while True:
            # Find the nearest feature in feats_remaining
            candidates = [(f, current_feat.geometry().distance(f.geometry())) 
                          for f in feats_remaining]
            if not candidates:
                break
            nearest_feat, distance = min(candidates, key=lambda t: t[1])
            if distance <= max_distance:
                chain.append(nearest_feat)
                feats_remaining.remove(nearest_feat)
                current_feat = nearest_feat
            else:
                # No nearby point within the threshold; end this chain.
                break
        chains.append(chain)
    
    # Create line features for each chain: connect consecutive points in each chain.
    for chain in chains:
        if len(chain) < 2:
            continue  # Skip chains with only one point
        for i in range(len(chain) - 1):
            p1 = chain[i].geometry().asPoint()
            p2 = chain[i + 1].geometry().asPoint()
            line_geom = QgsGeometry.fromPolylineXY([p1, p2])
            line_feat = QgsFeature(fields)
            line_feat.setGeometry(line_geom)
            # Optionally, transfer attributes from the first point of the pair.
            line_feat.setAttributes(chain[i].attributes())
            provider.addFeature(line_feat)

# Finalize and add the new line layer to the project.
line_layer.updateExtents()
QgsProject.instance().addMapLayer(line_layer)


