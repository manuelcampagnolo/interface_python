from shapely.geometry import Point
from shapely.ops import nearest_points
import geopandas as gpd
import pandas as pd

# def convert_3763_XY_into_urban_closest_vertex(X,Y, urban_path):
#     # Create a Point geometry from the input coordinates
#     urb = gpd.read_file(urban_path)
#     target_point = Point(X, Y)
#     # Initialize variables to store the closest point and distance
#     closest_vertex = None
#     min_distance = float('inf')
#     # Iterate through all polygons in the GeoDataFrame
#     for polygon in urb.geometry:
#         # Get the nearest point on the polygon boundary to the target point
#         nearest_point_on_polygon = nearest_points(polygon.boundary, target_point)[0]
#         # Calculate the distance between the target point and this nearest point
#         distance = target_point.distance(nearest_point_on_polygon)
#         # Update the closest vertex if this point is closer
#         if distance < min_distance:
#             min_distance = distance
#             closest_vertex = nearest_point_on_polygon
#     # Convert the closest vertex to a dictionary with keys 'X' and 'Y'
#     result = {'X': closest_vertex.x, 'Y': closest_vertex.y}
#     # convert to data frame
#     x0y0 = pd.DataFrame([result])
#     return x0y0

from shapely.geometry import Point
import geopandas as gpd
import pandas as pd

def convert_3763_XY_into_urban_closest_vertex(X, Y, urban_path):
    # Load the urban geometries
    urb = gpd.read_file(urban_path)

    # Create a Point object from input coordinates
    target_point = Point(X, Y)

    # Create a GeoSeries of polygon boundaries
    boundaries = urb.geometry.boundary

    # Compute distances to the target point
    distances = boundaries.distance(target_point)

    # Get the geometry with the minimum distance
    min_idx = distances.idxmin()

    # Compute the nearest point on the boundary
    nearest_pt = nearest_points(boundaries[min_idx], target_point)[0]

    # Return as a DataFrame
    return pd.DataFrame([{'X': nearest_pt.x, 'Y': nearest_pt.y}])
