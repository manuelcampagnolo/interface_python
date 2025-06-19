##############################################
    #    Main Library       #  
##############################################
import pandas as pd

##############################################
    #    Main Function       #  
##############################################
def extract_vertices(geodf):
    """
    Extract levels and coordinates from a GeoDataFrame.
    
    MULTIPOLYGON:
    L1 identifies main rings or holes ==> outer polygon or hole
    L2 specifies the ring ID within a particular polygon of the multipolygon ==> part
    L3 distinguishes between different multipolygons ==> feature
    
    The output is a pandas dataFrame
    """
    data = []
    for feature_index, geom in enumerate(geodf.geometry, start=1):  # Iterate over each geometry in the GeoDataFrame
        if geom is None:
            continue
        if geom.geom_type == 'MultiPolygon':  # Iterate over each Polygon in a MultiPolygon
            for part_index, poly in enumerate(geom.geoms, start=1):
                # Extract the exterior ring 
                exterior_coords = poly.exterior.coords
                for ring_index, (x, y) in enumerate(exterior_coords, start=1):
                    data.append([x, y, 1, part_index, feature_index])
                # Iterate over interior rings (holes)
                for interior_index, interior in enumerate(poly.interiors, start=1):
                    interior_coords = interior.coords
                    for ring_index, (x, y) in enumerate(interior_coords, start=1):
                        data.append([x, y, interior_index+1, part_index, feature_index])
        elif geom.geom_type == 'Polygon': 
            # Extract exterior ring of a single Polygon
            exterior_coords = geom.exterior.coords
            for ring_index, (x, y) in enumerate(exterior_coords, start=1):
                data.append([x, y, 1, 1, feature_index])
            # Extract interior rings (holes)
            for interior_index, interior in enumerate(geom.interiors, start=1):
                interior_coords = interior.coords
                for ring_index, (x, y) in enumerate(interior_coords, start=1):
                    data.append([x, y, interior_index+1, 1, feature_index])
        else:
            print(f"Unsupported geometry type: {geom.geom_type}")
    
    df = pd.DataFrame(data, columns=['x', 'y', 'L1', 'L2', 'L3'])
    return df

# # def extract_vertices(geodf):
# #     """
# #     Extract levels and coordinates from a GeoDataFrame with 0-based indexing.
    
# #     MULTIPOLYGON:
# #     L1 identifies main rings or holes ==> outer polygon or hole (0 for exterior, 1, 2... for holes)
# #     L2 specifies the ring ID within a particular polygon of the multipolygon ==> part (0-based)
# #     L3 distinguishes between different multipolygons ==> feature (0-based)
    
# #     The output is a pandas DataFrame.
# #     """
# #     data = []
# #     for feature_index, geom in enumerate(geodf.geometry):  # L3 starts from 0
# #         if geom is None:
# #             continue
# #         if geom.geom_type == 'MultiPolygon':
# #             for part_index, poly in enumerate(geom.geoms):  # L2 starts from 0
# #                 # Exterior ring: L1 = 0
# #                 for ring_index, (x, y) in enumerate(poly.exterior.coords):
# #                     data.append([x, y, 0, part_index, feature_index])
# #                 # Interior rings: L1 = 1, 2, ...
# #                 for interior_index, interior in enumerate(poly.interiors):
# #                     for ring_index, (x, y) in enumerate(interior.coords):
# #                         data.append([x, y, interior_index + 1, part_index, feature_index])
# #         elif geom.geom_type == 'Polygon':
# #             # Exterior ring: L1 = 0, L2 = 0
# #             for ring_index, (x, y) in enumerate(geom.exterior.coords):
# #                 data.append([x, y, 0, 0, feature_index])
# #             # Interior rings
# #             for interior_index, interior in enumerate(geom.interiors):
# #                 for ring_index, (x, y) in enumerate(interior.coords):
# #                     data.append([x, y, interior_index + 1, 0, feature_index])
# #         else:
# #             print(f"Unsupported geometry type: {geom.geom_type}")
    
# #     df = pd.DataFrame(data, columns=['x', 'y', 'L1', 'L2', 'L3'])
# #     return df

# def extract_vertices(geodf):
#     """
#     Extract levels and coordinates from a GeoDataFrame with 0-based indexing.

#     MULTIPOLYGON:
#     L1 identifies rings (0 = exterior, 1+ = interior)
#     L2 specifies the part (polygon) ID within a MultiPolygon
#     L3 identifies the feature index from the original GeoDataFrame

#     Returns:
#     pd.DataFrame: with columns ['x', 'y', 'L1', 'L2', 'L3']
#     """
#     data = []

#     for feature_index, geom in enumerate(geodf.geometry):
#         if geom is None or geom.is_empty:
#             continue

#         if geom.geom_type == 'MultiPolygon':
#             for part_index, poly in enumerate(geom.geoms):
#                 if poly.is_empty:
#                     continue
#                 # Exterior ring
#                 for x, y in poly.exterior.coords:
#                     data.append([x, y, 0, part_index, feature_index])
#                 # Interior rings
#                 for interior_index, interior in enumerate(poly.interiors):
#                     for x, y in interior.coords:
#                         data.append([x, y, interior_index + 1, part_index, feature_index])

#         elif geom.geom_type == 'Polygon':
#             if geom.is_empty:
#                 continue
#             # Exterior ring
#             for x, y in geom.exterior.coords:
#                 data.append([x, y, 0, 0, feature_index])
#             # Interior rings
#             for interior_index, interior in enumerate(geom.interiors):
#                 for x, y in interior.coords:
#                     data.append([x, y, interior_index + 1, 0, feature_index])
#         else:
#             print(f"Unsupported geometry type at index {feature_index}: {geom.geom_type}")
    
#     if not data:
#         raise ValueError("No valid vertices were extracted. Check the input geometries.")

#     df = pd.DataFrame(data, columns=['x', 'y', 'L1', 'L2', 'L3'])
#     return df
