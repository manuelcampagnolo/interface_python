##############################################
    #    Libraries       #  
##############################################
import sys
import numpy as np
import pandas as pd
bigN = 10**6 
from shapely.geometry import LineString
##############################################
    #    Functions       #  
##############################################
def decision(D, limiar, limiartheta, xV, yV, xF, yF, xW, yW, xWW, yWW, xWWW, yWWW, verbose=False, log_file="decision_table.csv"):
    """
    Determines if a flammable point is protected based on distances, angles, and whether
    the flammable-to-urban line intersects ANY urban edge.
    """

    # Distance calculations
    d2VF = (xV - xF)**2 + (yV - yF)**2
    d2VW = (xV - xW)**2 + (yV - yW)**2
    d2WF = (xW - xF)**2 + (yW - yF)**2

    sqrt_d2VW = np.sqrt(d2VW)
    sqrt_d2WF = np.sqrt(d2WF)
    sqrt_d2VF = np.sqrt(d2VF)

    # Angle calculation
    dot = (xW - xV) * (xF - xV) + (yW - yV) * (yF - yV)
    thetaV = np.full_like(dot, limiartheta)

    cond = (
        ~np.isnan(xV) & ~np.isnan(xF) & ~np.isnan(xW) &
        (d2VW > 0) & (d2VF > 0) & (dot**2 <= d2VW * d2VF)
    )
    thetaV[cond] = np.degrees(np.arccos(dot[cond] / np.sqrt(d2VW[cond] * d2VF[cond])))

    # Triangle perimeter condition
    Q = 0.1
    perimeter = sqrt_d2VW + sqrt_d2WF + sqrt_d2VF
    Q_condition = (
        (sqrt_d2VW > perimeter * Q) &
        (sqrt_d2WF > perimeter * Q) &
        (sqrt_d2VF > perimeter * Q)
    )

    # Build all FV segments
    fv_segments = [LineString([(xF[i], yF[i]), (xV[i], yV[i])]) for i in range(len(xF))]

    # Build all urban edges once
    urban_segments = [LineString([(xW[i], yW[i]), (xWW[i], yWW[i])]) for i in range(len(xW))]


    protedge_any = np.array([
        any(fv.intersects(urb) for urb in urban_segments)
        for fv in fv_segments
    ])

    protedge_prev = np.array([
        LineString([(xW[i], yW[i]), (xWW[i], yWW[i])]).intersects(LineString([(xWWW[i], yWWW[i]), (xF[i], yF[i])]))
        for i in range(len(xF))
    ])

    # Logical conditions
    
    condition_artifact = (xF == bigN) & (yF == bigN)
    condition_intersects = protedge_any | protedge_prev
    condition_valid_distances = ~np.isnan(d2VF) & ~np.isnan(d2VW) & ~np.isnan(d2WF)
    condition_positive_distances = (d2VW > 0) & (d2VF > 0) & (d2WF > 0)
    condition_threshold = (sqrt_d2VW + sqrt_d2WF < limiar * sqrt_d2VF)
    condition_angle = (thetaV < limiartheta)
    condition_closer_urban = (d2VF >= d2VW) & (d2VF >= d2WF)
    condition_triangle = Q_condition
    condition_outside_region = (d2VF > 2 * D**2)

    # Final protection decision
    result = (
        condition_intersects |condition_artifact|
        (
            ~condition_outside_region &
            condition_valid_distances &
            condition_positive_distances &
            condition_threshold &
            condition_angle &
            condition_closer_urban &
            condition_triangle
        )
    )

    # Log result
    table = pd.DataFrame({
        "xF": xF, "yF": yF,
        "xV": xV, "yV": yV,
        "xW": xW, "yW": yW,
        "xWW": xWW, "yWW": yWW,
        "xWWW": xWWW, "yWWW": yWWW,
        "Artifact Condition": condition_artifact,
        "Intersects Urban Edge": condition_intersects,
        "Valid Distances": condition_valid_distances,
        "Positive Distances": condition_positive_distances,
        "Threshold Condition": condition_threshold,
        "Angle Condition": condition_angle,
        "Closer Urban Condition": condition_closer_urban,
        "Triangle Condition": condition_triangle,
        "Outside Region Condition": condition_outside_region,
        "Final Decision": result
    })

    table.to_csv(log_file, index=False)

    if verbose:
        print(f"Results logged to {log_file}")

    return result


