from shapely.ops import nearest_points


def internal_points_to_coast(point, coasts):
    """Move points that are inside land to the nearest point on the coast"""
    in_land = point.within(coasts.geometry)

    if any(in_land):
        land_idx = in_land[in_land].index[0]
        new_pt, _ = nearest_points(coasts.loc[land_idx, 'geometry'].exterior, point)
        return new_pt
    else:
        return point
