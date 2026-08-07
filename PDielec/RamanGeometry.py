#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
"""Shared geometry helpers for layered Raman calculations."""


def resolve_collection_angle(incident_angle, collection_angle, collection_side):
    """Return the signed detector angle in the same units as the inputs.

    ``collection_angle=None`` selects the automatic convention: an antiparallel
    (retro-backscattered) ray on the superstrate side and a collinear forward
    ray on the substrate side. An explicit signed angle is returned unchanged,
    so a positive angle equal to the incident angle retains the historical
    specular-reflection geometry.
    """
    if collection_side not in {"superstrate", "substrate"}:
        msg = f"collection_side must be 'superstrate' or 'substrate', got {collection_side!r}"
        raise ValueError(msg)
    if collection_angle is not None:
        return float(collection_angle)
    incident_angle = float(incident_angle)
    return -incident_angle if collection_side == "superstrate" else incident_angle
