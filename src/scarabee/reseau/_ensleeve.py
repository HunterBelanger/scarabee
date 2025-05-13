from .._scarabee import (
    CrossSection,
    EmptyCell,
    Cartesian2D,
)
from typing import List, Tuple


def _ensleeve_quarter(
    pin_geom: Cartesian2D, pitch: float, gap_width: float, gap_xs: CrossSection
) -> Tuple[Cartesian2D, List[int]]:
    #  ---|
    #   PG|

    # First make the gap on the right hand side
    right_gap_Ny = round(
        (pin_geom.dy + gap_width - (0.5 * pitch + gap_width)) / (0.5 * pitch)
    )
    right_gap_dy = (
        pin_geom.dy + gap_width - (0.5 * pitch + gap_width)
    ) / right_gap_Ny  # Should be exactly 0.5*pitch !
    right_gap = Cartesian2D(
        [gap_width], right_gap_Ny * [right_gap_dy] + [0.5 * pitch + gap_width]
    )
    right_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width)]
    for i in range(right_gap_Ny):
        right_gap_tiles.append(EmptyCell(gap_xs, gap_width, right_gap_dy))
    right_gap.set_tiles(right_gap_tiles)

    # Now make the gap on the top
    top_gap_Nx = round(pin_geom.dx / (0.5 * pitch))
    top_gap_dx = pin_geom.dx / top_gap_Nx  # Should be exactly 0.5*pitch !
    top_gap = Cartesian2D(top_gap_Nx * [top_gap_dx], [gap_width])
    top_gap_tiles = []
    for i in range(top_gap_Nx):
        top_gap_tiles.append(EmptyCell(gap_xs, top_gap_dx, gap_width))
    top_gap.set_tiles(top_gap_tiles)

    # Make the temp geometry which is the pin_geom with the top gap
    temp_geom = Cartesian2D([pin_geom.dx], [pin_geom.dy, gap_width])
    temp_geom.set_tiles([top_gap, pin_geom])

    # Make the sleeved geometry
    sleeved_geom = Cartesian2D([pin_geom.dx, gap_width], [pin_geom.dy + gap_width])
    sleeved_geom.set_tiles([temp_geom, right_gap])

    # Get the FSR IDs for all the gap cells
    gap_fsr_ids = []
    for cell in right_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())
    for cell in top_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())

    return sleeved_geom, gap_fsr_ids


def _ensleeve_half(
    pin_geom: Cartesian2D, pitch: float, gap_width: float, gap_xs: CrossSection
) -> Tuple[Cartesian2D, List[int]]:
    #  |--|
    #  |PG|

    # First make the gap on the right hand side
    right_gap_Ny = round(
        (pin_geom.dy + gap_width - (0.5 * pitch + gap_width)) / (0.5 * pitch)
    )
    right_gap_dy = (
        pin_geom.dy + gap_width - (0.5 * pitch + gap_width)
    ) / right_gap_Ny  # Should be exactly 0.5*pitch !
    right_gap = Cartesian2D(
        [gap_width], right_gap_Ny * [right_gap_dy] + [0.5 * pitch + gap_width]
    )
    right_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width)]
    for i in range(right_gap_Ny):
        right_gap_tiles.append(EmptyCell(gap_xs, gap_width, right_gap_dy))
    right_gap.set_tiles(right_gap_tiles)

    # Make the gap on the left hand side
    left_gap_Ny = round(
        (pin_geom.dy + gap_width - (0.5 * pitch + gap_width)) / (0.5 * pitch)
    )
    left_gap_dy = (
        pin_geom.dy + gap_width - (0.5 * pitch + gap_width)
    ) / left_gap_Ny  # Should be exactly 0.5*pitch !
    left_gap = Cartesian2D(
        [gap_width], left_gap_Ny * [left_gap_dy] + [0.5 * pitch + gap_width]
    )
    left_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width)]
    for i in range(left_gap_Ny):
        left_gap_tiles.append(EmptyCell(gap_xs, gap_width, left_gap_dy))
    left_gap.set_tiles(left_gap_tiles)

    # Now make the gap on the top
    top_gap_Nx = round(pin_geom.dx / (0.5 * pitch))
    top_gap_dx = pin_geom.dx / top_gap_Nx  # Should be exactly 0.5*pitch !
    top_gap = Cartesian2D(top_gap_Nx * [top_gap_dx], [gap_width])
    top_gap_tiles = []
    for i in range(top_gap_Nx):
        top_gap_tiles.append(EmptyCell(gap_xs, top_gap_dx, gap_width))
    top_gap.set_tiles(top_gap_tiles)

    # Make the temp geometry which is the pin_geom with the top gap
    temp_geom = Cartesian2D([pin_geom.dx], [pin_geom.dy, gap_width])
    temp_geom.set_tiles([top_gap, pin_geom])

    # Make the sleeved geometry
    sleeved_geom = Cartesian2D(
        [gap_width, pin_geom.dx, gap_width], [pin_geom.dy + gap_width]
    )
    sleeved_geom.set_tiles([left_gap, temp_geom, right_gap])

    # Get the FSR IDs for all the gap cells
    gap_fsr_ids = []
    for cell in left_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())
    for cell in right_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())
    for cell in top_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())

    return sleeved_geom, gap_fsr_ids


def _ensleeve_full(
    pin_geom: Cartesian2D, pitch: float, gap_width: float, gap_xs: CrossSection
) -> Tuple[Cartesian2D, List[int]]:
    #  |--|
    #  |PG|
    #  |--|

    # First make the gap on the right hand side
    right_gap_Ny = round(
        (pin_geom.dy + 2.0 * gap_width - 2.0 * (0.5 * pitch + gap_width))
        / (0.5 * pitch)
    )
    right_gap_dy = (
        pin_geom.dy + 2.0 * gap_width - 2.0 * (0.5 * pitch + gap_width)
    ) / right_gap_Ny  # Should be exactly 0.5*pitch !
    right_gap = Cartesian2D(
        [gap_width],
        [0.5 * pitch + gap_width]
        + right_gap_Ny * [right_gap_dy]
        + [0.5 * pitch + gap_width],
    )
    right_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width)]
    for i in range(right_gap_Ny):
        right_gap_tiles.append(EmptyCell(gap_xs, gap_width, right_gap_dy))
    right_gap_tiles.append(EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width))
    right_gap.set_tiles(right_gap_tiles)

    # Make the gap on the left hand side
    left_gap_Ny = round(
        (pin_geom.dy + 2.0 * gap_width - 2.0 * (0.5 * pitch + gap_width))
        / (0.5 * pitch)
    )
    left_gap_dy = (
        pin_geom.dy + 2.0 * gap_width - 2.0 * (0.5 * pitch + gap_width)
    ) / left_gap_Ny  # Should be exactly 0.5*pitch !
    left_gap = Cartesian2D(
        [gap_width],
        [0.5 * pitch + gap_width]
        + left_gap_Ny * [left_gap_dy]
        + [0.5 * pitch + gap_width],
    )
    left_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width)]
    for i in range(left_gap_Ny):
        left_gap_tiles.append(EmptyCell(gap_xs, gap_width, left_gap_dy))
    left_gap_tiles.append(EmptyCell(gap_xs, gap_width, 0.5 * pitch + gap_width))
    left_gap.set_tiles(left_gap_tiles)

    # Now make the gap on the top
    top_gap_Nx = round(pin_geom.dx / (0.5 * pitch))
    top_gap_dx = pin_geom.dx / top_gap_Nx  # Should be exactly 0.5*pitch !
    top_gap = Cartesian2D(top_gap_Nx * [top_gap_dx], [gap_width])
    top_gap_tiles = []
    for i in range(top_gap_Nx):
        top_gap_tiles.append(EmptyCell(gap_xs, top_gap_dx, gap_width))
    top_gap.set_tiles(top_gap_tiles)

    # Now make the gap on the bottom
    bot_gap_Nx = round(pin_geom.dx / (0.5 * pitch))
    bot_gap_dx = pin_geom.dx / bot_gap_Nx  # Should be exactly 0.5*pitch !
    bot_gap = Cartesian2D(bot_gap_Nx * [bot_gap_dx], [gap_width])
    bot_gap_tiles = []
    for i in range(bot_gap_Nx):
        bot_gap_tiles.append(EmptyCell(gap_xs, bot_gap_dx, gap_width))
    bot_gap.set_tiles(bot_gap_tiles)

    # Make the temp geometry which is the pin_geom with the top gap
    temp_geom = Cartesian2D([pin_geom.dx], [gap_width, pin_geom.dy, gap_width])
    temp_geom.set_tiles([top_gap, pin_geom, bot_gap])

    # Make the sleeved geometry
    sleeved_geom = Cartesian2D(
        [gap_width, pin_geom.dx, gap_width], [pin_geom.dy + gap_width]
    )
    sleeved_geom.set_tiles([left_gap, temp_geom, right_gap])

    # Get the FSR IDs for all the gap cells
    gap_fsr_ids = []
    for cell in left_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())
    for cell in right_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())
    for cell in top_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())
    for cell in bot_gap_tiles:
        gap_fsr_ids.append(cell.get_all_fsr_ids())

    return sleeved_geom, gap_fsr_ids
