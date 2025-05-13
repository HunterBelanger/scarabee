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

    # First make the gap on the right hand side going from bottom to top
    LY = pin_geom.dy + gap_width
    dy = []
    right_gap_tiles = []
    while LY > 0.:
        ly = 0.5*pitch
        break_after = False
        if LY > 0.5*pitch and LY < pitch:
            ly = LY
            break_after = True
        LY -= ly
        dy.append(ly)
        right_gap_tiles.append(EmptyCell(gap_xs, gap_width, ly))
        if break_after:
            break
    right_gap_tiles.reverse()
    right_gap = Cartesian2D([gap_width], dy)
    right_gap.set_tiles(right_gap_tiles)
    
    # Make top gap 
    LX = pin_geom.dx
    dx = []
    top_gap_tiles = []
    while LX > 0.:
        lx = 0.5*pitch
        break_after = False
        if LX > 0.5*pitch and LX < pitch:
            lx = LX
            break_after = True
        LX -= lx
        dx.append(lx)
        top_gap_tiles.append(EmptyCell(gap_xs, lx, gap_width))
        if break_after:
            break
    top_gap = Cartesian2D(dx, [gap_width])
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
    
    # First make the gap on the right hand side going from bottom to top
    LY = pin_geom.dy + gap_width
    dy = []
    right_gap_tiles = []
    while LY > 0.:
        ly = 0.5*pitch
        break_after = False
        if LY > 0.5*pitch and LY < pitch:
            ly = LY
            break_after = True
        LY -= ly
        dy.append(ly)
        right_gap_tiles.append(EmptyCell(gap_xs, gap_width, ly))
        if break_after:
            break
    right_gap_tiles.reverse()
    right_gap = Cartesian2D([gap_width], dy)
    right_gap.set_tiles(right_gap_tiles)

    # Make the gap on the left hand side going from bottom to top
    LY = pin_geom.dy + gap_width
    dy = []
    left_gap_tiles = []
    while LY > 0.:
        ly = 0.5*pitch
        break_after = False
        if LY > 0.5*pitch and LY < pitch:
            ly = LY
            break_after = True
        LY -= ly
        dy.append(ly)
        left_gap_tiles.append(EmptyCell(gap_xs, gap_width, ly))
        if break_after:
            break
    left_gap_tiles.reverse()
    left_gap = Cartesian2D([gap_width], dy)
    left_gap.set_tiles(left_gap_tiles)

    # Make top gap 
    LX = pin_geom.dx
    dx = []
    top_gap_tiles = []
    while LX > 0.:
        lx = 0.5*pitch
        break_after = False
        if LX > 0.5*pitch and LX < pitch:
            lx = LX
            break_after = True
        LX -= lx
        dx.append(lx)
        top_gap_tiles.append(EmptyCell(gap_xs, lx, gap_width))
        if break_after:
            break
    top_gap = Cartesian2D(dx, [gap_width])
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

    # First make the gap on the right hand side going from bottom to top
    LY = pin_geom.dy + 2.*gap_width
    extra_y = 0.5*(LY - 0.5*pitch*round(LY / (0.5*pitch)))
    dy = [0.5*pitch + extra_y]
    right_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5*pitch+extra_y)]
    LY -= 0.5*pitch + extra_y
    while LY > 0.:
        ly = 0.5*pitch
        break_after = False
        if LY > 0.5*pitch and LY < pitch:
            ly = LY
            break_after = True
        LY -= ly
        dy.append(ly)
        right_gap_tiles.append(EmptyCell(gap_xs, gap_width, ly))
        if break_after:
            break
    right_gap_tiles.reverse()
    right_gap = Cartesian2D([gap_width], dy)
    right_gap.set_tiles(right_gap_tiles)

    # Make the gap on the left hand side going from bottom to top
    LY = pin_geom.dy + 2.*gap_width
    LY = pin_geom.dy + 2.*gap_width
    extra_y = 0.5*(LY - 0.5*pitch*round(LY / (0.5*pitch)))
    dy = [0.5*pitch + extra_y]
    left_gap_tiles = [EmptyCell(gap_xs, gap_width, 0.5*pitch+extra_y)]
    LY -= 0.5*pitch + extra_y
    while LY > 0.:
        ly = 0.5*pitch
        break_after = False
        if LY > 0.5*pitch and LY < pitch:
            ly = LY
            break_after = True
        LY -= ly
        dy.append(ly)
        left_gap_tiles.append(EmptyCell(gap_xs, gap_width, ly))
        if break_after:
            break
    left_gap_tiles.reverse()
    left_gap = Cartesian2D([gap_width], dy)
    left_gap.set_tiles(left_gap_tiles)

    # Make top gap 
    LX = pin_geom.dx
    dx = []
    top_gap_tiles = []
    while LX > 0.:
        lx = 0.5*pitch
        break_after = False
        if LX > 0.5*pitch and LX < pitch:
            lx = LX
            break_after = True
        LX -= lx
        dx.append(lx)
        top_gap_tiles.append(EmptyCell(gap_xs, lx, gap_width))
        if break_after:
            break
    top_gap = Cartesian2D(dx, [gap_width])
    top_gap.set_tiles(top_gap_tiles)

    # Make bottom gap 
    LX = pin_geom.dx
    dx = []
    bot_gap_tiles = []
    while LX > 0.:
        lx = 0.5*pitch
        break_after = False
        if LX > 0.5*pitch and LX < pitch:
            lx = LX
            break_after = True
        LX -= lx
        dx.append(lx)
        top_gap_tiles.append(EmptyCell(gap_xs, lx, gap_width))
        if break_after:
            break
    bot_gap = Cartesian2D(dx, [gap_width])
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
