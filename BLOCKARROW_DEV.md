# Block Arrow Development for ggplot2 Volcano Plots

## Overview
This document tracks the development of block arrow visualization for ggplot2-based volcano plots, following the pattern established by `ggh4x::geom_rectmargin()`.

## Completed Steps

### Step 1: Create geom_blockmargin base functionality ✓
- **File**: `R/jam-gg-blockarrow-geom.R`
- **Description**: Implemented a custom ggplot2 geom that draws rectangular blocks in plot margins
- **Key Features**:
  - Inherits from `ggplot2::GeomRect`
  - Supports placement on four sides: top, right, bottom, left (via `sides` parameter)
  - Uses `grid::unit()` for relative sizing (default 3% of plot size)
  - Respects coordinate transformations
  - Uses grid graphics for rendering

**Testing**: Confirmed rectangles appear in margins at bottom and top

### Step 2: Add text label support ✓
- **Description**: Extended `GeomBlockMargin` to display text labels
- **Implementation Details**:
  - Added `label` and `label_colour` aesthetic mappings
  - Labels position automatically based on which margin side(s) they're on
  - Rotation applied for vertical margins (left/right)
  - Centered alignment within margin rectangle
  - Customizable font size and color

**Testing**: Confirmed labels appear correctly in margin rectangles with proper alignment and rotation

## Next Steps (Pending)

### Step 3: Draw block arrow shapes
**Approach**:
- Modify `draw_panel()` to create arrow-head polygons instead of simple rectangles
- Use `grid::polygonGrob()` to construct arrow shapes
- Arrow configuration:
  - **Horizontal arrows** (top/bottom margins): pointing left/right
  - **Vertical arrows** (left/right margins): pointing up/down
  - Customizable arrow head size via parameters

**Implementation notes**:
- Create helper function to compute arrow polygon coordinates
- Consider parameters: `arrow_width`, `arrow_head_length` for fine-tuning appearance

### Step 4: Integrate into ggvolcano_plot()
**Approach**:
- Add parameters to `ggvolcano_plot()`:
  - `blockarrow = TRUE` (enable/disable)
  - `blockarrow_colors` (named color vector: "hit", "up", "down")
  - `blockarrow_length` (margin size)
  - `blockarrow_label_cex` (text size)
  
- Create data frame with hit counts and boundaries for each category:
  - **Top margin**: Two arrows for up/down regulated hits above significance threshold
  - **Right margin**: One arrow for total significant hits above fold-change threshold
  
- Add layers to plot with proper coordinate mapping

### Step 5: Move to jamma package
**Approach**:
- Copy `GeomBlockMargin` and `geom_blockmargin()` to jamma namespace
- Maintain compatibility with existing `blockArrowMargin()` function from base R
- Document in roxygen with proper parameters and aesthetics
- Add to NAMESPACE exports

---

## Technical Notes

### Grid Graphics Reference
The geom uses `grid` package functions:
- `grid::rectGrob()` - rectangles
- `grid::polygonGrob()` - arbitrary polygons (for block arrows)
- `grid::textGrob()` - text
- `grid::gList()` - collect multiple graphical objects
- `grid::gTree()` - create a drawable tree structure

### Coordinate System Considerations
- Uses `native` coordinates for data-space positioning
- Uses `npc` (normalized parent coordinates) for margin offsets
- The `coord$transform()` call handles coordinate system transformations
- Clipping with `clip = "off"` is automatic within margin regions

### Color and Styling
- Full support for `fill`, `colour`, `alpha`, `linewidth`, `linetype` aesthetics
- Uses `ggplot2::alpha()` to apply transparency
- Labels support custom color via `label_colour` aesthetic

---

## Key Design Decisions

1. **Pattern**: Following `ggh4x::geom_rectmargin()` ensures consistency with established ggplot2 extension practices

2. **Margin Placement**: Rectangles/arrows drawn in margins (outside plot area but inside plot+margin space) rather than completely outside, improving integration with facets

3. **Native Coordinates**: Using native data coordinates for x/y positioning ensures arrows scale with data limits

4. **Grid Graphics**: Using grid directly (rather than higher-level geoms) provides fine control over placement and shape

---

## Files Modified/Created

- ✓ `R/jam-gg-blockarrow-geom.R` - New geom implementation
- ✓ `R/jam-gg-volcano-plot.R` - Parameters added (blockarrow, blockarrow_colors)
- 🚧 Ready for Step 3: Block arrow polygon implementation
- 🚧 Ready for Step 4: ggvolcano_plot integration
- 🚧 Ready for Step 5: Move to jamma package

---

## Testing Checklist

- [x] Basic rectangle margin rendering
- [x] Multiple rectangles in same plot
- [x] Text labels in margins with proper rotation
- [x] Different sides (t, b, l, r)
- [ ] Arrow head shapes
- [ ] Faceted plots with block arrows
- [ ] Integration with actual volcano plot data
- [ ] Color and alpha styling
- [ ] Label color customization
