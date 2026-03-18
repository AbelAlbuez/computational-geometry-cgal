"""
VTK Interactive Polygon Drawer
================================
Click points on the canvas to create a polygon.

Controls:
  Left Click       - Add a point
  Right Click      - Close & finalize the polygon
  Middle Click     - Remove the last point (undo)
  'c' key          - Clear all points and start over
  'q' / ESC        - Quit
"""

import vtk


# colours 
POINT_COLOR       = (1.0, 0.3, 0.1)   # red-orange dots
LINE_COLOR        = (0.2, 0.8, 0.2)   # green outline
FILL_COLOR        = (0.2, 0.6, 1.0)   # blue fill
CLOSED_LINE_COLOR = (1.0, 0.9, 0.1)   # yellow when closed
BG_COLOR          = (0.12, 0.12, 0.18)

POINT_RADIUS = 6    # screen pixels (glyph size)
LINE_WIDTH   = 2.0


#  helpers 

def make_sphere_source(radius=POINT_RADIUS):
    s = vtk.vtkSphereSource()
    s.SetRadius(radius)
    s.SetPhiResolution(12)
    s.SetThetaResolution(12)
    return s


class PolygonDrawer:
    """Manages geometry and actors for an interactive polygon."""

    def __init__(self, renderer):
        self.ren = renderer
        self.points_3d: list[tuple[float, float, float]] = []
        self.closed = False

        #  point glyphs 
        self.pt_pts   = vtk.vtkPoints()
        self.pt_pd    = vtk.vtkPolyData()
        self.pt_pd.SetPoints(self.pt_pts)

        sphere = make_sphere_source()
        self.pt_glyph = vtk.vtkGlyph3D()
        self.pt_glyph.SetInputData(self.pt_pd)
        self.pt_glyph.SetSourceConnection(sphere.GetOutputPort())
        self.pt_glyph.SetScaleModeToDataScalingOff()

        pt_mapper = vtk.vtkPolyDataMapper()
        pt_mapper.SetInputConnection(self.pt_glyph.GetOutputPort())

        self.pt_actor = vtk.vtkActor()
        self.pt_actor.SetMapper(pt_mapper)
        self.pt_actor.GetProperty().SetColor(*POINT_COLOR)
        self.ren.AddActor(self.pt_actor)

        #  outline polyline 
        self.line_pts  = vtk.vtkPoints()
        self.line_cells = vtk.vtkCellArray()
        self.line_pd   = vtk.vtkPolyData()
        self.line_pd.SetPoints(self.line_pts)
        self.line_pd.SetLines(self.line_cells)

        line_mapper = vtk.vtkPolyDataMapper()
        line_mapper.SetInputData(self.line_pd)

        self.line_actor = vtk.vtkActor()
        self.line_actor.SetMapper(line_mapper)
        self.line_actor.GetProperty().SetColor(*LINE_COLOR)
        self.line_actor.GetProperty().SetLineWidth(LINE_WIDTH)
        self.ren.AddActor(self.line_actor)

        #  filled polygon 
        self.fill_actor = None

    #  public API 

    def add_point(self, x, y, z=0.0):
        if self.closed:
            return
        self.points_3d.append((x, y, z))
        self._rebuild()

    def undo_last(self):
        if self.points_3d and not self.closed:
            self.points_3d.pop()
            self._rebuild()

    def close(self):
        if len(self.points_3d) < 3 or self.closed:
            return
        self.closed = True
        self._rebuild()
        self._build_fill()

    def clear(self):
        self.points_3d.clear()
        self.closed = False
        if self.fill_actor:
            self.ren.RemoveActor(self.fill_actor)
            self.fill_actor = None
        self._rebuild()

    #  private 

    def _rebuild(self):
        """Refresh point glyphs and the polyline."""
        # points
        self.pt_pts.Reset()
        for p in self.points_3d:
            self.pt_pts.InsertNextPoint(*p)
        self.pt_pts.Modified()
        self.pt_pd.Modified()

        # polyline
        self.line_pts.Reset()
        self.line_cells.Reset()
        n = len(self.points_3d)
        if n >= 2:
            for p in self.points_3d:
                self.line_pts.InsertNextPoint(*p)
            if self.closed:
                # close the loop
                ids = list(range(n)) + [0]
                self.line_cells.InsertNextCell(n + 1, ids)
                self.line_actor.GetProperty().SetColor(*CLOSED_LINE_COLOR)
            else:
                self.line_cells.InsertNextCell(n, list(range(n)))
                self.line_actor.GetProperty().SetColor(*LINE_COLOR)
        self.line_pts.Modified()
        self.line_cells.Modified()
        self.line_pd.Modified()

    def _build_fill(self):
        """Create a filled polygon actor."""
        if self.fill_actor:
            self.ren.RemoveActor(self.fill_actor)

        pts = vtk.vtkPoints()
        for p in self.points_3d:
            pts.InsertNextPoint(*p)

        poly = vtk.vtkPolygon()
        poly.GetPointIds().SetNumberOfIds(len(self.points_3d))
        for i in range(len(self.points_3d)):
            poly.GetPointIds().SetId(i, i)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(poly)

        pd = vtk.vtkPolyData()
        pd.SetPoints(pts)
        pd.SetPolys(cells)

        # triangulate for proper rendering
        tri = vtk.vtkTriangleFilter()
        tri.SetInputData(pd)
        tri.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tri.GetOutputPort())

        self.fill_actor = vtk.vtkActor()
        self.fill_actor.SetMapper(mapper)
        self.fill_actor.GetProperty().SetColor(*FILL_COLOR)
        self.fill_actor.GetProperty().SetOpacity(0.35)
        self.fill_actor.SetPosition(0, 0, -1)  # behind outline/dots
        self.ren.AddActor(self.fill_actor)


#  interaction style 

class PolygonInteractorStyle(vtk.vtkInteractorStyleImage):
    """Custom mouse/key handler."""

    def __init__(self, drawer: PolygonDrawer, renderer, status_actor):
        super().__init__()
        self.drawer       = drawer
        self.ren          = renderer
        self.status_actor = status_actor
        self._update_status()

        self.AddObserver("LeftButtonPressEvent",   self._on_left)
        self.AddObserver("RightButtonPressEvent",  self._on_right)
        self.AddObserver("MiddleButtonPressEvent", self._on_middle)
        self.AddObserver("KeyPressEvent",          self._on_key)

    #  event callbacks 

    def _on_left(self, obj, event):
        pos = self.GetInteractor().GetEventPosition()
        world = self._display_to_world(pos[0], pos[1])
        self.drawer.add_point(*world)
        self._update_status()
        self.GetInteractor().GetRenderWindow().Render()

    def _on_right(self, obj, event):
        self.drawer.close()
        self._update_status()
        self.GetInteractor().GetRenderWindow().Render()

    def _on_middle(self, obj, event):
        self.drawer.undo_last()
        self._update_status()
        self.GetInteractor().GetRenderWindow().Render()

    def _on_key(self, obj, event):
        key = self.GetInteractor().GetKeySym().lower()
        if key == 'c':
            self.drawer.clear()
            self._update_status()
            self.GetInteractor().GetRenderWindow().Render()
        elif key in ('q', 'escape'):
            self.GetInteractor().GetRenderWindow().Finalize()
            self.GetInteractor().TerminateApp()
        elif key == 's':
            save_to_obj(self.drawer, "polygon.obj")

    #  helpers 

    def _display_to_world(self, x, y):
        """Convert 2-D screen coords → 3-D world coords (z=0 plane)."""
        self.ren.SetDisplayPoint(x, y, 0)
        self.ren.DisplayToWorld()
        wx, wy, wz, ww = self.ren.GetWorldPoint()
        if ww != 0:
            wx, wy, wz = wx / ww, wy / ww, wz / ww
        return wx, wy, 0.0

    def _update_status(self):
        n = len(self.drawer.points_3d)
        if self.drawer.closed:
            msg = f"Polygon closed  ({n} vertices)  |  'c' = clear  |  'q' = quit"
        else:
            msg = (f"Points: {n}  |  Left-click = add  |  Right-click = close  "
                   f"|  Middle-click = undo  |  'c' = clear  |  'q' = quit")
        self.status_actor.SetInput(msg)

def save_to_obj(drawer: PolygonDrawer, filename="polygon.obj"):
    """Save the current polygon to a Wavefront .obj file."""
    if len(drawer.points_3d) < 3:
        print("Not enough points to save (need at least 3).")
        return
    
    with open(filename, 'w') as f:
        f.write("# VTK Polygon Drawer - Exported Polygon\n")
        f.write(f"# Vertices: {len(drawer.points_3d)}\n\n")
        
        # Write vertices
        for x, y, z in drawer.points_3d:
            f.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        
        # Write face (1-indexed)
        indices = " ".join(str(i + 1) for i in range(len(drawer.points_3d)))
        f.write(f"\nf {indices}\n")
    
    print(f"Polygon saved to '{filename}'")

#  main 

def main():
    # renderer
    ren = vtk.vtkRenderer()
    ren.SetBackground(*BG_COLOR)

    # render window
    rw = vtk.vtkRenderWindow()
    rw.SetSize(900, 700)
    rw.SetWindowName("VTK Polygon Drawer")
    rw.AddRenderer(ren)

    # interactor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(rw)

    #  status text 
    status = vtk.vtkTextActor()
    status.GetTextProperty().SetFontSize(14)
    status.GetTextProperty().SetColor(0.9, 0.9, 0.9)
    status.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    status.SetPosition(0.01, 0.01)
    ren.AddActor2D(status)

    #  title 
    title = vtk.vtkTextActor()
    title.SetInput("VTK Polygon Drawer")
    title.GetTextProperty().SetFontSize(20)
    title.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    title.GetTextProperty().BoldOn()
    title.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    title.SetPosition(0.35, 0.95)
    ren.AddActor2D(title)

    #  camera (2-D orthographic) 
    cam = ren.GetActiveCamera()
    cam.ParallelProjectionOn()
    cam.SetPosition(0, 0, 100)
    cam.SetFocalPoint(0, 0, 0)
    cam.SetViewUp(0, 1, 0)
    cam.SetParallelScale(300)

    #  drawer & interaction style 
    drawer = PolygonDrawer(ren)
    style  = PolygonInteractorStyle(drawer, ren, status)
    iren.SetInteractorStyle(style)

    #  grid (optional visual aid) 
    _add_grid(ren)

    rw.Render()
    iren.Start()


def _add_grid(ren, extent=300, step=50):
    """Light background grid to help with spatial orientation."""
    pts   = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    pid   = 0
    for x in range(-extent, extent + 1, step):
        pts.InsertNextPoint(x, -extent, -2)
        pts.InsertNextPoint(x,  extent, -2)
        lines.InsertNextCell(2, [pid, pid + 1])
        pid += 2
    for y in range(-extent, extent + 1, step):
        pts.InsertNextPoint(-extent, y, -2)
        pts.InsertNextPoint( extent, y, -2)
        lines.InsertNextCell(2, [pid, pid + 1])
        pid += 2

    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    pd.SetLines(lines)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(pd)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.25, 0.25, 0.35)
    actor.GetProperty().SetLineWidth(0.5)
    ren.AddActor(actor)


if __name__ == "__main__":
    main()
