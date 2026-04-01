// =========================================================================
// Taller 2 — Grafo Dual: visualizador VTK (ventana + GIF)
// Reemplaza visualizer.py: mismas capas, colores y controles.
// =========================================================================

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkGlyph3D.h>
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleImage.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkWindowToImageFilter.h>

namespace
{
  namespace fs = std::filesystem;

  std::string
  resolve_output_base( const char* argv0 )
  {
    // -- Intenta resolver relativo al ejecutable (build/../output)
    fs::path exe = fs::weakly_canonical( fs::path( argv0 ) );
    fs::path candidate = exe.parent_path( ) / ".." / "output";
    if( fs::exists( candidate.parent_path( ) ) )
      return( fs::weakly_canonical( candidate ).string( ) );
    // -- Fallback: directorio de trabajo actual / output
    return( ( fs::current_path( ) / "output" ).string( ) );
  }

  using Vec2  = std::pair< double, double >;
  using Face  = std::vector< int >;
  using Edge  = std::pair< int, int >;

  // --- mismos valores que visualizer.py / polygon_drawer.py ----------------
  const double POINT_COLOR[3]       = { 1.0, 0.3, 0.1 };
  const double CLOSED_LINE_COLOR[3] = { 1.0, 0.9, 0.1 };
  const double BG_COLOR[3]          = { 0.12, 0.12, 0.18 };
  const double LINE_WIDTH           = 2.0;

  const double C_TRIANG[3]   = { 0.36, 0.72, 0.96 };
  const double C_DUAL_INT[3] = { 0.29, 0.80, 0.42 };
  const double C_DUAL_EXT[3] = { 0.88, 0.32, 0.32 };
  const double C_BARY[3]     = { 0.29, 0.80, 0.42 };
  const double C_PINF[3]     = { 0.88, 0.32, 0.32 };

  const double Z_POLY   = 0.00;
  const double Z_TRIANG = 0.10;
  const double Z_DUAL   = 0.20;
  const double Z_NODES  = 0.30;

  const double LW_TRIANG = 1.5;
  const double LW_DUAL   = 2.2;
  const double R_BARY    = 7.0;
  const double R_PINF    = 9.0;

  struct ObjData
  {
    std::vector< Vec2 > verts;
    std::vector< Edge > edges;
    std::vector< Face > faces;
  };

  void
  load_from_obj( const char* path, ObjData& out )
  {
    std::ifstream ifs( path );
    if( !ifs )
    {
      std::cerr << "No se pudo abrir: " << path << std::endl;
      std::exit( EXIT_FAILURE );
    }
    std::string line;
    while( std::getline( ifs, line ) )
    {
      if( line.empty( ) )
        continue;
      if( line[ 0 ] == '#' || line[ 0 ] == 'o' )
        continue;
      if( line.size( ) < 3 )
        continue;
      if( line[ 0 ] == 'v' && line[ 1 ] == ' ' )
      {
        std::istringstream iss( line.substr( 2 ) );
        double x, y;
        iss >> x >> y;
        out.verts.push_back( { x, y } );
      }
      else if( line[ 0 ] == 'l' && line[ 1 ] == ' ' )
      {
        std::istringstream iss( line.substr( 2 ) );
        int a, b;
        iss >> a >> b;
        out.edges.push_back( { a - 1, b - 1 } );
      }
      else if( line[ 0 ] == 'f' && line[ 1 ] == ' ' )
      {
        Face f;
        std::istringstream iss( line.substr( 2 ) );
        int id;
        while( iss >> id )
          f.push_back( id - 1 );
        out.faces.push_back( std::move( f ) );
      }
    }
  }

  void
  parse_dual(
    const std::vector< Vec2 >& dual_verts,
    const std::vector< Edge >& dual_edges,
    std::vector< Vec2 >& bary_pts,
    std::vector< Vec2 >& pinf_pt,
    std::vector< Edge >& int_edges,
    std::vector< Edge >& ext_edges,
    int& pinf_idx
    )
  {
    pinf_idx = static_cast< int >( dual_verts.size( ) ) - 1;
    bary_pts.assign( dual_verts.begin( ), dual_verts.begin( ) + pinf_idx );
    pinf_pt = { dual_verts[ pinf_idx ] };
    int_edges.clear( );
    ext_edges.clear( );
    for( const auto& e : dual_edges )
    {
      if( e.first != pinf_idx && e.second != pinf_idx )
        int_edges.push_back( e );
      if( e.first == pinf_idx || e.second == pinf_idx )
        ext_edges.push_back( e );
    }
  }

  std::vector< Edge >
  closed_edges( const Face& face )
  {
    std::vector< Edge > r;
    const int n = static_cast< int >( face.size( ) );
    for( int k = 0; k < n; ++k )
      r.push_back( { face[ k ], face[ ( k + 1 ) % n ] } );
    return( r );
  }

  std::vector< Edge >
  tri_unique_edges( const std::vector< Face >& faces )
  {
    std::vector< Edge > result;
    std::set< std::pair< int, int > > seen;
    for( const auto& tri : faces )
    {
      const int n = static_cast< int >( tri.size( ) );
      for( int k = 0; k < n; ++k )
      {
        int a = tri[ k ];
        int b = tri[ ( k + 1 ) % n ];
        auto k2 = std::make_pair( std::min( a, b ), std::max( a, b ) );
        if( seen.insert( k2 ).second )
          result.push_back( { a, b } );
      }
    }
    return( result );
  }

  vtkSmartPointer< vtkSphereSource >
  make_sphere_source( double radius )
  {
    vtkNew< vtkSphereSource > s;
    s->SetRadius( radius );
    s->SetPhiResolution( 6 );
    s->SetThetaResolution( 8 );
    s->Update( );
    return( s );
  }

  void
  add_grid( vtkRenderer* ren, int extent = 300, int step = 50 )
  {
    vtkNew< vtkPoints > pts;
    vtkNew< vtkCellArray > lines;
    vtkIdType pid = 0;
    for( int x = -extent; x <= extent; x += step )
    {
      pts->InsertNextPoint( x, -extent, -2 );
      pts->InsertNextPoint( x, extent, -2 );
      vtkNew< vtkIdList > ids;
      ids->InsertNextId( pid );
      ids->InsertNextId( pid + 1 );
      lines->InsertNextCell( ids );
      pid += 2;
    }
    for( int y = -extent; y <= extent; y += step )
    {
      pts->InsertNextPoint( -extent, y, -2 );
      pts->InsertNextPoint( extent, y, -2 );
      vtkNew< vtkIdList > ids;
      ids->InsertNextId( pid );
      ids->InsertNextId( pid + 1 );
      lines->InsertNextCell( ids );
      pid += 2;
    }
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetLines( lines );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( 0.25, 0.25, 0.35 );
    actor->GetProperty( )->SetLineWidth( 0.5 );
    ren->AddActor( actor );
  }

  vtkSmartPointer< vtkActor >
  make_glyph_actor(
    const std::vector< Vec2 >& pts_xy,
    const double color[3],
    double radius,
    double z
    )
  {
    vtkNew< vtkPoints > points;
    for( const auto& p : pts_xy )
      points->InsertNextPoint( p.first, p.second, z );
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( points );
    vtkNew< vtkGlyph3D > glyph;
    glyph->SetInputData( pd );
    glyph->SetSourceConnection( make_sphere_source( radius )->GetOutputPort( ) );
    glyph->SetScaleModeToDataScalingOff( );
    glyph->Update( );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputConnection( glyph->GetOutputPort( ) );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    return( actor );
  }

  vtkSmartPointer< vtkActor >
  make_lines_actor(
    const std::vector< Vec2 >& pts_xy,
    const std::vector< Edge >& edge_pairs,
    const double color[3],
    double line_width,
    double z
    )
  {
    vtkNew< vtkPoints > points;
    for( const auto& p : pts_xy )
      points->InsertNextPoint( p.first, p.second, z );
    vtkNew< vtkCellArray > cells;
    for( const auto& e : edge_pairs )
    {
      vtkNew< vtkIdList > ids;
      ids->InsertNextId( e.first );
      ids->InsertNextId( e.second );
      cells->InsertNextCell( ids );
    }
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( points );
    pd->SetLines( cells );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    actor->GetProperty( )->SetLineWidth( static_cast< float >( line_width ) );
    return( actor );
  }

  struct ViewerCtx
  {
    std::map< std::string, std::vector< vtkSmartPointer< vtkActor > > > layers;
    vtkRenderWindow* rw { nullptr };
  };

  void
  keypress_callback( vtkObject* caller, unsigned long, void* clientData, void* )
  {
    auto* style = vtkInteractorStyle::SafeDownCast( caller );
    auto* ctx   = static_cast< ViewerCtx* >( clientData );
    if( !style || !ctx || !ctx->rw )
      return;
    auto* iren = style->GetInteractor( );
    if( !iren )
      return;

    std::string key = iren->GetKeySym( );
    std::transform( key.begin( ), key.end( ), key.begin( ), []( unsigned char c ) {
      return( static_cast< char >( std::tolower( c ) ) );
    } );

    auto toggle = [&]( const char* name ) {
      for( auto& a : ctx->layers[ name ] )
        a->SetVisibility( !a->GetVisibility( ) );
    };

    if( key == "1" )
      toggle( "poly" );
    else if( key == "2" )
      toggle( "triang" );
    else if( key == "3" )
      toggle( "dual" );
    else if( key == "s" )
    {
      vtkNew< vtkWindowToImageFilter > w2i;
      w2i->SetInput( ctx->rw );
      w2i->Update( );
      vtkNew< vtkPNGWriter > writer;
      writer->SetFileName( "screenshot_dual.png" );
      writer->SetInputConnection( w2i->GetOutputPort( ) );
      writer->Write( );
      std::cout << "Screenshot guardado: screenshot_dual.png" << std::endl;
    }
    else if( key == "q" || key == "escape" )
    {
      ctx->rw->Finalize( );
      iren->TerminateApp( );
      return;
    }

    ctx->rw->Render( );
  }

  Face
  poly_contour_from_obj( const ObjData& poly )
  {
    if( !poly.faces.empty( ) )
      return( poly.faces[ 0 ] );
    Face c;
    for( int i = 0; i < static_cast< int >( poly.verts.size( ) ); ++i )
      c.push_back( i );
    return( c );
  }

  void
  setup_renderer_common(
    vtkRenderer* ren,
    vtkCamera* cam
    )
  {
    ren->SetBackground( BG_COLOR[ 0 ], BG_COLOR[ 1 ], BG_COLOR[ 2 ] );
    cam->ParallelProjectionOn( );
    cam->SetPosition( 0, 0, 100 );
    cam->SetFocalPoint( 0, 0, 0 );
    cam->SetViewUp( 0, 1, 0 );
  }

  void
  add_hud( vtkRenderer* ren )
  {
    vtkNew< vtkTextActor > status;
    status->SetInput(
      "1: polygon  |  2: triangulation  |  3: dual  |  s: screenshot  |  q: quit"
    );
    status->GetTextProperty( )->SetFontSize( 14 );
    status->GetTextProperty( )->SetColor( 0.9, 0.9, 0.9 );
    status->GetPositionCoordinate( )->SetCoordinateSystemToNormalizedViewport( );
    status->SetPosition( 0.01, 0.01 );
    ren->AddViewProp( status );

    vtkNew< vtkTextActor > title;
    title->SetInput( "Taller 2 - Grafo Dual de un Poligono Simple" );
    title->GetTextProperty( )->SetFontSize( 18 );
    title->GetTextProperty( )->SetColor( 1.0, 1.0, 1.0 );
    title->GetTextProperty( )->BoldOn( );
    title->GetPositionCoordinate( )->SetCoordinateSystemToNormalizedViewport( );
    title->SetPosition( 0.22, 0.95 );
    ren->AddViewProp( title );
  }

  void
  run_interactive(
    const char* poly_path,
    const char* tri_path,
    const char* dual_path
    )
  {
    ObjData poly, tri, dual;
    load_from_obj( poly_path, poly );
    load_from_obj( tri_path, tri );
    load_from_obj( dual_path, dual );

    Face poly_contour = poly_contour_from_obj( poly );
    std::vector< Vec2 > bary_pts, pinf_pt;
    std::vector< Edge > int_edges, ext_edges;
    int pinf_idx = 0;
    parse_dual( dual.verts, dual.edges, bary_pts, pinf_pt, int_edges, ext_edges, pinf_idx );

    vtkNew< vtkRenderer > ren;
    vtkCamera* cam = ren->GetActiveCamera( );
    setup_renderer_common( ren, cam );

    auto act_poly_line = make_lines_actor(
      poly.verts, closed_edges( poly_contour ), CLOSED_LINE_COLOR, LINE_WIDTH, Z_POLY
    );
    auto act_poly_verts = make_glyph_actor( poly.verts, POINT_COLOR, 6.0, Z_POLY );
    auto act_triang     = make_lines_actor(
      tri.verts, tri_unique_edges( tri.faces ), C_TRIANG, LW_TRIANG, Z_TRIANG
    );
    auto act_dual_int = make_lines_actor( dual.verts, int_edges, C_DUAL_INT, LW_DUAL, Z_DUAL );
    auto act_dual_ext = make_lines_actor( dual.verts, ext_edges, C_DUAL_EXT, LW_DUAL, Z_DUAL );
    auto act_bary     = make_glyph_actor( bary_pts, C_BARY, R_BARY, Z_NODES );
    auto act_pinf     = make_glyph_actor( pinf_pt, C_PINF, R_PINF, Z_NODES );

    ViewerCtx ctx;
    ctx.layers[ "poly" ]   = { act_poly_line, act_poly_verts };
    ctx.layers[ "triang" ] = { act_triang };
    ctx.layers[ "dual" ]   = { act_dual_int, act_dual_ext, act_bary, act_pinf };

    for( const auto& kv : ctx.layers )
      for( const auto& a : kv.second )
        ren->AddActor( a );

    add_hud( ren );
    add_grid( ren );

    ren->ResetCamera( );
    cam->SetParallelScale( cam->GetParallelScale( ) * 1.15 );

    vtkNew< vtkRenderWindow > rw;
    rw->SetSize( 1100, 750 );
    rw->SetWindowName( "Taller 2 - Grafo Dual" );
    rw->SetMultiSamples( 0 );
    rw->AddRenderer( ren );

    ctx.rw = rw;

    vtkNew< vtkRenderWindowInteractor > iren;
    iren->SetRenderWindow( rw );

    vtkNew< vtkInteractorStyleImage > style;
    vtkNew< vtkCallbackCommand > cb;
    cb->SetCallback( keypress_callback );
    cb->SetClientData( &ctx );
    style->AddObserver( vtkCommand::KeyPressEvent, cb );
    iren->SetInteractorStyle( style );

    rw->Render( );
    iren->Start( );
  }

  void
  set_actors_visible(
    const std::vector< vtkSmartPointer< vtkActor > >& actors,
    int vis
    )
  {
    for( auto& a : actors )
      a->SetVisibility( vis );
  }

  void
  capture_frame( vtkRenderWindow* rw, const char* filename )
  {
    rw->Render( );
    vtkNew< vtkWindowToImageFilter > w2i;
    w2i->SetInput( rw );
    w2i->ReadFrontBufferOff( );
    w2i->Update( );
    vtkNew< vtkPNGWriter > writer;
    writer->SetFileName( filename );
    writer->SetInputConnection( w2i->GetOutputPort( ) );
    writer->Write( );
  }

  std::string
  next_result_dir( const char* base )
  {
    fs::path root( base );
    int max_n = 0;
    if( fs::exists( root ) )
    {
      for( const auto& entry : fs::directory_iterator( root ) )
      {
        if( !entry.is_directory( ) )
          continue;
        std::string name = entry.path( ).filename( ).string( );
        if( name.size( ) == 10 && name.substr( 0, 7 ) == "result-" &&
            std::all_of(
              name.begin( ) + 7,
              name.end( ),
              []( unsigned char c ) {
                return( std::isdigit( c ) != 0 );
              }
              ) )
        {
          int n = std::stoi( name.substr( 7 ) );
          if( n > max_n )
            max_n = n;
        }
      }
    }
    char buf[ 32 ];
    std::snprintf( buf, sizeof( buf ), "result-%03d", max_n + 1 );
    fs::path result = root / buf;
    fs::create_directories( result );
    return( result.string( ) );
  }

  void
  run_output(
    const char* argv0,
    const char* poly_path,
    const char* tri_path,
    const char* dual_path,
    bool make_gif
    )
  {
    ObjData poly, tri, dual;
    load_from_obj( poly_path, poly );
    load_from_obj( tri_path, tri );
    load_from_obj( dual_path, dual );

    Face poly_contour = poly_contour_from_obj( poly );
    std::vector< Vec2 > bary_pts, pinf_pt;
    std::vector< Edge > int_edges, ext_edges;
    int pinf_idx = 0;
    parse_dual( dual.verts, dual.edges, bary_pts, pinf_pt, int_edges, ext_edges, pinf_idx );

    vtkNew< vtkRenderer > ren;
    vtkCamera* cam = ren->GetActiveCamera( );
    setup_renderer_common( ren, cam );

    auto act_poly_line = make_lines_actor(
      poly.verts, closed_edges( poly_contour ), CLOSED_LINE_COLOR, LINE_WIDTH, Z_POLY
    );
    auto act_poly_verts = make_glyph_actor( poly.verts, POINT_COLOR, 6.0, Z_POLY );
    auto act_triang     = make_lines_actor(
      tri.verts, tri_unique_edges( tri.faces ), C_TRIANG, LW_TRIANG, Z_TRIANG
    );
    auto act_dual_int = make_lines_actor( dual.verts, int_edges, C_DUAL_INT, LW_DUAL, Z_DUAL );
    auto act_dual_ext = make_lines_actor( dual.verts, ext_edges, C_DUAL_EXT, LW_DUAL, Z_DUAL );
    auto act_bary     = make_glyph_actor( bary_pts, C_BARY, R_BARY, Z_NODES );
    auto act_pinf     = make_glyph_actor( pinf_pt, C_PINF, R_PINF, Z_NODES );

    std::vector< vtkSmartPointer< vtkActor > > poly_actors   = { act_poly_line, act_poly_verts };
    std::vector< vtkSmartPointer< vtkActor > > tri_actors    = { act_triang };
    std::vector< vtkSmartPointer< vtkActor > > dual_actors   = {
      act_dual_int, act_dual_ext, act_bary, act_pinf
    };

    for( auto& a : poly_actors )
      ren->AddActor( a );
    for( auto& a : tri_actors )
      ren->AddActor( a );
    for( auto& a : dual_actors )
      ren->AddActor( a );

    add_grid( ren );

    ren->ResetCamera( );
    cam->SetParallelScale( cam->GetParallelScale( ) * 1.15 );

    vtkNew< vtkRenderWindow > rw;
    rw->SetSize( 1100, 750 );
    rw->SetOffScreenRendering( 1 );
    rw->AddRenderer( ren );

    auto hide_all_dual = [&]() {
      set_actors_visible( dual_actors, 0 );
    };

    std::string output_base = resolve_output_base( argv0 );
    std::string out_dir = next_result_dir( output_base.c_str( ) );
    std::cout << "Output: " << out_dir << "\n";

    auto out_png = [&]( const char* leaf ) {
      return( ( fs::path( out_dir ) / leaf ).string( ) );
    };

    set_actors_visible( poly_actors, 1 );
    set_actors_visible( tri_actors, 0 );
    hide_all_dual( );
    {
      std::string p = out_png( "00_polygon.png" );
      capture_frame( rw, p.c_str( ) );
      std::cout << "  " << fs::path( p ).filename( ).string( ) << "\n";
    }

    set_actors_visible( tri_actors, 1 );
    {
      std::string p = out_png( "01_triangulation.png" );
      capture_frame( rw, p.c_str( ) );
      std::cout << "  " << fs::path( p ).filename( ).string( ) << "\n";
    }

    act_dual_int->SetVisibility( 0 );
    act_dual_ext->SetVisibility( 0 );
    act_bary->SetVisibility( 1 );
    act_pinf->SetVisibility( 1 );
    {
      std::string p = out_png( "02_barycenters.png" );
      capture_frame( rw, p.c_str( ) );
      std::cout << "  " << fs::path( p ).filename( ).string( ) << "\n";
    }

    act_dual_int->SetVisibility( 1 );
    act_dual_ext->SetVisibility( 0 );
    {
      std::string p = out_png( "03_internal_edges.png" );
      capture_frame( rw, p.c_str( ) );
      std::cout << "  " << fs::path( p ).filename( ).string( ) << "\n";
    }

    act_dual_ext->SetVisibility( 1 );
    {
      std::string p = out_png( "04_dual_complete.png" );
      capture_frame( rw, p.c_str( ) );
      std::cout << "  " << fs::path( p ).filename( ).string( ) << "\n";
    }

    if( make_gif )
    {
      fs::path concat_path = fs::path( out_dir ) / "ffmpeg_concat.txt";
      {
        std::ofstream cf( concat_path );
        cf << "ffconcat version 1.0\n";
        const char* leaves[] = {
          "00_polygon.png",
          "01_triangulation.png",
          "02_barycenters.png",
          "03_internal_edges.png",
          "04_dual_complete.png",
        };
        auto write_file_line = [&]( const std::string& ap ) {
          cf << "file '";
          for( char c : ap )
          {
            if( c == '\'' )
              cf << "'\\''";
            else
              cf << c;
          }
          cf << "'\n";
        };
        for( const char* leaf : leaves )
        {
          std::string ap = fs::absolute( fs::path( out_dir ) / leaf ).string( );
          write_file_line( ap );
          cf << "duration 0.5\n";
        }
        write_file_line(
          fs::absolute( fs::path( out_dir ) / leaves[ 4 ] ).string( )
        );
      }

      fs::path gif_path = fs::path( out_dir ) / "dual_graph.gif";
      std::string cmd =
        std::string( "ffmpeg -y -loglevel error -f concat -safe 0 -i \"" ) +
        concat_path.string( ) +
        "\" -vf \"fps=10,scale=1100:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];"
        "[s1][p]paletteuse\" \"" +
        gif_path.string( ) + "\"";
      int ret = std::system( cmd.c_str( ) );
      if( ret != 0 )
      {
        std::cerr << "ffmpeg fallo (codigo " << ret << "). PNGs conservados en "
                  << out_dir << "\n";
      }
      else
        std::cout << "GIF: " << gif_path.string( ) << "\n";
    }
  }

} // namespace

int
main( int argc, char** argv )
{
  if( argc < 4 || argc > 5 )
  {
    std::cerr << "Uso: " << argv[ 0 ]
              << "  polygon.obj  triangulation.obj  dual.obj  [--gif | --interactive]\n"
              << "  (por defecto genera PNGs en output/result-NNN/)\n";
    return( EXIT_FAILURE );
  }

  bool interactive = false;
  bool make_gif    = false;

  if( argc == 5 )
  {
    std::string flag( argv[ 4 ] );
    if( flag == "--interactive" )
      interactive = true;
    else if( flag == "--gif" )
      make_gif = true;
    else
    {
      std::cerr << "Flag desconocido. Uso: [--gif | --interactive]\n";
      return( EXIT_FAILURE );
    }
  }

  if( interactive )
    run_interactive( argv[ 1 ], argv[ 2 ], argv[ 3 ] );
  else
    run_output( argv[ 0 ], argv[ 1 ], argv[ 2 ], argv[ 3 ], make_gif );

  return( EXIT_SUCCESS );
}

// eof - visualizer.cxx
