// =========================================================================
// VTK off-screen: 6 pasos del algoritmo de concavidades + GIF (ffmpeg).
// =========================================================================

#include "concavidades_viz.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkGlyph3D.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <map>
#include <set>
#include <vtkTriangle.h>
#include <vtkWindowToImageFilter.h>

namespace
{
  namespace fs = std::filesystem;

  using Vec2 = std::pair< double, double >;
  using Edge = std::pair< int, int >;

  std::string
  resolve_output_base( const char* argv0 )
  {
    fs::path exe = fs::weakly_canonical( fs::path( argv0 ) );
    fs::path candidate = exe.parent_path( ) / ".." / "output";
    if( fs::exists( candidate.parent_path( ) ) )
      return( fs::weakly_canonical( candidate ).string( ) );
    return( ( fs::current_path( ) / "output" ).string( ) );
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
        if(
          name.size( ) == 10 && name.substr( 0, 7 ) == "result-" &&
          std::all_of(
            name.begin( ) + 7,
            name.end( ),
            []( unsigned char c ) {
              return( std::isdigit( c ) != 0 );
            }
            )
          )
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
    vtkNew< vtkSphereSource > sph;
    sph->SetRadius( radius );
    sph->SetPhiResolution( 6 );
    sph->SetThetaResolution( 8 );
    sph->Update( );
    vtkNew< vtkGlyph3D > glyph;
    glyph->SetInputData( pd );
    glyph->SetSourceConnection( sph->GetOutputPort( ) );
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
  make_polyline_actor(
    const std::vector< Vec2 >& ring,
    const double color[3],
    double width,
    double z,
    bool closed,
    double opacity = 1.0
    )
  {
    vtkNew< vtkPoints > pts;
    for( const auto& p : ring )
      pts->InsertNextPoint( p.first, p.second, z );
    vtkNew< vtkCellArray > lines;
    if( closed && ring.size( ) >= 2 )
    {
      vtkNew< vtkIdList > poly;
      for( vtkIdType i = 0; i < static_cast< vtkIdType >( ring.size( ) ); ++i )
        poly->InsertNextId( i );
      poly->InsertNextId( 0 );
      lines->InsertNextCell( poly );
    }
    else
    {
      for( vtkIdType i = 0; i + 1 < static_cast< vtkIdType >( ring.size( ) ); ++i )
      {
        vtkNew< vtkIdList > seg;
        seg->InsertNextId( i );
        seg->InsertNextId( i + 1 );
        lines->InsertNextCell( seg );
      }
    }
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetLines( lines );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    actor->GetProperty( )->SetLineWidth( static_cast< float >( width ) );
    actor->GetProperty( )->SetOpacity( opacity );
    return( actor );
  }

  vtkSmartPointer< vtkActor >
  make_segment_actor(
    const Vec2& a,
    const Vec2& b,
    const double color[3],
    double width,
    double z,
    double opacity = 1.0
    )
  {
    std::vector< Vec2 > v = { a, b };
    return( make_polyline_actor( v, color, width, z, false, opacity ) );
  }

  std::map< Edge, int >
  edge_triangle_count( const ConcavidadVizInput& in )
  {
    std::map< Edge, int > cnt;
    auto bump = [&]( int a, int b ) {
      if( a == b )
        return;
      Edge e = { std::min( a, b ), std::max( a, b ) };
      cnt[ e ]++;
    };
    for( const auto& t : in.mesh_tris )
    {
      int a = static_cast< int >( t[ 0 ] );
      int b = static_cast< int >( t[ 1 ] );
      int c = static_cast< int >( t[ 2 ] );
      bump( a, b );
      bump( b, c );
      bump( a, c );
    }
    return( cnt );
  }

  vtkSmartPointer< vtkActor >
  make_triangulation_diagonals_actor(
    const ConcavidadVizInput& in,
    const double color[3],
    double width,
    double z
    )
  {
    std::map< Edge, int > cnt = edge_triangle_count( in );
    vtkNew< vtkPoints > pts;
    vtkNew< vtkCellArray > lines;
    vtkIdType nid = 0;
    for( const auto& kv : cnt )
    {
      if( kv.second != 2 )
        continue;
      const Edge& e = kv.first;
      Vec2 pa = in.mesh_points[ static_cast< std::size_t >( e.first ) ];
      Vec2 pb = in.mesh_points[ static_cast< std::size_t >( e.second ) ];
      pts->InsertNextPoint( pa.first, pa.second, z );
      pts->InsertNextPoint( pb.first, pb.second, z );
      vtkNew< vtkIdList > seg;
      seg->InsertNextId( nid );
      seg->InsertNextId( nid + 1 );
      lines->InsertNextCell( seg );
      nid += 2;
    }
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetLines( lines );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    actor->GetProperty( )->SetLineWidth( static_cast< float >( width ) );
    return( actor );
  }

  vtkSmartPointer< vtkActor >
  make_dashed_polyline_actor(
    const std::vector< Vec2 >& ring,
    const double color[3],
    double width,
    double z,
    bool closed,
    double opacity,
    int dash_segments,
    int gap_segments
    )
  {
    if( ring.size( ) < 2 )
      return( nullptr );
    auto edge_point = [&]( std::size_t ia, std::size_t ib, double t ) -> Vec2 {
      const Vec2& pa = ring[ ia ];
      const Vec2& pb = ring[ ib ];
      return(
        Vec2(
          pa.first + t * ( pb.first - pa.first ),
          pa.second + t * ( pb.second - pa.second )
          )
        );
    };
    vtkNew< vtkPoints > pts;
    vtkNew< vtkCellArray > lines;
    vtkIdType nid = 0;
    const std::size_t nseg = closed ? ring.size( ) : ring.size( ) - 1;
    for( std::size_t s = 0; s < nseg; ++s )
    {
      const std::size_t ia = s;
      const std::size_t ib = closed ? ( ( s + 1 ) % ring.size( ) ) : ( s + 1 );
      const Vec2& pa = ring[ ia ];
      const Vec2& pb = ring[ ib ];
      double len = std::hypot( pb.first - pa.first, pb.second - pa.second );
      if( len < 1e-18 )
        continue;
      const int total = dash_segments + gap_segments;
      if( total <= 0 )
        continue;
      const int steps = static_cast< int >( std::ceil( len / 2.5 ) );
      for( int k = 0; k < steps; ++k )
      {
        const int phase = k % total;
        if( phase >= dash_segments )
          continue;
        const double t0 = static_cast< double >( k ) / static_cast< double >( steps );
        const double t1 = static_cast< double >( k + 1 ) / static_cast< double >( steps );
        Vec2 q0 = edge_point( ia, ib, t0 );
        Vec2 q1 = edge_point( ia, ib, t1 );
        pts->InsertNextPoint( q0.first, q0.second, z );
        pts->InsertNextPoint( q1.first, q1.second, z );
        vtkNew< vtkIdList > seg;
        seg->InsertNextId( nid );
        seg->InsertNextId( nid + 1 );
        lines->InsertNextCell( seg );
        nid += 2;
      } // end for k
    } // end for s
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetLines( lines );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    actor->GetProperty( )->SetLineWidth( static_cast< float >( width ) );
    actor->GetProperty( )->SetOpacity( opacity );
    return( actor );
  }

  vtkSmartPointer< vtkActor >
  make_dual_lines(
    const ConcavidadVizInput& in,
    const std::vector< std::pair< std::size_t, std::size_t > >& edges,
    const double color[3],
    double width,
    double z,
    bool to_pinf,
    double opacity = 1.0
    )
  {
    vtkNew< vtkPoints > pts;
    vtkNew< vtkCellArray > lines;
    vtkIdType nid = 0;
    const std::size_t ntri = in.barycenters.size( );
    for( const auto& e : edges )
    {
      std::size_t i = e.first;
      std::size_t j = e.second;
      Vec2 pa = in.barycenters[ i ];
      Vec2 pb = to_pinf && j >= ntri ? in.p_inf : in.barycenters[ j ];
      vtkIdType ia = pts->InsertNextPoint( pa.first, pa.second, z );
      vtkIdType ib = pts->InsertNextPoint( pb.first, pb.second, z );
      (void)ia;
      (void)ib;
      vtkNew< vtkIdList > seg;
      seg->InsertNextId( nid );
      seg->InsertNextId( nid + 1 );
      lines->InsertNextCell( seg );
      nid += 2;
    }
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetLines( lines );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    actor->GetProperty( )->SetLineWidth( static_cast< float >( width ) );
    actor->GetProperty( )->SetOpacity( opacity );
    return( actor );
  }

  vtkSmartPointer< vtkActor >
  make_dual_ext_dashed(
    const ConcavidadVizInput& in,
    const double color[3],
    double width,
    double z,
    double opacity,
    int dash_segments,
    int gap_segments
    )
  {
    vtkNew< vtkPoints > pts;
    vtkNew< vtkCellArray > lines;
    vtkIdType nid = 0;
    const std::size_t ntri = in.barycenters.size( );
    for( const auto& e : in.dual_ext_edges )
    {
      std::size_t i = e.first;
      std::size_t j = e.second;
      Vec2 pa = in.barycenters[ i ];
      Vec2 pb = j >= ntri ? in.p_inf : in.barycenters[ j ];
      double len = std::hypot( pb.first - pa.first, pb.second - pa.second );
      if( len < 1e-18 )
        continue;
      const int total = dash_segments + gap_segments;
      if( total <= 0 )
        continue;
      const int steps = static_cast< int >( std::ceil( len / 2.5 ) );
      for( int k = 0; k < steps; ++k )
      {
        const int phase = k % total;
        if( phase >= dash_segments )
          continue;
        const double t0 = static_cast< double >( k ) / static_cast< double >( steps );
        const double t1 = static_cast< double >( k + 1 ) / static_cast< double >( steps );
        double x0 = pa.first + t0 * ( pb.first - pa.first );
        double y0 = pa.second + t0 * ( pb.second - pa.second );
        double x1 = pa.first + t1 * ( pb.first - pa.first );
        double y1 = pa.second + t1 * ( pb.second - pa.second );
        pts->InsertNextPoint( x0, y0, z );
        pts->InsertNextPoint( x1, y1, z );
        vtkNew< vtkIdList > seg;
        seg->InsertNextId( nid );
        seg->InsertNextId( nid + 1 );
        lines->InsertNextCell( seg );
        nid += 2;
      } // end for k
    } // end for e
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetLines( lines );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    actor->GetProperty( )->SetColor( color[ 0 ], color[ 1 ], color[ 2 ] );
    actor->GetProperty( )->SetLineWidth( static_cast< float >( width ) );
    actor->GetProperty( )->SetOpacity( opacity );
    return( actor );
  }

  void
  fit_cam_parallel_barycenters_only(
    vtkRenderer* ren,
    const ConcavidadVizInput& in
    )
  {
    if( in.barycenters.empty( ) )
    {
      ren->ResetCamera( );
      ren->ResetCameraClippingRange( );
      vtkCamera* cam = ren->GetActiveCamera( );
      cam->ParallelProjectionOn( );
      cam->SetParallelScale( cam->GetParallelScale( ) * 1.15 );
      return;
    }
    double xmin = in.barycenters[ 0 ].first;
    double xmax = xmin;
    double ymin = in.barycenters[ 0 ].second;
    double ymax = ymin;
    for( const auto& p : in.barycenters )
    {
      xmin = std::min( xmin, p.first );
      xmax = std::max( xmax, p.first );
      ymin = std::min( ymin, p.second );
      ymax = std::max( ymax, p.second );
    } // end for
    const double cx   = 0.5 * ( xmin + xmax );
    const double cy   = 0.5 * ( ymin + ymax );
    const double dx   = xmax - xmin;
    const double dy   = ymax - ymin;
    double span = std::max( dx, dy );
    if( span < 1e-12 )
      span = 1.0;
    vtkCamera* cam = ren->GetActiveCamera( );
    cam->ParallelProjectionOn( );
    cam->SetFocalPoint( cx, cy, 0.0 );
    cam->SetPosition( cx, cy, 1.0 );
    cam->SetViewUp( 0, 1, 0 );
    cam->SetParallelScale( 0.5 * span * 1.15 );
    ren->ResetCameraClippingRange( );
  }

  void
  add_title( vtkRenderer* ren, const char* text, double y = 0.92 )
  {
    vtkNew< vtkTextActor > ta;
    ta->SetInput( text );
    ta->GetTextProperty( )->SetFontSize( 16 );
    ta->GetTextProperty( )->SetColor( 1, 1, 1 );
    ta->GetTextProperty( )->SetJustificationToCentered( );
    ta->GetTextProperty( )->SetVerticalJustificationToTop( );
    ta->GetPositionCoordinate( )->SetCoordinateSystemToNormalizedViewport( );
    ta->SetPosition( 0.5, y );
    ren->AddViewProp( ta );
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

  static const double C_BG[3]     = { 0.11, 0.11, 0.16 };
  static const double C_CLOUD[3]  = { 1.0, 0.35, 0.12 };
  static const double C_CH[3]     = { 0.25, 0.55, 0.95 };
  static const double C_ADJ[3]    = { 1.0, 0.92, 0.35 };
  static const double C_TRI[3]    = { 0.25, 0.85, 0.45 };
  static const double C_DUAL_I[3] = { 0.2, 0.75, 0.35 };
  static const double C_DUAL_E[3] = { 0.95, 0.25, 0.25 };
  static const double C_TNODE[3]  = { 0.3, 0.5, 1.0 };
  static const double C_BNODE[3]  = { 1.0, 0.55, 0.15 };
  static const double C_PINF[3]   = { 0.95, 0.2, 0.2 };

  static const double PALETTE[][3] = {
    { 0.9, 0.3, 0.35 },
    { 0.35, 0.75, 0.95 },
    { 0.55, 0.9, 0.35 },
    { 0.85, 0.45, 0.9 },
    { 0.95, 0.85, 0.25 },
    { 0.5, 0.85, 0.85 },
    { 0.95, 0.5, 0.35 },
    { 0.45, 0.45, 0.95 },
  };

} // namespace

void
run_concavidades_viz( const char* argv0, const ConcavidadVizInput& in )
{
  vtkNew< vtkRenderer > ren;
  ren->SetBackground( C_BG[ 0 ], C_BG[ 1 ], C_BG[ 2 ] );
  vtkCamera* cam = ren->GetActiveCamera( );
  cam->ParallelProjectionOn( );

  auto act_cloud = make_glyph_actor( in.cloud, C_CLOUD, 6.0, 0.35 );
  auto act_ch    = make_polyline_actor( in.ch_ccw, C_CH, 2.5, 0.05, true, 1.0 );
  auto act_ch_faint =
    make_polyline_actor( in.ch_ccw, C_CH, 2.0, 0.05, true, 0.35 );
  auto act_ch_step2_dashed = make_dashed_polyline_actor(
    in.ch_ccw, C_CH, 2.0, 0.05, true, 0.45, 2, 2
    );
  auto act_adj = make_polyline_actor( in.adj_ccw, C_ADJ, 2.8, 0.12, true, 1.0 );
  auto act_adj_faint =
    make_polyline_actor( in.adj_ccw, C_ADJ, 2.4, 0.10, true, 0.4 );
  auto act_tri_green = make_triangulation_diagonals_actor(
    in, C_TRI, 1.4, 0.18
    );
  auto act_dual_int =
    make_dual_lines( in, in.dual_int_edges, C_DUAL_I, 2.0, 0.22, false, 1.0 );
  auto act_dual_ext = make_dual_ext_dashed(
    in, C_DUAL_E, 1.6, 0.22, 0.9, 2, 2
    );

  std::vector< vtkSmartPointer< vtkActor > > pocket_actors;
  const int npal = static_cast< int >( sizeof( PALETTE ) / sizeof( PALETTE[ 0 ] ) );
  for( std::size_t ti = 0; ti < in.mesh_tris.size( ); ++ti )
  {
    int comp = in.pocket_component[ ti ];
    if( comp < 0 )
      continue;
    const auto& t = in.mesh_tris[ ti ];
    const Vec2& p0 = in.mesh_points[ t[ 0 ] ];
    const Vec2& p1 = in.mesh_points[ t[ 1 ] ];
    const Vec2& p2 = in.mesh_points[ t[ 2 ] ];
    vtkNew< vtkPoints > pts;
    pts->InsertNextPoint( p0.first, p0.second, 0.14 );
    pts->InsertNextPoint( p1.first, p1.second, 0.14 );
    pts->InsertNextPoint( p2.first, p2.second, 0.14 );
    vtkNew< vtkTriangle > tri;
    tri->GetPointIds( )->SetId( 0, 0 );
    tri->GetPointIds( )->SetId( 1, 1 );
    tri->GetPointIds( )->SetId( 2, 2 );
    vtkNew< vtkCellArray > cells;
    cells->InsertNextCell( tri );
    vtkNew< vtkPolyData > pd;
    pd->SetPoints( pts );
    pd->SetPolys( cells );
    vtkNew< vtkPolyDataMapper > mapper;
    mapper->SetInputData( pd );
    vtkNew< vtkActor > actor;
    actor->SetMapper( mapper );
    const double* col = PALETTE[ comp % npal ];
    actor->GetProperty( )->SetColor( col[ 0 ], col[ 1 ], col[ 2 ] );
    actor->GetProperty( )->SetOpacity( 0.55 );
    pocket_actors.push_back( actor );
  }

  std::vector< std::pair< double, double > > bary_T, bary_B;
  for( std::size_t i = 0; i < in.barycenters.size( ); ++i )
  {
    if( in.is_bolsillo[ i ] )
      bary_B.push_back( in.barycenters[ i ] );
    else
      bary_T.push_back( in.barycenters[ i ] );
  }
  auto act_bary_T = make_glyph_actor( bary_T, C_TNODE, 5.5, 0.28 );
  auto act_bary_B = make_glyph_actor( bary_B, C_BNODE, 6.5, 0.28 );
  std::vector< Vec2 > pinf_v = { in.p_inf };
  auto act_pinf = make_glyph_actor( pinf_v, C_PINF, 7.0, 0.28 );

  std::string texto_resultado;
  {
    std::ostringstream oss;
    oss << "Resultado: " << in.num_concavities << " concavidad(es)";
    texto_resultado = oss.str( );
  }
  vtkNew< vtkTextActor > big_result;
  big_result->SetInput( texto_resultado.c_str( ) );
  big_result->GetTextProperty( )->SetFontSize( 22 );
  big_result->GetTextProperty( )->SetBold( true );
  big_result->GetTextProperty( )->SetColor( 1, 1, 0.85 );
  big_result->GetTextProperty( )->SetJustificationToCentered( );
  big_result->GetPositionCoordinate( )->SetCoordinateSystemToNormalizedViewport( );
  big_result->SetPosition( 0.5, 0.42 );

  vtkNew< vtkRenderWindow > rw;
  rw->SetSize( 1100, 750 );
  rw->SetOffScreenRendering( 1 );
  rw->AddRenderer( ren );

  auto clear = [&]() {
    ren->RemoveAllViewProps( );
  };

  // -- Igual que visualizer.cxx: solo geometria 3D antes de ResetCamera.
  //    Los vtkTextActor en viewport corrompen los bounds si se anaden antes.
  auto fit_cam = [&]() {
    ren->ResetCamera( );
    ren->ResetCameraClippingRange( );
    cam->SetParallelScale( cam->GetParallelScale( ) * 1.15 );
  };

  std::string output_base = resolve_output_base( argv0 );
  std::string out_dir     = next_result_dir( output_base.c_str( ) );
  std::cout << "Output: " << out_dir << "\n";

  auto out_png = [&]( const char* leaf ) -> std::string {
    return( ( fs::path( out_dir ) / leaf ).string( ) );
  };

  auto shot = [&]( const char* path ) {
    capture_frame( rw, path );
    std::cout << "  " << fs::path( path ).filename( ).string( ) << "\n";
  };

  // -- 00 entrada
  clear( );
  ren->AddActor( act_cloud );
  fit_cam( );
  add_title( ren, "Entrada: nube de puntos 2D" );
  shot( out_png( "00_entrada.png" ).c_str( ) );

  // -- 01 CH
  clear( );
  ren->AddActor( act_cloud );
  ren->AddActor( act_ch );
  fit_cam( );
  add_title( ren, "Paso 1: casco convexo (CH)" );
  shot( out_png( "01_casco_convexo.png" ).c_str( ) );

  // -- 02 ajustado
  clear( );
  ren->AddActor( act_cloud );
  ren->AddActor( act_ch_step2_dashed );
  ren->AddActor( act_adj );
  fit_cam( );
  add_title( ren, "Paso 2: polígono ajustado (CCW)" );
  shot( out_png( "02_poligono_ajustado.png" ).c_str( ) );

  // -- 03 triangulacion
  clear( );
  ren->AddActor( act_cloud );
  ren->AddActor( act_adj );
  ren->AddActor( act_tri_green );
  fit_cam( );
  add_title( ren, "Paso 3: triangulación (n-2 triángulos)" );
  shot( out_png( "03_triangulacion.png" ).c_str( ) );

  // -- 04 dual (camara solo baricentros; P puede quedar fuera del encuadre)
  clear( );
  ren->AddActor( act_cloud );
  ren->AddActor( act_adj_faint );
  ren->AddActor( act_dual_int );
  ren->AddActor( act_dual_ext );
  ren->AddActor( act_bary_T );
  ren->AddActor( act_bary_B );
  ren->AddActor( act_pinf );
  fit_cam_parallel_barycenters_only( ren, in );
  add_title( ren, "Paso 4: grafo dual (T=interior, B=bolsillo)" );
  shot( out_png( "04_grafo_dual.png" ).c_str( ) );

  // -- 05 resultado
  clear( );
  ren->AddActor( act_cloud );
  ren->AddActor( act_ch_faint );
  ren->AddActor( act_adj );
  for( auto& a : pocket_actors )
    ren->AddActor( a );
  fit_cam( );
  add_title( ren, "Paso 5: concavidades detectadas", 0.88 );
  ren->AddViewProp( big_result );
  shot( out_png( "05_resultado.png" ).c_str( ) );

  // -- GIF ffmpeg 1s por frame
  fs::path concat_path = fs::path( out_dir ) / "ffmpeg_concat.txt";
  {
    std::ofstream cf( concat_path );
    cf << "ffconcat version 1.0\n";
    const char* leaves[] = {
      "00_entrada.png",
      "01_casco_convexo.png",
      "02_poligono_ajustado.png",
      "03_triangulacion.png",
      "04_grafo_dual.png",
      "05_resultado.png",
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
      cf << "duration 1\n";
    }
    write_file_line(
      fs::absolute( fs::path( out_dir ) / leaves[ 5 ] ).string( )
      );
  }

  fs::path gif_path = fs::path( out_dir ) / "concavidades.gif";
  std::string cmd =
    std::string( "ffmpeg -y -loglevel error -f concat -safe 0 -i \"" ) +
    concat_path.string( ) +
    "\" -vf \"fps=1,scale=1100:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];"
    "[s1][p]paletteuse\" \"" +
    gif_path.string( ) + "\"";
  int ret = std::system( cmd.c_str( ) );
  if( ret != 0 )
    std::cerr << "ffmpeg fallo (codigo " << ret << "). PNGs en " << out_dir
              << "\n";
  else
    std::cout << "GIF: " << gif_path.string( ) << "\n";
}
