#ifndef __pujCGAL__SegmentsIntersection__hxx__
#define __pujCGAL__SegmentsIntersection__hxx__

#include <algorithm>
#include <iterator>

#include <CGAL/intersections.h>

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
void pujCGAL::SegmentsIntersection::
BruteForce( TSegmentsIt sB, TSegmentsIt sE, TPointsIt pIt )
{
  using TPoint = std::iter_value_t< typename TPointsIt::container_type >;

  for( auto i = sB; i != sE; ++i )
    for( auto j = i; j != sE; ++j )
      if( auto r = CGAL::intersection( *i, *j ) )
        if( const TPoint* p = std::get_if< TPoint >( &*r ) )
          *pIt = *p;
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
void pujCGAL::SegmentsIntersection::
BentleyOttmann( TSegmentsIt sB, TSegmentsIt sE, TPointsIt pIt )
{
  using THelpers
    =
    pujCGAL::SegmentsIntersection::
    BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >;
  using TReal    = typename THelpers::TReal;
  using TPoint   = typename THelpers::TPoint;
  using TSegment = typename THelpers::TSegment;
  using TQueue   = typename THelpers::TQueue;
  using TTree    = typename THelpers::TTree;

  // Set global variables
  THelpers::set_input_range( sB, sE );

  // Initialize event queue
  TQueue Q;
  for( auto sIt = sB; sIt != sE; ++sIt )
    THelpers::push_init_event( Q, sIt );

  // Prepare sweep status tree
  TTree T;

  // Main loop
  while( Q.size( ) > 0 )
  {
    auto e = THelpers::pop_event( Q );

    THelpers::move_sweep_line( e.P.x( ) );

    if( e.T == THelpers::EventType::LEFT )
    {
      auto j = T.insert( e.S1 ).first;
      if( j != T.begin( ) )
        THelpers::schedule( *( std::prev( j ) ), *j, Q );

      auto k = std::next( j );
      if( k != T.end( ) )
        THelpers::schedule( *j, *k, Q );
    }
    else if( e.T == THelpers::EventType::RIGHT )
    {
      auto j = T.find( e.S1 );
      if( j != T.end( ) )
      {
        auto i = ( j == T.begin( ) )? T.end( ): std::prev( j );
        auto k = std::next( j );
        if( i != T.end( ) && k != T.end( ) )
          THelpers::schedule( *i, *k, Q );
        T.erase( j );
      } // end if
    }
    else // if( e.T == THelpers::EventType::INTERSECTION )
    {
      // Keep intersection
      *pIt = e.P;

      // in T, swap the order of S1 and S2
      auto i1 = T.find( e.S1 );
      auto i2 = T.find( e.S2 );
      if( i1 != T.end( ) && i2 != T.end( ) )
      {
        T.erase( i1 );
        T.erase( i2 );
        THelpers::move_sweep_line( e.P.x( ) + 1e-9 );
        T.insert( e.S1 );
        T.insert( e.S2 );

        // now re-check neighbors around both segments and
        // schedule intersections with new neighbors
        for( auto i: { T.find( e.S1 ), T.find( e.S2 ) } )
        {
          if( i != T.end( ) )
          {
            if( i != T.begin( ) )
              THelpers::schedule( *( std::prev( i ) ), *i, Q );

            auto j = std::next( i );
            if( j != T.end( ) )
              THelpers::schedule( *i, *j, Q );
          } // end if
        } // end for
      } // end if
    } // end if
  } // end while
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
void pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
set_input_range( TSegmentsIt b, TSegmentsIt e )
{
  Self::SegmentsBegin = b;
  Self::SegmentsEnd = e;
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
void pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
move_sweep_line( const TReal& x )
{
  Self::SweepX = x;
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
typename pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
TSegment
pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
normal( TSegmentsIt i )
{
  if( CGAL::lexicographically_xy_larger( i->source( ), i->target( ) ) )
    return( i->opposite( ) );
  else
    return( *i );
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
typename pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
TReal
pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
y( const TSegment& s )
{
  static const typename TKernel::Direction_2 d( TReal( 0 ), TReal( 1 ) );
  typename TKernel::Line_2 l( TPoint( Self::SweepX, TReal( 0 ) ), d );

  if( auto r = CGAL::intersection( s.supporting_line( ), l ) )
    if( const TPoint* p = std::get_if< TPoint >( &*r ) )
      return( p->y( ) );

  return( s.source( ).y( ) );
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
bool pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::SEvent::
operator<( const SEvent& other ) const
{
  auto c = CGAL::compare_xy( this->P, other.P );
  if( c == CGAL::SMALLER )
    return( true );
  else if( c == CGAL::LARGER )
    return( false );
  else // if( c == CGAL::EQUAL )
    return( other.T < this->T );
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
void pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
push_init_event( TQueue& Q, TSegmentsIt i )
{
  TSegment s = Self::normal( i );

  Q.push_back(
    SEvent( { s.source( ), Self::EventType::LEFT, i, Self::SegmentsEnd } )
    );
  std::push_heap( Q.begin( ), Q.end( ) );

  Q.push_back(
    SEvent( { s.target( ), Self::EventType::RIGHT, i, Self::SegmentsEnd } )
    );
  std::push_heap( Q.begin( ), Q.end( ) );
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
typename pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
SEvent
pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
pop_event( TQueue& Q )
{
  SEvent e = Q.front( );
  std::pop_heap( Q.begin( ), Q.end( ) );
  Q.pop_back( );
  return( e );
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
bool pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::TTreeCmp::
operator()( TSegmentsIt i, TSegmentsIt j ) const
{
  if( i == j )
    return( false );

  TSegment s1 = Self::normal( i );
  TSegment s2 = Self::normal( j );

  TPoint y1( TReal( 0 ), Self::y( s1 ) );
  TPoint y2( TReal( 0 ), Self::y( s2 ) );

  auto c = CGAL::compare_y( y1, y2 );
  if( c == CGAL::EQUAL )
  {
    c = CGAL::compare_xy( s1.source( ), s2.source( ) );
    if( c == CGAL::SMALLER )
      return( true );
    else if( c == CGAL::LARGER )
      return( false );
    else // if( c == CGAL::EQUAL )
    {
      auto di = std::distance( Self::SegmentsBegin, i );
      auto dj = std::distance( Self::SegmentsBegin, j );
      return( di < dj );
    } // end if
  }
  else
    return( c == CGAL::SMALLER );
}

// -------------------------------------------------------------------------
template< class TSegmentsIt, class TPointsIt >
void pujCGAL::SegmentsIntersection::
BentleyOttmann_Helpers< TSegmentsIt, TPointsIt >::
schedule( TSegmentsIt i, TSegmentsIt j, TQueue& Q )
{
  if( i == j )
    return;

  if( auto r = CGAL::intersection( Self::normal( i ), Self::normal( j ) ) )
  {
    if( const TPoint* p = std::get_if< TPoint >( &*r ) )
    {
      if(
        CGAL::compare_x( *p, TPoint( Self::SweepX, TReal( 0 ) ) )
        !=
        CGAL::SMALLER
        )
      {
        Q.push_back( SEvent( { *p, EventType::INTERSECTION, i, j } ) );
        std::push_heap( Q.begin( ), Q.end( ) );
      } // end if
    } // end if
  } // end if
}

#endif // __pujCGAL__SegmentsIntersection__hxx__

// eof - SegmentsIntersection.hxx
