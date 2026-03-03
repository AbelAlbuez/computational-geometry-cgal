#ifndef __pujCGAL__SegmentsIntersection__h__
#define __pujCGAL__SegmentsIntersection__h__

#include <set>
#include <vector>

namespace pujCGAL
{
  namespace SegmentsIntersection
  {
    template< class TSegmentsIt, class TPointsIt >
    void BruteForce( TSegmentsIt sB, TSegmentsIt sE, TPointsIt pIt );

    template< class TSegmentsIt, class TPointsIt >
    void BentleyOttmann( TSegmentsIt sB, TSegmentsIt sE, TPointsIt pIt );

    /**
     */
    template< class TSegmentsIt, class TPointsIt >
    class BentleyOttmann_Helpers
    {
    public:
      using Self     = BentleyOttmann_Helpers;
      using TSegment = std::iter_value_t< TSegmentsIt >;
      using TKernel  = typename TSegment::R;
      using TReal    = typename TKernel::RT;
      using TPoint   = typename TKernel::Point_2;

    public:

      static inline TSegmentsIt SegmentsBegin;
      static inline TSegmentsIt SegmentsEnd;
      static inline TReal       SweepX;

      static void set_input_range( TSegmentsIt b, TSegmentsIt e );
      static void move_sweep_line( const TReal& x );

      static TSegment normal( TSegmentsIt i );
      static TReal y( const TSegment& s );

      // Event queue types and algorithms
      enum EventType { LEFT = 0, RIGHT = 0, INTERSECTION = 2 };
      struct SEvent
      {
        TPoint P;
        EventType T;
        TSegmentsIt S1, S2;

        bool operator<( const SEvent& other ) const;
      };
      using TQueue = std::vector< SEvent >;

      static void push_init_event( TQueue& Q, TSegmentsIt i );
      static SEvent pop_event( TQueue& Q );

      struct TTreeCmp
      {
        bool operator()( TSegmentsIt i, TSegmentsIt j ) const;
      };
      using TTree = std::set< TSegmentsIt, TTreeCmp >;

      static void schedule( TSegmentsIt i, TSegmentsIt j, TQueue& Q );
    };
  } // end namespace
} // end namespace

#include <pujCGAL/SegmentsIntersection.hxx>

#endif // __pujCGAL__SegmentsIntersection__h__

// eof - SegmentsIntersection.h
