#ifndef __pujCGAL__Algoritmos__h__
#define __pujCGAL__Algoritmos__h__

#include <cstddef>
#include <string>

namespace pujCGAL
{
  // =========================================================================
  // Estructura de resultados
  // =========================================================================
  /**
   * Resultado de un algoritmo de simplificacion.
   *
   * Reune metricas de tamano (vertices, triangulos), de tiempo y de
   * calidad (errores L_inf, L2/RMSE, promedio, y angulos). Sirve como
   * formato comun para comparar algoritmos heterogeneos.
   */
  struct ResultadoSimplificacion
  {
    std::string nombre_algoritmo;            ///< etiqueta legible
    std::size_t vertices_antes      { 0 };   ///< pixeles totales del heightmap
    std::size_t vertices_despues    { 0 };   ///< vertices finales en la malla
    std::size_t triangulos_despues  { 0 };   ///< caras finitas finales
    double      reduccion_porcentaje { 0.0 };///< 100 * (1 - despues/antes)
    double      tiempo_milisegundos { 0.0 }; ///< tiempo total del algoritmo

    // Errores de aproximacion (mallapproximada vs heightmap original)
    double      error_L_infinito   { 0.0 };  ///< max |H(p) - T_interp(p)|
    double      error_L2_RMSE      { 0.0 };  ///< sqrt(sum(err^2) / N)
    double      error_promedio     { 0.0 };  ///< sum(|err|) / N

    // Calidad geometrica de la triangulacion
    double      angulo_minimo_grados { 0.0 };
    double      angulo_maximo_grados { 0.0 };
  };

  // =========================================================================
  // Algoritmo 1: Submuestreo uniforme + Delaunay
  // =========================================================================
  /**
   * Construye una triangulacion de Delaunay tomando 1 de cada @p paso
   * pixeles del heightmap. No usa metrica de error: es el baseline trivial.
   *
   * @tparam TDelaunay  CGAL::Delaunay_triangulation_2 con vertice que
   *                    expone info() para almacenar la altura.
   * @tparam THeightmap pujCGAL::Heightmap<TReal>.
   *
   * @param[out] T     triangulacion destino. DEBE entrar vacia: la
   *                   funcion la construye desde cero.
   * @param[in]  hm    mapa de alturas fuente.
   * @param[in]  paso  espaciado en pixeles (>= 1). Aproximadamente
   *                   conserva (W * H) / paso^2 vertices.
   *
   * @return metricas basicas (sin errores ni angulos; calcular_metricas
   *         las completa).
   */
  template< class TDelaunay, class THeightmap >
  ResultadoSimplificacion simplificar_submuestreo_uniforme(
      TDelaunay& T,
      const THeightmap& hm,
      int paso
      );

  // =========================================================================
  // Algoritmo 2: Decimacion con error L1-promedio (algoritmo del taller)
  // =========================================================================
  /**
   * Inserta todos los pixeles del heightmap en la triangulacion y luego
   * elimina iterativamente el vertice de menor "error de planitud":
   *
   *     err(v) = (1/n) * sum_q | H(v) - H(q) |
   *
   * donde q recorre los vecinos de v en la estrella DCEL. Si err(v) <
   * @p epsilon, v es redundante y se elimina. CGAL re-triangula el hueco
   * automaticamente. Despues se recalculan errores en la vecindad de
   * orden @p orden_k mediante face circulator. Lazy deletion en heap.
   *
   * Equivale al algoritmo del enunciado original del taller, extraido
   * aqui como funcion para poder compararlo con otros.
   *
   * @param[in,out] T        triangulacion destino. DEBE entrar vacia: la
   *                         funcion inserta primero todos los pixeles y
   *                         luego decima. Esto incluye el costo de
   *                         construccion en el tiempo medido.
   * @param[in]     hm       mapa de alturas fuente.
   * @param[in]     epsilon  umbral de error de planitud para eliminacion.
   * @param[in]     orden_k  numero de anillos de vecindad para recalculo
   *                         de error tras cada eliminacion (>= 1).
   *
   * @return metricas basicas.
   */
  template< class TDelaunay, class THeightmap >
  ResultadoSimplificacion simplificar_decimacion_L1(
      TDelaunay& T,
      const THeightmap& hm,
      double epsilon,
      int orden_k
      );

  // =========================================================================
  // Algoritmo 3: Insercion golosa de Garland-Heckbert (1995)
  // =========================================================================
  /**
   * Aproxima el heightmap de forma "coarse-to-fine":
   *
   *   1. Inicializa la triangulacion con las 4 esquinas del grid.
   *   2. Asigna cada pixel restante al triangulo que lo contiene.
   *   3. Para cada cara calcula el residuo vertical maximo:
   *         err(f) = max_{p en f} | H(p) - T_interp(p) |
   *      donde T_interp es la interpolacion baricentrica.
   *   4. Saca del max-heap la cara con mayor error e inserta en la
   *      triangulacion el pixel responsable. CGAL re-triangula la
   *      cavidad de Delaunay (Bowyer-Watson) automaticamente.
   *   5. Re-bucketing: los pixeles huerfanos de las caras invalidadas
   *      se reasignan a las nuevas caras incidentes al vertice insertado.
   *      Solo se recalcula el error de las caras nuevas.
   *   6. Termina cuando max_f err(f) <= @p epsilon_max o se alcanza
   *      @p max_vertices vertices.
   *
   * Garantia L_inf vertical bajo @p epsilon_max. Metodo canonico para
   * mapas de altura sobre grilla regular (Garland & Heckbert,
   * "Fast Polygonal Approximation of Terrains and Height Fields",
   * CMU-CS-95-181, 1995).
   *
   * @param[out] T             triangulacion destino. DEBE entrar vacia.
   * @param[in]  hm            mapa de alturas fuente.
   * @param[in]  epsilon_max   umbral L_inf vertical (criterio de parada).
   * @param[in]  max_vertices  cota dura de vertices (0 = sin limite).
   *
   * @return metricas basicas.
   */
  template< class TDelaunay, class THeightmap >
  ResultadoSimplificacion simplificar_insercion_greedy(
      TDelaunay& T,
      const THeightmap& hm,
      double epsilon_max,
      std::size_t max_vertices = 0
      );

  // =========================================================================
  // Calculo de metricas de calidad
  // =========================================================================
  /**
   * Recorre todos los pixeles del heightmap y, para cada uno, localiza
   * el triangulo que lo contiene en @p T y compara la interpolacion
   * baricentrica con la altura real. Reporta L_inf, RMSE y promedio.
   * Adicionalmente recorre las caras finitas y reporta angulos minimo
   * y maximo de la malla.
   *
   * Requiere que @p T tenga info() = altura para cada vertice (los
   * algoritmos de este modulo lo dejan asi al terminar).
   *
   * @param[in,out] resultado  estructura a poblar (los campos de tamano
   *                           y tiempo no se tocan).
   * @param[in]     T          triangulacion a evaluar.
   * @param[in]     hm         heightmap de referencia.
   */
  template< class TDelaunay, class THeightmap >
  void calcular_metricas(
      ResultadoSimplificacion& resultado,
      const TDelaunay& T,
      const THeightmap& hm
      );

} // end namespace pujCGAL

#include <pujCGAL/Algoritmos.hxx>

#endif // __pujCGAL__Algoritmos__h__

// eof - Algoritmos.h
