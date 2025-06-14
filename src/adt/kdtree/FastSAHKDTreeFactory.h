#pragma once

#include <SAHKDTreeFactory.h>

/**
 * @author Alberto M. Esmoris Pena
 * @version 1.0
 *
 * @brief Class providing building methods for k-dimensional trees with a
 *  fast strategy to approximate Surface Area Heuristic (SAH)
 *
 * The fast strategy to approximate Surface Area Heuristic (SAH) is based on
 *  min-max histograms. Now, let \f$H_*\f$ and \f$H^*\f$ be the min and max
 *  histograms respectively with bins
 *  \f$\left\{[a_1, b_1), \ldots, [a_{m-1}, b_{m-1}), [a_m, b_m]\right\}\f$,
 *  where \f$a_i\f$ is the start point of \f$i\f$-th bin and \f$b_i\f$ is the
 *  end point of \f$i\f$-th bin.
 *  The min histogram is populated from minimum vertices of axis
 *  aligned bounding boxes containing each primitive, while the max histogram
 *  does the same with maximum vertices. For the sake of convenience, let
 *  \f$b_0 = a_1\f$ so the optimum split position can be approximated as
 *  \f$b_i\f$ with \f$i \in [0, m]\f$ such that:
 *
 * \f[
 *  \mathrm{argmin}_{i} \;\;\;\;
 *      \frac{i}{m} \left|\left\{x \in H_*: x < b_i\right\}\right| +
 *      \left(1-\frac{i}{m}\right)
 *          \left|\left\{x \in H^*: x \geq b_i\right\}\right|
 * \f]
 *
 * Thus, by means of modifying how the loss function is computed, sorting is
 *  no longer required and can be replaced by a cheap approximation based on
 *  min-max method.
 *
 * <i>The approach here proposed is based on "Highly parallel fast KD-tree
 *  construction for interactive ray tracing of dynamic scenes" by
 *  Maxim Shevtsov, Alexei Soupikov and Alexander Kapustin</i>
 *
 * @see SAHKDTreeFactory::defineSplit
 */

class FastSAHKDTreeFactory : public SAHKDTreeFactory
{
  using SAHKDTreeFactory::ci;
  using SAHKDTreeFactory::cl;
  using SAHKDTreeFactory::co;
  using SAHKDTreeFactory::lossNodes;

private:
  // ***  SERIALIZATION  *** //
  // *********************** //
  friend class boost::serialization::access;
  /**
   * @brief Serialize a fast surface area heuristic KDTree factory to a
   *  stream of bytes
   * @tparam Archive Type of rendering
   * @param ar Specific rendering for the stream of bytes
   * @param version Version number for the fast surface area heuristic KDTRee
   *  factory
   */
  template<typename Archive>
  void serialize(Archive& ar, unsigned int const version)
  {
    boost::serialization::void_cast_register<FastSAHKDTreeFactory,
                                             SAHKDTreeFactory>();
    ar& boost::serialization::base_object<SAHKDTreeFactory>(*this);
  }

public:
  // ***  CONSTRUCTION / DESTRUCTION  *** //
  // ************************************ //
  /**
   * @brief Fast surface area heuristic KDTree factory default constructor
   * @param lossNodes How many bins use for the min-max approximation of
   *  loss function
   * @see SAHKDTreeFactory::SAHKDTreeFactory
   * @see FastSAHKDTreeFactory
   */
  FastSAHKDTreeFactory(size_t const lossNodes = 32,
                       double const ci = 1,
                       double const cl = 1,
                       double const co = 1);
  ~FastSAHKDTreeFactory() override = default;

  // ***  CLONE  *** //
  // *************** //
  /**
   * @see KDTreeFactory::clone
   */
  KDTreeFactory* clone() const override;
  /**
   * @brief Assign attributes from FastSAHKDTreeFactory to its clone
   */
  void _clone(KDTreeFactory* kdtf) const override;

  // ***  SAH UTILS  *** //
  // ******************* //
  /**
   * @brief Find the best split position using an min-max like approximation
   *  for the loss function of Surface Area Heuristic (SAH)
   * @return Approximate loss of best split position. The position itself is
   *  already stored in given node
   * @see SAHKDTreeFactory::findSplitPositionBySAH
   * @see SAHKDTreeFactory::defineSplit
   */
  double findSplitPositionBySAH(
    KDTreeNode* node,
    std::vector<Primitive*>& primitives) const override;

  /**
   * @brief The recipe for finding split position by Fast-SAH algorithm. It
   *  is meant to be used by the FastSAHKDTreeFactory::findSplitPositionBySAH
   *  but also by any alternative implementation which shares the same recipe
   *  (global logic) but changes the way some parts are computed. For
   *  instance, it is used to handle geometry-level parallelization.
   *
   * @param f_recount Function to count min and max vertices forward and
   *  backward from binning
   *
   * @see FastSAHKDTreeFactory::findSplitPositionBySAH
   * @see FastSAHKDTreeFactory::GEOM_findSplitPositionBySAH
   * @see MultiThreadKDTreeFactory
   *
   * @return Loss of best split position. The position itself is already
   *  stored in given node.
   */
  virtual double findSplitPositionByFastSAHRecipe(
    KDTreeNode* node,
    std::vector<Primitive*>& primitives,
    std::function<void(std::vector<Primitive*>& primitives,
                       int const splitAxis,
                       double const minp,
                       double const deltap,
                       std::size_t const lossNodes,
                       std::size_t const lossCases,
                       std::vector<std::size_t>& cForward,
                       std::vector<std::size_t>& cBackward)> f_recount) const;
};
