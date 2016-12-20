//
// Created by iskakoff on 20/12/16.
//

#ifndef GFTOOLS_INDEX_GRID_HPP
#define GFTOOLS_INDEX_GRID_HPP

#include "grid_base.hpp"
#include "enum_grid.hpp"

namespace gftools {
/**
 * @brief this class implements integer index grid. One of the purpose is to represent the spin of Green's function
 *
 * @author iskakoff
 */
  class index_grid : public enum_grid {
  public:
    index_grid(int N): enum_grid(0, N){};
  };
}

#endif //GFTOOLS_INDEX_GRID_HPP
