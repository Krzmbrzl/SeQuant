#ifndef SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP
#define SEQUANT_EXTERNAL_INTERFACE_UTILS_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr_fwd.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>

#include <functional>
#include <map>

class IndexSpaceMeta {
public:
	IndexSpaceMeta() = default;

	std::size_t getSize(const sequant::IndexSpace &space) const;
	std::size_t getSize(const sequant::Index &index) const;

	std::function< std::size_t(const sequant::IndexSpace &) > getSizeProxy() const;
	std::function< std::size_t(const sequant::Index &) > getIndexSizeProxy() const;

private:
	std::map< sequant::IndexSpace, std::size_t > m_sizes;
};

sequant::container::svector< sequant::container::svector< sequant::Index > >
	getExternalIndexPairs(const sequant::ExprPtr &expression);

#endif
