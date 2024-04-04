#include "utils.hpp"

#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/utility/string.hpp>

#include <functional>

using namespace sequant;

std::size_t IndexSpaceMeta::getSize(const IndexSpace &space) const {
	auto iter = m_sizes.find(space);

	if (iter == m_sizes.end()) {
		throw std::runtime_error("No known size for index space with label '" + toUtf8(to_wstring(space.type())) + "'");
	}

	return iter->second;
}

std::size_t IndexSpaceMeta::getSize(const Index &index) const {
	return getSize(index.space());
}

std::function< std::size_t(const IndexSpace &) > IndexSpaceMeta::getSizeProxy() const {
	return [this](const IndexSpace &space) { return this->getSize(space); };
}

std::function< std::size_t(const Index &) > IndexSpaceMeta::getIndexSizeProxy() const {
	return [this](const Index &index) { return this->getSize(index); };
}


container::svector< container::svector< Index > > getExternalIndexPairs(const ExprPtr &expression) {
	// TODO
	return {};
}
