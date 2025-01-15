#ifndef SEQUANT_VERTEX_PAINTER_H
#define SEQUANT_VERTEX_PAINTER_H

#include <SeQuant/core/abstract_tensor.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/tensor_network_v2.hpp>

#include <utility>
#include <variant>

namespace sequant {

using ProtoBundle =
    std::decay_t<decltype(std::declval<const Index &>().proto_indices())>;

struct BraGroup {
  explicit BraGroup(std::size_t id) : id(id) {}

  std::size_t id;
};
struct KetGroup {
  explicit KetGroup(std::size_t id) : id(id) {}

  std::size_t id;
};
struct AuxGroup {
  explicit AuxGroup(std::size_t id) : id(id) {}

  std::size_t id;
};
struct ParticleGroup {
  explicit ParticleGroup(std::size_t id) : id(id) {}

  std::size_t id;
};

class VertexPainter {
 public:
  using Color = TensorNetworkV2::Graph::VertexColor;
  using VertexData =
      std::variant<const AbstractTensor *, Index, const ProtoBundle *, BraGroup,
                   KetGroup, AuxGroup, ParticleGroup>;
  using ColorMap = container::map<Color, VertexData>;

  VertexPainter(const TensorNetworkV2::NamedIndexSet &named_indices);

  const ColorMap &used_colors() const;

  Color operator()(const AbstractTensor &tensor);
  Color operator()(const BraGroup &group);
  Color operator()(const KetGroup &group);
  Color operator()(const AuxGroup &group);
  Color operator()(const ParticleGroup &group);
  Color operator()(const Index &idx);
  Color operator()(const ProtoBundle &bundle);

 private:
  ColorMap used_colors_;
  const TensorNetworkV2::NamedIndexSet &named_indices_;

  Color to_color(std::size_t color) const;

  template <typename T>
  Color ensure_uniqueness(Color color, const T &val) {
    auto it = used_colors_.find(color);
    while (it != used_colors_.end() && !may_have_same_color(it->second, val)) {
      // Color collision: val was computed to have the same color
      // as another object, but these objects do not compare equal (for
      // the purpose of color assigning).
      // -> Need to modify color until conflict is resolved.
      color++;
      it = used_colors_.find(color);
    }

    if (it == used_colors_.end()) {
      // We have not yet seen this color before -> add it to cache
      if constexpr (std::is_same_v<T, AbstractTensor> ||
                    std::is_same_v<T, ProtoBundle>) {
        used_colors_[color] = &val;
      } else {
        used_colors_[color] = val;
      }
    }

    return color;
  }

  bool may_have_same_color(const VertexData &data,
                           const AbstractTensor &tensor);
  bool may_have_same_color(const VertexData &data, const BraGroup &group);
  bool may_have_same_color(const VertexData &data, const KetGroup &group);
  bool may_have_same_color(const VertexData &data, const AuxGroup &group);
  bool may_have_same_color(const VertexData &data, const ParticleGroup &group);
  bool may_have_same_color(const VertexData &data, const Index &idx);
  bool may_have_same_color(const VertexData &data, const ProtoBundle &bundle);
};

}  // namespace sequant

#endif
